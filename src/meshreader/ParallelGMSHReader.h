// SPDX-FileCopyrightText: 2020 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_
#define PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_

#include "GMSHBuilder.h"
#include "third_party/MPITraits.h"
#include "utils/logger.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "helper/Distributor.h"
#include "helper/MPIConvenience.h"

namespace puml {

template <typename P, std::size_t OrderP> class ParallelGMSHReader {
  public:
  constexpr static std::size_t Order = OrderP;
  constexpr static std::size_t Dim = 3;
  using bc_t = std::array<int, Dim + 1u>;
  constexpr static std::size_t Facet2Nodes[Dim + 1][Dim] = {
      {1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  constexpr static int BoundaryConditionOffset = 100;

  explicit ParallelGMSHReader(MPI_Comm comm = MPI_COMM_WORLD) : comm_(comm) {}

  void open(char const* meshFile) {
    int rank;
    MPI_Comm_rank(comm_, &rank);

    if (rank == 0) {
      P parser(&builder_);
      bool ok = parser.parseFile(meshFile);
      if (!ok) {
        logError() << meshFile << std::endl << parser.getErrorMessage();
      }
      builder_.postprocess();
      convertBoundaryConditions();

      nVertices_ = builder_.vertices.size();
      nElements_ = builder_.elements.size();
      hasIdentify_ = builder_.identify.empty() ? 0 : 1;
    }

    MPI_Bcast(&nVertices_, 1, tndm::mpi_type_t<decltype(nVertices_)>(), 0, comm_);
    MPI_Bcast(&nElements_, 1, tndm::mpi_type_t<decltype(nElements_)>(), 0, comm_);
    MPI_Bcast(&hasIdentify_, 1, tndm::mpi_type_t<decltype(hasIdentify_)>(), 0, comm_);
  }

  [[nodiscard]] std::size_t nVertices() const { return nVertices_; }
  [[nodiscard]] std::size_t nElements() const { return nElements_; }
  void readElements(std::size_t* elements) const {
    static_assert(sizeof(typename GMSHBuilder<Dim, Order>::element_t) ==
                  nodeCount(Dim, Order) * sizeof(std::size_t));
    scatter(builder_.elements.data()->data(), elements, nElements(), nodeCount(Dim, Order));
  }
  void readVertices(double* vertices) const {
    static_assert(sizeof(typename GMSHBuilder<Dim, Order>::vertex_t) == Dim * sizeof(double));
    scatter(builder_.vertices.data()->data(), vertices, nVertices(), Dim);
  }
  void readBoundaries(int* boundaries) const {
    static_assert(sizeof(bc_t) == (Dim + 1) * sizeof(int));
    scatter(bcs_.data()->data(), boundaries, nElements(), Dim + 1);
  }
  void readGroups(int* groups) const { scatter(builder_.groups.data(), groups, nElements(), 1); }

  constexpr static bool SupportsIdentify = true;
  bool hasIdentify() const { return hasIdentify_ != 0; }
  void readIdentify(std::size_t* vertices) const {
    scatter(builder_.identify.data(), vertices, nVertices(), 1);
  }

  private:
  /**
   * GMSH stores boundary conditions on a surface mesh whereas SeisSol expects
   * boundary conditions to be stored per element. In the following we convert
   * the boundary condition representation by matching the 4 faces of an element
   * with the surface mesh.
   */
  void convertBoundaryConditions() {
    using Face = std::array<std::size_t, Dim>;
    struct BoundaryFace {
      Face vertices;
      int bc;
    };
    const auto sortedFace = [](Face face) {
      std::sort(face.begin(), face.end());
      return face;
    };

    // the surface mesh as sorted vertex triples; the vertices of a facet come first
    std::vector<std::uint8_t> onBoundary(builder_.vertices.size(), 0);
    std::vector<BoundaryFace> boundaryFaces(builder_.facets.size());
    for (std::size_t fctNo = 0; fctNo < builder_.facets.size(); ++fctNo) {
      Face face{};
      for (std::size_t i = 0; i < Dim; ++i) {
        face[i] = builder_.facets[fctNo][i];
        onBoundary[face[i]] = 1;
      }
      boundaryFaces[fctNo] = {sortedFace(face), builder_.bcs[fctNo]};
    }
    const auto byVertices = [](const BoundaryFace& a, const BoundaryFace& b) {
      return a.vertices < b.vertices;
    };
    std::sort(boundaryFaces.begin(), boundaryFaces.end(), byVertices);

    const auto nElements = builder_.elements.size();
    bcs_.assign(nElements, bc_t{});
    for (std::size_t elNo = 0; elNo < nElements; ++elNo) {
      const auto& element = builder_.elements[elNo];
      // only faces with all vertices on the surface mesh can be part of it
      std::array<bool, Dim + 1> vertexOnBoundary{};
      std::size_t verticesOnBoundary = 0;
      for (std::size_t i = 0; i < Dim + 1; ++i) {
        vertexOnBoundary[i] = onBoundary[element[i]] != 0;
        verticesOnBoundary += vertexOnBoundary[i] ? 1 : 0;
      }
      if (verticesOnBoundary < Dim) {
        continue;
      }
      for (std::size_t localFctNo = 0; localFctNo < Dim + 1; ++localFctNo) {
        Face face{};
        bool candidate = true;
        for (std::size_t i = 0; i < Dim; ++i) {
          const auto localVertex = Facet2Nodes[localFctNo][i];
          candidate = candidate && vertexOnBoundary[localVertex];
          face[i] = element[localVertex];
        }
        if (!candidate) {
          continue;
        }
        const BoundaryFace key{sortedFace(face), 0};
        const auto [first, last] =
            std::equal_range(boundaryFaces.begin(), boundaryFaces.end(), key, byVertices);
        if (last - first > 1) {
          logError() << "A face of an element exists multiple times in the surface mesh.";
        }
        if (first != last) {
          bcs_[elNo][localFctNo] = adjustBoundaryCondition(first->bc);
        }
      }
    }
  }

  /**
   * Boundary conditions in SeisSol used to start at 100, e.g. 101 = free
   * surface. In the hdf5 format one starts counting from 0, e.g. 1 = free
   * surface. In order to be compatible with legacy gmsh scripts we subtract
   * 100 here if the boundary condition is larger or equal 100.
   */
  [[nodiscard]] int adjustBoundaryCondition(int bc) const {
    if (bc >= BoundaryConditionOffset) {
      bc -= BoundaryConditionOffset;
    }
    return bc;
  }

  template <typename T>
  void scatter(T const* sendbuf, T* recvbuf, std::size_t numElements,
               std::size_t numPerElement) const {
    int rank;
    int procs;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &procs);

    auto sendcounts = std::vector<std::size_t>(procs);
    auto displs = std::vector<std::size_t>(procs + 1);
    displs[0] = 0;
    for (int i = 0; i < procs; ++i) {
      sendcounts[i] = numPerElement * getChunksize(numElements, i, procs);
      displs[i + 1] = displs[i] + sendcounts[i];
    }

    auto recvcount = sendcounts[rank];
    largeScatterv(sendbuf, sendcounts.data(), displs.data(), tndm::mpi_type_t<T>(), recvbuf,
                  recvcount, tndm::mpi_type_t<T>(), 0, comm_);
  }

  MPI_Comm comm_;
  GMSHBuilder<Dim, Order> builder_;
  std::vector<bc_t> bcs_;
  std::size_t nVertices_ = 0;
  std::size_t nElements_ = 0;
  // only rank 0 parses the file, so all ranks need this flag for the collective reads and writes
  int hasIdentify_ = 0;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_
