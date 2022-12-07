#ifndef PARALLELGMSHREADER_20201014_H
#define PARALLELGMSHREADER_20201014_H

#include "GMSHBuilder.h"
#include "third_party/MPITraits.h"
#include "utils/logger.h"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace puml {

enum class gmsh_version { v2, v4 };

template <gmsh_version V> class ParallelGMSHReader {
  public:
  constexpr static std::size_t Dim = 3;
  using bc_t = std::array<int, Dim + 1u>;
  constexpr static std::size_t Facet2Nodes[Dim + 1][Dim] = {
      {1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  constexpr static int BoundaryConditionOffset = 100;

  explicit ParallelGMSHReader(MPI_Comm comm = MPI_COMM_WORLD) : comm_(comm) {}

  void open(char const* meshFile);
  [[nodiscard]] unsigned int nVertices() const { return nVertices_; }
  [[nodiscard]] unsigned int nElements() const { return nElements_; }
  void readElements(int* elements) const {
    static_assert(sizeof(GMSHBuilder<Dim>::element_t) == (Dim + 1) * sizeof(int));
    scatter(builder_.elements.data()->data(), elements, nElements(), Dim + 1);
  }
  void readVertices(double* vertices) const {
    static_assert(sizeof(GMSHBuilder<Dim>::vertex_t) == Dim * sizeof(double));
    scatter(builder_.vertices.data()->data(), vertices, nVertices(), Dim);
  }
  void readBoundaries(int* boundaries) const {
    static_assert(sizeof(bc_t) == (Dim + 1) * sizeof(int));
    scatter(bcs_.data()->data(), boundaries, nElements(), Dim + 1);
  }
  void readGroups(int* groups) const { scatter(builder_.groups.data(), groups, nElements(), 1); }

  private:
  /**
   * GMSH stores boundary conditions on a surface mesh whereas SeisSol expects
   * boundary conditions to be stored per element. In the following we convert
   * the boundary condition representation by matching the 4 faces of an element
   * with the surface mesh.
   */
  void convertBoundaryConditions() {
    const auto nVerts = builder_.vertices.size();
    const auto nFacets = builder_.facets.size();
    const auto nElements = builder_.elements.size();

    bcs_.clear();
    bcs_.resize(nElements);

    constexpr auto nFacetsPerTet = Dim + 1u;
    constexpr auto nVertsPerTri = Dim;

    // vertex to triangle map
    std::vector<std::vector<std::size_t>> vertex2facets(nVerts);
    for (std::size_t fctNo = 0; fctNo < nFacets; ++fctNo) {
      for (auto const& vtxNo : builder_.facets[fctNo]) {
        vertex2facets[vtxNo].push_back(fctNo);
      }
    }
    for (auto& v2f : vertex2facets) {
      std::sort(v2f.begin(), v2f.end());
    }

    for (unsigned elNo = 0; elNo < nElements; ++elNo) {
      for (unsigned localFctNo = 0; localFctNo < nFacetsPerTet; ++localFctNo) {
        unsigned nodes[nVertsPerTri];
        for (unsigned localNo = 0; localNo < nVertsPerTri; ++localNo) {
          const auto localNoElement = Facet2Nodes[localFctNo][localNo];
          nodes[localNo] = builder_.elements[elNo][localNoElement];
        }
        std::vector<unsigned> intersect[nVertsPerTri - 1];
        std::set_intersection(vertex2facets[nodes[0]].begin(), vertex2facets[nodes[0]].end(),
                              vertex2facets[nodes[1]].begin(), vertex2facets[nodes[1]].end(),
                              std::back_inserter(intersect[0]));
        for (unsigned node = 2; node < nVertsPerTri; ++node) {
          std::set_intersection(intersect[node - 2].begin(), intersect[node - 2].end(),
                                vertex2facets[nodes[node]].begin(),
                                vertex2facets[nodes[node]].end(),
                                std::back_inserter(intersect[node - 1]));
        }
        if (!intersect[nVertsPerTri - 2].empty()) {
          if (intersect[nVertsPerTri - 2].size() > 1) {
            logError() << "A face of an element exists multiple times in the "
                          "surface mesh.";
          }
          const auto fctNo = intersect[nVertsPerTri - 2][0];
          bcs_[elNo][localFctNo] = adjustBoundaryCondition(builder_.bcs[fctNo]);
        } else {
          bcs_[elNo][localFctNo] = 0;
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

    const std::size_t chunkSize = (numElements + procs - 1) / procs;
    auto sendcounts = std::vector<int>(procs);
    auto displs = std::vector<int>(procs + 1);
    std::fill(sendcounts.begin(), sendcounts.end(), numPerElement * chunkSize);
    sendcounts.back() = numPerElement * (numElements - (procs - 1u) * chunkSize);
    displs[0] = 0;
    for (int i = 0; i < procs; ++i) {
      displs[i + 1] = displs[i] + sendcounts[i];
    }

    auto recvcount = sendcounts[rank];
    MPI_Scatterv(sendbuf, sendcounts.data(), displs.data(), tndm::mpi_type_t<T>(), recvbuf,
                 recvcount, tndm::mpi_type_t<T>(), 0, comm_);
  }

  MPI_Comm comm_;
  GMSHBuilder<Dim> builder_;
  std::vector<bc_t> bcs_;
  unsigned int nVertices_ = 0;
  unsigned int nElements_ = 0;
};

} // namespace puml

#endif // PARALLELGMSHREADER_20201014_H
