// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "InsphereCalculator.h"
#include "third_party/MPITraits.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <mpi.h>
#include <unordered_map>

#include "MPIConvenience.h"

CellVertices::CellVertices(const std::vector<puml::CellType>& cellTypes,
                           const std::vector<std::uint64_t>& connectivity,
                           const std::vector<double>& geometry, MPI_Comm comm)
    : cellTypes(cellTypes), connectivity(connectivity), geometry(geometry) {
  if (!cellTypes.empty() &&
      std::any_of(cellTypes.begin(), cellTypes.end(),
                  [&](puml::CellType type) { return type != cellTypes.front(); })) {
    firstVertex.resize(cellTypes.size() + 1, 0);
    for (std::size_t cell = 0; cell < cellTypes.size(); ++cell) {
      firstVertex[cell + 1] = firstVertex[cell] + puml::shapeOf(cellTypes[cell]).vertexCount;
    }
  }

  int commsize;
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);

  MPI_Datatype vertexType;

  MPI_Type_contiguous(3, MPI_DOUBLE, &vertexType);
  MPI_Type_commit(&vertexType);

  std::vector<std::size_t> outrequests(commsize);
  std::vector<std::size_t> inrequests(commsize);

  std::size_t localVertices = geometry.size() / 3;

  vertexDist.resize(commsize + 1);

  MPI_Allgather(&localVertices, 1, tndm::mpi_type_t<std::size_t>(), vertexDist.data() + 1, 1,
                tndm::mpi_type_t<std::size_t>(), comm);

  for (std::size_t i = 0; i < commsize; ++i) {
    vertexDist[i + 1] += vertexDist[i];
  }

  outidxmap.resize(commsize);

  for (std::size_t node = 0; node < connectivity.size(); ++node) {
    const auto vertex = connectivity[node];
    auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
    auto position = std::distance(vertexDist.begin(), itPosition) - 1;
    auto localVertex = vertex - vertexDist[position];

    // only transfer each vertex coordinate once, and only if it's not already on the same rank
    if (position != commrank &&
        outidxmap[position].find(localVertex) == outidxmap[position].end()) {
      outidxmap[position][localVertex] = outrequests[position];
      ++outrequests[position];
    }
  }

  MPI_Alltoall(outrequests.data(), 1, tndm::mpi_type_t<std::size_t>(), inrequests.data(), 1,
               tndm::mpi_type_t<std::size_t>(), comm);

  outdisp.resize(commsize + 1);
  std::vector<std::size_t> indisp(commsize + 1);
  for (std::size_t i = 1; i < commsize + 1; ++i) {
    outdisp[i] = outdisp[i - 1] + outrequests[i - 1];
    indisp[i] = indisp[i - 1] + inrequests[i - 1];
  }

  std::vector<std::size_t> inidx(indisp[commsize]);

  {
    std::vector<std::size_t> outidx(outdisp[commsize]);

    std::vector<std::size_t> counter(commsize);

    for (std::size_t i = 0; i < commsize; ++i) {
      for (const auto& [localVertex, j] : outidxmap[i]) {
        outidx[j + outdisp[i]] = localVertex;
      }
    }

    // transfer indices

    sparseAlltoallv(outidx.data(), outrequests.data(), outdisp.data(),
                    tndm::mpi_type_t<std::size_t>(), inidx.data(), inrequests.data(), indisp.data(),
                    tndm::mpi_type_t<std::size_t>(), comm);
  }

  outvertices.resize(3 * outdisp[commsize]);

  {
    std::vector<double> invertices(3 * indisp[commsize]);

    for (std::size_t i = 0; i < inidx.size(); ++i) {
      for (int j = 0; j < 3; ++j) {
        invertices[3 * i + j] = geometry[3 * inidx[i] + j];
      }
    }

    // transfer vertices

    sparseAlltoallv(invertices.data(), inrequests.data(), indisp.data(), vertexType,
                    outvertices.data(), outrequests.data(), outdisp.data(), vertexType, comm);
  }

  MPI_Type_free(&vertexType);
}

std::array<std::array<double, 3>, puml::MaxCellVertices>
CellVertices::operator()(std::size_t cell) const {
  std::array<std::array<double, 3>, puml::MaxCellVertices> vertices{};
  const auto count = puml::shapeOf(cellTypes[cell]).vertexCount;
  const auto first = firstVertex.empty() ? cell * count : firstVertex[cell];
  for (std::size_t j = 0; j < count; ++j) {
    auto vertex = connectivity[first + j];
    auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
    auto position = std::distance(vertexDist.begin(), itPosition) - 1;
    auto localVertex = vertex - vertexDist[position];
    if (position == commrank) {
      vertices[j][0] = geometry[localVertex * 3 + 0];
      vertices[j][1] = geometry[localVertex * 3 + 1];
      vertices[j][2] = geometry[localVertex * 3 + 2];
    } else {
      auto transferidx = outidxmap[position].at(localVertex) + outdisp[position];
      vertices[j][0] = outvertices[transferidx * 3 + 0];
      vertices[j][1] = outvertices[transferidx * 3 + 1];
      vertices[j][2] = outvertices[transferidx * 3 + 2];
    }
  }
  return vertices;
}

std::vector<double> calculateInsphere(const CellVertices& cells) {
  std::vector<double> inspheres(cells.numCells());

  for (std::size_t i = 0; i < cells.numCells(); ++i) {
    const auto& shape = puml::shapeOf(cells.type(i));
    const auto vertices = cells(i);
    std::array<double, 3> center{};
    for (std::size_t v = 0; v < shape.vertexCount; ++v) {
      for (int d = 0; d < 3; ++d) {
        center[d] += vertices[v][d] / static_cast<double>(shape.vertexCount);
      }
    }

    // the faces split into triangles; the volume adds up the pyramids from the centre to them
    double volume = 0;
    double area = 0;
    for (std::size_t f = 0; f < shape.faceCount; ++f) {
      const auto& face = shape.faceVertices[f];
      const auto& a = vertices[face[0]];
      for (std::size_t t = 1; t + 1 < shape.faceVertexCount[f]; ++t) {
        const auto& b = vertices[face[t]];
        const auto& c = vertices[face[t + 1]];
        const std::array<double, 3> u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const std::array<double, 3> w{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const std::array<double, 3> normal{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                           u[0] * w[1] - u[1] * w[0]};
        area +=
            std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) / 2;
        volume += (normal[0] * (a[0] - center[0]) + normal[1] * (a[1] - center[1]) +
                   normal[2] * (a[2] - center[2])) /
                  6;
      }
    }

    inspheres[i] = 3 * std::abs(volume) / area;
  }

  return inspheres;
}
