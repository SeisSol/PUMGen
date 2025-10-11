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

std::vector<double> calculateInsphere(const std::vector<std::size_t>& connectivity,
                                      const std::vector<double>& geometry, MPI_Comm comm) {
  int commsize;
  int commrank;

  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);

  MPI_Datatype vertexType;

  MPI_Type_contiguous(3, MPI_DOUBLE, &vertexType);
  MPI_Type_commit(&vertexType);

  std::vector<std::size_t> outrequests(commsize);
  std::vector<std::size_t> inrequests(commsize);

  std::size_t localVertices = geometry.size() / 3;

  std::vector<std::size_t> vertexDist(commsize + 1);

  MPI_Allgather(&localVertices, 1, tndm::mpi_type_t<std::size_t>(), vertexDist.data() + 1, 1,
                tndm::mpi_type_t<std::size_t>(), comm);

  for (std::size_t i = 0; i < commsize; ++i) {
    vertexDist[i + 1] += vertexDist[i];
  }

  std::vector<std::unordered_map<std::size_t, std::size_t>> outidxmap(commsize);

  for (const auto& vertex : connectivity) {
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

  std::vector<std::size_t> outdisp(commsize + 1);
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

  std::vector<double> outvertices(3 * connectivity.size());

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

  std::vector<std::size_t> counter(commsize);
  std::vector<double> inspheres(connectivity.size() / 4);

  for (std::size_t i = 0; i < connectivity.size() / 4; ++i) {
    std::array<std::array<double, 3>, 4> vertices;
    for (int j = 0; j < 4; ++j) {
      auto vertex = connectivity[i * 4 + j];
      auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
      auto position = std::distance(vertexDist.begin(), itPosition) - 1;
      auto localVertex = vertex - vertexDist[position];
      if (position == commrank) {
        vertices[j][0] = geometry[localVertex * 3 + 0];
        vertices[j][1] = geometry[localVertex * 3 + 1];
        vertices[j][2] = geometry[localVertex * 3 + 2];
      } else {
        auto transferidx = outidxmap[position][localVertex] + outdisp[position];
        vertices[j][0] = outvertices[transferidx * 3 + 0];
        vertices[j][1] = outvertices[transferidx * 3 + 1];
        vertices[j][2] = outvertices[transferidx * 3 + 2];
      }
    }
    double a11 = vertices[1][0] - vertices[0][0];
    double a12 = vertices[1][1] - vertices[0][1];
    double a13 = vertices[1][2] - vertices[0][2];
    double a21 = vertices[2][0] - vertices[0][0];
    double a22 = vertices[2][1] - vertices[0][1];
    double a23 = vertices[2][2] - vertices[0][2];
    double a31 = vertices[3][0] - vertices[0][0];
    double a32 = vertices[3][1] - vertices[0][1];
    double a33 = vertices[3][2] - vertices[0][2];
    double b11 = vertices[2][0] - vertices[1][0];
    double b12 = vertices[2][1] - vertices[1][1];
    double b13 = vertices[2][2] - vertices[1][2];
    double b21 = vertices[3][0] - vertices[1][0];
    double b22 = vertices[3][1] - vertices[1][1];
    double b23 = vertices[3][2] - vertices[1][2];
    double det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 -
                 a12 * a21 * a33 - a11 * a23 * a32;
    double gram = std::abs(det);

    double a1a2x1 = a12 * a23 - a13 * a22;
    double a1a2x2 = a13 * a21 - a11 * a23;
    double a1a2x3 = a11 * a22 - a12 * a21;
    double a1a3x1 = a12 * a33 - a13 * a32;
    double a1a3x2 = a13 * a31 - a11 * a33;
    double a1a3x3 = a11 * a32 - a12 * a31;
    double a2a3x1 = a22 * a33 - a23 * a32;
    double a2a3x2 = a23 * a31 - a21 * a33;
    double a2a3x3 = a21 * a32 - a22 * a31;
    double b1b2x1 = b12 * b23 - b13 * b22;
    double b1b2x2 = b13 * b21 - b11 * b23;
    double b1b2x3 = b11 * b22 - b12 * b21;

    double faces = std::sqrt(a1a2x1 * a1a2x1 + a1a2x2 * a1a2x2 + a1a2x3 * a1a2x3) +
                   std::sqrt(a1a3x1 * a1a3x1 + a1a3x2 * a1a3x2 + a1a3x3 * a1a3x3) +
                   std::sqrt(a2a3x1 * a2a3x1 + a2a3x2 * a2a3x2 + a2a3x3 * a2a3x3) +
                   std::sqrt(b1b2x1 * b1b2x1 + b1b2x2 * b1b2x2 + b1b2x3 * b1b2x3);

    inspheres[i] = gram / faces;
  }

  MPI_Type_free(&vertexType);

  return inspheres;
}
