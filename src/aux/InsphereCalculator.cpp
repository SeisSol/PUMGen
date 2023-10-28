#include "InsphereCalculator.h"
#include "third_party/MPITraits.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <mpi.h>

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

  std::vector<int> outrequests(commsize);
  std::vector<int> inrequests(commsize);

  std::size_t localVertices = geometry.size() / 3;

  std::vector<std::size_t> vertexDist(commsize + 1);

  MPI_Allgather(&localVertices, 1, tndm::mpi_type_t<std::size_t>(), vertexDist.data() + 1, 1,
                tndm::mpi_type_t<std::size_t>(), comm);

  for (std::size_t i = 0; i < commsize; ++i) {
    vertexDist[i + 1] += vertexDist[i];
  }

  for (const auto& vertex : connectivity) {
    auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
    auto position = std::distance(vertexDist.begin(), itPosition) - 1;
    ++outrequests[position];
  }

  MPI_Alltoall(outrequests.data(), 1, MPI_INT, inrequests.data(), 1, MPI_INT, comm);

  std::vector<std::size_t> outdisp(commsize + 1);
  std::vector<std::size_t> indisp(commsize + 1);
  for (std::size_t i = 1; i < commsize + 1; ++i) {
    outdisp[i] = outdisp[i - 1] + outrequests[i - 1];
    indisp[i] = indisp[i - 1] + inrequests[i - 1];
  }

  std::vector<std::size_t> inidx(indisp[commsize]);

  {
    std::vector<std::size_t> outidx(connectivity.size());

    std::vector<std::size_t> counter(commsize);

    for (std::size_t i = 0; i < connectivity.size(); ++i) {
      auto vertex = connectivity[i];
      auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
      auto position = std::distance(vertexDist.begin(), itPosition) - 1;
      outidx[counter[position] + outdisp[position]] = vertex - vertexDist[position];
      ++counter[position];
    }

    // transfer indices

    sparseAlltoallv(outidx.data(), outrequests.data(), outdisp.data(),
                    tndm::mpi_type_t<std::size_t>(), inidx.data(), inrequests.data(), indisp.data(),
                    tndm::mpi_type_t<std::size_t>(), comm);
  }

  // this is a bit inefficient, since we don't look for duplicates... TODO: maybe improve
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
    std::array<std::size_t, 4> vertexIds{0, 0, 0, 0};
    for (int j = 0; j < 4; ++j) {
      auto vertex = connectivity[i * 4 + j];
      auto itPosition = std::upper_bound(vertexDist.begin(), vertexDist.end(), vertex);
      auto position = std::distance(vertexDist.begin(), itPosition) - 1;
      vertexIds[j] = counter[position] + outdisp[position];
      ++counter[position];
    }
    double a11 = outvertices[3 * vertexIds[1] + 0] - outvertices[3 * vertexIds[0] + 0];
    double a12 = outvertices[3 * vertexIds[1] + 1] - outvertices[3 * vertexIds[0] + 1];
    double a13 = outvertices[3 * vertexIds[1] + 2] - outvertices[3 * vertexIds[0] + 2];
    double a21 = outvertices[3 * vertexIds[2] + 0] - outvertices[3 * vertexIds[0] + 0];
    double a22 = outvertices[3 * vertexIds[2] + 1] - outvertices[3 * vertexIds[0] + 1];
    double a23 = outvertices[3 * vertexIds[2] + 2] - outvertices[3 * vertexIds[0] + 2];
    double a31 = outvertices[3 * vertexIds[3] + 0] - outvertices[3 * vertexIds[0] + 0];
    double a32 = outvertices[3 * vertexIds[3] + 1] - outvertices[3 * vertexIds[0] + 1];
    double a33 = outvertices[3 * vertexIds[3] + 2] - outvertices[3 * vertexIds[0] + 2];
    double b11 = outvertices[3 * vertexIds[2] + 0] - outvertices[3 * vertexIds[1] + 0];
    double b12 = outvertices[3 * vertexIds[2] + 1] - outvertices[3 * vertexIds[1] + 1];
    double b13 = outvertices[3 * vertexIds[2] + 2] - outvertices[3 * vertexIds[1] + 2];
    double b21 = outvertices[3 * vertexIds[3] + 0] - outvertices[3 * vertexIds[1] + 0];
    double b22 = outvertices[3 * vertexIds[3] + 1] - outvertices[3 * vertexIds[1] + 1];
    double b23 = outvertices[3 * vertexIds[3] + 2] - outvertices[3 * vertexIds[1] + 2];
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
