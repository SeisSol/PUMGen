/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2017 Technical University of Munich
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @remark This class is taken from XdmfWriter
 * (https://github.com/TUM-I5/XdmfWriter)
 */

#ifndef PARALLEL_VERTEX_FILTER_H
#define PARALLEL_VERTEX_FILTER_H

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>
#include <array>

#include "utils/logger.h"

#include "third_party/MPITraits.h"

/**
 * Filters duplicate vertices in parallel
 */
class ParallelVertexFilter {
  private:
  /**
   * Compares 3D-vertex indices according to the vertices
   */
  class IndexedVertexComparator {
private:
    const std::vector<double>& m_vertices;

public:
    IndexedVertexComparator(const std::vector<double>& vertices) : m_vertices(vertices) {}

    bool operator()(std::size_t i, std::size_t j) {
      i *= 3;
      j *= 3;

      return (m_vertices[i] < m_vertices[j]) ||
             (m_vertices[i] == m_vertices[j] && m_vertices[i + 1] < m_vertices[j + 1]) ||
             (m_vertices[i] == m_vertices[j] && m_vertices[i + 1] == m_vertices[j + 1] &&
              m_vertices[i + 2] < m_vertices[j + 2]);
    }
  };

  private:
  /** The communicator we use */
  MPI_Comm m_comm;

  /** Our rank */
  int m_rank;

  /** #Processes */
  int m_numProcs;

  /** Global id after filtering */
  std::vector<std::size_t> m_globalIds;

  /** Number of local vertices after filtering */
  std::size_t m_numLocalVertices;

  /** Local vertices after filtering */
  std::vector<double> m_localVertices;

  public:
  ParallelVertexFilter(MPI_Comm comm = MPI_COMM_WORLD)
      : m_comm(comm), m_numLocalVertices(0) {
    MPI_Comm_rank(comm, &m_rank);
    MPI_Comm_size(comm, &m_numProcs);

    if (vertexType == MPI_DATATYPE_NULL) {
      MPI_Type_contiguous(3, MPI_DOUBLE, &vertexType);
      MPI_Type_commit(&vertexType);
    }
  }

  virtual ~ParallelVertexFilter() {
  }

  /**
   * @param vertices Vertices that should be filtered, must have the size
   * <code>numVertices * 3</code>
   */
  void filter(std::size_t numVertices, const std::vector<double>& vertices) {
    // Chop the last 4 bits to avoid numerical errors
    std::vector<double> roundVertices (numVertices * 3);
    removeRoundError(vertices, roundVertices);

    // Create indices and sort them locally
    std::vector<std::size_t> sortIndices(numVertices);
    createSortedIndices(roundVertices, sortIndices);

    // Select BUCKETS_PER_RANK-1 splitter elements
    std::array<double, BUCKETS_PER_RANK - 1> localSplitters;
#if 0 // Use omp only if we create a larger amount of buckets
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
#endif
    for (int i = 0; i < BUCKETS_PER_RANK - 1; i++) {
      unsigned long vrtxIndex = static_cast<unsigned long>(i) *
                                static_cast<unsigned long>(numVertices) /
                                static_cast<unsigned long>(BUCKETS_PER_RANK - 1);
      assert(vrtxIndex < numVertices);

      localSplitters[i] = roundVertices[sortIndices[vrtxIndex] * 3];
    }

    // Collect all splitter elements on rank 0
    std::vector<double> allSplitters;

    if (m_rank == 0) {
      allSplitters.resize(m_numProcs * (BUCKETS_PER_RANK - 1));
    }

    MPI_Gather(localSplitters.data(), BUCKETS_PER_RANK - 1, MPI_DOUBLE, allSplitters.data(), BUCKETS_PER_RANK - 1,
               MPI_DOUBLE, 0, m_comm);

    // Sort splitter elements
    if (m_rank == 0) {
      std::sort(allSplitters.begin(), allSplitters.end());
    }

    // Distribute splitter to all processes
    std::vector<double> splitters (m_numProcs - 1);

    if (m_rank == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < m_numProcs - 1; i++) {
        unsigned long spltIndex = (i + 1) * (BUCKETS_PER_RANK - 1);
        assert(spltIndex < static_cast<unsigned int>(m_numProcs * (BUCKETS_PER_RANK - 1)));

        splitters[i] = allSplitters[spltIndex];
      }
    }

    MPI_Bcast(splitters.data(), m_numProcs - 1, MPI_DOUBLE, 0, m_comm);

    // Determine the bucket for each vertex
    std::vector<std::size_t> bucket (numVertices);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < numVertices; i++) {
      auto ub = std::upper_bound(splitters.begin(), splitters.end(), roundVertices[i * 3]);

      bucket[i] = std::distance(splitters.begin(), ub);
    }

    // Determine the (local and total) bucket size
    std::vector<int> bucketSize (m_numProcs);
    for (std::size_t i = 0; i < numVertices; i++) {
      ++bucketSize[bucket[i]];
    }

    // Tell all processes what we are going to send them
    std::vector<int> recvSize(m_numProcs);

    MPI_Alltoall(bucketSize.data(), 1, MPI_INT, recvSize.data(), 1, MPI_INT, m_comm);

    std::size_t numSortVertices = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : numSortVertices)
#endif
    for (int i = 0; i < m_numProcs; i++) {
      numSortVertices += recvSize[i];
    }

    // Create sorted send buffer
    std::vector<double> sendVertices (3 * numVertices);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < numVertices; i++) {
      std::copy_n(vertices.begin() + sortIndices[i] * 3, 3, sendVertices.begin() + i * 3);
    }

    // Allocate buffer for the vertices and exchange them
    std::vector<double> sortVertices (3 * numSortVertices);

    std::vector<int> sDispls (m_numProcs);
    std::vector<int> rDispls (m_numProcs);
    sDispls[0] = 0;
    rDispls[0] = 0;
    for (int i = 1; i < m_numProcs; i++) {
      sDispls[i] = sDispls[i - 1] + bucketSize[i - 1];
      rDispls[i] = rDispls[i - 1] + recvSize[i - 1];
    }
    MPI_Alltoallv(sendVertices.data(), bucketSize.data(), sDispls.data(), vertexType, sortVertices.data(), recvSize.data(), rDispls.data(),
                  vertexType, m_comm);

    // Chop the last 4 bits to avoid numerical errors
    roundVertices.resize(numSortVertices * 3);
    removeRoundError(sortVertices, roundVertices);

    // Create indices and sort them (such that the vertices are sorted)
    std::vector<std::size_t> sortSortIndices (numSortVertices);
    createSortedIndices(roundVertices, sortSortIndices);

    // Initialize the global ids we send back to the other processors
    std::vector<std::size_t> gids (numSortVertices);

    if (numSortVertices > 0) {
      gids[sortSortIndices[0]] = 0;
      for (std::size_t i = 1; i < numSortVertices; i++) {
        if (equals(&sortVertices[sortSortIndices[i - 1] * 3],
                   &sortVertices[sortSortIndices[i] * 3])) {
          gids[sortSortIndices[i]] = gids[sortSortIndices[i - 1]];
                   }
        else {
          gids[sortSortIndices[i]] = gids[sortSortIndices[i - 1]] + 1;
        }
      }
    }

    // Create the local vertices list
    if (numSortVertices > 0) {
      m_numLocalVertices = gids[sortSortIndices[numSortVertices - 1]] + 1;
    }
    else {
      m_numLocalVertices = 0;
    }

    m_localVertices.resize(m_numLocalVertices * 3);
    for (std::size_t i = 0; i < numSortVertices; i++) {
      std::copy_n(sortVertices.begin() + i * 3, 3, m_localVertices.begin() + gids[i] * 3);
    }

    // Get the vertices offset
    std::size_t offset = m_numLocalVertices;
    MPI_Scan(MPI_IN_PLACE, &offset, 1, tndm::mpi_type_t<std::size_t>(), MPI_SUM, m_comm);
    offset -= m_numLocalVertices;

    // Add offset to the global ids
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < numSortVertices; i++) {
      gids[i] += offset;
    }

    // Send result back
    std::vector<std::size_t> globalIds (numVertices);
    MPI_Alltoallv(gids.data(), recvSize.data(), rDispls.data(), tndm::mpi_type_t<std::size_t>(), globalIds.data(), bucketSize.data(), sDispls.data(),
                  tndm::mpi_type_t<std::size_t>(), m_comm);

    // Assign the global ids to the correct vertices
    m_globalIds.resize(numVertices);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < numVertices; i++) {
      m_globalIds[sortIndices[i]] = globalIds[i];
    }
  }

  /**
   * @return The list of the global identifiers after filtering
   */
  const std::vector<std::size_t>& globalIds() const { return m_globalIds; }

  /**
   * @return Number of vertices this process is responsible for after filtering
   */
  std::size_t numLocalVertices() const { return m_numLocalVertices; }

  /**
   * @return The list of vertices this process is responsible for after
   * filtering
   */
  const std::vector<double>& localVertices() const { return m_localVertices; }

  private:
  /**
   * Removes round errors of double values by setting the last 4 bits
   * (of the significand) to zero.
   *
   * @warning Only works if <code>value</code> ist not nan or infinity
   * @todo This should work for arbitrary precision
   */
  static double removeRoundError(double value) {
    static const uint64_t mask = ~0xF;

    union FloatUnion {
      double f;
      uint64_t bits;
    };

    FloatUnion result;
    result.f = value;

    result.bits &= mask;

    return result.f;
  }

  /**
   * Removes the round errors using {@link removeRoundError(double)}
   *
   * @param values The list of floating point values
   * @param count Number of values
   * @param[out] roundValues The list of rounded values
   *  (the caller is responsible for allocating the memory)
   */
  static void removeRoundError(const std::vector<double>& values, std::vector<double>& roundValues) {
    assert(values.size() == roundValues.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < roundValues.size(); i++) {
      roundValues[i] = removeRoundError(values[i]);
    }
  }

  /**
   * Creates the list of sorted indices for the vertices.
   * The caller is responsible for allocating the memory.
   */
  static void createSortedIndices(const std::vector<double>& vertices,
                                  std::vector<std::size_t>& sortedIndices) {

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < sortedIndices.size(); i++) {
      sortedIndices[i] = i;
    }

    IndexedVertexComparator comparator(vertices);
    std::sort(sortedIndices.begin(), sortedIndices.end(), comparator);
  }

  /**
   * Compares to vertices for equality
   * Assumes that the rounding errors are removed.
   */
  static bool equals(const double* vertexA, const double* vertexB) {
    return vertexA[0] == vertexB[0] && vertexA[1] == vertexB[1] && vertexA[2] == vertexB[2];
  }

  /** MPI data type consisting of three doubles */
  static MPI_Datatype vertexType;

  /** The total buckets we create is <code>BUCKETS_PER_RANK * numProcs</code> */
  constexpr static int BUCKETS_PER_RANK = 8;
};

#endif // PARALLEL_VERTEX_FILTER_H
