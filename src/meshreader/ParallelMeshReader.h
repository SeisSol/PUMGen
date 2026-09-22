// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_MESHREADER_PARALLELMESHREADER_H_
#define PUMGEN_SRC_MESHREADER_PARALLELMESHREADER_H_

#include <mpi.h>

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

#include "helper/Distributor.h"
#include "third_party/MPITraits.h"
#include "utils/logger.h"

template <class R> class ParallelMeshReader {
  private:
  // Some variables that are required by all processes
  std::size_t m_nVertices;
  std::size_t m_nElements;
  std::size_t m_nBoundaries;

  protected:
  const MPI_Comm m_comm;

  int m_rank;
  int m_nProcs;

  R m_serialReader;

  public:
  constexpr static std::size_t Order = 1;
  constexpr static std::size_t Dim = 3;
  constexpr static bool SupportsIdentify = false;
  ParallelMeshReader(MPI_Comm comm = MPI_COMM_WORLD)
      : m_nVertices(0), m_nElements(0), m_nBoundaries(0), m_comm(comm) {
    init();
  }

  /**
   * @param meshFile Only required by rank 0
   */
  ParallelMeshReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : m_nVertices(0), m_nElements(0), m_nBoundaries(0), m_comm(comm) {
    init();
    open(meshFile);
  }

  virtual ~ParallelMeshReader() {}

  void open(const char* meshFile) {
    std::array<std::size_t, 3> vars;

    if (m_rank == 0) {
      m_serialReader.open(meshFile);

      vars[0] = m_serialReader.nVertices();
      vars[1] = m_serialReader.nElements();
      vars[2] = m_serialReader.nBoundaries();
    }

    MPI_Bcast(vars.data(), 3, tndm::mpi_type_t<std::size_t>(), 0, m_comm);

    m_nVertices = vars[0];
    m_nElements = vars[1];
    m_nBoundaries = vars[2];
  }

  std::size_t nVertices() const { return m_nVertices; }

  std::size_t nElements() const { return m_nElements; }

  /**
   * @return Number of boundary faces
   */
  std::size_t nBoundaries() const { return m_nBoundaries; }

  /**
   * Reads all vertices. Each process gets the chunk given by getChunksize/getChunksum.
   *
   * This is a collective operation.
   *
   * @todo Only 3 dimensional meshes are supported
   */
  void readVertices(double* vertices) {
    distributeChunks(m_nVertices, 3, vertices, "vertices",
                     [&](std::size_t start, std::size_t count, double* buffer) {
                       m_serialReader.readVertices(start, count, buffer);
                     });
  }

  /**
   * Reads all elements. Each process gets the chunk given by getChunksize/getChunksum.
   *
   * This is a collective operation.
   *
   * @todo Only tetrahedral meshes are supported
   */
  virtual void readElements(std::size_t* elements) {
    distributeChunks(m_nElements, 4, elements, "elements",
                     [&](std::size_t start, std::size_t count, std::size_t* buffer) {
                       m_serialReader.readElements(start, count, buffer);
                     });
  }

  protected:
  static int messageSize(std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      logError() << "MPI message too large:" << count;
    }
    return static_cast<int>(count);
  }

  private:
  /**
   * Reads the items on rank 0 and sends each process its chunk.
   *
   * @param read Callable (start, count, buffer) reading <code>count</code> items starting at
   * <code>start</code>
   */
  template <typename T, typename ReadF>
  void distributeChunks(std::size_t total, std::size_t valuesPerItem, T* local, const char* name,
                        ReadF&& read) {
    const std::size_t localCount = getChunksize(total, m_rank, m_nProcs);

    if (m_rank != 0) {
      MPI_Recv(local, messageSize(localCount * valuesPerItem), tndm::mpi_type_t<T>(), 0, 0, m_comm,
               MPI_STATUS_IGNORE);
      return;
    }

    // two send buffers, so that reading the next chunk overlaps with sending the previous one
    const std::size_t maxCount = getChunksize(total, 0, m_nProcs);
    std::array<std::vector<T>, 2> buffers;
    std::array<MPI_Request, 2> requests = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    for (int rank = 1; rank < m_nProcs; ++rank) {
      const auto slot = rank % 2;
      MPI_Wait(&requests[slot], MPI_STATUS_IGNORE);
      buffers[slot].resize(maxCount * valuesPerItem);

      const std::size_t count = getChunksize(total, rank, m_nProcs);
      logInfo() << "Reading" << name << "for rank" << rank << "of" << m_nProcs;
      read(getChunksum(total, rank, m_nProcs), count, buffers[slot].data());
      MPI_Isend(buffers[slot].data(), messageSize(count * valuesPerItem), tndm::mpi_type_t<T>(),
                rank, 0, m_comm, &requests[slot]);
    }

    logInfo() << "Reading" << name << "for rank 0 of" << m_nProcs;
    read(0, localCount, local);
    MPI_Waitall(2, requests.data(), MPI_STATUSES_IGNORE);
  }

  /**
   * Initialize some parameters
   */
  void init() {
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_nProcs);
  }
};

#endif // PUMGEN_SRC_MESHREADER_PARALLELMESHREADER_H_
