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
 */

#ifndef PARALLEL_MESH_READER_H
#define PARALLEL_MESH_READER_H

#include <mpi.h>

#include <array>
#include <vector>

#include "third_party/MPITraits.h"

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
   * Reads all vertices
   * Each process gets <code>nVertices() + processes - 1) / processes</code>
   * elements, except for the last, which gets the remaining elements
   *
   * This is a collective operation
   *
   * @todo Only 3 dimensional meshes are supported
   */
  void readVertices(double* vertices) {
    std::size_t chunkSize = (m_nVertices + m_nProcs - 1) / m_nProcs;

    if (m_rank == 0) {
      // Allocate second buffer so we can read and send in parallel
      std::vector<double> tempVertices(chunkSize * 3);
      double* vertices2 = tempVertices.data();
      if (m_nProcs % 2 == 0)
        // swap once so we have the correct buffer at the end
        swap(vertices, vertices2);

      MPI_Request request = MPI_REQUEST_NULL;
      for (int i = 1; i < m_nProcs - 1; i++) {
        logInfo() << "Reading vertices part" << i << "of" << m_nProcs;
        m_serialReader.readVertices(i * chunkSize, chunkSize, vertices);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(vertices, chunkSize * 3, MPI_DOUBLE, i, 0, m_comm, &request);
        swap(vertices, vertices2);
      }

      if (m_nProcs > 1) {
        // Read last one
        const std::size_t lastChunkSize = m_nVertices - (m_nProcs - 1) * chunkSize;
        logInfo() << "Reading vertices part" << (m_nProcs - 1) << "of" << m_nProcs;
        m_serialReader.readVertices((m_nProcs - 1) * chunkSize, lastChunkSize, vertices);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(vertices, lastChunkSize * 3, MPI_DOUBLE, m_nProcs - 1, 0, m_comm, &request);
        swap(vertices, vertices2);
      }

      // Finally read the first part
      logInfo() << "Reading vertices part" << m_nProcs << "of" << m_nProcs;
      m_serialReader.readVertices(0, chunkSize, vertices);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    } else {
      if (m_rank == m_nProcs - 1)
        chunkSize = m_nVertices - (m_nProcs - 1) * chunkSize;

      MPI_Recv(vertices, chunkSize * 3, MPI_DOUBLE, 0, 0, m_comm, MPI_STATUS_IGNORE);
    }
  }

  /**
   * Reads all elements.
   * Each process gets <code>nElements() + processes - 1) / processes</code>
   * elements, except for the last, which gets the remaining elements
   *
   * This is a collective operation.
   *
   * @todo Only tetrahedral meshes are supported
   */
  virtual void readElements(std::size_t* elements) {
    std::size_t chunkSize = (m_nElements + m_nProcs - 1) / m_nProcs;

    if (m_rank == 0) {
      // Allocate second buffer so we can read and send in parallel
      std::vector<std::size_t> tempElements(chunkSize * 4);
      std::size_t* elements2 = tempElements.data();
      if (m_nProcs % 2 == 0)
        // swap once so we have the correct buffer at the end
        swap(elements, elements2);

      MPI_Request request = MPI_REQUEST_NULL;

      for (int i = 1; i < m_nProcs - 1; i++) {
        logInfo() << "Reading elements part" << i << "of" << m_nProcs;
        m_serialReader.readElements(i * chunkSize, chunkSize, elements);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(elements, chunkSize * 4, tndm::mpi_type_t<std::size_t>(), i, 0, m_comm, &request);
        swap(elements, elements2);
      }

      if (m_nProcs > 1) {
        // Read last one
        const std::size_t lastChunkSize = m_nElements - (m_nProcs - 1) * chunkSize;
        logInfo() << "Reading elements part" << (m_nProcs - 1) << "of" << m_nProcs;
        m_serialReader.readElements((m_nProcs - 1) * chunkSize, lastChunkSize, elements);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(elements, lastChunkSize * 4, tndm::mpi_type_t<std::size_t>(), m_nProcs - 1, 0, m_comm, &request);
        swap(elements, elements2);
      }

      // Finally read the first part
      logInfo() << "Reading elements part" << m_nProcs << "of" << m_nProcs;
      m_serialReader.readElements(0, chunkSize, elements);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    } else {
      if (m_rank == m_nProcs - 1)
        chunkSize = m_nElements - (m_nProcs - 1) * chunkSize;

      MPI_Recv(elements, chunkSize * 4, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);
    }
  }

  private:
  /**
   * Initialize some parameters
   */
  void init() {
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_nProcs);
  }

  protected:
  /**
   * Swaps two pointers
   */
  template <typename T> static void swap(T*& p1, T*& p2) {
    T* tmp = p1;
    p1 = p2;
    p2 = tmp;
  }
};

#endif // PARALLEL_MESH_READER_H
