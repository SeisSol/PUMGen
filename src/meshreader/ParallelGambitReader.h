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

#ifndef PARALLEL_GAMBIT_READER
#define PARALLEL_GAMBIT_READER

#include <mpi.h>

#include "GambitReader.h"
#include "ParallelMeshReader.h"
#include "helper/Distributor.h"

#include "third_party/MPITraits.h"

namespace puml {

class ParallelGambitReader : public ParallelMeshReader<GambitReader> {
  public:
  ParallelGambitReader(MPI_Comm comm = MPI_COMM_WORLD) : ParallelMeshReader<GambitReader>(comm) {}

  ParallelGambitReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : ParallelMeshReader<GambitReader>(meshFile, comm) {}

  /**
   * Reads all group numbers.
   * Each process gets <code>nElements() + processes - 1) / processes</code>
   * group numbers, except for the last, which gets the remaining
   *
   * This is a collective operation.
   */
  void readGroups(int* groups) {
    std::size_t chunkSize = getChunksize(nElements(), m_rank, m_nProcs);

    if (m_rank == 0) {
      std::size_t maxChunkSize = chunkSize;
      std::vector<ElementGroup> map(maxChunkSize);

      std::vector<std::vector<int>> aggregator(m_nProcs - 1);
      std::vector<std::size_t> sizes(m_nProcs - 1);
      std::vector<MPI_Request> requests((m_nProcs - 1) * 2, MPI_REQUEST_NULL);

      std::size_t position = 0;

      for (int i = 0; i < m_nProcs; i++) {
        logInfo() << "Reading group information part" << (i + 1) << "of" << m_nProcs;

        chunkSize = getChunksize(nElements(), i, m_nProcs);

        m_serialReader.readGroups(position, chunkSize, map.data());
        position += chunkSize;

        // Wait for all sending from last iteration
        MPI_Waitall((m_nProcs - 1) * 2, requests.data(), MPI_STATUSES_IGNORE);
        for (int j = 0; j < m_nProcs - 1; j++) {
          aggregator[j].clear();
        }

        // Sort group numbers into the corresponding aggregator
        for (std::size_t j = 0; j < chunkSize; j++) {
          if (map[j].element < maxChunkSize)
            // Local element
            groups[map[j].element] = map[j].group;
          else {
            // Element for another processor
            // Serials the struct to make sending easier
            unsigned int proc = map[j].element / maxChunkSize;
            aggregator[proc - 1].push_back(map[j].element % maxChunkSize);
            aggregator[proc - 1].push_back(map[j].group);
          }
        }

        // Send send aggregated mapping
        for (int j = 0; j < m_nProcs - 1; j++) {
          if (aggregator[j].empty()) {
            continue;
          }

          sizes[j] = aggregator[j].size() / 2; // element id + group number
          MPI_Isend(&sizes[j], 1, tndm::mpi_type_t<std::size_t>(), j + 1, 0, m_comm,
                    &requests[j * 2]);
          MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_INT, j + 1, 0, m_comm,
                    &requests[j * 2 + 1]);
        }
      }

      MPI_Waitall((m_nProcs - 1) * 2, requests.data(), MPI_STATUSES_IGNORE);
    } else {
      // Allocate enough memory
      std::vector<int> buf(chunkSize * 2);

      unsigned int recieved = 0;

      while (recieved < chunkSize) {
        std::size_t size;
        MPI_Recv(&size, 1, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);
        MPI_Recv(buf.data(), size * 2, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);

        for (std::size_t i = 0; i < size * 2; i += 2) {
          groups[buf[i]] = buf[i + 1];
        }

        recieved += size;
      }
    }
  }

  /**
   * Reads all boundaries.
   * Each process gets the boundaries for <code>nElements() + processes - 1) /
   * processes</code> elements. Boundaries not specified in the mesh are not
   * modified.
   *
   * This is a collective operation
   *
   * @todo Only tetrahedral meshes are supported
   */
  void readBoundaries(int* boundaries) {
    std::size_t chunkSize = getChunksize(nElements(), m_rank, m_nProcs);

    if (m_rank == 0) {
      std::size_t maxChunkSize = chunkSize;
      std::vector<GambitBoundaryFace> faces(maxChunkSize);

      std::vector<std::vector<int>> aggregator(m_nProcs - 1);
      std::vector<std::size_t> sizes(m_nProcs - 1);
      std::vector<MPI_Request> requests((m_nProcs - 1) * 2, MPI_REQUEST_NULL);

      std::size_t nChunks = (nBoundaries() + chunkSize - 1) / chunkSize;
      std::size_t position = 0;
      for (std::size_t i = 0; i < nChunks; i++) {
        logInfo() << "Reading boundary conditions part" << (i + 1) << "of" << nChunks;

        chunkSize = getChunksize(nBoundaries(), i, m_nProcs);

        m_serialReader.readBoundaries(position, chunkSize, faces.data());
        position += chunkSize;

        // Wait for all sending from last iteration
        MPI_Waitall((m_nProcs - 1) * 2, requests.data(), MPI_STATUSES_IGNORE);
        for (int j = 0; j < m_nProcs - 1; j++) {
          aggregator[j].clear();
        }

        // Sort boundary conditions into the corresponding aggregator
        for (std::size_t j = 0; j < chunkSize; j++) {
          if (faces[j].element < maxChunkSize)
            // Local element
            boundaries[faces[j].element * 4 + faces[j].face] = faces[j].type;
          else {
            // Face for another processor
            // Serials the struct to make sending easier
            unsigned int proc = faces[j].element / maxChunkSize;
            aggregator[proc - 1].push_back((faces[j].element % maxChunkSize) * 4 + faces[j].face);
            aggregator[proc - 1].push_back(faces[j].type);
          }
        }

        // Send send aggregated values
        for (int j = 0; j < m_nProcs - 1; j++) {
          if (aggregator[j].empty()) {
            continue;
          }

          sizes[j] = aggregator[j].size() / 2; // element id + face type
          MPI_Isend(&sizes[j], 1, tndm::mpi_type_t<std::size_t>(), j + 1, 0, m_comm,
                    &requests[j * 2]);
          MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_INT, j + 1, 0, m_comm,
                    &requests[j * 2 + 1]);
        }
      }

      MPI_Waitall((m_nProcs - 1) * 2, requests.data(), MPI_STATUSES_IGNORE);

      // Send finish signal to all other processes
      std::fill(sizes.begin(), sizes.end(), 0);
      for (int i = 0; i < m_nProcs - 1; i++) {
        MPI_Isend(&sizes[i], 1, tndm::mpi_type_t<std::size_t>(), i + 1, 0, m_comm, &requests[i]);
      }
      MPI_Waitall(m_nProcs - 1, requests.data(), MPI_STATUSES_IGNORE);
    } else {
      // Allocate enough memory
      std::vector<int> buf(chunkSize * 4);

      while (true) {
        std::size_t size;
        MPI_Recv(&size, 1, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);
        if (size == 0) {
          // Finished
          break;
        }

        MPI_Recv(buf.data(), size * 2, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);

        for (std::size_t i = 0; i < size * 2; i += 2) {
          boundaries[buf[i]] = buf[i + 1];
        }
      }
    }
  }
};

} // namespace puml

#endif // PARALLEL_GAMBIT_READER
