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

#ifndef PARALLEL_FIDAP_READER
#define PARALLEL_FIDAP_READER

#include <mpi.h>

#include <PCU.h>

#include "FidapReader.h"
#include "ParallelMeshReader.h"
#include "third_party/MPITraits.h"

class ParallelFidapReader : public ParallelMeshReader<FidapReader> {
  private:
  /** Face map required for boundaries */
  std::map<FaceVertex, FaceElement> m_faceMap;

  public:
  ParallelFidapReader(MPI_Comm comm = MPI_COMM_WORLD) : ParallelMeshReader<FidapReader>(comm) {}

  ParallelFidapReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : ParallelMeshReader<FidapReader>(meshFile, comm) {}

  void readElements(std::size_t* elements) {
    ParallelMeshReader<FidapReader>::readElements(elements);

    // Create the face map
    std::size_t chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;
    const std::size_t maxChunkSize = chunkSize;
    if (m_rank == m_nProcs - 1) {
      chunkSize = nElements() - (m_nProcs - 1) * chunkSize;
    }

    FidapReader::createFaceMap(elements, chunkSize, m_faceMap);
    for (auto& face : m_faceMap) {
      face.second.element += m_rank * maxChunkSize;
    }
  }

  /**
   * @copydoc ParallelGambitReader::readGroups
   */
  void readGroups(int* groups) {
    std::size_t chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;

    if (m_rank == 0) {
      // Allocate second buffer so we can read and send in parallel
      std::vector<int> tempGroups(chunkSize);
      int* groups2 = tempGroups.data();
      if (m_nProcs % 2 == 0)
        // swap once so we have the correct buffer at the end
        swap(groups, groups2);

      MPI_Request request = MPI_REQUEST_NULL;

      for (int i = 1; i < m_nProcs - 1; i++) {
        logInfo() << "Reading group information part" << i << "of" << m_nProcs;
        m_serialReader.readGroups(i * chunkSize, chunkSize, groups);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(groups, chunkSize, MPI_INT, i, 0, m_comm, &request);
        swap(groups, groups2);
      }

      if (m_nProcs > 1) {
        // Read last one
        const std::size_t lastChunkSize = nElements() - (m_nProcs - 1) * chunkSize;
        logInfo() << "Reading group information part" << (m_nProcs - 1) << "of" << m_nProcs;
        m_serialReader.readGroups((m_nProcs - 1) * chunkSize, lastChunkSize, groups);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(groups, lastChunkSize, MPI_INT, m_nProcs - 1, 0, m_comm, &request);
        swap(groups, groups2);
      }

      // Finally read the first part
      logInfo() << "Reading group information part" << m_nProcs << "of" << m_nProcs;
      m_serialReader.readGroups(0, chunkSize, groups);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    } else {
      if (m_rank == m_nProcs - 1)
        chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

      MPI_Recv(groups, chunkSize, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);
    }
  }

  /**
   * @copydoc ParallelGambitReader::readBoundaries
   */
  void readBoundaries(int* boundaries) {
    // Create the distributed face map
    int vertexChunk = (nVertices() + m_nProcs - 1) / m_nProcs;
    PCU_Comm_Begin();
    for (std::map<FaceVertex, FaceElement>::iterator i = m_faceMap.begin(); i != m_faceMap.end();) {
      int proc = i->first.vertices[1] / vertexChunk;
      if (proc == m_rank)
        i++;
      else {
        PCU_COMM_PACK(proc, i->first);
        PCU_COMM_PACK(proc, i->second);
        m_faceMap.erase(i++);
      }
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      FaceVertex v;
      PCU_COMM_UNPACK(v);
      FaceElement e;
      PCU_COMM_UNPACK(e);
      m_faceMap[v] = e;
    }

    // Distribute the faces
    std::vector<FaceElement> facePos;
    std::vector<int> faceType;

    std::size_t chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;
    const std::size_t maxChunkSize = chunkSize;

    if (m_rank == 0) {
      std::vector<FidapBoundaryFace> faces(maxChunkSize);

      std::vector<std::vector<int>> aggregator(m_nProcs - 1);
      std::vector<std::size_t> sizes(m_nProcs - 1);
      std::vector<MPI_Request> requests((m_nProcs - 1) * 2, MPI_REQUEST_NULL);

      const std::size_t nChunks = (nBoundaries() + chunkSize - 1) / chunkSize;
      for (std::size_t i = 0; i < nChunks; i++) {
        logInfo() << "Reading boundary conditions part" << (i + 1) << "of" << nChunks;

        if (i == nChunks - 1) {
          chunkSize = nBoundaries() - (nChunks - 1) * chunkSize;
        }

        m_serialReader.readBoundaries(i * maxChunkSize, chunkSize, faces.data());

        // Wait for all sending from last iteration
        MPI_Waitall((m_nProcs - 1) * 2, requests.data(), MPI_STATUSES_IGNORE);
        for (int j = 0; j < m_nProcs - 1; j++) {
          aggregator[j].clear();
        }

        // Sort boundary conditions into the corresponding aggregator
        for (std::size_t j = 0; j < chunkSize; j++) {
          FaceVertex v(faces[j].vertices);

          if (v.vertices[1] < vertexChunk) {
            facePos.push_back(m_faceMap.at(v));
            faceType.push_back(faces[j].type);
          } else {
            // Face for another processor
            // Serials the struct to make sending easier
            int proc = v.vertices[1] / vertexChunk;
            aggregator[proc - 1].insert(aggregator[proc - 1].end(), v.vertices.begin(), v.vertices.end());
            aggregator[proc - 1].push_back(faces[j].type);
          }
        }

        // Send send aggregated values
        for (int j = 0; j < m_nProcs - 1; j++) {
          if (aggregator[j].empty())
            continue;

          sizes[j] = aggregator[j].size() / 4; // 3 vertices + face type
          MPI_Isend(&sizes[j], 1, tndm::mpi_type_t<std::size_t>(), j + 1, 0, m_comm, &requests[j * 2]);
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
      if (m_rank == m_nProcs - 1)
        chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

      // Allocate enough memory
      std::vector<std::size_t> buf(chunkSize * 4);

      while (true) {
        std::size_t size;
        MPI_Recv(&size, 1, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);
        if (size == 0)
          // Finished
          break;

        MPI_Recv(buf.data(), size * 4, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);

        for (std::size_t i = 0; i < size * 4; i += 4) {
          FaceVertex v(buf[i], buf[i+1], buf[i+2]);
          facePos.push_back(m_faceMap.at(v));
          faceType.push_back(buf[i + 3]);
        }
      }
    }

    // Distribute the faces to the final ranks
    PCU_Comm_Begin();
    for (std::size_t i = 0; i < facePos.size(); i++) {
      if (facePos[i].element / maxChunkSize == static_cast<unsigned int>(m_rank))
        boundaries[(facePos[i].element % maxChunkSize) * 4 + facePos[i].side] = faceType[i];
      else {
        PCU_COMM_PACK(facePos[i].element / maxChunkSize, facePos[i]);
        PCU_COMM_PACK(facePos[i].element / maxChunkSize, faceType[i]);
      }
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      FaceElement e;
      PCU_COMM_UNPACK(e);
      int type;
      PCU_COMM_UNPACK(type);
      boundaries[(e.element % maxChunkSize) * 4 + e.side] = type;
    }
  }
};

#endif // PARALLEL_FIDAP_READER
