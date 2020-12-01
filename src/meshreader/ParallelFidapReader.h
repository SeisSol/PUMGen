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

class ParallelFidapReader : public ParallelMeshReader<FidapReader> {
  private:
  /** Face map required for boundaries */
  std::map<FaceVertex, FaceElement> m_faceMap;

  public:
  ParallelFidapReader(MPI_Comm comm = MPI_COMM_WORLD) : ParallelMeshReader<FidapReader>(comm) {}

  ParallelFidapReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : ParallelMeshReader<FidapReader>(meshFile, comm) {}

  void readElements(int* elements) {
    ParallelMeshReader<FidapReader>::readElements(elements);

    // Create the face map
    unsigned int chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;
    unsigned int maxChunkSize = chunkSize;
    if (m_rank == m_nProcs - 1)
      chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

    FidapReader::createFaceMap(elements, chunkSize, m_faceMap);
    for (std::map<FaceVertex, FaceElement>::iterator i = m_faceMap.begin(); i != m_faceMap.end();
         i++)
      i->second.element += m_rank * maxChunkSize;
  }

  /**
   * @copydoc ParallelGambitReader::readGroups
   */
  void readGroups(int* groups) {
    unsigned int chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;

    if (m_rank == 0) {
      // Allocate second buffer so we can read and send in parallel
      int* groups2 = new int[chunkSize];
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
        unsigned int lastChunkSize = nElements() - (m_nProcs - 1) * chunkSize;
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

      delete[] groups2;
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

    unsigned int chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;
    unsigned int maxChunkSize = chunkSize;

    if (m_rank == 0) {
      FidapBoundaryFace* faces = new FidapBoundaryFace[maxChunkSize];

      std::vector<int>* aggregator = new std::vector<int>[m_nProcs - 1];
      unsigned int* sizes = new unsigned int[m_nProcs - 1];
      MPI_Request* requests = new MPI_Request[(m_nProcs - 1) * 2];
      for (int i = 0; i < m_nProcs - 1; i++) {
        requests[i * 2] = MPI_REQUEST_NULL;
        requests[i * 2 + 1] = MPI_REQUEST_NULL;
      }

      unsigned int nChunks = (nBoundaries() + chunkSize - 1) / chunkSize;
      for (unsigned int i = 0; i < nChunks; i++) {
        logInfo() << "Reading boundary conditions part" << (i + 1) << "of" << nChunks;

        if (i == nChunks - 1)
          chunkSize = nBoundaries() - (nChunks - 1) * chunkSize;

        m_serialReader.readBoundaries(i * maxChunkSize, chunkSize, faces);

        // Wait for all sending from last iteration
        MPI_Waitall((m_nProcs - 1) * 2, requests, MPI_STATUSES_IGNORE);
        for (int j = 0; j < m_nProcs - 1; j++)
          aggregator[j].clear();

        // Sort boundary conditions into the corresponding aggregator
        for (unsigned int j = 0; j < chunkSize; j++) {
          FaceVertex v(faces[j].vertices);

          if (v.vertices[1] < vertexChunk) {
            facePos.push_back(m_faceMap.at(v));
            faceType.push_back(faces[j].type);
          } else {
            // Face for another processor
            // Serials the struct to make sending easier
            unsigned int proc = v.vertices[1] / vertexChunk;
            aggregator[proc - 1].insert(aggregator[proc - 1].end(), v.vertices, v.vertices + 3);
            aggregator[proc - 1].push_back(faces[j].type);
          }
        }

        // Send send aggregated values
        for (int j = 0; j < m_nProcs - 1; j++) {
          if (aggregator[j].empty())
            continue;

          sizes[j] = aggregator[j].size() / 4; // 3 vertices + face type
          MPI_Isend(&sizes[j], 1, MPI_UNSIGNED, j + 1, 0, m_comm, &requests[j * 2]);
          MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_INT, j + 1, 0, m_comm,
                    &requests[j * 2 + 1]);
        }
      }

      MPI_Waitall((m_nProcs - 1) * 2, requests, MPI_STATUSES_IGNORE);

      delete[] faces;
      delete[] aggregator;

      // Send finish signal to all other processes
      memset(sizes, 0, (m_nProcs - 1) * sizeof(unsigned int));
      for (int i = 0; i < m_nProcs - 1; i++)
        MPI_Isend(&sizes[i], 1, MPI_UNSIGNED, i + 1, 0, m_comm, &requests[i]);
      MPI_Waitall(m_nProcs - 1, requests, MPI_STATUSES_IGNORE);

      delete[] sizes;
      delete[] requests;
    } else {
      if (m_rank == m_nProcs - 1)
        chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

      // Allocate enough memory
      int* buf = new int[chunkSize * 4];

      while (true) {
        unsigned int size;
        MPI_Recv(&size, 1, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
        if (size == 0)
          // Finished
          break;

        MPI_Recv(buf, size * 4, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);

        for (unsigned int i = 0; i < size * 4; i += 4) {
          FaceVertex v(&buf[i]);
          facePos.push_back(m_faceMap.at(v));
          faceType.push_back(buf[i + 3]);
        }
      }

      delete[] buf;
    }

    // Distribute the faces to the final ranks
    PCU_Comm_Begin();
    for (unsigned int i = 0; i < facePos.size(); i++) {
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
