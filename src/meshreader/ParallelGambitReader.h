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
        unsigned int chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;

        if (m_rank == 0) {
            unsigned int maxChunkSize = chunkSize;
            ElementGroup* map = new ElementGroup[maxChunkSize];

            std::vector<int>* aggregator = new std::vector<int>[m_nProcs - 1];
            unsigned int* sizes = new unsigned int[m_nProcs - 1];
            MPI_Request* requests = new MPI_Request[(m_nProcs - 1) * 2];
            for (int i = 0; i < m_nProcs - 1; i++) {
                requests[i * 2] = MPI_REQUEST_NULL;
                requests[i * 2 + 1] = MPI_REQUEST_NULL;
            }

            for (int i = 0; i < m_nProcs; i++) {
                logInfo() << "Reading group information part" << (i + 1) << "of" << m_nProcs;

                if (i == m_nProcs - 1)
                    chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

                m_serialReader.readGroups(i * maxChunkSize, chunkSize, map);

                // Wait for all sending from last iteration
                MPI_Waitall((m_nProcs - 1) * 2, requests, MPI_STATUSES_IGNORE);
                for (int j = 0; j < m_nProcs - 1; j++)
                    aggregator[j].clear();

                // Sort group numbers into the corresponding aggregator
                for (unsigned int j = 0; j < chunkSize; j++) {
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
                    if (aggregator[j].empty())
                        continue;

                    sizes[j] = aggregator[j].size() / 2; // element id + group number
                    MPI_Isend(&sizes[j], 1, MPI_UNSIGNED, j + 1, 0, m_comm, &requests[j * 2]);
                    MPI_Isend(&aggregator[j][0], aggregator[j].size(), MPI_INT, j + 1, 0, m_comm,
                              &requests[j * 2 + 1]);
                }
            }

            MPI_Waitall((m_nProcs - 1) * 2, requests, MPI_STATUSES_IGNORE);

            delete[] map;
            delete[] aggregator;
            delete[] sizes;
            delete[] requests;
        } else {
            if (m_rank == m_nProcs - 1)
                chunkSize = nElements() - (m_nProcs - 1) * chunkSize;

            // Allocate enough memory
            unsigned int* buf = new unsigned int[chunkSize * 2];

            unsigned int recieved = 0;

            while (recieved < chunkSize) {
                unsigned int size;
                MPI_Recv(&size, 1, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
                MPI_Recv(buf, size * 2, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);

                for (unsigned int i = 0; i < size * 2; i += 2)
                    groups[buf[i]] = buf[i + 1];

                recieved += size;
            }

            delete[] buf;
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
        unsigned int chunkSize = (nElements() + m_nProcs - 1) / m_nProcs;

        if (m_rank == 0) {
            unsigned int maxChunkSize = chunkSize;
            GambitBoundaryFace* faces = new GambitBoundaryFace[maxChunkSize];

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
                    if (faces[j].element < maxChunkSize)
                        // Local element
                        boundaries[faces[j].element * 4 + faces[j].face] = faces[j].type;
                    else {
                        // Face for another processor
                        // Serials the struct to make sending easier
                        unsigned int proc = faces[j].element / maxChunkSize;
                        aggregator[proc - 1].push_back((faces[j].element % maxChunkSize) * 4 +
                                                       faces[j].face);
                        aggregator[proc - 1].push_back(faces[j].type);
                    }
                }

                // Send send aggregated values
                for (int j = 0; j < m_nProcs - 1; j++) {
                    if (aggregator[j].empty())
                        continue;

                    sizes[j] = aggregator[j].size() / 2; // element id + face type
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
            int* buf = new int[chunkSize * 2];

            while (true) {
                unsigned int size;
                MPI_Recv(&size, 1, MPI_UNSIGNED, 0, 0, m_comm, MPI_STATUS_IGNORE);
                if (size == 0)
                    // Finished
                    break;

                MPI_Recv(buf, size * 2, MPI_INT, 0, 0, m_comm, MPI_STATUS_IGNORE);

                for (unsigned int i = 0; i < size * 2; i += 2)
                    boundaries[buf[i]] = buf[i + 1];
            }

            delete[] buf;
        }
    }
};

} // namespace puml

#endif // PARALLEL_GAMBIT_READER
