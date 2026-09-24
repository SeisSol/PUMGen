// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_MESHREADER_PARALLELGAMBITREADER_H_
#define PUMGEN_SRC_MESHREADER_PARALLELGAMBITREADER_H_

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include "GambitReader.h"
#include "ParallelMeshReader.h"
#include "helper/Distributor.h"

#include "third_party/MPITraits.h"
#include "utils/logger.h"

namespace puml {

class ParallelGambitReader : public ParallelMeshReader<GambitReader> {
  public:
  ParallelGambitReader(MPI_Comm comm = MPI_COMM_WORLD) : ParallelMeshReader<GambitReader>(comm) {}

  ParallelGambitReader(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : ParallelMeshReader<GambitReader>(meshFile, comm) {}

  /**
   * Reads all group numbers. Each process gets the group numbers of its elements (see
   * ParallelMeshReader::readElements).
   *
   * This is a collective operation.
   */
  void readGroups(int* groups) {
    std::size_t numEntries = 0;
    if (m_rank == 0) {
      numEntries = m_serialReader.nGroupedElements();
      if (numEntries < nElements()) {
        logWarning() << nElements() - numEntries << "elements belong to no group and get group 0";
      } else if (numEntries > nElements()) {
        logWarning() << "The groups list" << numEntries - nElements()
                     << "elements more than the mesh has; an element listed in several groups "
                        "gets the last of them";
      }
    }
    scatterToElementOwners<ElementGroup>(
        numEntries, 1, groups, "group information",
        [&](std::size_t start, std::size_t count, ElementGroup* entries) {
          m_serialReader.readGroups(start, count, entries);
        },
        [](const ElementGroup& entry) {
          return std::make_tuple(entry.element, std::size_t{0}, entry.group);
        });
  }

  /**
   * Reads all boundaries. Each process gets the boundaries of the faces of its elements (see
   * ParallelMeshReader::readElements). Boundaries not specified in the mesh are not modified.
   *
   * This is a collective operation.
   *
   * @todo Only tetrahedral meshes are supported
   */
  void readBoundaries(int* boundaries) {
    scatterToElementOwners<GambitBoundaryFace>(
        nBoundaries(), 4, boundaries, "boundary conditions",
        [&](std::size_t start, std::size_t count, GambitBoundaryFace* entries) {
          m_serialReader.readBoundaries(start, count, entries);
        },
        [](const GambitBoundaryFace& entry) {
          return std::make_tuple(entry.element, static_cast<std::size_t>(entry.face), entry.type);
        });
  }

  private:
  /**
   * Reads per-element entries in chunks on rank 0 and stores the value of each entry at
   * <code>local[localElement * slotsPerElement + slot]</code> on the process owning the element.
   *
   * @param read Callable (start, count, entries) reading <code>count</code> entries starting at
   * <code>start</code>
   * @param decode Callable returning the tuple (element, slot, value) for an entry
   */
  template <typename Entry, typename ReadF, typename DecodeF>
  void scatterToElementOwners(std::size_t numEntries, std::size_t slotsPerElement, int* local,
                              const char* name, ReadF&& read, DecodeF&& decode) {
    const auto messageType = tndm::mpi_type_t<std::int64_t>();

    if (m_rank != 0) {
      std::vector<std::int64_t> buffer;
      while (true) {
        std::size_t size = 0;
        MPI_Recv(&size, 1, tndm::mpi_type_t<std::size_t>(), 0, 0, m_comm, MPI_STATUS_IGNORE);
        if (size == 0) {
          break;
        }
        buffer.resize(2 * size);
        MPI_Recv(buffer.data(), messageSize(buffer.size()), messageType, 0, 0, m_comm,
                 MPI_STATUS_IGNORE);
        for (std::size_t i = 0; i < size; ++i) {
          local[buffer[2 * i]] = static_cast<int>(buffer[2 * i + 1]);
        }
      }
      return;
    }

    const std::size_t chunkSize = std::max(getChunksize(nElements(), 0, m_nProcs), std::size_t{1});
    std::vector<Entry> entries(std::min(chunkSize, numEntries));
    std::vector<std::vector<std::int64_t>> outgoing(m_nProcs);
    std::vector<std::size_t> sizes(m_nProcs, 0);
    std::vector<MPI_Request> requests(2 * m_nProcs, MPI_REQUEST_NULL);

    for (std::size_t start = 0; start < numEntries; start += chunkSize) {
      const std::size_t count = std::min(chunkSize, numEntries - start);
      logInfo() << "Reading" << name << start << "to" << start + count << "of" << numEntries;
      read(start, count, entries.data());

      // the send buffers of the previous chunk are reused
      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      for (auto& buffer : outgoing) {
        buffer.clear();
      }

      for (std::size_t i = 0; i < count; ++i) {
        const auto [element, slot, value] = decode(entries[i]);
        if (element >= nElements() || slot >= slotsPerElement) {
          logError() << "Invalid" << name << "entry: element" << element << "slot" << slot;
        }
        const int owner = getChunkOwner(nElements(), element, m_nProcs);
        const std::size_t index =
            (element - getChunksum(nElements(), owner, m_nProcs)) * slotsPerElement + slot;
        if (owner == 0) {
          local[index] = value;
        } else {
          outgoing[owner].push_back(static_cast<std::int64_t>(index));
          outgoing[owner].push_back(value);
        }
      }

      for (int rank = 1; rank < m_nProcs; ++rank) {
        if (outgoing[rank].empty()) {
          continue;
        }
        sizes[rank] = outgoing[rank].size() / 2;
        MPI_Isend(&sizes[rank], 1, tndm::mpi_type_t<std::size_t>(), rank, 0, m_comm,
                  &requests[2 * rank]);
        MPI_Isend(outgoing[rank].data(), messageSize(outgoing[rank].size()), messageType, rank, 0,
                  m_comm, &requests[2 * rank + 1]);
      }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // an empty message ends the transfer
    std::fill(sizes.begin(), sizes.end(), 0);
    for (int rank = 1; rank < m_nProcs; ++rank) {
      MPI_Isend(&sizes[rank], 1, tndm::mpi_type_t<std::size_t>(), rank, 0, m_comm, &requests[rank]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  }
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_PARALLELGAMBITREADER_H_
