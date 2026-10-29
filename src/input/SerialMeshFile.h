// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_SERIALMESHFILE_H_
#define PUMGEN_SRC_INPUT_SERIALMESHFILE_H_

#include <mpi.h>

#include <cstddef>
#include <vector>

#include "MeshData.h"
#include "helper/Distributor.h"
#include "utils/logger.h"

/**
 * Read a mesh from a serial file. A reader either hands over the part of the mesh of this rank
 * (ProvidesLocalMesh), or it fills the arrays of a tetrahedral mesh.
 */
template <typename T> class SerialMeshFile : public FullStorageMeshData {
  public:
  explicit SerialMeshFile(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : m_meshReader(comm) {
    int rank = 0;
    int processes = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &processes);

    m_meshReader.open(meshFile);
    if constexpr (T::ProvidesLocalMesh) {
      take(m_meshReader.read());
    } else {
      readTetrahedra(rank, processes);
    }
  }

  private:
  T m_meshReader;

  void readTetrahedra(int rank, int processes) {
    const std::size_t nLocalVertices = getChunksize(m_meshReader.nVertices(), rank, processes);
    const std::size_t nLocalElements = getChunksize(m_meshReader.nElements(), rank, processes);

    setup(nLocalElements, nLocalVertices);

    logInfo() << "Read vertex coordinates";
    m_meshReader.readVertices(geometryData.data());

    logInfo() << "Read cell vertices";
    m_meshReader.readElements(connectivityData.data());

    logInfo() << "Read cell groups";
    m_meshReader.readGroups(groupData.data());

    logInfo() << "Read boundary conditions";
    constexpr std::size_t Faces = 4;
    std::vector<int> faces(nLocalElements * Faces);
    m_meshReader.readBoundaries(faces.data());
    for (std::size_t i = 0; i < nLocalElements; ++i) {
      for (std::size_t j = 0; j < Faces; ++j) {
        setBoundary(i, static_cast<int>(j), faces[Faces * i + j]);
      }
    }
  }
};

#endif // PUMGEN_SRC_INPUT_SERIALMESHFILE_H_
