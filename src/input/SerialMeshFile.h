// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_SERIALMESHFILE_H_
#define PUMGEN_SRC_INPUT_SERIALMESHFILE_H_

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include "helper/Distributor.h"
#include <cstddef>
#include <vector>

#include "MeshData.h"
#include "utils/logger.h"

/**
 * Read a mesh from a serial file
 */
template <typename T> class SerialMeshFile : public FullStorageMeshData {
  public:
  virtual ~SerialMeshFile() = default;

  std::size_t vertexSize() const override { return T::Dim; }
  std::size_t cellSize() const override { return puml::nodeCount(T::Dim, T::Order); }

  private:
#ifdef PARALLEL
  MPI_Comm m_comm;
#endif // PARALLEL

  int m_rank;
  int m_nProcs;

  T m_meshReader;

  public:
#ifdef PARALLEL
  SerialMeshFile(const char* meshFile, int boundarySize, MPI_Comm comm = MPI_COMM_WORLD)
      : FullStorageMeshData(boundarySize), m_comm(comm), m_meshReader(comm) {
    init();
    open(meshFile);
  }
#else  // PARALLEL
  SerialMeshFile(const char* meshFile, int boundarySize) : FullStorageMeshData(boundarySize) {
    init();
    open(meshFile);
  }
#endif // PARALLEL

  private:
  /**
   * Sets some parameters (called from the constructor)
   */
  void init() {
#ifdef PARALLEL
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_nProcs);
#else  // PARALLLEL
    m_rank = 0;
    m_nProcs = 1;
#endif // PARALLEL
  }

  void open(const char* meshFile) {
    m_meshReader.open(meshFile);

    const std::size_t nVertices = m_meshReader.nVertices();
    const std::size_t nElements = m_meshReader.nElements();
    const std::size_t nLocalVertices = getChunksize(nVertices, m_rank, m_nProcs);
    const std::size_t nLocalElements = getChunksize(nElements, m_rank, m_nProcs);

    bool identify = false;
    if constexpr (T::SupportsIdentify) {
      identify = m_meshReader.hasIdentify();
    }

    setup(nLocalElements, nLocalVertices, identify);

    logInfo() << "Read vertex coordinates";
    m_meshReader.readVertices(geometryData.data());

    logInfo() << "Read cell vertices";
    m_meshReader.readElements(connectivityData.data());

    logInfo() << "Read cell groups";
    m_meshReader.readGroups(groupData.data());

    logInfo() << "Read boundary conditions";
    std::vector<int> preBoundaryData(nLocalElements * (vertexSize() + 1));
    m_meshReader.readBoundaries(preBoundaryData.data());
    for (std::size_t i = 0; i < nLocalElements; ++i) {
      for (int j = 0; j < (vertexSize() + 1); ++j) {
        setBoundary(i, j, preBoundaryData[(vertexSize() + 1) * i + j]);
      }
    }

    if constexpr (T::SupportsIdentify) {
      if (identify) {
        m_meshReader.readIdentify(identifyData.data());
      }
    }
  }
};

#endif // PUMGEN_SRC_INPUT_SERIALMESHFILE_H_
