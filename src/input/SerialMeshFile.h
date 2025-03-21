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

#ifndef SERIAL_MESH_FILE_H
#define SERIAL_MESH_FILE_H

#ifdef PARALLEL
#include <mpi.h>
#endif // PARALLEL

#include "aux/Distributor.h"
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

    setup(nLocalElements, nLocalVertices);

    logInfo(m_rank) << "Read vertex coordinates";
    m_meshReader.readVertices(geometryData.data());

    logInfo(m_rank) << "Read cell vertices";
    m_meshReader.readElements(connectivityData.data());

    logInfo(m_rank) << "Read cell groups";
    m_meshReader.readGroups(groupData.data());

    logInfo(m_rank) << "Read boundary conditions";
    std::vector<int> preBoundaryData(nLocalElements * (vertexSize() + 1));
    m_meshReader.readBoundaries(preBoundaryData.data());
    for (std::size_t i = 0; i < nLocalElements; ++i) {
      for (int j = 0; j < (vertexSize() + 1); ++j) {
        setBoundary(i, j, preBoundaryData[(vertexSize() + 1) * i + j]);
      }
    }
  }
};

#endif // SERIAL_MESH_FILE_H
