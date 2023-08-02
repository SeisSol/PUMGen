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

#include "ApfConvertWrapper.h"
#include <PCU.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <vector>

#include "MeshInput.h"

/**
 * Read a mesh from a serial file
 */
template <typename T> class SerialMeshFile : public MeshInput {
  private:
#ifdef PARALLEL
  MPI_Comm m_comm;
#endif // PARALLEL

  int m_rank;
  int m_nProcs;

  T m_meshReader;

  public:
#ifdef PARALLEL
  SerialMeshFile(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
      : m_comm(comm), m_meshReader(comm) {
    init();
    open(meshFile);
  }
#else  // PARALLEL
  SerialMeshFile(const char* meshFile) {
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

    unsigned int nVertices = m_meshReader.nVertices();
    unsigned int nElements = m_meshReader.nElements();
    unsigned int nLocalVertices = (nVertices + m_nProcs - 1) / m_nProcs;
    unsigned int nLocalElements = (nElements + m_nProcs - 1) / m_nProcs;
    if (m_rank == m_nProcs - 1) {
      nLocalVertices = nVertices - (m_nProcs - 1) * nLocalVertices;
      nLocalElements = nElements - (m_nProcs - 1) * nLocalElements;
    }

    gmi_register_null();
    gmi_model* model = gmi_load(".null");
    m_mesh = apf::makeEmptyMdsMesh(model, 3, false);

    // Create elements
    apf::GlobalToVert vertMap;
    {
      std::vector<ElementID> elements;
      {
        std::vector<int> elementsInt(nLocalElements * 4);
        m_meshReader.readElements(elementsInt.data());
        elements = std::vector<ElementID>(elementsInt.begin(), elementsInt.end());
      }
      logInfo(m_rank) << "Create APF connectivity";
      apf::construct(m_mesh, elements.data(), nLocalElements, apf::Mesh::TET, vertMap);
    }

    apf::alignMdsRemotes(m_mesh);
    apf::deriveMdsModel(m_mesh);

    // Set vertices
    {
      std::vector<double> vertices(nLocalVertices * 3);
      m_meshReader.readVertices(vertices.data());
      logInfo(m_rank) << "Set coordinates in APF";
      apf::setCoords(m_mesh, vertices.data(), nLocalVertices, vertMap);
    }

    // Set boundaries
    apf::MeshTag* boundaryTag = m_mesh->createIntTag("boundary condition", 1);
    {
      std::vector<int> boundaries(nLocalElements * 4);
      m_meshReader.readBoundaries(boundaries.data());
      apf::MeshIterator* it = m_mesh->begin(3);
      unsigned int i = 0;
      while (apf::MeshEntity* element = m_mesh->iterate(it)) {
        apf::Adjacent adjacent;
        m_mesh->getAdjacent(element, 2, adjacent);

        for (unsigned int j = 0; j < 4; j++) {
          if (!boundaries[i * 4 + j]) {
            continue;
          }

          m_mesh->setIntTag(adjacent[j], boundaryTag, &boundaries[i * 4 + j]);
        }

        i++;
      }
      m_mesh->end(it);
    }

    // Set groups
    apf::MeshTag* groupTag = m_mesh->createIntTag("group", 1);
    {
      std::vector<int> groups(nLocalElements);
      m_meshReader.readGroups(groups.data());
      apf::MeshIterator* it = m_mesh->begin(3);
      unsigned int i = 0;
      while (apf::MeshEntity* element = m_mesh->iterate(it)) {
        m_mesh->setIntTag(element, groupTag, &groups[i]);
        i++;
      }
      m_mesh->end(it);
    }
  }
};

#endif // SERIAL_MESH_FILE_H
