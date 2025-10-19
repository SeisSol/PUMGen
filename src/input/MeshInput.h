// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_MESHINPUT_H_
#define PUMGEN_SRC_INPUT_MESHINPUT_H_
#include "utils/logger.h"
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <functional>

#ifndef MESH_INTPUT_H
#define MESH_INTPUT_H

/**
 * Interface for mesh input
 */
class ApfMeshInput : public FullStorageMeshData {
  private:
  template <typename F> static void iterate(apf::Mesh* mesh, int dim, F&& function) {
    apf::MeshIterator* it = mesh->begin(dim);
    std::size_t index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      std::invoke(function, index, element);
      ++index;
    }
    mesh->end(it);
  }

  protected:
  apf::Mesh2* m_mesh = nullptr;

  ApfMeshInput(int boundarySize) : FullStorageMeshData(boundarySize) {}

  public:
  apf::Mesh2* getMesh() { return m_mesh; }

  virtual ~ApfMeshInput() {
    delete m_mesh;
    m_mesh = nullptr;
  }

  void generate() {
    if (alignMdsMatches(m_mesh)) {
      logWarning() << "Fixed misaligned matches";
    }
    m_mesh->verify();
    apf::GlobalNumbering* vertexNum = apf::makeGlobal(apf::numberOwnedNodes(m_mesh, "vertices"));
    apf::synchronize(vertexNum);

    apf::Sharing* sharing = apf::getSharing(m_mesh);

    setup(apf::countOwned(m_mesh, 3), apf::countOwned(m_mesh, 0));

    // connectivity
    iterate(m_mesh, 3, [&](auto index, auto* element) {
      apf::NewArray<long> vn;
      apf::getElementNumbers(vertexNum, element, vn);

      for (int i = 0; i < 4; i++) {
        connectivityData[index * 4 + i] = vn[i];
      }
    });

    std::size_t vertexIndex = 0;

    // geometry
    iterate(m_mesh, 0, [&](auto index, auto* element) {
      if (!sharing->isOwned(element)) {
        return;
      }

      apf::Vector3 point;
      m_mesh->getPoint(element, 0, point);
      double geometry[3];
      point.toArray(geometry);

      for (int i = 0; i < 3; ++i) {
        geometryData[vertexIndex * 3 + i] = geometry[i];
      }
      ++vertexIndex;
    });

    // groups
    apf::MeshTag* groupTag = m_mesh->findTag("group");
    iterate(m_mesh, 3, [&](auto index, auto* element) {
      int group;
      m_mesh->getIntTag(element, groupTag, &group);
      groupData[index] = group;
    });

    // boundary
    apf::MeshTag* boundaryTag = m_mesh->findTag("boundary condition");
    iterate(m_mesh, 3, [&](auto index, auto* element) {
      apf::Downward faces;
      m_mesh->getDownward(element, 2, faces);

      for (int i = 0; i < 4; i++) {
        if (m_mesh->hasTag(faces[i], boundaryTag)) {
          int b;
          m_mesh->getIntTag(faces[i], boundaryTag, &b);

          setBoundary(index, i, b);
        }
      }
    });

    delete sharing;
  }
};

#endif // PUMGEN_SRC_INPUT_MESHINPUT_H_
