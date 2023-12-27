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

#include "MeshData.h"
#include <apfMesh2.h>

#ifndef MESH_INTPUT_H
#define MESH_INTPUT_H

/**
 * Interface for mesh input
 */
class ApfMeshInput {
  protected:
  apf::Mesh2* m_mesh;

  public:
  apf::Mesh2* getMesh() { return m_mesh; }

  MeshData* generate() {
    FullStorageMeshData* data = new FullStorageMeshData();

    return data;
  }
};

#endif // MESH_INPUT_H
