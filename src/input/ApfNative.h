// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef APF_NATIVE_H
#define APF_NATIVE_H

#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include "utils/logger.h"

#include "MeshInput.h"

class ApfNative : public ApfMeshInput {
  public:
  ApfNative(const char* mesh, int boundarySize, const char* model = 0L)
      : ApfMeshInput(boundarySize) {
    if (model)
      gmi_register_mesh();
    else {
      gmi_register_null();
      model = ".null";
    }
    m_mesh = apf::loadMdsMesh(model, mesh);
  }
};

#endif // APF_NATIVE_H
