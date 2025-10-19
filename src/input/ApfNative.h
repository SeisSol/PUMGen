// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_APFNATIVE_H_
#define PUMGEN_SRC_INPUT_APFNATIVE_H_

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

#endif // PUMGEN_SRC_INPUT_APFNATIVE_H_
