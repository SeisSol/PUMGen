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

#ifndef APF_NATIVE_H
#define APF_NATIVE_H

#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>

#include "utils/logger.h"

#include "MeshInput.h"

class ApfNative : public MeshInput {
public:
  ApfNative(const char *mesh, const char *model = 0L) {
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
