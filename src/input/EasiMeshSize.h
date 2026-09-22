// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_INPUT_EASIMESHSIZE_H_
#define PUMGEN_SRC_INPUT_EASIMESHSIZE_H_

#include "sizing/VelocityAwareMeshSize.h"
#include <MeshTypes.h>
#include <SimModel.h>
#include <array>
#include <memory>
#include <unordered_map>

/**
 * Velocity-aware mesh size for SimModSuite, which asks for one point at a time.
 */
class EasiMeshSize {
  private:
  std::shared_ptr<VelocityAwareMeshSize> meshSize;
  pGModel simModel;
  std::unordered_map<pGRegion, int> groupMap;

  int findGroup(std::array<double, 3> point);

  public:
  EasiMeshSize(const VelocityAwareRefinementSettings& refinementSettings, pGModel simModel,
               std::unordered_map<pGRegion, int> groupMap);

  double getMeshSize(const std::array<double, 3>& point);
};

#endif // PUMGEN_SRC_INPUT_EASIMESHSIZE_H_
