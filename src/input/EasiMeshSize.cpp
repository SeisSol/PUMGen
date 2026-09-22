// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "EasiMeshSize.h"

#include <cassert>
#include <limits>
#include <utility>

EasiMeshSize::EasiMeshSize(const VelocityAwareRefinementSettings& refinementSettings,
                           pGModel simModel, std::unordered_map<pGRegion, int> groupMap)
    : meshSize(std::make_shared<VelocityAwareMeshSize>(refinementSettings)), simModel(simModel),
      groupMap(std::move(groupMap)) {}

int EasiMeshSize::findGroup(std::array<double, 3> point) {
  // GR_containsPoint can be expensive for large geometry,
  // therefore we bypass it the simple case of one region
  if (groupMap.size() == 1) {
    return 1;
  }
  GRIter regionIt = GM_regionIter(simModel);
  while (pGRegion region = GRIter_next(regionIt)) {
    // Note: Does not work on discrete regions!
    if (GR_containsPoint(region, point.data()) > 0) {
      assert(groupMap.count(region) > 0);
      return groupMap[region];
    }
  }
  return -1;
}

double EasiMeshSize::getMeshSize(const std::array<double, 3>& point) {
  constexpr double DefaultMeshSize = std::numeric_limits<double>::max();
  const auto [targetedFrequency, bypassFindRegionAndUseGroup] =
      meshSize->targetedFrequencyAndGroup(point);
  if (targetedFrequency == 0.0) {
    return DefaultMeshSize;
  }

  const int group =
      bypassFindRegionAndUseGroup != 0 ? bypassFindRegionAndUseGroup : findGroup(point);
  // This means the point is slightly outside the geometry
  if (group <= 0) {
    return DefaultMeshSize;
  }

  return meshSize->meshSizes({point}, {group}, {targetedFrequency})[0];
}
