// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "EasiMeshSize.h"
#include <unordered_map>
#include <utils/logger.h>

#include <tuple>
#include <utility>

EasiMeshSize::EasiMeshSize() : parser(nullptr), model(nullptr) {}

EasiMeshSize::EasiMeshSize(VelocityAwareRefinementSettings refinementSettings, pGModel simModel,
                           std::unordered_map<pGRegion, int> groupMap)
    : refinementSettings(refinementSettings), simModel(simModel), parser(new easi::YAMLParser(3)),
      groupMap(std::move(groupMap)) {
  model = parser->parse(refinementSettings.getEasiFileName());
}

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

std::tuple<const double, const int>
EasiMeshSize::getTargetedFrequencyAndRegion(const std::array<double, 3>& point) {
  const auto& refinementRegions = refinementSettings.getRefinementRegions();
  double targetedFrequency = 0.0;
  int bypassFindRegionAndUseGroup = 0;
  for (const auto& refinementCube : refinementRegions) {
    auto& cuboid = refinementCube.cuboid;
    bool isInCuboid = true;

    double u0 = (point[0] - cuboid.center[0]) * cuboid.cosSinRotationZ[0] +
                (point[1] - cuboid.center[1]) * cuboid.cosSinRotationZ[1];
    isInCuboid &= std::abs(u0) <= cuboid.halfSize[0];
    double u1 = (point[0] - cuboid.center[0]) * -cuboid.cosSinRotationZ[1] +
                (point[1] - cuboid.center[1]) * cuboid.cosSinRotationZ[0];
    isInCuboid &= std::abs(u1) <= cuboid.halfSize[1];
    double u2 = (point[2] - cuboid.center[2]);
    isInCuboid &= std::abs(u2) <= cuboid.halfSize[2];

    if (isInCuboid) {
      if (refinementCube.targetedFrequency >= targetedFrequency) {
        bypassFindRegionAndUseGroup = refinementCube.bypassFindRegionAndUseGroup;
      }
      targetedFrequency = std::max(targetedFrequency, refinementCube.targetedFrequency);
    }
  }
  return std::make_tuple(targetedFrequency, bypassFindRegionAndUseGroup);
}

double EasiMeshSize::getMeshSize(const std::array<double, 3>& point) {
  if (!model) {
    logError() << "Model was not parsed correctly";
  }

  constexpr double defaultMeshSize = std::numeric_limits<double>::max();
  const auto [targetedFrequency, bypassFindRegionAndUseGroup] =
      getTargetedFrequencyAndRegion(point);

  if (targetedFrequency == 0.0) {
    return defaultMeshSize;
  }

  easi::Query query(1, 3);
  query.x(0, 0) = point[0];
  query.x(0, 1) = point[1];
  query.x(0, 2) = point[2];
  query.group(0) = bypassFindRegionAndUseGroup ? bypassFindRegionAndUseGroup : findGroup(point);

  // This means the point is slightly outside the geometry
  if (query.group(0) <= 0) {
    return defaultMeshSize;
  }

  // Need array for easi interface
  auto materials = std::array<ElasticMaterial, 1>();
  auto& material = materials[0];
  auto adapter = easi::ArrayOfStructsAdapter<ElasticMaterial>(materials.data());
  adapter.addBindingPoint("lambda", &ElasticMaterial::lambda);
  adapter.addBindingPoint("mu", &ElasticMaterial::mu);
  adapter.addBindingPoint("rho", &ElasticMaterial::rho);

  model->evaluate(query, adapter);

  double waveSpeed = 0.0;
  if (material.mu < 10e-14) {
    // acoustic
    waveSpeed = std::sqrt((material.lambda + 2 * material.mu) / material.rho);
  } else {
    // elastic
    waveSpeed = std::sqrt(material.mu / material.rho);
  }

  const auto waveLength = waveSpeed / targetedFrequency;
  return waveLength / refinementSettings.getElementsPerWaveLength();
}
