#include "EasiMeshSize.h"
#include <utils/logger.h>

#include <tuple>
#include <utility>

EasiMeshSize::EasiMeshSize() : parser(nullptr), model(nullptr), query(easi::Query(1, 3)) {}

EasiMeshSize::EasiMeshSize(VelocityAwareRefinementSettings refinementSettings, pGModel simModel,
                           std::unordered_map<pGRegion, int> groupMap)
    : refinementSettings(refinementSettings), simModel(simModel), parser(new easi::YAMLParser(3)),
      groupMap(std::move(groupMap)), query(easi::Query(1, 3)) {
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
EasiMeshSize::getTargetedFrequencyAndRegion(std::array<double, 3> point) {
  const auto& refinementRegions = refinementSettings.getRefinementRegions();
  double targetedFrequency = 0.0;
  int bypassFindRegionAndUseGroup = 0;
  for (const auto& refinementCube : refinementRegions) {
    auto& cuboid = refinementCube.cuboid;
    bool isInCuboid = true;
    for (int i = 0; i < 3; ++i) {
      const auto minCoord = cuboid.center[i] - cuboid.halfSize[i];
      const auto maxCoord = cuboid.center[i] + cuboid.halfSize[i];
      isInCuboid &= minCoord <= point[i] && point[i] <= maxCoord;
    }
    if (isInCuboid) {
      if (refinementCube.targetedFrequency >= targetedFrequency) {
        bypassFindRegionAndUseGroup = refinementCube.bypassFindRegionAndUseGroup;
      }
      targetedFrequency = std::max(targetedFrequency, refinementCube.targetedFrequency);
    }
  }
  return std::make_tuple(targetedFrequency, bypassFindRegionAndUseGroup);
}

double EasiMeshSize::getMeshSize(std::array<double, 3> point) {
  if (!model) {
    logError() << "Model was not parsed correctly";
  }

  constexpr double defaultMeshSize = std::numeric_limits<double>::max();
  auto freqRegion = getTargetedFrequencyAndRegion(point);
  const double targetedFrequency = std::get<0>(freqRegion);
  const int bypassFindRegionAndUseGroup = std::get<1>(freqRegion);
  if (targetedFrequency == 0.0) {
    return defaultMeshSize;
  }

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
