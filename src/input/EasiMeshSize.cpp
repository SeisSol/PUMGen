#include "EasiMeshSize.h"
#include <easi/YAMLParser.h>
#include <utils/logger.h>


EasiMeshSize::EasiMeshSize() : parser(nullptr), model(nullptr) {
  std::cout << "EasiMeshSize() called" << std::endl;
}

EasiMeshSize::EasiMeshSize(std::string easiFileName, double targetedFrequency,
                           double elementsPerWaveLength, pGModel simModel,
                           std::unordered_map<pGRegion, int> groupMap)
    :
    easiFileName(easiFileName), targetedFrequency(targetedFrequency),
    elementsPerWaveLength(elementsPerWaveLength),
    simModel(simModel),
    parser(new easi::YAMLParser(3)),
    groupMap(groupMap) {
  std::cout << "EasiMeshSize(easiFileName) called" << std::endl;
  model = parser->parse(easiFileName);
}

int EasiMeshSize::findGroup(std::array<double, 3> point) {
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

double EasiMeshSize::getMeshSize(std::array<double, 3> point) {
  if (!model) {
    logError() << "Model was not parsed correctly";
  }
  auto query = easi::Query(1,3);
  query.x(0,0) = point[0];
  query.x(0,1) = point[1];
  query.x(0,2) = point[2];
  query.group(0) = findGroup(point);
  assert(query.group(0) > 0);

  // Need array for easi interface
  auto materials = std::array<ElasticMaterial,1>();
  auto& material = materials[0];
  auto adapter = easi::ArrayOfStructsAdapter<ElasticMaterial>(materials.data());
  adapter.addBindingPoint("lambda", &ElasticMaterial::lambda);
  adapter.addBindingPoint("mu", &ElasticMaterial::mu);
  adapter.addBindingPoint("rho", &ElasticMaterial::rho);

  model->evaluate(query, adapter);

  const auto pWaveSpeed = std::sqrt(
      (material.lambda + 2 * material.mu ) / material.rho
      );

  const auto waveLength = pWaveSpeed / targetedFrequency;
  return waveLength / elementsPerWaveLength;
}
EasiMeshSize::~EasiMeshSize() {
  std::cout << "Destructor of EasiMeshSize called" << std::endl;
  //delete parser;
  //delete model;
}

