#include "EasiMeshSize.h"
#include "easi/YAMLParser.h"
#include <utils/logger.h>


EasiMeshSize::EasiMeshSize() : parser(nullptr), model(nullptr) {
  std::cout << "EasiMeshSize() called" << std::endl;
}

EasiMeshSize::EasiMeshSize(std::string easiFileName) :
    parser(new easi::YAMLParser(3)), easiFileName(easiFileName){
  std::cout << "EasiMeshSize(easiFileName) called" << std::endl;
  model = parser->parse(easiFileName);
}
double EasiMeshSize::getMeshSize(std::array<double, 3> point) {
  if (!model) {
    logError() << "Model was not parsed correctly";
  }
  auto query = easi::Query(1,3);
  query.x(0,0) = point[0];
  query.x(0,1) = point[1];
  query.x(0,2) = point[2];
  query.group(0) = 0;

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

  std::cout << pWaveSpeed << std::endl;
  // TODO: What to do about groups?
  if(point[2] > 0 ) {
    return 1.0;
  } else {
    return 0.1;
  }
}
EasiMeshSize::~EasiMeshSize() {
  std::cout << "Destructor of EasiMeshSize called" << std::endl;
  delete parser;
  delete model;
}
