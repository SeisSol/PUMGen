#ifndef PUMGEN_EASIMESHSIZE_H
#define PUMGEN_EASIMESHSIZE_H

#include <easi/YAMLParser.h>
#include <easi/Component.h>
#include <MeshTypes.h>
#include <SimModel.h>
#include <array>
#include <memory>
#include <string>

struct ElasticMaterial {
  double lambda, mu, rho;
};

class EasiMeshSize {
private:
  std::string easiFileName;
  double targetedFrequency;
  double elementsPerWaveLength;
  easi::YAMLParser* parser;
  easi::Component* model; // Unique ptr to model leads to segfault
  pGModel simModel;
  std::unordered_map<pGRegion, int> groupMap;

  int findGroup(std::array<double, 3> point);

  public:
  EasiMeshSize();;
  EasiMeshSize(std::string easiFileName, double targetedFrequency, double elementsPerWaveLength,
               pGModel simModel, std::unordered_map<pGRegion, int> groupMap);

  double getMeshSize(std::array<double, 3> point);

  ~EasiMeshSize();
};

#endif // PUMGEN_EASIMESHSIZE_H
