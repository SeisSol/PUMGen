#ifndef PUMGEN_EASIMESHSIZE_H
#define PUMGEN_EASIMESHSIZE_H

#include "MeshAttributes.h"
#include <MeshTypes.h>
#include <SimModel.h>
#include <array>
#include <easi/Component.h>
#include <easi/YAMLParser.h>
#include <memory>
#include <string>

struct ElasticMaterial {
  double lambda, mu, rho;
};

class EasiMeshSize {
  private:
  VelocityAwareRefinementSettings refinementSettings;

  easi::YAMLParser* parser;
  easi::Component* model; // Unique ptr to model leads to segfault
  pGModel simModel;
  std::unordered_map<pGRegion, int> groupMap;

  int findGroup(const std::array<double, 3>& point);

  std::tuple<const double, const int>
  getTargetedFrequencyAndRegion(const std::array<double, 3>& point);

  public:
  EasiMeshSize();
  ;
  EasiMeshSize(VelocityAwareRefinementSettings refinementSettings, pGModel simModel,
               std::unordered_map<pGRegion, int> groupMap);

  double getMeshSize(const std::array<double, 3>& point);
};

#endif // PUMGEN_EASIMESHSIZE_H
