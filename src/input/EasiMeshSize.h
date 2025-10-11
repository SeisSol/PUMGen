// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_INPUT_EASIMESHSIZE_H_
#define PUMGEN_SRC_INPUT_EASIMESHSIZE_H_

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

  int findGroup(std::array<double, 3> point);

  std::tuple<const double, const int>
  getTargetedFrequencyAndRegion(const std::array<double, 3>& point);

  public:
  EasiMeshSize();
  ;
  EasiMeshSize(VelocityAwareRefinementSettings refinementSettings, pGModel simModel,
               std::unordered_map<pGRegion, int> groupMap);

  double getMeshSize(const std::array<double, 3>& point);
};

#endif // PUMGEN_SRC_INPUT_EASIMESHSIZE_H_
