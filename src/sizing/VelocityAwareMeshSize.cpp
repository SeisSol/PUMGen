// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "VelocityAwareMeshSize.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <easi/Component.h>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>

#include "utils/logger.h"

namespace {
struct ElasticMaterial {
  double lambda;
  double mu;
  double rho;
};

// the material model is evaluated in batches of this many points
constexpr std::size_t QueryBatch = std::size_t{1} << 20;
} // namespace

VelocityAwareMeshSize::VelocityAwareMeshSize(VelocityAwareRefinementSettings settings)
    : settings(std::move(settings)), parser(std::make_unique<easi::YAMLParser>(3)) {
  model = parser->parse(this->settings.getEasiFileName());
  if (model == nullptr) {
    logError() << "Could not parse the easi file" << this->settings.getEasiFileName();
  }
}

VelocityAwareMeshSize::~VelocityAwareMeshSize() { delete model; }

std::pair<double, int>
VelocityAwareMeshSize::targetedFrequencyAndGroup(const std::array<double, 3>& point) const {
  double targetedFrequency = 0.0;
  int bypassFindRegionAndUseGroup = 0;
  for (const auto& refinementCube : settings.getRefinementRegions()) {
    if (refinementCube.cuboid.contains(point)) {
      if (refinementCube.targetedFrequency >= targetedFrequency) {
        bypassFindRegionAndUseGroup = refinementCube.bypassFindRegionAndUseGroup;
      }
      targetedFrequency = std::max(targetedFrequency, refinementCube.targetedFrequency);
    }
  }
  return {targetedFrequency, bypassFindRegionAndUseGroup};
}

std::vector<double>
VelocityAwareMeshSize::meshSizes(const std::vector<std::array<double, 3>>& points,
                                 const std::vector<int>& groups,
                                 const std::vector<double>& frequencies) const {
  std::vector<double> sizes(points.size(), std::numeric_limits<double>::max());

  std::vector<std::size_t> evaluated;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (frequencies[i] > 0.0 && groups[i] > 0) {
      evaluated.push_back(i);
    }
  }

  std::vector<ElasticMaterial> materials;
  for (std::size_t first = 0; first < evaluated.size(); first += QueryBatch) {
    const std::size_t count = std::min(QueryBatch, evaluated.size() - first);
    easi::Query query(count, 3);
    for (std::size_t j = 0; j < count; ++j) {
      const auto i = evaluated[first + j];
      for (int k = 0; k < 3; ++k) {
        query.x(j, k) = points[i][k];
      }
      query.group(j) = groups[i];
    }

    materials.assign(count, ElasticMaterial{});
    easi::ArrayOfStructsAdapter<ElasticMaterial> adapter(materials.data());
    adapter.addBindingPoint("lambda", &ElasticMaterial::lambda);
    adapter.addBindingPoint("mu", &ElasticMaterial::mu);
    adapter.addBindingPoint("rho", &ElasticMaterial::rho);
    model->evaluate(query, adapter);

    for (std::size_t j = 0; j < count; ++j) {
      const auto& material = materials[j];
      // the shear wave speed, or the P wave speed in acoustic material
      const double waveSpeed = material.mu < 10e-14
                                   ? std::sqrt((material.lambda + 2 * material.mu) / material.rho)
                                   : std::sqrt(material.mu / material.rho);
      const auto i = evaluated[first + j];
      const double waveLength = waveSpeed / frequencies[i];
      sizes[i] = waveLength / settings.getElementsPerWaveLength();
    }
  }
  return sizes;
}
