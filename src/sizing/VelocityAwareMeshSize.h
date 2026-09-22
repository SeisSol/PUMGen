// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_SIZING_VELOCITYAWAREMESHSIZE_H_
#define PUMGEN_SRC_SIZING_VELOCITYAWAREMESHSIZE_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "VelocityAwareSettings.h"

namespace easi {
class Component;
class YAMLParser;
} // namespace easi

/**
 * Mesh sizes from the wave speeds of an easi material model: the wavelength of the targeted
 * frequency divided by the number of elements per wavelength.
 */
class VelocityAwareMeshSize {
  public:
  explicit VelocityAwareMeshSize(VelocityAwareRefinementSettings settings);
  ~VelocityAwareMeshSize();
  VelocityAwareMeshSize(const VelocityAwareMeshSize&) = delete;
  VelocityAwareMeshSize& operator=(const VelocityAwareMeshSize&) = delete;

  /**
   * The highest targeted frequency of the refinement regions containing the point (0 if there is
   * none) and the group fixed by that region (0 if the region does not fix one).
   */
  [[nodiscard]] std::pair<double, int>
  targetedFrequencyAndGroup(const std::array<double, 3>& point) const;

  /**
   * The mesh sizes at the points for the given frequencies, with the material of the given
   * groups; points with a frequency of 0 or a group <= 0 get std::numeric_limits<double>::max().
   * The material model is evaluated for all points at once.
   */
  [[nodiscard]] std::vector<double> meshSizes(const std::vector<std::array<double, 3>>& points,
                                              const std::vector<int>& groups,
                                              const std::vector<double>& frequencies) const;

  [[nodiscard]] const VelocityAwareRefinementSettings& getSettings() const { return settings; }

  private:
  VelocityAwareRefinementSettings settings;
  std::unique_ptr<easi::YAMLParser> parser;
  easi::Component* model = nullptr;
};

#endif // PUMGEN_SRC_SIZING_VELOCITYAWAREMESHSIZE_H_
