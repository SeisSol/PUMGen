// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_SIZING_VELOCITYAWARESETTINGS_H_
#define PUMGEN_SRC_SIZING_VELOCITYAWARESETTINGS_H_

#include <array>
#include <string>
#include <vector>

namespace tinyxml2 {
class XMLDocument;
} // namespace tinyxml2

/**
 * A cuboid, rotated around the z axis.
 */
struct SimpleCuboid {
  std::array<double, 3> center;
  std::array<double, 3> halfSize;
  std::array<double, 2> cosSinRotationZ;
  double rotationZ;

  /**
   * Whether the point lies in the cuboid, enlarged by margin in each direction.
   */
  [[nodiscard]] bool contains(const std::array<double, 3>& point, double margin = 0) const;

  /**
   * The axis-aligned bounding box (minimum and maximum corner) of the enlarged cuboid.
   */
  [[nodiscard]] std::array<std::array<double, 3>, 2> boundingBox(double margin = 0) const;
};

struct VelocityRefinementCube {
  VelocityRefinementCube(SimpleCuboid cuboid, double targetedFrequency,
                         int bypassFindRegionAndUseGroup)
      : cuboid(cuboid), targetedFrequency(targetedFrequency),
        bypassFindRegionAndUseGroup(bypassFindRegionAndUseGroup) {}

  SimpleCuboid cuboid;
  double targetedFrequency;
  int bypassFindRegionAndUseGroup;
};

class VelocityAwareRefinementSettings {
  public:
  VelocityAwareRefinementSettings() = default;
  VelocityAwareRefinementSettings(double elementsPerWaveLength, std::string easiFileName);

  void addRefinementRegion(SimpleCuboid cuboid, double targetedFrequency,
                           int bypassFindRegionAndUseGroup);

  [[nodiscard]] bool isVelocityAwareRefinementOn() const;

  [[nodiscard]] const std::string& getEasiFileName() const;

  [[nodiscard]] double getElementsPerWaveLength() const;

  [[nodiscard]] const std::vector<VelocityRefinementCube>& getRefinementRegions() const;

  private:
  double elementsPerWaveLength{};
  std::string easiFileName;
  std::vector<VelocityRefinementCube> refinementRegions{};
};

/**
 * Reads the VelocityAwareMeshing element of a mesh attributes document; the settings are empty
 * (isVelocityAwareRefinementOn() is false) if the document has none.
 */
VelocityAwareRefinementSettings readVelocityAwareSettings(const tinyxml2::XMLDocument& doc);

#endif // PUMGEN_SRC_SIZING_VELOCITYAWARESETTINGS_H_
