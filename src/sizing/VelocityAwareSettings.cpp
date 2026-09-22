// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "VelocityAwareSettings.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

#include "tinyxml2/tinyxml2.h"
#include "utils/logger.h"

namespace {
const double ToRadians = std::acos(-1) / 180.0;
} // namespace

bool SimpleCuboid::contains(const std::array<double, 3>& point, double margin) const {
  const double u0 =
      (point[0] - center[0]) * cosSinRotationZ[0] + (point[1] - center[1]) * cosSinRotationZ[1];
  const double u1 =
      (point[0] - center[0]) * -cosSinRotationZ[1] + (point[1] - center[1]) * cosSinRotationZ[0];
  const double u2 = point[2] - center[2];
  return std::abs(u0) <= halfSize[0] + margin && std::abs(u1) <= halfSize[1] + margin &&
         std::abs(u2) <= halfSize[2] + margin;
}

std::array<std::array<double, 3>, 2> SimpleCuboid::boundingBox(double margin) const {
  const double hx = halfSize[0] + margin;
  const double hy = halfSize[1] + margin;
  const double hz = halfSize[2] + margin;
  const double c = std::abs(cosSinRotationZ[0]);
  const double s = std::abs(cosSinRotationZ[1]);
  const double ex = hx * c + hy * s;
  const double ey = hx * s + hy * c;
  return {{{center[0] - ex, center[1] - ey, center[2] - hz},
           {center[0] + ex, center[1] + ey, center[2] + hz}}};
}

VelocityAwareRefinementSettings::VelocityAwareRefinementSettings(double elementsPerWaveLength,
                                                                 std::string easiFileName)
    : elementsPerWaveLength(elementsPerWaveLength), easiFileName(std::move(easiFileName)) {}

void VelocityAwareRefinementSettings::addRefinementRegion(SimpleCuboid cuboid,
                                                          double targetedFrequency,
                                                          int bypassFindRegionAndUseGroup) {
  refinementRegions.emplace_back(cuboid, targetedFrequency, bypassFindRegionAndUseGroup);
}

bool VelocityAwareRefinementSettings::isVelocityAwareRefinementOn() const {
  return !refinementRegions.empty();
}

const std::string& VelocityAwareRefinementSettings::getEasiFileName() const { return easiFileName; }

double VelocityAwareRefinementSettings::getElementsPerWaveLength() const {
  return elementsPerWaveLength;
}

const std::vector<VelocityRefinementCube>&
VelocityAwareRefinementSettings::getRefinementRegions() const {
  return refinementRegions;
}

VelocityAwareRefinementSettings readVelocityAwareSettings(const tinyxml2::XMLDocument& doc) {
  VelocityAwareRefinementSettings settings;
  int numChilds = 0;
  const auto name = "VelocityAwareMeshing";
  for (auto velocityAwareMeshingElement = doc.FirstChildElement(name); velocityAwareMeshingElement;
       velocityAwareMeshingElement = velocityAwareMeshingElement->NextSiblingElement(name)) {
    const auto easiFileName = velocityAwareMeshingElement->Attribute("easiFile");
    const auto elementsPerWaveLength =
        std::stof(velocityAwareMeshingElement->Attribute("elementsPerWaveLength"));
    settings = VelocityAwareRefinementSettings(elementsPerWaveLength, easiFileName);
    constexpr auto cuboidName = "VelocityRefinementCuboid";
    logInfo() << "Activating velocity aware meshing, using" << elementsPerWaveLength
              << "elements per wavelength and easi file" << easiFileName;
    for (auto child = velocityAwareMeshingElement->FirstChildElement(cuboidName); child;
         child = child->NextSiblingElement(cuboidName)) {

      const char* attr = child->Attribute("rotationZAnticlockwiseFromX");
      double rotationZAnticlockwiseFromX = 0.0;
      if (attr != nullptr) {
        rotationZAnticlockwiseFromX = std::stof(attr);
      }
      attr = child->Attribute("bypassFindRegionAndUseGroup");
      int bypassFindRegionAndUseGroup = 0;
      if (attr != nullptr) {
        bypassFindRegionAndUseGroup = std::stoi(attr);
      }
      auto cuboid = SimpleCuboid{{
                                     std::stof(child->Attribute("centerX")),
                                     std::stof(child->Attribute("centerY")),
                                     std::stof(child->Attribute("centerZ")),
                                 },
                                 {
                                     std::stof(child->Attribute("halfSizeX")),
                                     std::stof(child->Attribute("halfSizeY")),
                                     std::stof(child->Attribute("halfSizeZ")),
                                 },
                                 {std::cos(rotationZAnticlockwiseFromX * ToRadians),
                                  std::sin(rotationZAnticlockwiseFromX * ToRadians)},
                                 rotationZAnticlockwiseFromX};
      const auto targetedFrequency = std::stof(child->Attribute("frequency"));

      settings.addRefinementRegion(cuboid, targetedFrequency, bypassFindRegionAndUseGroup);
      logInfo() << "Adding velocity aware refinement region targeting" << targetedFrequency
                << "Hz, centered at x =" << cuboid.center[0] << "y=" << cuboid.center[1]
                << "z=" << cuboid.center[2] << "with half sizes"
                << "x =" << cuboid.halfSize[0] << "y =" << cuboid.halfSize[1]
                << "z =" << cuboid.halfSize[2];
      if (std::abs(cuboid.rotationZ) > 0.0) {
        logInfo() << "rotated around z axis by " << cuboid.rotationZ
                  << "degree(s) counterclockwise from x axis.";
      }
      if (bypassFindRegionAndUseGroup) {
        logInfo() << "bypass findRegion and use group =" << bypassFindRegionAndUseGroup;
      }
    }
    if (!settings.isVelocityAwareRefinementOn()) {
      logWarning() << "Activated velocity aware meshing but did not specify any refinement region!";
    }
    ++numChilds;
  }
  if (numChilds > 1) {
    logError() << "Multiple definitions of velocityAwareMeshing";
  }
  return settings;
}
