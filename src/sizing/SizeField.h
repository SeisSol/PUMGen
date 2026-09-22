// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_SIZING_SIZEFIELD_H_
#define PUMGEN_SRC_SIZING_SIZEFIELD_H_

#include <array>
#include <cstddef>
#include <vector>

#include "VelocityAwareMeshSize.h"

/**
 * A regular grid for gmsh's Structured field: the value of node (i, j, k) is stored at
 * (i * count[1] + j) * count[2] + k.
 */
struct StructuredGrid {
  std::array<double, 3> origin;
  double spacing;
  std::array<std::size_t, 3> count;

  [[nodiscard]] std::array<double, 3> node(std::size_t i, std::size_t j, std::size_t k) const {
    return {origin[0] + static_cast<double>(i) * spacing,
            origin[1] + static_cast<double>(j) * spacing,
            origin[2] + static_cast<double>(k) * spacing};
  }
};

/**
 * The value of grid nodes which do not constrain the mesh size.
 */
constexpr double UnconstrainedSize = 1e22;

/**
 * The grid covering a cuboid with a margin of two grid spacings, so that all grid cells which
 * intersect the cuboid have their corners within the margin.
 */
StructuredGrid gridForCuboid(const SimpleCuboid& cuboid, double spacing);

/**
 * The mesh sizes of the grid nodes with x index in [first, last) for one refinement cuboid.
 *
 * The nodes within the margin of the cuboid get the mesh size of the material of the given group;
 * the others get UnconstrainedSize. Each node then takes the minimum of its 3x3x3 neighbourhood,
 * so that the trilinear interpolation in a grid cell never exceeds the smallest size sampled at the
 * cell's corners.
 */
std::vector<double> cuboidMeshSizes(const VelocityAwareMeshSize& meshSize,
                                    const VelocityRefinementCube& region,
                                    const StructuredGrid& grid, std::size_t first,
                                    std::size_t last);

#endif // PUMGEN_SRC_SIZING_SIZEFIELD_H_
