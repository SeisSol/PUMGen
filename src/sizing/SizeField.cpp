// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "SizeField.h"

#include <algorithm>
#include <cmath>

StructuredGrid gridForCuboid(const SimpleCuboid& cuboid, double spacing) {
  const auto box = cuboid.boundingBox(2 * spacing);
  StructuredGrid grid{box[0], spacing, {}};
  for (int d = 0; d < 3; ++d) {
    grid.count[d] = static_cast<std::size_t>(std::ceil((box[1][d] - box[0][d]) / spacing)) + 1;
  }
  return grid;
}

std::vector<double> cuboidMeshSizes(const VelocityAwareMeshSize& meshSize,
                                    const VelocityRefinementCube& region,
                                    const StructuredGrid& grid, std::size_t first,
                                    std::size_t last) {
  const std::size_t ny = grid.count[1];
  const std::size_t nz = grid.count[2];
  const std::size_t layer = ny * nz;
  // one layer more on each side for the neighbourhood minimum
  const std::size_t haloFirst = first > 0 ? first - 1 : first;
  const std::size_t haloLast = std::min(last + 1, grid.count[0]);
  const std::size_t numLayers = haloLast - haloFirst;

  std::vector<double> sizes(numLayers * layer, UnconstrainedSize);
  std::vector<std::array<double, 3>> points;
  std::vector<std::size_t> indices;
  for (std::size_t i = haloFirst; i < haloLast; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        const auto point = grid.node(i, j, k);
        if (region.cuboid.contains(point, 2 * grid.spacing)) {
          points.push_back(point);
          indices.push_back(((i - haloFirst) * ny + j) * nz + k);
        }
      }
    }
  }
  const std::vector<int> groups(points.size(), region.bypassFindRegionAndUseGroup);
  const std::vector<double> frequencies(points.size(), region.targetedFrequency);
  const auto values = meshSize.meshSizes(points, groups, frequencies);
  for (std::size_t p = 0; p < points.size(); ++p) {
    sizes[indices[p]] = std::min(values[p], UnconstrainedSize);
  }

  // the minimum over the 3x3x3 neighbourhood, one direction after the other
  std::vector<double> previous;
  const auto minimumAlong = [&](std::size_t count, std::size_t stride) {
    previous = sizes;
    for (std::size_t index = 0; index < sizes.size(); ++index) {
      const std::size_t position = (index / stride) % count;
      if (position > 0) {
        sizes[index] = std::min(sizes[index], previous[index - stride]);
      }
      if (position + 1 < count) {
        sizes[index] = std::min(sizes[index], previous[index + stride]);
      }
    }
  };
  minimumAlong(nz, 1);
  minimumAlong(ny, nz);
  minimumAlong(numLayers, layer);

  // drop the halo layers
  const std::size_t offset = (first - haloFirst) * layer;
  return {sizes.begin() + static_cast<std::ptrdiff_t>(offset),
          sizes.begin() + static_cast<std::ptrdiff_t>(offset + (last - first) * layer)};
}
