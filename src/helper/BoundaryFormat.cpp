// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "BoundaryFormat.h"

#include <algorithm>

#include "utils/logger.h"

std::vector<std::uint64_t> packBoundaries(const std::vector<int>& faces,
                                          const std::vector<puml::CellType>& cellTypes,
                                          int bitsPerFace) {
  constexpr std::size_t PackedFaces = 4;
  const std::uint64_t largest = (std::uint64_t{1} << bitsPerFace) - 1;
  std::vector<std::uint64_t> packed(cellTypes.size(), 0);
  std::size_t first = 0;
  for (std::size_t cell = 0; cell < cellTypes.size(); ++cell) {
    const auto faceCount = puml::shapeOf(cellTypes[cell]).faceCount;
    for (std::size_t face = 0; face < std::min(faceCount, PackedFaces); ++face) {
      const int value = faces[first + face];
      if (value < 0 || static_cast<std::uint64_t>(value) > largest) {
        logError() << "Cannot handle boundary condition" << value;
      }
      packed[cell] |= static_cast<std::uint64_t>(value) << (face * bitsPerFace);
    }
    first += faceCount;
  }
  return packed;
}

std::vector<std::int32_t> boundariesPerFace(const std::vector<int>& faces,
                                            const std::vector<puml::CellType>& cellTypes,
                                            std::size_t facesPerCell) {
  std::vector<std::int32_t> values(cellTypes.size() * facesPerCell, 0);
  std::size_t first = 0;
  for (std::size_t cell = 0; cell < cellTypes.size(); ++cell) {
    const auto faceCount = puml::shapeOf(cellTypes[cell]).faceCount;
    std::copy_n(&faces[first], std::min(faceCount, facesPerCell), &values[cell * facesPerCell]);
    first += faceCount;
  }
  return values;
}
