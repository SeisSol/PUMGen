// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <cstdint>
#include <vector>

#include "helper/BoundaryFormat.h"
#include "mesh/CellType.h"

using puml::CellType;

TEST_CASE("Boundary conditions are packed with 8 bits per face (i32)") {
  CHECK(packBoundaries({1, 3, 5, 255}, {CellType::Tetrahedron}, 8)[0] == 0xFF050301ULL);
}

TEST_CASE("Boundary conditions are packed with 16 bits per face (i64), up to the sign bit") {
  CHECK(packBoundaries({1, 0, 0x8000, 0xFFFF}, {CellType::Tetrahedron}, 16)[0] ==
        0xFFFF800000000001ULL);
}

TEST_CASE("Boundary conditions are stored per face, as many as a cell has") {
  const std::vector<int> faces{100, 101, 102, 103, 104, 105, 1, 2, 3, 4};
  const auto values = boundariesPerFace(faces, {CellType::Hexahedron, CellType::Tetrahedron}, 6);
  REQUIRE(values.size() == 12);
  for (int face = 0; face < 6; ++face) {
    CHECK(values[face] == 100 + face);
  }
  CHECK(values[6] == 1);
  CHECK(values[9] == 4);
  CHECK(values[10] == 0);
  CHECK(values[11] == 0);
}
