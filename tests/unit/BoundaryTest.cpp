// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <cstdint>
#include <cstring>

#include "input/MeshData.h"

namespace {
class OneCell : public FullStorageMeshData {
  public:
  explicit OneCell(int boundarySize) : FullStorageMeshData(boundarySize) { setup(1, 4); }
  using FullStorageMeshData::setBoundary;
};

std::uint64_t bitsOf(std::int64_t value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}
} // namespace

TEST_CASE("Boundary conditions are packed with 8 bits per face (i32)") {
  OneCell mesh(8);
  mesh.setBoundary(0, 0, 1);
  mesh.setBoundary(0, 1, 3);
  mesh.setBoundary(0, 2, 5);
  mesh.setBoundary(0, 3, 255);
  CHECK(bitsOf(mesh.boundary()[0]) == 0xFF050301ULL);
}

TEST_CASE("Boundary conditions are packed with 16 bits per face (i64), up to the sign bit") {
  OneCell mesh(16);
  mesh.setBoundary(0, 0, 1);
  mesh.setBoundary(0, 2, 0x8000);
  mesh.setBoundary(0, 3, 0xFFFF);
  CHECK(bitsOf(mesh.boundary()[0]) == 0xFFFF800000000001ULL);
}

TEST_CASE("Boundary conditions are stored per face (i32x4)") {
  OneCell mesh(-1);
  for (int face = 0; face < 4; ++face) {
    mesh.setBoundary(0, face, 100 + face);
  }
  for (int face = 0; face < 4; ++face) {
    CHECK(mesh.boundary()[face] == 100 + face);
  }
}
