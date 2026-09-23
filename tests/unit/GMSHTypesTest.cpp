// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <cstddef>

#include "meshreader/GMSHParser.h"

namespace {
std::size_t nodesOf(std::size_t type) { return puml::GMSHParser::NumNodes[type - 1]; }
} // namespace

TEST_CASE("The node counts of the MSH element types are the ones of gmsh") {
  // a sample over the whole table, as gmsh 4.15.2 reports them
  CHECK(nodesOf(4) == 4);     // tetrahedron
  CHECK(nodesOf(11) == 10);   // tetrahedron of order 2
  CHECK(nodesOf(12) == 27);   // hexahedron of order 2
  CHECK(nodesOf(13) == 18);   // prism of order 2
  CHECK(nodesOf(14) == 14);   // pyramid of order 2
  CHECK(nodesOf(75) == 286);  // tetrahedron of order 10
  CHECK(nodesOf(76) == 0);    // no element type
  CHECK(nodesOf(79) == 34);   // incomplete tetrahedron of order 4
  CHECK(nodesOf(84) == 1);    // line of order 0
  CHECK(nodesOf(92) == 64);   // hexahedron of order 3
  CHECK(nodesOf(98) == 1000); // hexahedron of order 9
  CHECK(nodesOf(118) == 30);  // pyramid of order 3
  CHECK(nodesOf(131) == 69);  // incomplete pyramid of order 6
}
