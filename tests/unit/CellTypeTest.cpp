// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <array>
#include <cstddef>
#include <set>
#include <utility>

#include "mesh/CellType.h"
#include "meshreader/GMSHParser.h"

using puml::CellType;

namespace {
// the linear reference cells of gmsh
std::array<std::array<double, 3>, 8> referenceVertices(CellType type) {
  switch (type) {
  case CellType::Hexahedron:
    return {{{-1, -1, -1},
             {1, -1, -1},
             {1, 1, -1},
             {-1, 1, -1},
             {-1, -1, 1},
             {1, -1, 1},
             {1, 1, 1},
             {-1, 1, 1}}};
  case CellType::Wedge:
    return {{{0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}}};
  case CellType::Pyramid:
    return {{{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, 0, 1}}};
  case CellType::Tetrahedron:
  default:
    return {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  }
}
} // namespace

TEST_CASE("The faces and edges of every cell kind fit together") {
  for (const auto& shape : puml::CellShapes) {
    CAPTURE(shape.name);
    CHECK(&puml::shapeOf(shape.type) == &shape);
    CHECK(puml::isCellType(static_cast<std::uint8_t>(shape.type)));
    // Euler's formula for a polyhedron
    CHECK(shape.vertexCount + shape.faceCount == shape.edgeCount + 2);

    std::set<std::pair<std::size_t, std::size_t>> edges;
    for (std::size_t e = 0; e < shape.edgeCount; ++e) {
      const auto [a, b] = shape.edgeVertices[e];
      edges.emplace(std::min(a, b), std::max(a, b));
    }
    REQUIRE(edges.size() == shape.edgeCount);

    // every face is a cycle of edges, and every edge lies in two faces
    std::array<int, puml::MaxCellEdges> incident{};
    std::size_t corners = 0;
    for (std::size_t f = 0; f < shape.faceCount; ++f) {
      const auto n = shape.faceVertexCount[f];
      corners += n;
      for (std::size_t i = 0; i < n; ++i) {
        const auto a = shape.faceVertices[f][i];
        const auto b = shape.faceVertices[f][(i + 1) % n];
        const auto edge = edges.find({std::min(a, b), std::max(a, b)});
        REQUIRE(edge != edges.end());
        ++incident[static_cast<std::size_t>(std::distance(edges.begin(), edge))];
      }
    }
    CHECK(corners == 2 * shape.edgeCount);
    for (std::size_t e = 0; e < shape.edgeCount; ++e) {
      CHECK(incident[e] == 2);
    }

    // the faces point outwards on the reference cell of gmsh
    const auto vertices = referenceVertices(shape.type);
    std::array<double, 3> center{};
    for (std::size_t v = 0; v < shape.vertexCount; ++v) {
      for (int d = 0; d < 3; ++d) {
        center[d] += vertices[v][d] / static_cast<double>(shape.vertexCount);
      }
    }
    for (std::size_t f = 0; f < shape.faceCount; ++f) {
      const auto& a = vertices[shape.faceVertices[f][0]];
      const auto& b = vertices[shape.faceVertices[f][1]];
      const auto& c = vertices[shape.faceVertices[f][2]];
      const std::array<double, 3> u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
      const std::array<double, 3> w{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
      const std::array<double, 3> normal{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                         u[0] * w[1] - u[1] * w[0]};
      const double outwards = normal[0] * (a[0] - center[0]) + normal[1] * (a[1] - center[1]) +
                              normal[2] * (a[2] - center[2]);
      CAPTURE(f);
      CHECK(outwards > 0);
    }
  }
}

TEST_CASE("The MSH element types are classified with their node counts") {
  for (long type = 1; type <= static_cast<long>(puml::GMSHParser::NumTypes); ++type) {
    CAPTURE(type);
    const auto element = puml::gmshElementType(type);
    const auto nodes = puml::GMSHParser::NumNodes[type - 1];
    if (element.kind == puml::GmshElementType::Kind::Cell && element.complete) {
      CHECK(nodes == puml::lagrangeNodeCount(element.cellType, element.order));
    }
    if (element.kind == puml::GmshElementType::Kind::Face) {
      CHECK(nodes >= element.faceVertices);
    }
  }
  CHECK(puml::gmshElementType(4).cellType == CellType::Tetrahedron);
  CHECK(puml::gmshElementType(92).order == 3);
  CHECK(puml::gmshElementType(92).cellType == CellType::Hexahedron);
  CHECK(puml::gmshElementType(14).cellType == CellType::Pyramid);
  CHECK(puml::gmshElementType(13).cellType == CellType::Wedge);
  CHECK_FALSE(puml::gmshElementType(17).complete);
  CHECK(puml::gmshElementType(10).faceVertices == 4);
  CHECK(puml::gmshElementType(1).kind == puml::GmshElementType::Kind::Other);
  CHECK(puml::gmshElementType(15).kind == puml::GmshElementType::Kind::Other);
}
