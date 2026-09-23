// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESH_CELLTYPE_H_
#define PUMGEN_SRC_MESH_CELLTYPE_H_

#include <array>
#include <cstddef>
#include <cstdint>

namespace puml {

/**
 * The kind of a cell. The values are the cell types of VTK, which the Types array of a VTKHDF
 * file holds and which PUML reads.
 */
enum class CellType : std::uint8_t {
  Tetrahedron = 10,
  Hexahedron = 12,
  Wedge = 13,
  Pyramid = 14,
};

constexpr std::size_t MaxCellVertices = 8;
constexpr std::size_t MaxCellFaces = 6;
constexpr std::size_t MaxCellEdges = 12;
constexpr std::size_t MaxFaceVertices = 4;

/**
 * How a cell of one kind is built from its vertices.
 *
 * The vertices are in the order of gmsh, which is the one of VTK for these linear cells. The faces
 * are numbered as PUML numbers them, and every face lists its vertices such that its normal points
 * out of a positively oriented cell.
 */
struct CellShape {
  CellType type;
  /** As written to the cell-type attribute */
  const char* name;
  /** The XDMF TopologyType of a mesh of cells of this kind only */
  const char* xdmfName;
  std::size_t vertexCount;
  std::size_t faceCount;
  std::array<std::size_t, MaxCellFaces> faceVertexCount;
  std::array<std::array<std::size_t, MaxFaceVertices>, MaxCellFaces> faceVertices;
  std::size_t edgeCount;
  std::array<std::array<std::size_t, 2>, MaxCellEdges> edgeVertices;
};

inline constexpr std::array<CellShape, 4> CellShapes = {{
    {CellType::Tetrahedron,
     "tetrahedron",
     "Tetrahedron",
     4,
     4,
     {3, 3, 3, 3, 0, 0},
     {{{1, 0, 2, 0}, {0, 1, 3, 0}, {1, 2, 3, 0}, {2, 0, 3, 0}, {}, {}}},
     6,
     {{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}}},
    {CellType::Hexahedron,
     "hexahedron",
     "Hexahedron",
     8,
     6,
     {4, 4, 4, 4, 4, 4},
     {{{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}}},
     12,
     {{{0, 1},
       {1, 2},
       {2, 3},
       {3, 0},
       {4, 5},
       {5, 6},
       {6, 7},
       {7, 4},
       {0, 4},
       {1, 5},
       {2, 6},
       {3, 7}}}},
    {CellType::Wedge,
     "wedge",
     "Wedge",
     6,
     5,
     {3, 3, 4, 4, 4, 0},
     {{{0, 2, 1, 0}, {3, 4, 5, 0}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {}}},
     9,
     {{{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}}}},
    {CellType::Pyramid,
     "pyramid",
     "Pyramid",
     5,
     5,
     {4, 3, 3, 3, 3, 0},
     {{{0, 3, 2, 1}, {0, 1, 4, 0}, {1, 2, 4, 0}, {2, 3, 4, 0}, {3, 0, 4, 0}, {}}},
     8,
     {{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}}}},
}};

/**
 * Whether the value is one of the cell kinds.
 */
constexpr bool isCellType(std::uint8_t value) {
  for (const auto& shape : CellShapes) {
    if (static_cast<std::uint8_t>(shape.type) == value) {
      return true;
    }
  }
  return false;
}

constexpr const CellShape& shapeOf(CellType type) {
  switch (type) {
  case CellType::Hexahedron:
    return CellShapes[1];
  case CellType::Wedge:
    return CellShapes[2];
  case CellType::Pyramid:
    return CellShapes[3];
  case CellType::Tetrahedron:
  default:
    return CellShapes[0];
  }
}

/**
 * The number of nodes of a complete Lagrange cell of the given kind and order.
 */
constexpr std::size_t lagrangeNodeCount(CellType type, std::size_t order) {
  const std::size_t p = order;
  switch (type) {
  case CellType::Hexahedron:
    return (p + 1) * (p + 1) * (p + 1);
  case CellType::Wedge:
    return (p + 1) * (p + 1) * (p + 2) / 2;
  case CellType::Pyramid:
    return (p + 1) * (p + 2) * (2 * p + 3) / 6;
  case CellType::Tetrahedron:
  default:
    return (p + 1) * (p + 2) * (p + 3) / 6;
  }
}

/**
 * What an element type of the MSH format describes.
 */
struct GmshElementType {
  enum class Kind { Cell, Face, Other };
  Kind kind = Kind::Other;
  /** The kind of a cell */
  CellType cellType = CellType::Tetrahedron;
  /** The number of vertices of a face: 3 or 4 */
  std::size_t faceVertices = 0;
  std::size_t order = 0;
  /** Whether a cell holds all nodes of the Lagrange cell of its order */
  bool complete = false;
};

/**
 * Classifies an element type of the MSH format. Cells are the three-dimensional types, faces the
 * triangles and quadrilaterals of any order, whose first nodes are their vertices; all other types
 * (points, lines, and the elements of a single node) are of the kind Other.
 */
GmshElementType gmshElementType(long type);

} // namespace puml

#endif // PUMGEN_SRC_MESH_CELLTYPE_H_
