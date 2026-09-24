// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "CellType.h"

#include <algorithm>
#include <initializer_list>

namespace puml {

namespace {
struct CellTypes {
  CellType type;
  // the MSH types of the complete cells, from order 1 on
  std::initializer_list<long> complete;
  // the MSH types of the other cells of this kind (incomplete, and of order 0)
  std::initializer_list<long> incomplete;
};

const CellTypes MshCells[] = {
    {CellType::Tetrahedron,
     {4, 11, 29, 30, 31, 71, 72, 73, 74, 75},
     {32, 33, 79, 80, 81, 82, 83, 87}},
    {CellType::Hexahedron,
     {5, 12, 92, 93, 94, 95, 96, 97, 98},
     {17, 88, 99, 100, 101, 102, 103, 104, 105}},
    {CellType::Wedge,
     {6, 13, 90, 91, 106, 107, 108, 109, 110},
     {18, 89, 111, 112, 113, 114, 115, 116, 117}},
    {CellType::Pyramid,
     {7, 14, 118, 119, 120, 121, 122, 123, 124},
     {19, 125, 126, 127, 128, 129, 130, 131, 132}},
};

// the orders of the triangles and quadrilaterals; incomplete faces are given with the order of
// their edges, their first nodes are the vertices all the same
const std::pair<long, std::size_t> MshTriangles[] = {
    {2, 1},  {9, 2},  {20, 3}, {21, 3},  {22, 4}, {23, 4}, {24, 5}, {25, 5}, {42, 6},
    {43, 7}, {44, 8}, {45, 9}, {46, 10}, {52, 6}, {53, 7}, {54, 8}, {55, 9}, {56, 10}};
const std::pair<long, std::size_t> MshQuadrilaterals[] = {
    {3, 1},  {10, 2}, {16, 2}, {36, 3},  {37, 4}, {38, 5}, {39, 3}, {40, 4}, {41, 5}, {47, 6},
    {48, 7}, {49, 8}, {50, 9}, {51, 10}, {57, 6}, {58, 7}, {59, 8}, {60, 9}, {61, 10}};
} // namespace

GmshElementType gmshElementType(long type) {
  GmshElementType result;
  for (const auto& cells : MshCells) {
    std::size_t order = 1;
    for (const long complete : cells.complete) {
      if (complete == type) {
        result.kind = GmshElementType::Kind::Cell;
        result.cellType = cells.type;
        result.order = order;
        result.complete = true;
        return result;
      }
      ++order;
    }
    if (std::find(cells.incomplete.begin(), cells.incomplete.end(), type) !=
        cells.incomplete.end()) {
      result.kind = GmshElementType::Kind::Cell;
      result.cellType = cells.type;
      return result;
    }
  }
  for (const auto& [msh, order] : MshTriangles) {
    if (msh == type) {
      result.kind = GmshElementType::Kind::Face;
      result.faceVertices = 3;
      result.order = order;
      return result;
    }
  }
  for (const auto& [msh, order] : MshQuadrilaterals) {
    if (msh == type) {
      result.kind = GmshElementType::Kind::Face;
      result.faceVertices = 4;
      result.order = order;
      return result;
    }
  }
  return result;
}

} // namespace puml
