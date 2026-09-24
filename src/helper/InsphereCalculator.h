// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
#define PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <unordered_map>
#include <vector>

#include "mesh/CellType.h"

/**
 * The vertex coordinates of the local cells; the coordinates of vertices owned by other ranks are
 * fetched once. Assumes a contiguous distribution of all vertices over all processes.
 */
class CellVertices {
  public:
  CellVertices(const std::vector<puml::CellType>& cellTypes,
               const std::vector<std::uint64_t>& connectivity, const std::vector<double>& geometry,
               MPI_Comm comm);

  [[nodiscard]] std::size_t numCells() const { return cellTypes.size(); }

  [[nodiscard]] puml::CellType type(std::size_t cell) const { return cellTypes[cell]; }

  /**
   * The coordinates of the vertices of a cell, as many as its kind has.
   */
  [[nodiscard]] std::array<std::array<double, 3>, puml::MaxCellVertices>
  operator()(std::size_t cell) const;

  private:
  const std::vector<puml::CellType>& cellTypes;
  const std::vector<std::uint64_t>& connectivity;
  const std::vector<double>& geometry;
  /** Where the vertices of every cell begin; empty for cells of one kind */
  std::vector<std::size_t> firstVertex;
  int commrank = 0;
  std::vector<std::size_t> vertexDist;
  std::vector<std::unordered_map<std::size_t, std::size_t>> outidxmap;
  std::vector<std::size_t> outdisp;
  std::vector<double> outvertices;
};

/**
 * The insphere radius of each cell, as three times its volume over its surface area; for a cell
 * with a sphere touching all its faces, as a tetrahedron, that is the radius of that sphere.
 */
std::vector<double> calculateInsphere(const CellVertices& cells);

#endif // PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
