// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
#define PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_

#include <array>
#include <cstddef>
#include <mpi.h>
#include <unordered_map>
#include <vector>

/**
 * The vertex coordinates of the local cells; the coordinates of vertices owned by other ranks are
 * fetched once. Assumes a contiguous distribution of all vertices over all processes. Each cell
 * has cellSize nodes, the first four of which are its vertices.
 */
class CellVertices {
  public:
  static constexpr std::size_t VerticesPerCell = 4;

  CellVertices(const std::vector<std::size_t>& connectivity, const std::vector<double>& geometry,
               std::size_t cellSize, MPI_Comm comm);

  [[nodiscard]] std::size_t numCells() const { return connectivity.size() / cellSize; }

  [[nodiscard]] std::array<std::array<double, 3>, VerticesPerCell>
  operator()(std::size_t cell) const;

  private:
  const std::vector<std::size_t>& connectivity;
  const std::vector<double>& geometry;
  std::size_t cellSize;
  int commrank = 0;
  std::vector<std::size_t> vertexDist;
  std::vector<std::unordered_map<std::size_t, std::size_t>> outidxmap;
  std::vector<std::size_t> outdisp;
  std::vector<double> outvertices;
};

/**
 * The insphere radius of each tetrahedron.
 */
std::vector<double> calculateInsphere(const CellVertices& cells);

#endif // PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
