// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_NETCDFPARTITION_H_
#define PUMGEN_SRC_INPUT_NETCDFPARTITION_H_

#include <algorithm>
#include <cstddef>
#include <vector>

/**
 * Describes one partition (required for reading netCDF meshes)
 */
class Partition {
  private:
  std::size_t m_nElements = 0;
  std::size_t m_nVertices = 0;

  std::vector<int> m_elements;
  std::vector<double> m_vertices;
  std::vector<int> m_boundaries;
  std::vector<int> m_groups;

  public:
  void setElemSize(std::size_t nElements) {
    if (m_nElements != 0)
      return;

    m_nElements = nElements;

    m_elements.resize(nElements * 4);
    m_boundaries.resize(nElements * 4);
    // group 0 unless the file gives one
    m_groups.assign(nElements, 0);
  }

  void setVrtxSize(std::size_t nVertices) {
    if (m_nVertices != 0)
      return;

    m_nVertices = nVertices;

    m_vertices.resize(nVertices * 3);
  }

  void convertBoundary() {
    int ncBoundaries[4];

    for (std::size_t i = 0; i < m_nElements * 4; i += 4) {
      std::copy_n(&m_boundaries[i], 4, ncBoundaries);
      for (unsigned int j = 0; j < 4; j++)
        m_boundaries[i + j] = ncBoundaries[INTERNAL2EX_ORDER[j]];
    }
  }

  std::size_t nElements() const { return m_nElements; }

  std::size_t nVertices() const { return m_nVertices; }

  int* elements() { return m_elements.data(); }

  double* vertices() { return m_vertices.data(); }

  int* boundaries() { return m_boundaries.data(); }

  int* groups() { return m_groups.data(); }

  private:
  const static int INTERNAL2EX_ORDER[4];
};

#endif // PUMGEN_SRC_INPUT_NETCDFPARTITION_H_
