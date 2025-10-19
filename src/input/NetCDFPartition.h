// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_NETCDFPARTITION_H_
#define PUMGEN_SRC_INPUT_NETCDFPARTITION_H_

#include <cstring>

/**
 * Describes one partition (required for reading netCDF meshes)
 */
class Partition {
  private:
  std::size_t m_nElements;
  std::size_t m_nVertices;

  int* m_elements;
  double* m_vertices;
  int* m_boundaries;
  int* m_groups;

  public:
  Partition()
      : m_nElements(0), m_nVertices(0), m_elements(0L), m_vertices(0L), m_boundaries(0L),
        m_groups(0L) {}

  ~Partition() {
    delete[] m_elements;
    delete[] m_vertices;
    delete[] m_boundaries;
    delete[] m_groups;
  }

  void setElemSize(std::size_t nElements) {
    if (m_nElements != 0)
      return;

    m_nElements = nElements;

    m_elements = new int[nElements * 4];
    m_boundaries = new int[nElements * 4];
    m_groups = new int[nElements];
    // Set default value for groups
    memset(m_groups, 0, nElements * sizeof(int));
  }

  void setVrtxSize(std::size_t nVertices) {
    if (m_nVertices != 0)
      return;

    m_nVertices = nVertices;

    m_vertices = new double[nVertices * 3];
  }

  void convertBoundary() {
    int ncBoundaries[4];

    for (unsigned int i = 0; i < m_nElements * 4; i += 4) {
      memcpy(ncBoundaries, &m_boundaries[i], 4 * sizeof(int));
      for (unsigned int j = 0; j < 4; j++)
        m_boundaries[i + j] = ncBoundaries[INTERNAL2EX_ORDER[j]];
    }
  }

  std::size_t nElements() const { return m_nElements; }

  std::size_t nVertices() const { return m_nVertices; }

  int* elements() { return m_elements; }

  double* vertices() { return m_vertices; }

  int* boundaries() { return m_boundaries; }

  int* groups() { return m_groups; }

  private:
  const static int INTERNAL2EX_ORDER[4];
};

#endif // PUMGEN_SRC_INPUT_NETCDFPARTITION_H_
