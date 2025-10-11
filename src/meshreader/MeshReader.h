// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef MESH_READER_H
#define MESH_READER_H

#include <fstream>

#include "utils/logger.h"

class MeshReader {
  protected:
  struct FileSection {
    /** Start of the section */
    size_t seekPosition;
    /** Number of elements (lines) in the section */
    std::size_t nLines;
    /** Line size */
    size_t lineSize;
  };

  std::ifstream m_mesh;

  public:
  MeshReader() {}

  MeshReader(const char* meshFile) { open(meshFile); }

  virtual ~MeshReader() {}

  virtual void open(const char* meshFile) {
    m_mesh.open(meshFile);

    // Files are readable?
    if (!m_mesh)
      logError() << "Could not open mesh file" << meshFile;
  }

  virtual std::size_t nVertices() const = 0;
  virtual std::size_t nElements() const = 0;
  /**
   * @return Number of boundary faces
   */
  virtual std::size_t nBoundaries() const = 0;

  /**
   * Reads vertices from start tp start+count from the file and stores them in
   * vertices
   *
   * @param vertices Buffer to store coordinates of each vertex. The caller is
   * responsible for allocating the buffer. The size of the buffer must be
   * count*dimensions
   */
  virtual void readVertices(std::size_t start, std::size_t count, double* vertices) = 0;

  /**
   * Read all vertices
   *
   * @see readVertices(size_t, size_t, double*)
   */
  void readVertices(double* vertices) {
    logInfo() << "Reading vertices";
    readVertices(0, nVertices(), vertices);
  }

  /**
   * Reads elements from start to start+count from the file and stores them in
   * elements
   *
   * @param elements Buffer to store vertices of each element. The caller is
   * responsible for allocating the buffer. The Size of the buffer must be
   * count*vertices_per_element.
   */
  virtual void readElements(std::size_t start, std::size_t count, std::size_t* elements) = 0;

  /**
   * Reads all elements
   *
   * @see readElements(size_t, size_t, size_t*)
   */
  void readElements(std::size_t* elements) {
    logInfo() << "Reading elements";
    readElements(0, nElements(), elements);
  }
};

#endif // MESH_READER_H
