/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2017 Technical University of Munich
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 */

#ifndef FIDAP_READER_H
#define FIDAP_READER_H

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "utils/stringutils.h"

#include "MeshReader.h"

struct FidapBoundaryFace {
  /** The vertices of the face */
  int vertices[3];
  /** The type of the boundary */
  int type;
};

/** Defines a face by three vertices */
struct FaceVertex {
  FaceVertex() {}

  FaceVertex(int v1, int v2, int v3) {
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    std::sort(vertices, vertices + 3);
  }

  FaceVertex(const int v[3]) {
    memcpy(vertices, v, sizeof(vertices));
    std::sort(vertices, vertices + 3);
  }

  bool operator<(const FaceVertex &other) const {
    return vertices[0] < other.vertices[0] ||
           (vertices[0] == other.vertices[0] &&
            vertices[1] < other.vertices[1]) ||
           (vertices[0] == other.vertices[0] &&
            vertices[1] == other.vertices[1] &&
            vertices[2] < other.vertices[2]);
  }

  int vertices[3];
};

/** Defines a face by element and face side */
struct FaceElement {
  FaceElement() {
    element = -1;
    side = -1;
  }

  FaceElement(int e, int s) {
    element = e;
    side = s;
  }

  int element;
  int side;
};

class FidapReader : public MeshReader {
private:
  struct ElementSection : FileSection {
    /** Group id for elements in this section */
    int group;
  };

  struct BoundarySection : FileSection {
    /** Type of the boundary */
    int type;
  };

  // Section descriptions
  FileSection m_vertices;
  std::vector<ElementSection> m_elements;
  std::vector<BoundarySection> m_boundaries;

public:
  void open(const char *meshFile) {
    MeshReader::open(meshFile);

    std::string line;

    // Read header information
    getline(m_mesh, line);
    utils::StringUtils::trim(line);
    if (line != FIDAP_FILE_ID)
      logError() << "Not an Fidap mesh file:" << meshFile;

    std::string name;
    getline(m_mesh, name); // Internal name

    getline(m_mesh, line); // Skip line:
                           // VERSION: x.y
    getline(m_mesh, line); // Date
    getline(m_mesh, line); // Empty line

    unsigned int nElements, nZones, t, dimensions;
    m_mesh >> m_vertices.nLines;
    m_mesh >> nElements;
    m_mesh >> nZones;
    m_mesh >> t; // I don't know what this variable contains
    m_mesh >> dimensions;
    if (dimensions != 3)
      logError() << "Fidap file does not contain a 3 dimensional mesh";
    getline(m_mesh, line); // Read newline

    getline(m_mesh, line); // Empty line
    getline(m_mesh, line); // Unknown line
    getline(m_mesh, line); // Empty line
    getline(m_mesh, line); // Unknown line
    getline(m_mesh, line); // Empty line
    getline(m_mesh, line); // Unknown line

    // Find seek positions
    // Vertices
    getline(m_mesh, line);
    if (line.find(NODAL_COORDINATES) == std::string::npos)
      logError() << "Invalid Fidap format: Coordinates expected, found" << line;
    m_vertices.seekPosition = m_mesh.tellg();
    getline(m_mesh, line);
    m_vertices.lineSize = line.length() + 1;
    m_mesh.seekg(m_vertices.seekPosition +
                 m_vertices.nLines * m_vertices.lineSize);

    // Boundary conditions (currently ignored)
    getline(m_mesh, line);
    if (line.find(BOUNDARY_CONDITIONS) == std::string::npos)
      logError() << "Invalid Fidap format: Boundary conditions expected, found"
                 << line;
    getline(m_mesh, line);

    // Element groups
    getline(m_mesh, line);
    if (line.find(ELEMENT_GROUPS) == std::string::npos)
      logError() << "Invalid Fidap format: Element groups expected, found"
                 << line;

    for (unsigned int i = 0; i < nZones; i++) {
      FileSection sec;

      // Group
      m_mesh >> name;
      if (name.find(ZONE_GROUP) == std::string::npos)
        logError() << "Invalid Fidap format: Expected group, found" << name;
      unsigned int id;
      m_mesh >> id;
      // Elements
      m_mesh >> name;
      if (name.find(ZONE_ELEMENTS) == std::string::npos)
        logError() << "Invalid Fidap format: Expected elements, found" << name;
      m_mesh >> sec.nLines;
      // Nodes
      m_mesh >> name;
      if (name.find(ZONE_NODES) == std::string::npos)
        logError() << "Invalid Fidap format: Expected nodes, found" << name;
      unsigned int nodeType;
      m_mesh >> nodeType;
      // Geometry
      m_mesh >> name;
      if (name.find(ZONE_GEOMETRY) == std::string::npos)
        logError() << "Invalid Fidap format: Expected geometry, found" << name;
      unsigned int geoType;
      m_mesh >> geoType;
      // Type
      m_mesh >> name;
      if (name.find(ZONE_TYPE) == std::string::npos)
        logError() << "Invalid Fidap format: Expected type, found" << name;
      unsigned int type;
      m_mesh >> type;
      getline(m_mesh, line); // Read newline

      // Entity name
      m_mesh >> name;
      m_mesh >> name;
      std::string entityName;
      getline(m_mesh, entityName);
      utils::StringUtils::trim(entityName);

      sec.seekPosition = m_mesh.tellg();
      getline(m_mesh, line);
      sec.lineSize = line.size() + 1;

      switch (nodeType) {
      case 3: {
        BoundarySection boundary;
        boundary.seekPosition = sec.seekPosition;
        boundary.nLines = sec.nLines;
        boundary.lineSize = sec.lineSize;

        boundary.type = convertType(entityName);

        m_boundaries.push_back(boundary);
        break;
      }
      case 4: {
        ElementSection element;
        element.seekPosition = sec.seekPosition;
        element.nLines = sec.nLines;
        element.lineSize = sec.lineSize;

        element.group = m_elements.size() + 1;

        m_elements.push_back(element);
        break;
      }
      default:
        logError() << "Unsupported node type" << nodeType << "found";
      }

      m_mesh.seekg(sec.seekPosition + sec.nLines * sec.lineSize);
    }
  }

  unsigned int nVertices() const { return m_vertices.nLines; }

  unsigned int nElements() const {
    unsigned int count = 0;
    for (std::vector<ElementSection>::const_iterator i = m_elements.begin();
         i != m_elements.end(); i++) {
      count += i->nLines;
    }

    return count;
  }

  unsigned int nBoundaries() const {
    unsigned int count = 0;
    for (std::vector<BoundarySection>::const_iterator i = m_boundaries.begin();
         i != m_boundaries.end(); i++) {
      count += i->nLines;
    }

    return count;
  }

  void readVertices(unsigned int start, unsigned int count, double *vertices) {
    m_mesh.clear();

    m_mesh.seekg(m_vertices.seekPosition + start * m_vertices.lineSize +
                 m_vertices.lineSize - 3 * COORDINATE_SIZE - 1);

    for (unsigned int i = 0; i < count; i++) {
      for (int j = 0; j < 3; j++)
        m_mesh >> vertices[i * 3 + j];

      m_mesh.seekg(m_vertices.lineSize - 3 * COORDINATE_SIZE,
                   std::fstream::cur);
    }
  }

  void readElements(unsigned int start, unsigned int count, int *elements) {
    m_mesh.clear();

    // Get the boundary, were we should start reading
    std::vector<ElementSection>::const_iterator section;
    for (section = m_elements.begin();
         section != m_elements.end() && section->nLines < start; section++)
      start -= section->nLines;

    m_mesh.seekg(section->seekPosition + start * section->lineSize +
                 section->lineSize - 4 * VERTEX_ID_SIZE - 1);

    for (unsigned int i = 0; i < count; i++) {
      if (start >= section->nLines) {
        // Are we starting a new section in this iteration?
        start = 0;
        section++;

        m_mesh.seekg(section->seekPosition + section->lineSize -
                     4 * VERTEX_ID_SIZE - 1);
      }

      for (unsigned int j = 0; j < 4; j++) {
        m_mesh >> elements[i * 4 + j];
        elements[i * 4 + j]--;
      }

      m_mesh.seekg(section->lineSize - 4 * VERTEX_ID_SIZE, std::fstream::cur);

      start++; // Line in the current section
    }
  }

  /**
   * Read group ids from start to start+count. The group ids are sorted
   * in the same order as the elements.
   */
  void readGroups(unsigned int start, unsigned int count, int *groups) {
    // Get the boundary, were we should start reading
    std::vector<ElementSection>::const_iterator section;
    for (section = m_elements.begin();
         section != m_elements.end() && section->nLines < start; section++)
      start -= section->nLines;

    for (unsigned int i = 0; i < count; i++) {
      if (start >= section->nLines) {
        // Are we starting a new section in this iteration?
        start = 0;
        section++;
      }

      groups[i] = section->group;

      start++; // "Line" in the current section
    }
  }

  /**
   * Read all groups
   */
  void readGroups(int *groups) { readGroups(0, nElements(), groups); }

  void readBoundaries(unsigned int start, unsigned int count,
                      FidapBoundaryFace *boundaries) {
    m_mesh.clear();

    // Get the boundary, were we should start reading
    std::vector<BoundarySection>::const_iterator section;
    for (section = m_boundaries.begin();
         section != m_boundaries.end() && section->nLines < start; section++)
      start -= section->nLines;

    m_mesh.seekg(section->seekPosition + start * section->lineSize +
                 section->lineSize - 3 * VERTEX_ID_SIZE - 1);

    for (unsigned int i = 0; i < count; i++) {
      if (start >= section->nLines) {
        // Are we starting a new section in this iteration?
        start = 0;
        section++;

        m_mesh.seekg(section->seekPosition + section->lineSize -
                     3 * VERTEX_ID_SIZE - 1);
      }

      for (unsigned int j = 0; j < 3; j++) {
        m_mesh >> boundaries[i].vertices[j];
        boundaries[i].vertices[j]--;
      }
      boundaries[i].type = section->type;

      // Seek to next position
      m_mesh.seekg(section->lineSize - 3 * VERTEX_ID_SIZE, std::fstream::cur);

      start++; // Line in the current section
    }
  }

  // TODO void readBoundaries(int* boundaries) {}

public:
  /**
   * Creates a map from vertex definition of a face to a element/side definition
   *
   * @param elements Elements (4 vertices per element)
   * @param n Number of elements
   * @param[out] map The map
   */
  static void createFaceMap(const int *elements, unsigned int n,
                            std::map<FaceVertex, FaceElement> &map) {
    for (unsigned int i = 0; i < n; i++) {
      FaceVertex v1(elements[i * 4], elements[i * 4 + 1], elements[i * 4 + 2]);
      FaceElement e1(i, 0);
      map[v1] = e1;

      FaceVertex v2(elements[i * 4], elements[i * 4 + 1], elements[i * 4 + 3]);
      FaceElement e2(i, 1);
      map[v2] = e2;

      FaceVertex v3(elements[i * 4 + 1], elements[i * 4 + 2],
                    elements[i * 4 + 3]);
      FaceElement e3(i, 2);
      map[v3] = e3;

      FaceVertex v4(elements[i * 4], elements[i * 4 + 2], elements[i * 4 + 3]);
      FaceElement e4(i, 3);
      map[v4] = e4;
    }
  }

private:
  static int convertType(const std::string &name) {
    std::istringstream ss(name.substr(2));

    int type;
    ss >> type;

    return type % 100;
  }

  /** Number of characters required to store a coordinate */
  static const size_t COORDINATE_SIZE = 20ul;
  /** Number of characters required to store a vertex id */
  static const size_t VERTEX_ID_SIZE = 8ul;

  static const char *FIDAP_FILE_ID;
  static const char *NODAL_COORDINATES;
  static const char *BOUNDARY_CONDITIONS;
  static const char *ELEMENT_GROUPS;

  static const char *ZONE_GROUP;
  static const char *ZONE_ELEMENTS;
  static const char *ZONE_NODES;
  static const char *ZONE_GEOMETRY;
  static const char *ZONE_TYPE;
};

#endif // FIDAP_READER_H
