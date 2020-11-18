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

#ifndef GAMBIT_READER_H
#define GAMBIT_READER_H

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

#include "utils/logger.h"
#include "utils/stringutils.h"

#include "MeshReader.h"

namespace puml {

/**
 * Describes the group of an element
 */
struct ElementGroup {
  /** The element */
  unsigned int element;
  /** The group for this element */
  int group;
};

/**
 * Describes a boundary face
 */
struct GambitBoundaryFace {
  /** The element the face belongs to */
  unsigned int element;
  /** The face of the element */
  unsigned int face;
  /** The type of the boundary */
  int type;
};

class GambitReader : public MeshReader {
  private:
  struct ElementSection : FileSection {
    /** Start of vertices in an element */
    size_t vertexStart;
    /** Size per vertex id */
    size_t vertexSize;
  };

  struct GroupSection : FileSection {
    /** Size per element id */
    size_t elementSize;
    /** Group number */
    int group;
  };

  struct BoundarySection : FileSection {
    /** Type of the boundary */
    int type;
    /** Size per element id (only fixed line length) */
    size_t elementSize;
    /** Size per element type (only fixed line length) */
    size_t elementTypeSize;
    /** Size per face id (only fixed line length) */
    size_t faceSize;
    /** Is line length variable? */
    bool variableLineLength;
  };

  // Section descriptions
  FileSection m_vertices;
  ElementSection m_elements;
  std::vector<GroupSection> m_groups;
  std::vector<BoundarySection> m_boundaries;

  public:
  void open(const char* meshFile) {
    MeshReader::open(meshFile);

    std::string line;

    // Read header information
    getline(m_mesh, line); // First line contains version
                           // we ignore this line for know
    line.clear();
    getline(m_mesh, line);
    utils::StringUtils::trim(line);
    if (line != GAMBIT_FILE_ID)
      logError() << "Not a Gambit mesh file:" << meshFile;

    std::string name;
    getline(m_mesh, name); // Internal name

    getline(m_mesh, line); // Skip line:
                           // PROGRAM: Gambit VERSION: x.y.z

    getline(m_mesh, line); // Date

    getline(m_mesh, line); // Skip problem size names

    unsigned int nGroups, nBoundaries, dimensions;
    m_mesh >> m_vertices.nLines;
    m_mesh >> m_elements.nLines;
    m_mesh >> nGroups;
    m_mesh >> nBoundaries;
    m_mesh >> dimensions;
    getline(m_mesh, line); // Skip rest of the line
    if (dimensions != 3)
      logError() << "Gambit file does not contain a 3 dimensional mesh";

    getline(m_mesh, line);
    utils::StringUtils::trim(line);
    if (line != ENDSECTION)
      logError() << "Invalid Gambit format: End of header expected, found" << line;

    // Find seek positions
    // Vertices
    getline(m_mesh, line);
    if (line.find(NODAL_COORDINATES) == std::string::npos)
      logError() << "Invalid Gambit format: Coordinates expected, found" << line;
    m_vertices.seekPosition = m_mesh.tellg();
    getline(m_mesh, line);
    m_vertices.lineSize = line.length() + 1;
    m_mesh.seekg(m_vertices.seekPosition + m_vertices.nLines * m_vertices.lineSize);
    getline(m_mesh, line);
    utils::StringUtils::rtrim(line); // remove \r
    if (line != ENDSECTION)
      logError() << "Invalid Gambit format: End of coordinates expected, found" << line;

    // Elements
    getline(m_mesh, line);
    if (line.find(ELEMENT_CELLS) == std::string::npos)
      logError() << "Invalid Gambit format: Elements expected, found" << line;
    m_elements.seekPosition = m_mesh.tellg();
    getline(m_mesh, line);
    m_elements.lineSize = line.length() + 1;
    m_mesh.seekg(m_elements.seekPosition + m_elements.nLines * m_elements.lineSize);

    std::istringstream ss(line);
    int tmp;
    ss >> tmp; // Element id
    ss >> tmp;
    ss >> tmp;
    ss.seekg(1, std::stringstream::cur); // White space
    m_elements.vertexStart = ss.tellg();
    // TODO works only for tetrahedrals
    utils::StringUtils::rtrim(line);
    if ((line.length() - m_elements.vertexStart) % 4 != 0)
      logError() << "Invalid Gambit format: Mesh does not seem to contain "
                    "tetrahedrals";
    m_elements.vertexSize = (line.length() - m_elements.vertexStart) / 4;

    getline(m_mesh, line);
    utils::StringUtils::rtrim(line); // remove \r
    if (line != ENDSECTION)
      logError() << "Invalid Gambit format: End of elements expected, found" << line;

    // Groups
    m_groups.resize(nGroups);
    for (unsigned int i = 0; i < nGroups; i++) {
      getline(m_mesh, line);
      if (line.find(ELEMENT_GROUP) == std::string::npos)
        logError() << "Invalid Gambit format: Group expected, found" << line;
      getline(m_mesh, line);

      // Group size and material
      std::string y;
      std::istringstream ss(line);
      ss >> y;
      ss >> m_groups[i].group;
      ss >> y;
      ss >> m_groups[i].nLines; // This is not the actual number of lines
                                // Because Gambit stores more than one element
                                // per line.

      getline(m_mesh, line);
      getline(m_mesh, line);
      m_groups[i].seekPosition = m_mesh.tellg();

      getline(m_mesh, line);
      m_groups[i].lineSize = line.size() + 1;

      utils::StringUtils::rtrim(line);
      if (line.length() % ELEMENTS_PER_LINE_GROUP != 0)
        logError() << "Invalid Gambit format: Mesh does not contain" << ELEMENTS_PER_LINE_GROUP
                   << "in one group line";
      m_groups[i].elementSize = line.length() / ELEMENTS_PER_LINE_GROUP;

      m_mesh.seekg(m_groups[i].seekPosition +
                   (m_groups[i].nLines / ELEMENTS_PER_LINE_GROUP) * m_groups[i].lineSize);

      if (m_groups[i].nLines % ELEMENTS_PER_LINE_GROUP != 0)
        // Last line
        m_mesh.seekg(m_groups[i].lineSize -
                         (ELEMENTS_PER_LINE_GROUP - m_groups[i].nLines % ELEMENTS_PER_LINE_GROUP) *
                             m_groups[i].elementSize,
                     std::ifstream::cur);

      getline(m_mesh, line);
      utils::StringUtils::rtrim(line); // remove \r
      if (line != ENDSECTION)
        logError() << "Invalid Gambit format: End of group expected, found" << line;
    }

    // Boundaries
    m_boundaries.resize(nBoundaries);
    for (unsigned int i = 0; i < nBoundaries; i++) {
      getline(m_mesh, line);
      if (line.find(BOUNDARY_CONDITIONS) == std::string::npos)
        logError() << "Invalid Gambit format: Boundaries expected, found" << line;
      m_boundaries[i].variableLineLength = false;

      getline(m_mesh, line);
      m_boundaries[i].seekPosition = m_mesh.tellg();

      // Boundary type and size
      int x;
      std::istringstream ss(line);
      ss >> m_boundaries[i].type;
      m_boundaries[i].type %= 100; // Fix boundary type TODO not sure where this should be placed
      ss >> x;
      ss >> m_boundaries[i].nLines;

      // Try boundary with fixed line length
      getline(m_mesh, line);
      m_boundaries[i].lineSize = line.size() + 1;

      // Get element size
      std::istringstream ssE(line);
      unsigned int element, type, face;
      ssE >> element;
      m_boundaries[i].elementSize = ssE.tellg();
      ssE >> type;
      m_boundaries[i].elementTypeSize =
          static_cast<size_t>(ssE.tellg()) - m_boundaries[i].elementSize;
      ssE >> face;
      if (ssE.tellg() == -1)
        m_boundaries[i].faceSize =
            line.length() - m_boundaries[i].elementSize - m_boundaries[i].elementTypeSize;
      else
        m_boundaries[i].faceSize = static_cast<size_t>(ssE.tellg()) - m_boundaries[i].elementSize -
                                   m_boundaries[i].elementTypeSize;

      m_mesh.seekg(m_boundaries[i].seekPosition +
                   m_boundaries[i].nLines * m_boundaries[i].lineSize);

      getline(m_mesh, line);
      utils::StringUtils::rtrim(line); // remove \r
      if (line != ENDSECTION) {
        logWarning() << "Gambit format does not seem to have a fixed boundary "
                        "line length. Trying with variable line length";

        m_boundaries[i].variableLineLength = true;
        m_boundaries[i].lineSize = 0; // Variable line size
        m_boundaries[i].elementSize = 0;

        m_mesh.clear(); // The fixed boundary try may seek beyond the end of the
                        // file
        m_mesh.seekg(m_boundaries[i].seekPosition);

        for (size_t j = 0; j < m_boundaries[i].nLines; j++)
          m_mesh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        getline(m_mesh, line);
        utils::StringUtils::rtrim(line); // remove \r
        if (line != ENDSECTION)
          logError() << "Invalid Gambit format: End of boundaries expected, found" << line;
      }
    }
  }

  unsigned int nVertices() const { return m_vertices.nLines; }

  unsigned int nElements() const { return m_elements.nLines; }

  unsigned int nBoundaries() const {
    unsigned int count = 0;
    for (std::vector<BoundarySection>::const_iterator i = m_boundaries.begin();
         i != m_boundaries.end(); i++) {
      count += i->nLines;
    }

    return count;
  }

  /**
   * @copydoc MeshReader::readVertices(unsigned int, unsigned int, double*)
   *
   * @todo Only 3 dimensional meshes are supported
   */
  void readVertices(unsigned int start, unsigned int count, double* vertices) {
    m_mesh.clear();

    m_mesh.seekg(m_vertices.seekPosition + start * m_vertices.lineSize + m_vertices.lineSize -
                 3 * COORDINATE_SIZE - 1);

    char* buf = new char[3 * COORDINATE_SIZE];

    for (unsigned int i = 0; i < count; i++) {
      m_mesh.read(buf, 3 * COORDINATE_SIZE);

      for (int j = 0; j < 3; j++) {
        std::istringstream ss(std::string(&buf[j * COORDINATE_SIZE], COORDINATE_SIZE));
        ss >> vertices[i * 3 + j];
      }

      m_mesh.seekg(m_vertices.lineSize - 3 * COORDINATE_SIZE, std::fstream::cur);
    }

    delete[] buf;
  }

  /**
   * @copydoc MeshReader::readElements(unsigned int, unsigned int, int*)
   *
   * @todo Only tetrahedral meshes are supported
   * @todo Support for varying coordinate/vertexid fields
   */
  void readElements(unsigned int start, unsigned int count, int* elements) {
    m_mesh.clear();

    m_mesh.seekg(m_elements.seekPosition + start * m_elements.lineSize + m_elements.vertexStart);

    char* buf = new char[4 * m_elements.vertexSize];

    for (unsigned int i = 0; i < count; i++) {
      m_mesh.read(buf, 4 * m_elements.vertexSize);

      for (int j = 0; j < 4; j++) {
        std::istringstream ss(std::string(&buf[j * m_elements.vertexSize], m_elements.vertexSize));
        ss >> elements[i * 4 + j];
        elements[i * 4 + j]--;
      }

      m_mesh.seekg(m_elements.lineSize - 4 * m_elements.vertexSize, std::fstream::cur);
    }

    delete[] buf;
  }

  /**
   * Reads group number for elements from start to start+count
   *
   * @param groups Buffer for storing the group numbers. The caller is
   * responsible for allocating the buffer.
   */
  void readGroups(unsigned int start, unsigned int count, ElementGroup* groups) {
    m_mesh.clear();

    // Get the group, were we should start reading
    std::vector<GroupSection>::const_iterator section;
    for (section = m_groups.begin(); section != m_groups.end() && section->nLines < start;
         section++)
      start -= section->nLines;

    m_mesh.seekg(section->seekPosition + (start / ELEMENTS_PER_LINE_GROUP) * section->lineSize +
                 (start % ELEMENTS_PER_LINE_GROUP) * section->elementSize);

    char* buf = new char[section->elementSize];

    for (unsigned int i = 0; i < count; i++) {
      m_mesh.read(buf, section->elementSize);

      std::istringstream ss(std::string(buf, section->elementSize));
      ss >> groups[i].element;
      groups[i].element--;
      groups[i].group = section->group;

      start++;
      // May need to jump from one section to the next
      if (start >= section->nLines) {
        start = 0;
        section++;

        m_mesh.seekg(section->seekPosition);
      } else if (start % ELEMENTS_PER_LINE_GROUP == 0)
        // Skip newline char at end of line
        m_mesh.seekg(section->lineSize - section->elementSize * ELEMENTS_PER_LINE_GROUP,
                     std::fstream::cur);
    }

    delete[] buf;
  }

  /**
   * Reads all groups numbers.
   * In contrast to readGroups(unsigned int, unsigned int, ElementGroup*) it
   * returns the group numbers sorted according to the elements.
   */
  void readGroups(int* groups) {
    logInfo() << "Reading group information";

    ElementGroup* map = new ElementGroup[nElements()];
    readGroups(0, nElements(), map);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned int i = 0; i < nElements(); i++)
      groups[map[i].element] = map[i].group;

    delete[] map;
  }

  /**
   * Reads boundaries from start to start+count from the file and stores them in
   * <code>boundaries</code>.
   *
   * @param elements Buffer to store boundaries. Each boundary consists of the
   * element, the face and the boundary type.
   */
  void readBoundaries(unsigned int start, unsigned int count, GambitBoundaryFace* boundaries) {
    m_mesh.clear();

    // Get the boundary, were we should start reading
    std::vector<BoundarySection>::const_iterator section;
    for (section = m_boundaries.begin(); section != m_boundaries.end() && section->nLines < start;
         section++)
      start -= section->nLines;

    char* buf;

    if (section->variableLineLength) {
      m_mesh.seekg(section->seekPosition);
      for (size_t i = 0; i < start; i++)
        m_mesh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      buf = 0L;
    } else {
      m_mesh.seekg(section->seekPosition + start * section->lineSize);

      buf = new char[section->elementSize + section->elementTypeSize + section->faceSize];
    }

    for (unsigned int i = 0; i < count; i++) {
      if (start >= section->nLines) {
        // Are we starting a new section in this iteration?
        start = 0;
        section++;

        m_mesh.seekg(section->seekPosition);

        delete[] buf;
        if (section->variableLineLength)
          buf = 0L;
        else
          buf = new char[section->elementSize + section->elementTypeSize + section->faceSize];
      }

      if (section->variableLineLength) {
        unsigned int elementType; // Ignored
        m_mesh >> boundaries[i].element;
        m_mesh >> elementType;
        m_mesh >> boundaries[i].face;
      } else {
        m_mesh.read(buf, section->elementSize + section->elementTypeSize + section->faceSize);

        std::istringstream ssE(std::string(buf, section->elementSize));
        ssE >> boundaries[i].element;

        std::istringstream ssF(
            std::string(&buf[section->elementSize + section->elementTypeSize], section->faceSize));
        ssF >> boundaries[i].face;

        // Seek to next position
        m_mesh.seekg(section->lineSize -
                         (section->elementSize + section->elementTypeSize + section->faceSize),
                     std::fstream::cur);
      }

      boundaries[i].element--;
      boundaries[i].face--;
      boundaries[i].type = section->type;

      start++; // Line in the current section
    }

    delete[] buf;
  }

  /**
   * Reads all boundaries.
   *
   * @param Boundary condition for all faces. The caller is responsible
   *  for allocation the memory (<code>nElements()*4</code>). Only the faces
   *  for which boundary conditions are available are modified.
   *
   * @todo Only tetrahedral meshes are supported
   */
  void readBoundaries(int* boundaries) {
    logInfo() << "Reading boundary conditions";

    unsigned int nBnds = nBoundaries();
    GambitBoundaryFace* faces = new GambitBoundaryFace[nBoundaries()];
    readBoundaries(0, nBnds, faces);

    for (unsigned int i = 0; i < nBnds; i++)
      boundaries[faces[i].element * 4 + faces[i].face] = faces[i].type;

    delete[] faces;
  }

  private:
  /** Number of character required to store a coordinate */
  static const size_t COORDINATE_SIZE = 20ul;
  /** Number of elements stored in one group line */
  static const size_t ELEMENTS_PER_LINE_GROUP = 10ul;

  static const char* GAMBIT_FILE_ID;
  static const char* ENDSECTION;
  static const char* NODAL_COORDINATES;
  static const char* ELEMENT_CELLS;
  static const char* ELEMENT_GROUP;
  static const char* BOUNDARY_CONDITIONS;
};

} // namespace puml

#endif // GAMBIT_READER_H
