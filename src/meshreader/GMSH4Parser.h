// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_
#define PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_

#include <cstddef>
#include <map>
#include <set>

#include "GMSHParser.h"

namespace puml {

/**
 * Parser for MSH 4.1 files.
 */
class GMSH4Parser : public GMSHParser {
  public:
  using GMSHParser::GMSHParser;

  private:
  void parse_() override;
  void parseEntities();
  void parseNodes();
  void parseElements();
  void parsePeriodic();

  // the values of the sections
  std::size_t readSize() { return expectSize(); }
  long readInt() { return expectInteger(); }
  double readDouble() { return expectNumber(); }

  /**
   * Reads a node tag and returns the index of the node.
   */
  std::size_t readNodeIndex();

  /**
   * The physical tag of the elements of an entity; 0 if the entity has no physical group.
   */
  long physicalTag(std::size_t dim, long entityTag);

  std::map<long, long> physicalSurfaceIds;
  std::map<long, long> physicalVolumeIds;
  // surface and volume entities with elements, but without a physical group
  std::set<long> unassignedSurfaces;
  std::set<long> unassignedVolumes;

  // the node tags form the contiguous range [firstNodeTag, firstNodeTag + numNodes)
  bool hasNodes = false;
  std::size_t firstNodeTag = 0;
  std::size_t numNodes = 0;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_
