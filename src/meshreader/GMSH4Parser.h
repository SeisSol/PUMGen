// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_
#define PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_

#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <vector>

#include "GMSHParser.h"

namespace puml {

/**
 * Parser for MSH 4.1 files, ASCII or binary.
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

  /**
   * Starts the data of a section; in binary files, the data begins after the line break.
   */
  void beginSectionData();

  // the values of the sections: text, or binary with sizeof(size_t) == dataSize and 4-byte ints
  std::size_t readSize();
  long readInt();
  double readDouble();
  // bulk reading of binary sizes and doubles
  void readSizes(std::size_t* values, std::size_t count);
  void readDoubles(double* values, std::size_t count);
  void readBinary(void* data, std::size_t bytes);

  /**
   * The offset of the next value, e.g. for error messages.
   */
  std::size_t position();

  /**
   * The index of the node with the given tag; the tag started at the given file offset.
   */
  std::size_t nodeIndex(std::size_t tag, std::size_t offset);

  bool binary = false;
  bool swapBytes = false;
  std::size_t dataSize = 8;
  std::vector<std::uint32_t> narrowSizes;

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
