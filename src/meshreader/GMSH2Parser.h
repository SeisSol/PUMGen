// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSH2PARSER_H_
#define PUMGEN_SRC_MESHREADER_GMSH2PARSER_H_

#include <cstddef>

#include "GMSHParser.h"

namespace puml {

/**
 * Parser for MSH 2 files (ASCII).
 */
class GMSH2Parser : public GMSHParser {
  public:
  using GMSHParser::GMSHParser;

  private:
  void parse_() override;
  void parseNodes();
  void parseElements();
  void parsePeriodic();

  /**
   * Reads a node tag and returns the index of the node.
   */
  std::size_t readNodeIndex();

  // the node tags are 1, ..., numNodes
  bool hasNodes = false;
  std::size_t numNodes = 0;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSH2PARSER_H_
