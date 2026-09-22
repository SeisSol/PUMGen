// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_
#define PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_

#include <cstdio>
#include <map>
#include <string>
#include <string_view>

#include "meshreader/GMSHBuilder.h"
#include "third_party/GMSHLexer.h"
#include "third_party/GMSHParser.h"

namespace puml {

class GMSH4Parser : public tndm::GMSHParser {
  public:
  using tndm::GMSHParser::GMSHParser;

  private:
  std::map<unsigned long, long> physicalSurfaceIds;
  std::map<unsigned long, long> physicalVolumeIds;

  // the node tags form the contiguous range [firstNodeTag, firstNodeTag + numNodes)
  bool hasNodes = false;
  std::size_t firstNodeTag = 0;
  std::size_t numNodes = 0;

  /**
   * Reads a node tag and returns the index of the node.
   */
  std::size_t expectNodeIndex() {
    const std::size_t tag = expectNonNegativeInt();
    if (tag < firstNodeTag || tag - firstNodeTag >= numNodes) {
      char buf[128];
      snprintf(buf, sizeof(buf), "Unknown node tag %zu", tag);
      return logErrorAnnotated<std::size_t>(buf);
    }
    return tag - firstNodeTag;
  }

  unsigned long expectNonNegativeInt() {
    if (curTok != tndm::GMSHToken::integer || lexer.getInteger() < 0) {
      return logErrorAnnotated<bool>("Expected non-negative integer");
    }
    return static_cast<unsigned long>(lexer.getInteger());
  }

  double expectNumber() {
    auto num = getNumber();
    if (!num) {
      return logErrorAnnotated<bool>("Expected number");
    }
    return num.value();
  }

  bool parseEntities();
  bool parseNodes();
  bool parseElements();
  bool parsePeriodic();
  virtual bool parse_() override;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSH4PARSER_H_
