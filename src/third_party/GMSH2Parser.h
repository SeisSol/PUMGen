// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_THIRD_PARTY_GMSH2PARSER_H_
#define PUMGEN_SRC_THIRD_PARTY_GMSH2PARSER_H_

#include "GMSHLexer.h"
#include "GMSHMeshBuilder.h"
#include "GMSHParser.h"

#include <array>
#include <optional>
#include <streambuf>
#include <string>
#include <string_view>

namespace tndm {

class GMSH2Parser : public GMSHParser {
  public:
  using GMSHParser::GMSHParser;

  private:
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

  // the node tags are 1, ..., numNodes
  bool hasNodes = false;
  std::size_t numNodes = 0;

  /**
   * Reads a node tag and returns the index of the node.
   */
  std::size_t expectNodeIndex() {
    if (curTok != GMSHToken::integer || lexer.getInteger() < 1 ||
        static_cast<std::size_t>(lexer.getInteger()) > numNodes) {
      return logErrorAnnotated<std::size_t>("Unknown node tag");
    }
    return static_cast<std::size_t>(lexer.getInteger()) - 1;
  }

  bool parseNodes();
  bool parseElements();
  bool parsePeriodic();
  virtual bool parse_() override;
};

}; // namespace tndm

#endif // PUMGEN_SRC_THIRD_PARTY_GMSH2PARSER_H_
