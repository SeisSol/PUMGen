#ifndef GMSH2PARSER_20200901_H
#define GMSH2PARSER_20200901_H

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

  bool parseNodes();
  bool parseElements();
  bool parsePeriodic();
  virtual bool parse_() override;
};

}; // namespace tndm

#endif // GMSH2PARSER_20200901_H
