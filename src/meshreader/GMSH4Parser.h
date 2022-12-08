#ifndef PUMGEN_GMSH4PARSER_H
#define PUMGEN_GMSH4PARSER_H

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
  virtual bool parse_() override;
};

} // namespace puml

#endif // PUMGEN_GMSH4PARSER_H
