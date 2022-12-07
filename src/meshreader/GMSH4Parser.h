#ifndef PUMGEN_GMSH4PARSER_H
#define PUMGEN_GMSH4PARSER_H

#include <string>
#include <string_view>
#include <map>

#include "third_party/GMSH2Lexer.h"
#include "meshreader/GMSHBuilder.h"
#include "meshreader/GMSHParser.h"

namespace puml {

class GMSH4Parser : public GMSHParser {
  public:
  explicit GMSH4Parser(puml::GMSHMeshBuilder* builder) : GMSHParser(builder) {
    lexer = new tndm::GMSH2Lexer();
  }

  private:
  std::map<long, long> physicalSurfaceIds;
  std::map<long, long> physicalVolumeIds;
  long expectInt() {
    if (curTok != puml::GMSHToken::integer || lexer->getInteger() < 0) {
      return logErrorAnnotated<bool>("Expected non-zero integer");
    }
    return lexer->getInteger();
  };
  double expectNumber() {
    auto num = getNumber();
    if (!num) {
      return logErrorAnnotated<bool>("Expected number");
    }
    return num.value();
  };
  bool parseEntities();
  bool parseNodes();
  bool parseElements();
  virtual bool parse_() override;
};

} // namespace puml

#endif // PUMGEN_GMSH4PARSER_H
