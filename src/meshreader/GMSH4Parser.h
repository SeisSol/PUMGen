#ifndef PUMGEN_GMSH4PARSER_H
#define PUMGEN_GMSH4PARSER_H

#include <string>
#include <string_view>

#include "GMSH4Lexer.h"
#include "meshreader/GMSHBuilder.h"
#include "meshreader/GMSHParser.h"

namespace puml {

class GMSH4Parser : public GMSHParser {
  public:
  explicit GMSH4Parser(puml::GMSHMeshBuilder* builder) : GMSHParser(builder) {
    lexer = new GMSH4Lexer();
  }

  private:
  double parseMeshFormat() { return 1; };
  bool parseNodes() { return 1; };
  bool parseElements() { return 1; };
  virtual bool parse_() override;
};

} // namespace puml

#endif // PUMGEN_GMSH4PARSER_H
