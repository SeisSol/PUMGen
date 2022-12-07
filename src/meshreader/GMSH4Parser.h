#ifndef PUMGEN_GMSH4PARSER_H
#define PUMGEN_GMSH4PARSER_H

#include <string>
#include <string_view>

#include "GMSH4Lexer.h"
#include "third_party/GMSHBuilder.h"

namespace puml {

class GMSH4Parser {
  public:
  GMSH4Parser(tndm::GMSHMeshBuilder* builder) : builder(builder) {}
  bool parseFile(std::string const& fileName);
  std::string_view getErrorMessage() const { return errorMsg; }

  private:
  template <typename T> T logError(std::string_view msg);
  template <typename T> T logErrorAnnotated(std::string_view msg);
  bool parse_() { return true; };
  GMSH4Lexer lexer;
  GMSHSourceLocation curLoc;
  std::string errorMsg;
  tndm::GMSHMeshBuilder* builder;
};

} // namespace puml

#endif // PUMGEN_GMSH4PARSER_H
