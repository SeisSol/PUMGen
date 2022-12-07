#ifndef PUMGEN_GMSH4PARSER_H
#define PUMGEN_GMSH4PARSER_H

#include "third_party/GMSHBuilder.h"
#include <string>

namespace puml {

class GMSH4Parser {
  public:
  GMSH4Parser(tndm::GMSHMeshBuilder* builder) : builder(builder) {}
  bool parseFile(std::string const& fileName);
  std::string_view getErrorMessage() const { return errorMsg; }

  private:
  tndm::GMSHMeshBuilder* builder;
  std::string errorMsg;
};

} // namespace puml

#endif // PUMGEN_GMSH4PARSER_H
