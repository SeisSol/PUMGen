#include "GMSH4Parser.h"

#include "utils/logger.h"

namespace puml {
bool GMSH4Parser::parseFile(std::string const& fileName) {
  logInfo() << "Read file" << fileName;
  return true;
}
} // namespace puml
