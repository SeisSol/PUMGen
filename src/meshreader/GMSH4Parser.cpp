#include "GMSH4Parser.h"

#include <cstdio>
#include <fstream>
#include <sstream>

#include "utils/logger.h"

namespace puml {

template <typename T> T GMSH4Parser::logErrorAnnotated(std::string_view msg) {
  std::stringstream ss;
  ss << "GMSH parser error in line " << curLoc.line << " in column " << curLoc.col << ":\n";
  ss << '\t' << msg << '\n';
  errorMsg += ss.str();
  return {};
}

template <typename T> T GMSH4Parser::logError(std::string_view msg) {
  errorMsg += "GMSH parser error:\n\t";
  errorMsg += msg;
  errorMsg += '\n';
  return {};
}

bool GMSH4Parser::parseFile(std::string const& fileName) {
  std::ifstream in(fileName);
  if (!in.is_open()) {
    return logError<bool>("Unable to open MSH file");
  }
  lexer.setIStream(&in);
  return parse_();
}
} // namespace puml
