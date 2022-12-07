#include "GMSH4Parser.h"

#include <cstdio>
#include <fstream>
#include <sstream>

#include "utils/logger.h"

namespace puml {
bool GMSH4Parser::parse_() {
  errorMsg.clear();
  getNextToken();

  double version = parseMeshFormat();
  if (version < 2.0 || version >= 3.0) {
    char buf[128];
    sprintf(buf, "Unsupported MSH version %.1lf", version);
    return logError<bool>(buf);
  }

  bool hasNodes = false;
  bool hasElements = false;

  while (curTok != GMSHToken::eof) {
    switch (curTok) {
    case GMSHToken::nodes:
      hasNodes = parseNodes();
      break;
    case GMSHToken::elements:
      hasElements = parseElements();
      break;
    default:
      getNextToken();
      break;
    }
  }

  return hasNodes && hasElements;
}

} // namespace puml
