#include "GMSHParser.h"

namespace puml {
double GMSHParser::parseMeshFormat() {
  if (curTok != puml::GMSHToken::mesh_format) {
    return logErrorAnnotated<double>("Expected $MeshFormat");
  }
  getNextToken();
  auto version = getNumber();
  if (!version) {
    return logErrorAnnotated<double>("Expected version number");
  }
  getNextToken();
  if (curTok != puml::GMSHToken::integer || lexer->getInteger() != 0) {
    return logErrorAnnotated<double>("Expected 0");
  }
  getNextToken(); // skip data-size
  getNextToken();
  if (curTok != puml::GMSHToken::end_mesh_format) {
    return logErrorAnnotated<double>("Expected $EndMeshFormat");
  }
  getNextToken();
  return *version;
}
};