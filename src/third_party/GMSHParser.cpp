// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSHParser.h"

namespace tndm {

double GMSHParser::parseMeshFormat() {
  if (curTok != GMSHToken::mesh_format) {
    return logErrorAnnotated<double>("Expected $MeshFormat");
  }
  getNextToken();
  auto version = getNumber();
  if (!version) {
    return logErrorAnnotated<double>("Expected version number");
  }
  getNextToken();
  if (curTok == GMSHToken::integer && lexer.getInteger() == 1) {
    return logErrorAnnotated<double>(
        "Binary MSH files are not supported; write the mesh as ASCII (Mesh.Binary = 0)");
  }
  if (curTok != GMSHToken::integer || lexer.getInteger() != 0) {
    return logErrorAnnotated<double>("Expected file type 0 (ASCII)");
  }
  getNextToken(); // skip data-size
  getNextToken();
  if (curTok != GMSHToken::end_mesh_format) {
    return logErrorAnnotated<double>("Expected $EndMeshFormat");
  }
  getNextToken();
  return *version;
}
} // namespace tndm
