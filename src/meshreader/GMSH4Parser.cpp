#include "GMSH4Parser.h"

#include <cassert>
#include <cstdio>

#include "utils/logger.h"

namespace puml {
bool GMSH4Parser::parse_() {
  errorMsg.clear();
  getNextToken();

  double version = parseMeshFormat();
  if (version < 4.0 || version >= 5.0) {
    char buf[128];
    sprintf(buf, "Unsupported MSH version %.1lf", version);
    return logError<bool>(buf);
  }

  bool hasEntities = false;
  bool hasNodes = false;
  bool hasElements = false;

  while (curTok != tndm::GMSHToken::eof) {
    switch (curTok) {
    case tndm::GMSHToken::entities:
      hasEntities = parseEntities();
      break;
    case tndm::GMSHToken::nodes:
      hasNodes = parseNodes();
      break;
    case tndm::GMSHToken::elements:
      hasElements = parseElements();
      break;
    default:
      getNextToken();
      break;
    }
  }

  return hasEntities && hasNodes && hasElements;
}

bool GMSH4Parser::parseEntities() {
  getNextToken();
  std::size_t numPoints = expectNonNegativeInt();
  getNextToken();
  std::size_t numCurves = expectNonNegativeInt();
  getNextToken();
  std::size_t numSurfaces = expectNonNegativeInt();
  getNextToken();
  std::size_t numVolumes = expectNonNegativeInt();

  // ignore points
  for (std::size_t i = 0; i < numPoints; ++i) {
    getNextToken();
    [[maybe_unused]] std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_1 = expectNumber();
    getNextToken();
    std::size_t numTags = expectNonNegativeInt();
    // ignore tags
    for (std::size_t tagIdx = 0; tagIdx < numTags; tagIdx++) {
      getNextToken();
    }
  }

  //  ignore curves
  for (std::size_t i = 0; i < numCurves; ++i) {
    getNextToken();
    [[maybe_unused]] std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] auto z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_2 = expectNumber();
    getNextToken();
    std::size_t numPhysicalTags = expectNonNegativeInt();
    getNextToken();
    // ignore tags
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; tagIdx++) {
      getNextToken();
    }
    std::size_t numBoundaryTags = expectNonNegativeInt();
    // ignore boundary
    for (std::size_t tagIdx = 0; tagIdx < numBoundaryTags; tagIdx++) {
      getNextToken();
    }
  }

  //  read surfaces
  for (std::size_t i = 0; i < numSurfaces; ++i) {
    getNextToken();
    std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_2 = expectNumber();
    getNextToken();
    std::size_t numPhysicalTags = expectNonNegativeInt();
    std::vector<std::size_t> tags;
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; ++tagIdx) {
      getNextToken();
      tags.push_back(expectNonNegativeInt());
    }
    if (!tags.empty()) {
      physicalSurfaceIds.insert_or_assign(id, tags.at(0));
    }
    getNextToken();
    std::size_t numBoundaryTags = expectNonNegativeInt();
    // ignore boundary
    for (std::size_t curveIdx = 0; curveIdx < numBoundaryTags; ++curveIdx) {
      getNextToken();
    }
  }

  //  read volumes
  for (std::size_t i = 0; i < numVolumes; ++i) {
    getNextToken();
    std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] double z_2 = expectNumber();
    getNextToken();
    std::size_t numPhysicalTags = expectNonNegativeInt();
    std::vector<std::size_t> tags;
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; ++tagIdx) {
      getNextToken();
      tags.push_back(expectNonNegativeInt());
    }
    physicalVolumeIds.insert_or_assign(id, tags.at(0));
    getNextToken();
    std::size_t numBoundaryTags = expectNonNegativeInt();
    // ignore boundary
    for (std::size_t curveIdx = 0; curveIdx < numBoundaryTags; ++curveIdx) {
      getNextToken();
    }
  }
  getNextToken();
  if (curTok != tndm::GMSHToken::end_entities) {
    return logErrorAnnotated<bool>("Expected $EndEntities");
  }
  getNextToken();
  return true;
}

bool GMSH4Parser::parseNodes() {
  getNextToken();
  std::size_t numBlocks = expectNonNegativeInt();
  getNextToken();
  std::size_t numVertices = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] std::size_t minNodeTag = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] std::size_t maxNodeTag = expectNonNegativeInt();
  builder->setNumVertices(numVertices);

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    [[maybe_unused]] std::size_t dim = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] std::size_t entityTag = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] std::size_t parametric = expectNonNegativeInt();
    getNextToken();
    std::size_t numVerticesInBlock = expectNonNegativeInt();
    std::vector<std::size_t> vertexIds;
    vertexIds.reserve(numVerticesInBlock);
    // first read vertex ids
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      getNextToken();
      vertexIds.push_back(expectNonNegativeInt() - 1);
    }
    // then read vertex data
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      std::array<double, 3> x{};
      for (std::size_t i = 0; i < 3; i++) {
        getNextToken();
        x[i] = expectNumber();
      }
      assert(vertexIds.at(vertexIdx) < numVertices);
      builder->setVertex(vertexIds.at(vertexIdx), x);
    }
  }
  getNextToken();
  if (curTok != tndm::GMSHToken::end_nodes) {
    return logErrorAnnotated<bool>("Expected $EndNodes");
  }
  getNextToken();
  return true;
}

bool GMSH4Parser::parseElements() {
  getNextToken();
  std::size_t numBlocks = expectNonNegativeInt();
  getNextToken();
  std::size_t numElements = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] std::size_t minElementTag = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] std::size_t maxElementTag = expectNonNegativeInt();
  builder->setNumElements(numElements);

  constexpr std::size_t MaxNodes = sizeof(NumNodes) / sizeof(std::size_t);
  std::array<long, MaxNodes> nodes{};
  builder->setNumElements(numElements);

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    [[maybe_unused]] std::size_t dim = expectNonNegativeInt();
    getNextToken();
    std::size_t entityTag = expectNonNegativeInt();
    getNextToken();
    std::size_t type = expectNonNegativeInt();
    getNextToken();
    std::size_t numElementsInBlock = expectNonNegativeInt();
    for (std::size_t elementIdx = 0; elementIdx < numElementsInBlock; ++elementIdx) {
      getNextToken();
      [[maybe_unused]] std::size_t id = expectNonNegativeInt();
      for (std::size_t nodeIdx = 0; nodeIdx < NumNodes[type - 1]; nodeIdx++) {
        getNextToken();
        nodes.at(nodeIdx) = expectNonNegativeInt() - 1;
      }
      std::size_t tag = (type == 2) ? physicalSurfaceIds[entityTag] : physicalVolumeIds[entityTag];
      builder->addElement(type, tag, nodes.data(), NumNodes[type - 1]);
    }
  }
  getNextToken();
  if (curTok != tndm::GMSHToken::end_elements) {
    return logErrorAnnotated<bool>("Expected $EndElements");
  }
  getNextToken();

  return true;
}

} // namespace puml
