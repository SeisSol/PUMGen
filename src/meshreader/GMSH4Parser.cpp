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

  while (curTok != GMSHToken::eof) {
    switch (curTok) {
    case GMSHToken::entities:
      hasEntities = parseEntities();
      break;
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

  return hasEntities && hasNodes && hasElements;
}

bool GMSH4Parser::parseEntities() {
  //logInfo() << "Read entities, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  std::size_t numPoints = expectInt();
  getNextToken();
  std::size_t numCurves = expectInt();
  getNextToken();
  std::size_t numSurfaces = expectInt();
  getNextToken();
  std::size_t numVolumes = expectInt();
  //logInfo() << numPoints << "points," << numCurves << "curves," << numSurfaces << "surfaces," << numVolumes << "volumes.";

  // ignore points
  for (std::size_t i = 0; i < numPoints; ++i) {
    getNextToken();
    long id = expectInt();
    getNextToken();
    double x_1 = getNumber().value();
    getNextToken();
    double y_1 = getNumber().value();
    getNextToken();
    double z_1 = getNumber().value();
    getNextToken();
    long numTags = expectInt();
    // ignore tags
    for (std::size_t tagIdx = 0; tagIdx < numTags; tagIdx++) {
      getNextToken();
    }
    //logInfo() << "Point" << id << ":" << x_1 << "," << y_1 << "," << z_1 << "with" << numTags << "tags.";
  }
  //logInfo() << "Ignored all points, currently at" << curLoc.line << "," << curLoc.col;
  // ignore curves
  for (std::size_t i = 0; i < numCurves; ++i) {
    getNextToken();
    long id = expectInt();
    getNextToken();
    double x_1 = getNumber().value();
    getNextToken();
    double y_1 = getNumber().value();
    getNextToken();
    auto z_1 = getNumber().value();
    getNextToken();
    double x_2 = getNumber().value();
    getNextToken();
    double y_2 = getNumber().value();
    getNextToken();
    double z_2 = getNumber().value();
    getNextToken();
    long numPhysicalTags = expectInt();
    getNextToken();
    // ignore tags
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; tagIdx++) {
      getNextToken();
    }
    long numBoundaryTags = expectInt();
    // ignore boundary
    for (std::size_t tagIdx = 0; tagIdx < numBoundaryTags; tagIdx++) {
      getNextToken();
    }
    //logInfo() << "Line" << id << ": from" << x_1 << "," << y_1 << "," << z_1 << "to" << x_2 << "," << y_2 << "," << z_2 << "with" << numPhysicalTags << "physical tags and" << numBoundaryTags <<"boundary tags.";
  }
  //logInfo() << "Ignored all curves, currently at" << curLoc.line << "," << curLoc.col;
  // read surfaces
  for (std::size_t i = 0; i < numSurfaces; ++i) {
    getNextToken();
    long id = expectInt();
    getNextToken();
    double x_1 = getNumber().value();
    getNextToken();
    double y_1 = getNumber().value();
    getNextToken();
    double z_1 = getNumber().value();
    getNextToken();
    double x_2 = getNumber().value();
    getNextToken();
    double y_2 = getNumber().value();
    getNextToken();
    double z_2 = getNumber().value();
    getNextToken();
    long numPhysicalTags = expectInt();
    std::vector<long> tags;
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; ++tagIdx) {
      getNextToken();
      tags.push_back(expectInt());
    }
    if(!tags.empty()) {
      physicalSurfaceIds.insert_or_assign(id, tags.at(0));
    }
    getNextToken();
    long numBoundaryTags = expectInt();
    // ignore boundary
    for (std::size_t curveIdx = 0; curveIdx < numBoundaryTags; ++curveIdx) {
      getNextToken();
    }
    //logInfo() << "Surface" << id << ": from" << x_1 << "," << y_1 << "," << z_1 << "to" << x_2 << "," << y_2 << "," << z_2 << "with" << numPhysicalTags << "physical tags and" << numBoundaryTags <<"boundary tags.";
  }
  //logInfo() << "Processed all surfaces, currently at" << curLoc.line << "," << curLoc.col;
  // read volumes
  for (std::size_t i = 0; i < numVolumes; ++i) {
    getNextToken();
    long id = expectInt();
    getNextToken();
    double x_1 = getNumber().value();
    getNextToken();
    double y_1 = getNumber().value();
    getNextToken();
    double z_1 = getNumber().value();
    getNextToken();
    double x_2 = getNumber().value();
    getNextToken();
    double y_2 = getNumber().value();
    getNextToken();
    double z_2 = getNumber().value();
    getNextToken();
    long numPhysicalTags = expectInt();
    std::vector<long> tags;
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; ++tagIdx) {
      getNextToken();
      tags.push_back(expectInt());
    }
    physicalVolumeIds.insert_or_assign(id, tags.at(0));
    getNextToken();
    long numBoundaryTags = expectInt();
    // ignore boundary
    for (std::size_t curveIdx = 0; curveIdx < numBoundaryTags; ++curveIdx) {
      getNextToken();
    }
    //logInfo() << "Volume" << id << ": from" << x_1 << "," << y_1 << "," << z_1 << "to" << x_2 << "," << y_2 << "," << z_2 << "with" << numPhysicalTags << "physical tags and" << numBoundaryTags <<"boundary tags.";
  }
  //logInfo() << "Processed all volumes, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  if (curTok != puml::GMSHToken::end_entities) {
    return logErrorAnnotated<bool>("Expected $EndEntities");
  }
  getNextToken();
  return true;
}

bool GMSH4Parser::parseNodes() {
  //logInfo() << "Read nodes, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  std::size_t numBlocks = expectInt();
  getNextToken();
  std::size_t numVertices = expectInt();
  getNextToken();
  std::size_t minNodeTag = expectInt();
  getNextToken();
  std::size_t maxNodeTag = expectInt();
  builder->setNumVertices(numVertices);
  //logInfo() << numVertices << "vertices in" << numBlocks << "blocks with ids" << minNodeTag << "to" << maxNodeTag;

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    std::size_t dim = expectInt();
    getNextToken();
    std::size_t entityTag = expectInt();
    getNextToken();
    std::size_t parametric = expectInt();
    getNextToken();
    std::size_t numVerticesInBlock = expectInt();
    //logInfo() << "Read block" << blockIdx << "with" << numVerticesInBlock
    //          << "vertices: dim =" << dim << ", entityTag =" << entityTag
    //          << "parametric =" << parametric;
    std::vector<long> vertexIds;
    vertexIds.reserve(numVerticesInBlock);
    // first read vertex ids
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      getNextToken();
      vertexIds.push_back(expectInt() - 1);
    }
    // then read vertex data
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      std::array<double, 3> x;
      for (std::size_t i = 0; i < 3; i++) {
        getNextToken();
        x[i] = expectNumber();
      }
      //logInfo() << "Vertex" << vertexIds.at(vertexIdx) << "at" << x[0] << "," << x[1] << ","
      //          << x[2];
      assert(vertexIds.at(vertexIdx) < numVertices);
      builder->setVertex(vertexIds.at(vertexIdx), x);
    }
  }
  //logInfo() << "Processed all nodes, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  if (curTok != puml::GMSHToken::end_nodes) {
    return logErrorAnnotated<bool>("Expected $EndNodes");
  }
  getNextToken();
  return true;
}

bool GMSH4Parser::parseElements() {
  //logInfo() << "Read elements, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  std::size_t numBlocks = expectInt();
  getNextToken();
  std::size_t numElements = expectInt();
  getNextToken();
  std::size_t minElementTag = expectInt();
  getNextToken();
  std::size_t maxElementTag = expectInt();
  builder->setNumElements(numElements);
  //logInfo() << numElements << "elements in" << numBlocks << "blocks with ids" << minElementTag << "to" << maxElementTag;

  constexpr std::size_t MaxNodes = sizeof(NumNodes) / sizeof(std::size_t);
  std::array<long, MaxNodes> nodes;
  builder->setNumElements(numElements);

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    std::size_t dim = expectInt();
    getNextToken();
    std::size_t entityTag = expectInt();
    getNextToken();
    std::size_t type = expectInt();
    getNextToken();
    std::size_t numElementsInBlock = expectInt();
    //logInfo() << "Read block" << blockIdx << "with" << numElementsInBlock
    //          << "vertices: dim =" << dim << ", entityTag =" << entityTag << "type =" << ElementTypes[type-1];
    for (std::size_t elementIdx = 0; elementIdx < numElementsInBlock; ++elementIdx) {
      getNextToken();
      long id = expectInt();
      for (std::size_t nodeIdx = 0; nodeIdx < NumNodes[type-1]; nodeIdx++) {
        getNextToken();
        nodes.at(nodeIdx) = expectInt() - 1;
      }
      long tag = (type == 2) ? physicalSurfaceIds[entityTag] : physicalVolumeIds[entityTag];
      //logInfo() << "Found" << ElementTypes[type-1] << id << "with physical tag" << tag;
      builder->addElement(type, tag, nodes.data(), NumNodes[type - 1]);
    }
  }
  //logInfo() << "Processed all elements, currently at" << curLoc.line << "," << curLoc.col;
  getNextToken();
  if (curTok != puml::GMSHToken::end_elements) {
    return logErrorAnnotated<bool>("Expected $EndElements");
  }
  getNextToken();

  return true;
}


} // namespace puml
