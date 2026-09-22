// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSH4Parser.h"

#include <cstdio>

#include "utils/logger.h"

namespace puml {
bool GMSH4Parser::parse_() {
  getNextToken();

  const double version = parseMeshFormat();
  if (version < 4.0 || version >= 5.0) {
    char buf[128];
    snprintf(buf, sizeof(buf), "Unsupported MSH version %.1lf", version);
    return logError<bool>(buf);
  }
  if (version < 4.1) {
    return logError<bool>(
        "MSH 4.0 is not supported; write the mesh as MSH 4.1 (Mesh.MshFileVersion = 4.1)");
  }

  bool hasEntities = false;
  bool hasElements = false;
  bool hasPeriodic = false;

  while (curTok != tndm::GMSHToken::eof) {
    switch (curTok) {
    case tndm::GMSHToken::entities:
      hasEntities = parseEntities();
      break;
    case tndm::GMSHToken::nodes:
      parseNodes();
      break;
    case tndm::GMSHToken::elements:
      hasElements = parseElements();
      break;
    case tndm::GMSHToken::periodic:
      hasPeriodic = parsePeriodic();
      break;
    case tndm::GMSHToken::partitioned_entities:
      return logErrorAnnotated<bool>(
          "Partitioned MSH files are not supported; save the mesh without partitions");
    default:
      getNextToken();
      break;
    }
  }

  if (!hasEntities) {
    return logError<bool>("Missing $Entities section");
  }
  if (!hasNodes) {
    return logError<bool>("Missing $Nodes section");
  }
  if (!hasElements) {
    return logError<bool>("Missing $Elements section");
  }
  if (!unassignedVolumes.empty() || !unassignedSurfaces.empty()) {
    logWarning() << "The elements of" << unassignedVolumes.size() << "volume(s) and"
                 << unassignedSurfaces.size()
                 << "surface(s) without a physical group get group or boundary condition 0";
  }
  return true;
}

long GMSH4Parser::physicalTag(std::size_t dim, unsigned long entityTag) {
  if (dim != 2 && dim != 3) {
    return 0;
  }
  const auto& physicalIds = dim == 3 ? physicalVolumeIds : physicalSurfaceIds;
  const auto it = physicalIds.find(entityTag);
  if (it == physicalIds.end()) {
    (dim == 3 ? unassignedVolumes : unassignedSurfaces).insert(entityTag);
    return 0;
  }
  return it->second;
}

bool GMSH4Parser::parseEntities() {
  getNextToken();
  const std::size_t numPoints = expectNonNegativeInt();
  getNextToken();
  const std::size_t numCurves = expectNonNegativeInt();
  getNextToken();
  const std::size_t numSurfaces = expectNonNegativeInt();
  getNextToken();
  const std::size_t numVolumes = expectNonNegativeInt();

  // ignore points
  for (std::size_t i = 0; i < numPoints; ++i) {
    getNextToken();
    [[maybe_unused]] const std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_1 = expectNumber();
    getNextToken();
    const std::size_t numTags = expectNonNegativeInt();
    // ignore tags
    for (std::size_t tagIdx = 0; tagIdx < numTags; tagIdx++) {
      getNextToken();
    }
  }

  //  ignore curves
  for (std::size_t i = 0; i < numCurves; ++i) {
    getNextToken();
    [[maybe_unused]] const std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_2 = expectNumber();
    getNextToken();
    const std::size_t numPhysicalTags = expectNonNegativeInt();
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
    const std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_2 = expectNumber();
    getNextToken();
    const std::size_t numPhysicalTags = expectNonNegativeInt();
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
    const std::size_t id = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const double x_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_1 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double x_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double y_2 = expectNumber();
    getNextToken();
    [[maybe_unused]] const double z_2 = expectNumber();
    getNextToken();
    const std::size_t numPhysicalTags = expectNonNegativeInt();
    std::vector<std::size_t> tags;
    for (std::size_t tagIdx = 0; tagIdx < numPhysicalTags; ++tagIdx) {
      getNextToken();
      tags.push_back(expectNonNegativeInt());
    }
    if (!tags.empty()) {
      physicalVolumeIds.insert_or_assign(id, tags.at(0));
    }
    getNextToken();
    const std::size_t numBoundaryTags = expectNonNegativeInt();
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
  const std::size_t numBlocks = expectNonNegativeInt();
  getNextToken();
  const std::size_t numVertices = expectNonNegativeInt();
  getNextToken();
  const std::size_t minNodeTag = expectNonNegativeInt();
  getNextToken();
  const std::size_t maxNodeTag = expectNonNegativeInt();
  if (numVertices > 0 &&
      (minNodeTag == 0 || maxNodeTag < minNodeTag || maxNodeTag - minNodeTag + 1 != numVertices)) {
    char buf[192];
    snprintf(buf, sizeof(buf),
             "Non-contiguous node tags are not supported (%zu nodes with tags from %zu to %zu)",
             numVertices, minNodeTag, maxNodeTag);
    return logErrorAnnotated<bool>(buf);
  }
  hasNodes = true;
  firstNodeTag = minNodeTag;
  numNodes = numVertices;
  builder->setNumVertices(numVertices);

  // with as many distinct tags as the tag range is long, every node is defined exactly once
  std::vector<bool> defined(numVertices, false);

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    const std::size_t dim = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const std::size_t entityTag = expectNonNegativeInt();
    getNextToken();
    const std::size_t parametric = expectNonNegativeInt();
    getNextToken();
    const std::size_t numVerticesInBlock = expectNonNegativeInt();
    std::vector<std::size_t> vertexIds;
    vertexIds.reserve(numVerticesInBlock);
    // first read vertex ids
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      getNextToken();
      const std::size_t index = expectNodeIndex();
      if (defined[index]) {
        char buf[128];
        snprintf(buf, sizeof(buf), "Duplicate node tag %zu", firstNodeTag + index);
        return logErrorAnnotated<bool>(buf);
      }
      defined[index] = true;
      vertexIds.push_back(index);
    }
    // then read vertex data; parametric nodes carry dim additional parametric coordinates
    const std::size_t numParametric = parametric != 0 ? dim : 0;
    for (std::size_t vertexIdx = 0; vertexIdx < numVerticesInBlock; ++vertexIdx) {
      std::array<double, 3> x{};
      for (std::size_t i = 0; i < 3; i++) {
        getNextToken();
        x[i] = expectNumber();
      }
      for (std::size_t i = 0; i < numParametric; ++i) {
        getNextToken();
        expectNumber();
      }
      builder->setVertex(vertexIds[vertexIdx], x);
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
  if (!hasNodes) {
    return logErrorAnnotated<bool>("Expected $Nodes before $Elements");
  }
  getNextToken();
  const std::size_t numBlocks = expectNonNegativeInt();
  getNextToken();
  const std::size_t numElements = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] const std::size_t minElementTag = expectNonNegativeInt();
  getNextToken();
  [[maybe_unused]] const std::size_t maxElementTag = expectNonNegativeInt();
  builder->setNumElements(numElements);

  std::array<long, MaxNodesPerElement> nodes{};

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; blockIdx++) {
    getNextToken();
    const std::size_t dim = expectNonNegativeInt();
    getNextToken();
    const std::size_t entityTag = expectNonNegativeInt();
    getNextToken();
    const std::size_t type = expectNonNegativeInt();
    if (type < 1 || type > NumTypes || NumNodes[type - 1] == 0) {
      char buf[128];
      snprintf(buf, sizeof(buf), "Unknown element type %zu", type);
      return logErrorAnnotated<bool>(buf);
    }
    getNextToken();
    const std::size_t numElementsInBlock = expectNonNegativeInt();
    const long tag = numElementsInBlock > 0 ? physicalTag(dim, entityTag) : 0;
    for (std::size_t elementIdx = 0; elementIdx < numElementsInBlock; ++elementIdx) {
      getNextToken();
      [[maybe_unused]] const std::size_t id = expectNonNegativeInt();
      for (std::size_t nodeIdx = 0; nodeIdx < NumNodes[type - 1]; nodeIdx++) {
        getNextToken();
        nodes[nodeIdx] = static_cast<long>(expectNodeIndex());
      }
      builder->addElement(static_cast<long>(type), tag, nodes.data(), NumNodes[type - 1]);
    }
  }
  getNextToken();
  if (curTok != tndm::GMSHToken::end_elements) {
    return logErrorAnnotated<bool>("Expected $EndElements");
  }
  getNextToken();

  return true;
}

bool GMSH4Parser::parsePeriodic() {
  if (!hasNodes) {
    return logErrorAnnotated<bool>("Expected $Nodes before $Periodic");
  }
  getNextToken();
  const auto numPeriodic = expectNonNegativeInt();

  // ignore everything but the node identification

  for (std::size_t blockIdx = 0; blockIdx < numPeriodic; blockIdx++) {
    getNextToken();
    [[maybe_unused]] const std::size_t entityDim = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const std::size_t entityId = expectNonNegativeInt();
    getNextToken();
    [[maybe_unused]] const std::size_t entityIdentifyId = expectNonNegativeInt();

    getNextToken();
    const std::size_t affineSize = expectNonNegativeInt();
    for (std::size_t i = 0; i < affineSize; ++i) {
      getNextToken();
      [[maybe_unused]] const double affineValue = expectNumber();
    }

    getNextToken();
    const std::size_t identifySize = expectNonNegativeInt();

    for (std::size_t i = 0; i < identifySize; ++i) {
      getNextToken();
      const std::size_t nodeId = expectNodeIndex();
      getNextToken();
      const std::size_t nodeIdentifyId = expectNodeIndex();
      builder->addVertexLink(nodeId, nodeIdentifyId);
    }
  }

  getNextToken();
  if (curTok != tndm::GMSHToken::end_periodic) {
    return logErrorAnnotated<bool>("Expected $EndPeriodic");
  }
  getNextToken();

  return true;
}

} // namespace puml
