// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSH4Parser.h"

#include <array>
#include <cstdio>
#include <string>
#include <string_view>
#include <vector>

#include "utils/logger.h"

namespace puml {

void GMSH4Parser::parse_() {
  const auto format = parseMeshFormatHeader();
  if (format.version < 4.0 || format.version >= 5.0) {
    char buf[128];
    snprintf(buf, sizeof(buf), "Unsupported MSH version %.1lf", format.version);
    failFile(buf);
  }
  if (format.version < 4.1) {
    failFile("MSH 4.0 is not supported; write the mesh as MSH 4.1 (Mesh.MshFileVersion = 4.1)");
  }
  if (format.fileType == 1) {
    failFile("Binary MSH files are not supported; write the mesh as ASCII (Mesh.Binary = 0)");
  }
  if (format.fileType != 0) {
    failFile("Expected file type 0 (ASCII)");
  }
  expectToken("$EndMeshFormat");

  bool hasEntities = false;
  bool hasElements = false;
  while (true) {
    const auto section = nextSection();
    if (section.empty()) {
      break;
    }
    if (section == "$Entities") {
      parseEntities();
      hasEntities = true;
    } else if (section == "$Nodes") {
      parseNodes();
    } else if (section == "$Elements") {
      parseElements();
      hasElements = true;
    } else if (section == "$Periodic") {
      parsePeriodic();
    } else if (section == "$PartitionedEntities") {
      failFile("Partitioned MSH files are not supported; save the mesh without partitions");
    } else {
      skipSection(section);
    }
  }

  if (!hasEntities) {
    failFile("Missing $Entities section");
  }
  if (!hasNodes) {
    failFile("Missing $Nodes section");
  }
  if (!hasElements) {
    failFile("Missing $Elements section");
  }
  if (!unassignedVolumes.empty() || !unassignedSurfaces.empty()) {
    logWarning() << "The elements of" << unassignedVolumes.size() << "volume(s) and"
                 << unassignedSurfaces.size()
                 << "surface(s) without a physical group get group or boundary condition 0";
  }
}

void GMSH4Parser::parseEntities() {
  const std::size_t numPoints = readSize();
  const std::size_t numCurves = readSize();
  const std::size_t numSurfaces = readSize();
  const std::size_t numVolumes = readSize();

  for (std::size_t i = 0; i < numPoints; ++i) {
    readInt();
    for (int k = 0; k < 3; ++k) {
      readDouble();
    }
    const std::size_t numPhysicalTags = readSize();
    for (std::size_t j = 0; j < numPhysicalTags; ++j) {
      readInt();
    }
  }

  const auto readEntities = [&](std::size_t count, std::map<long, long>* physicalIds) {
    for (std::size_t i = 0; i < count; ++i) {
      const long tag = readInt();
      for (int k = 0; k < 6; ++k) {
        readDouble();
      }
      const std::size_t numPhysicalTags = readSize();
      for (std::size_t j = 0; j < numPhysicalTags; ++j) {
        const long physical = readInt();
        if (j == 0 && physicalIds != nullptr) {
          physicalIds->insert_or_assign(tag, physical);
        }
      }
      const std::size_t numBoundingEntities = readSize();
      for (std::size_t j = 0; j < numBoundingEntities; ++j) {
        readInt();
      }
    }
  };
  readEntities(numCurves, nullptr);
  readEntities(numSurfaces, &physicalSurfaceIds);
  readEntities(numVolumes, &physicalVolumeIds);

  expectToken("$EndEntities");
}

std::size_t GMSH4Parser::readNodeIndex() {
  input->skipWhitespace();
  const auto offset = input->offset();
  const std::size_t tag = readSize();
  if (tag < firstNodeTag || tag - firstNodeTag >= numNodes) {
    failAt(offset, "Unknown node tag " + std::to_string(tag));
  }
  return tag - firstNodeTag;
}

void GMSH4Parser::parseNodes() {
  const std::size_t numBlocks = readSize();
  const std::size_t numVertices = readSize();
  const std::size_t minNodeTag = readSize();
  const std::size_t maxNodeTag = readSize();
  if (numVertices > 0 &&
      (minNodeTag == 0 || maxNodeTag < minNodeTag || maxNodeTag - minNodeTag + 1 != numVertices)) {
    char buf[192];
    snprintf(buf, sizeof(buf),
             "Non-contiguous node tags are not supported (%zu nodes with tags from %zu to %zu)",
             numVertices, minNodeTag, maxNodeTag);
    fail(buf);
  }
  hasNodes = true;
  firstNodeTag = minNodeTag;
  numNodes = numVertices;
  builder->setNumVertices(numVertices);

  // with as many distinct tags as the tag range is long, every node is defined exactly once
  std::vector<bool> defined(numVertices, false);
  std::vector<std::size_t> indices;

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    readInt(); // entity tag
    const long parametric = readInt();
    const std::size_t numVerticesInBlock = readSize();

    indices.resize(numVerticesInBlock);
    for (auto& index : indices) {
      input->skipWhitespace();
      const auto offset = input->offset();
      index = readNodeIndex();
      if (defined[index]) {
        failAt(offset, "Duplicate node tag " + std::to_string(firstNodeTag + index));
      }
      defined[index] = true;
    }

    // parametric nodes carry dim additional parametric coordinates
    const long numParametric = parametric != 0 ? dim : 0;
    for (const auto index : indices) {
      std::array<double, 3> x{};
      for (auto& coordinate : x) {
        coordinate = readDouble();
      }
      for (long i = 0; i < numParametric; ++i) {
        readDouble();
      }
      builder->setVertex(static_cast<long>(index), x);
    }
  }
  expectToken("$EndNodes");
}

void GMSH4Parser::parseElements() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Elements");
  }
  const std::size_t numBlocks = readSize();
  const std::size_t numElements = readSize();
  readSize(); // minimal element tag
  readSize(); // maximal element tag
  builder->setNumElements(numElements);

  std::array<long, MaxNodesPerElement> nodes{};

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    const long entityTag = readInt();
    input->skipWhitespace();
    const auto typeOffset = input->offset();
    const long type = readInt();
    if (type < 1 || type > static_cast<long>(NumTypes) || NumNodes[type - 1] == 0) {
      failAt(typeOffset, "Unknown element type " + std::to_string(type));
    }
    const std::size_t numNodesPerElement = NumNodes[type - 1];
    const std::size_t numElementsInBlock = readSize();
    const long tag = numElementsInBlock > 0 ? physicalTag(dim, entityTag) : 0;

    for (std::size_t elementIdx = 0; elementIdx < numElementsInBlock; ++elementIdx) {
      readSize(); // element tag
      for (std::size_t nodeIdx = 0; nodeIdx < numNodesPerElement; ++nodeIdx) {
        nodes[nodeIdx] = static_cast<long>(readNodeIndex());
      }
      builder->addElement(type, tag, nodes.data(), numNodesPerElement);
    }
  }
  expectToken("$EndElements");
}

void GMSH4Parser::parsePeriodic() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Periodic");
  }
  const std::size_t numPeriodic = readSize();

  // only the node identification is needed
  for (std::size_t blockIdx = 0; blockIdx < numPeriodic; ++blockIdx) {
    readInt(); // entity dimension
    readInt(); // entity tag
    readInt(); // master entity tag
    const std::size_t numAffine = readSize();
    for (std::size_t i = 0; i < numAffine; ++i) {
      readDouble();
    }
    const std::size_t numIdentified = readSize();
    for (std::size_t i = 0; i < numIdentified; ++i) {
      const std::size_t node = readNodeIndex();
      const std::size_t master = readNodeIndex();
      builder->addVertexLink(node, master);
    }
  }
  expectToken("$EndPeriodic");
}

long GMSH4Parser::physicalTag(std::size_t dim, long entityTag) {
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

} // namespace puml
