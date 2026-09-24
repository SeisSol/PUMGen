// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSH4Parser.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include "utils/logger.h"

namespace puml {

namespace {
template <typename T> T swapBytesOf(T value) {
  std::array<unsigned char, sizeof(T)> bytes{};
  std::memcpy(bytes.data(), &value, sizeof(T));
  std::reverse(bytes.begin(), bytes.end());
  std::memcpy(&value, bytes.data(), sizeof(T));
  return value;
}

// blocks of binary element data are read in pieces of this many elements
constexpr std::size_t ElementChunk = 1 << 16;
} // namespace

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
    binary = true;
    dataSize = format.dataSize;
    if (dataSize != 4 && dataSize != 8) {
      failFile("Unsupported data size " + std::to_string(dataSize) + " (expected 4 or 8)");
    }
    beginSectionData();
    std::int32_t one = 0;
    readBinary(&one, sizeof(one));
    if (swapBytesOf(one) == 1) {
      swapBytes = true;
    } else if (one != 1) {
      failFile("Invalid byte order mark in the binary MSH file");
    }
  } else if (format.fileType != 0) {
    failFile("Expected file type 0 (ASCII) or 1 (binary)");
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
    } else if (section == "$PhysicalNames") {
      parsePhysicalNames();
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

void GMSH4Parser::beginSectionData() {
  if (binary) {
    if (!input->readLineBreak()) {
      fail("Expected a line break before the binary data");
    }
    binaryData = true;
  }
}

void GMSH4Parser::readBinary(void* data, std::size_t bytes) {
  if (!input->readRaw(data, bytes)) {
    failFile("Unexpected end of file in binary data");
  }
}

std::size_t GMSH4Parser::readSize() {
  if (!binary) {
    return expectSize();
  }
  std::size_t value = 0;
  readSizes(&value, 1);
  return value;
}

long GMSH4Parser::readInt() {
  if (!binary) {
    return expectInteger();
  }
  std::int32_t value = 0;
  readBinary(&value, sizeof(value));
  return swapBytes ? swapBytesOf(value) : value;
}

double GMSH4Parser::readDouble() {
  if (!binary) {
    return expectNumber();
  }
  double value = 0;
  readDoubles(&value, 1);
  return value;
}

void GMSH4Parser::readSizes(std::size_t* values, std::size_t count) {
  if (!binary) {
    for (std::size_t i = 0; i < count; ++i) {
      values[i] = expectSize();
    }
    return;
  }
  if (dataSize == sizeof(std::uint64_t)) {
    static_assert(sizeof(std::size_t) == sizeof(std::uint64_t));
    readBinary(values, count * sizeof(std::uint64_t));
    if (swapBytes) {
      std::transform(values, values + count, values, swapBytesOf<std::size_t>);
    }
  } else {
    narrowSizes.resize(count);
    readBinary(narrowSizes.data(), count * sizeof(std::uint32_t));
    for (std::size_t i = 0; i < count; ++i) {
      values[i] = swapBytes ? swapBytesOf(narrowSizes[i]) : narrowSizes[i];
    }
  }
}

void GMSH4Parser::readDoubles(double* values, std::size_t count) {
  if (!binary) {
    for (std::size_t i = 0; i < count; ++i) {
      values[i] = expectNumber();
    }
    return;
  }
  readBinary(values, count * sizeof(double));
  if (swapBytes) {
    std::transform(values, values + count, values, swapBytesOf<double>);
  }
}

std::size_t GMSH4Parser::position() {
  if (!binary) {
    input->skipWhitespace();
  }
  return input->offset();
}

std::size_t GMSH4Parser::nodeIndex(std::size_t tag, std::size_t offset) {
  if (tag < firstNodeTag || tag - firstNodeTag >= numNodes) {
    failAt(offset, "Unknown node tag " + std::to_string(tag));
  }
  return tag - firstNodeTag;
}

void GMSH4Parser::parseEntities() {
  beginSectionData();
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

void GMSH4Parser::parseNodes() {
  beginSectionData();
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
  std::vector<double> coordinates;

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    readInt(); // entity tag
    const long parametric = readInt();
    const std::size_t numVerticesInBlock = readSize();

    indices.resize(numVerticesInBlock);
    const auto tagsOffset = input->offset();
    if (binary) {
      readSizes(indices.data(), numVerticesInBlock);
    }
    for (std::size_t i = 0; i < numVerticesInBlock; ++i) {
      auto offset = tagsOffset + i * dataSize;
      if (!binary) {
        offset = position();
        indices[i] = expectSize();
      }
      const std::size_t index = nodeIndex(indices[i], offset);
      if (defined[index]) {
        failAt(offset, "Duplicate node tag " + std::to_string(indices[i]));
      }
      defined[index] = true;
      indices[i] = index;
    }

    // parametric nodes carry dim additional parametric coordinates
    const std::size_t valuesPerNode = 3 + (parametric != 0 ? static_cast<std::size_t>(dim) : 0);
    coordinates.resize(numVerticesInBlock * valuesPerNode);
    readDoubles(coordinates.data(), coordinates.size());
    for (std::size_t i = 0; i < numVerticesInBlock; ++i) {
      const double* x = &coordinates[i * valuesPerNode];
      builder->setVertex(static_cast<long>(indices[i]), {x[0], x[1], x[2]});
    }
  }
  expectToken("$EndNodes");
}

void GMSH4Parser::parseElements() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Elements");
  }
  beginSectionData();
  const std::size_t numBlocks = readSize();
  const std::size_t numElements = readSize();
  readSize(); // minimal element tag
  readSize(); // maximal element tag
  builder->setNumElements(numElements);

  std::array<long, MaxNodesPerElement> nodes{};
  std::vector<std::size_t> values;

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    const long entityTag = readInt();
    const auto typeOffset = position();
    const long type = readInt();
    if (type < 1 || type > static_cast<long>(NumTypes) || NumNodes[type - 1] == 0) {
      failAt(typeOffset, "Unknown element type " + std::to_string(type));
    }
    const std::size_t numNodesPerElement = NumNodes[type - 1];
    const std::size_t numElementsInBlock = readSize();
    const long tag = numElementsInBlock > 0 ? physicalTag(dim, entityTag) : 0;

    // each element consists of its tag and its node tags
    const std::size_t valuesPerElement = 1 + numNodesPerElement;
    for (std::size_t first = 0; first < numElementsInBlock; first += ElementChunk) {
      const std::size_t count = std::min(ElementChunk, numElementsInBlock - first);
      const auto chunkOffset = input->offset();
      if (binary) {
        values.resize(count * valuesPerElement);
        readSizes(values.data(), values.size());
      }
      for (std::size_t element = 0; element < count; ++element) {
        if (!binary) {
          readSize(); // element tag
        }
        for (std::size_t node = 0; node < numNodesPerElement; ++node) {
          if (binary) {
            const std::size_t valueIdx = element * valuesPerElement + 1 + node;
            nodes[node] =
                static_cast<long>(nodeIndex(values[valueIdx], chunkOffset + valueIdx * dataSize));
          } else {
            const auto offset = position();
            nodes[node] = static_cast<long>(nodeIndex(expectSize(), offset));
          }
        }
        builder->addElement(type, tag, nodes.data(), numNodesPerElement);
      }
    }
  }
  expectToken("$EndElements");
}

void GMSH4Parser::parsePeriodic() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Periodic");
  }
  beginSectionData();
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
      const auto nodeOffset = position();
      const std::size_t node = nodeIndex(readSize(), nodeOffset);
      const auto masterOffset = position();
      const std::size_t master = nodeIndex(readSize(), masterOffset);
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
