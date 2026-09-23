// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "Msh4Index.h"

#include "mesh/CellType.h"

#include <cstdio>
#include <string>

namespace puml {

void Msh4Indexer::skipData(std::size_t bytes) {
  if (!input->seek(input->offset() + bytes)) {
    failFile("Unexpected end of file in binary data");
  }
}

void Msh4Indexer::parseNodes() {
  if (!binary) {
    failFile("Only binary MSH 4.1 files can be read in parallel");
  }
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
  index.dataSize = dataSize;
  index.swapBytes = swapBytes;
  index.firstNodeTag = minNodeTag;
  index.numNodes = numVertices;

  std::size_t numInBlocks = 0;
  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    readInt(); // entity tag
    const long parametric = readInt();
    const std::size_t count = readSize();
    Msh4NodeBlock block{};
    block.count = count;
    block.valuesPerNode = 3 + (parametric != 0 ? static_cast<std::uint64_t>(dim) : 0);
    block.tagsOffset = input->offset();
    block.coordinatesOffset = block.tagsOffset + count * dataSize;
    skipData(count * dataSize + count * block.valuesPerNode * sizeof(double));
    if (count > 0) {
      index.nodeBlocks.push_back(block);
    }
    numInBlocks += count;
  }
  if (numInBlocks != numVertices) {
    failFile("The node blocks contain " + std::to_string(numInBlocks) + " nodes instead of " +
             std::to_string(numVertices));
  }
  expectToken("$EndNodes");
}

void Msh4Indexer::parseElements() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Elements");
  }
  beginSectionData();
  const std::size_t numBlocks = readSize();
  readSize(); // number of elements
  readSize(); // minimal element tag
  readSize(); // maximal element tag

  for (std::size_t blockIdx = 0; blockIdx < numBlocks; ++blockIdx) {
    const long dim = readInt();
    const long entityTag = readInt();
    const auto typeOffset = position();
    const long type = readInt();
    if (type < 1 || type > static_cast<long>(NumTypes) || NumNodes[type - 1] == 0) {
      failAt(typeOffset, "Unknown element type " + std::to_string(type));
    }
    const std::size_t count = readSize();
    const Msh4ElementBlock block{count, input->offset(),
                                 count > 0 ? physicalTag(dim, entityTag) : 0, type,
                                 NumNodes[type - 1]};
    skipData(count * (1 + NumNodes[type - 1]) * dataSize);
    const auto element = gmshElementType(type);
    if (count > 0 && element.kind == GmshElementType::Kind::Cell) {
      if (!element.complete) {
        failAt(typeOffset, "Element type " + std::to_string(type) + " is a " +
                               shapeOf(element.cellType).name +
                               " without all nodes of a Lagrange cell, which is not supported");
      }
      index.highOrder = index.highOrder || element.order > 1;
      index.cellBlocks.push_back(block);
    } else if (count > 0 && element.kind == GmshElementType::Kind::Face) {
      index.facetBlocks.push_back(block);
    }
  }
  expectToken("$EndElements");
}

void Msh4Indexer::parsePeriodic() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Periodic");
  }
  beginSectionData();
  const std::size_t numPeriodic = readSize();
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
      const std::size_t node = readSize();
      nodeIndex(node, nodeOffset);
      const auto masterOffset = position();
      const std::size_t master = readSize();
      nodeIndex(master, masterOffset);
      index.periodic.emplace_back(node, master);
    }
  }
  expectToken("$EndPeriodic");
}

bool hasHighOrderCells(const std::string& fileName) {
  Msh4Indexer indexer;
  return indexer.parseFile(fileName) && indexer.getIndex().highOrder;
}

bool isBinaryMsh4(const std::string& fileName) {
  try {
    MshInput input(fileName);
    if (input.readToken() != "$MeshFormat") {
      return false;
    }
    const auto version = input.readReal();
    const auto fileType = input.readInteger<long>();
    return version && fileType && *version >= 4.1 && *version < 5.0 && *fileType == 1;
  } catch (const std::exception&) {
    return false;
  }
}

} // namespace puml
