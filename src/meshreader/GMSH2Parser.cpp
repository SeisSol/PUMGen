// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSH2Parser.h"

#include <array>
#include <cstdio>
#include <string>
#include <vector>

namespace puml {

void GMSH2Parser::parse_() {
  const auto format = parseMeshFormatHeader();
  if (format.version < 2.0 || format.version >= 3.0) {
    char buf[128];
    snprintf(buf, sizeof(buf), "Unsupported MSH version %.1lf", format.version);
    failFile(buf);
  }
  if (format.fileType != 0) {
    failFile("Binary MSH 2 files are not supported; write the mesh as ASCII (Mesh.Binary = 0)");
  }
  expectToken("$EndMeshFormat");

  bool hasElements = false;
  while (true) {
    const auto section = nextSection();
    if (section.empty()) {
      break;
    }
    if (section == "$Nodes") {
      parseNodes();
    } else if (section == "$Elements") {
      parseElements();
      hasElements = true;
    } else if (section == "$Periodic") {
      parsePeriodic();
    } else {
      skipSection(section);
    }
  }

  if (!hasNodes) {
    failFile("Missing $Nodes section");
  }
  if (!hasElements) {
    failFile("Missing $Elements section");
  }
}

std::size_t GMSH2Parser::readNodeIndex() {
  input->skipWhitespace();
  const auto offset = input->offset();
  const auto tag = input->readInteger<std::size_t>();
  if (!tag || *tag < 1 || *tag > numNodes) {
    failAt(offset, "Unknown node tag");
  }
  return *tag - 1;
}

void GMSH2Parser::parseNodes() {
  const std::size_t numVertices = expectSize();
  builder->setNumVertices(numVertices);
  hasNodes = true;
  numNodes = numVertices;

  // with numVertices distinct tags in [1, numVertices], every node is defined exactly once
  std::vector<bool> defined(numVertices, false);

  for (std::size_t i = 0; i < numVertices; ++i) {
    input->skipWhitespace();
    const auto offset = input->offset();
    const auto tag = input->readInteger<std::size_t>();
    if (!tag || *tag < 1 || *tag > numVertices) {
      failAt(offset, "Expected node-tag with 1 <= node-tag <= " + std::to_string(numVertices));
    }
    const std::size_t id = *tag - 1;
    if (defined[id]) {
      failAt(offset, "Duplicate node tag " + std::to_string(*tag));
    }
    defined[id] = true;

    std::array<double, 3> x{};
    for (auto& coordinate : x) {
      coordinate = expectNumber();
    }
    builder->setVertex(static_cast<long>(id), x);
  }
  expectToken("$EndNodes");
}

void GMSH2Parser::parseElements() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Elements");
  }
  const std::size_t numElements = expectSize();
  builder->setNumElements(numElements);

  std::array<long, MaxNodesPerElement> nodes{};

  for (std::size_t elementIdx = 0; elementIdx < numElements; ++elementIdx) {
    expectInteger(); // element tag
    input->skipWhitespace();
    const auto typeOffset = input->offset();
    const long type = expectInteger();
    if (type < 1 || type > static_cast<long>(NumTypes) || NumNodes[type - 1] == 0) {
      failAt(typeOffset, "Unknown element type " + std::to_string(type));
    }
    const std::size_t numTags = expectSize();
    // the physical tag; elements without tags get 0
    long tag = 0;
    for (std::size_t tagIdx = 0; tagIdx < numTags; ++tagIdx) {
      const long value = expectInteger();
      if (tagIdx == 0) {
        tag = value;
      }
    }
    const std::size_t numNodesPerElement = NumNodes[type - 1];
    for (std::size_t nodeIdx = 0; nodeIdx < numNodesPerElement; ++nodeIdx) {
      nodes[nodeIdx] = static_cast<long>(readNodeIndex());
    }
    builder->addElement(type, tag, nodes.data(), numNodesPerElement);
  }
  expectToken("$EndElements");
}

void GMSH2Parser::parsePeriodic() {
  if (!hasNodes) {
    fail("Expected $Nodes before $Periodic");
  }
  const std::size_t numPeriodic = expectSize();

  // only the node identification is needed
  for (std::size_t blockIdx = 0; blockIdx < numPeriodic; ++blockIdx) {
    expectSize();    // entity dimension
    expectInteger(); // entity tag
    expectInteger(); // master entity tag
    if (input->peek() == 'A') {
      // "Affine" and 16 values
      expectToken("Affine");
      for (int i = 0; i < 16; ++i) {
        expectNumber();
      }
    }
    const std::size_t numIdentified = expectSize();
    for (std::size_t i = 0; i < numIdentified; ++i) {
      const std::size_t node = readNodeIndex();
      const std::size_t master = readNodeIndex();
      builder->addVertexLink(node, master);
    }
  }
  expectToken("$EndPeriodic");
}

} // namespace puml
