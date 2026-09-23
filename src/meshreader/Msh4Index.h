// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_MSH4INDEX_H_
#define PUMGEN_SRC_MESHREADER_MSH4INDEX_H_

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "GMSH4Parser.h"

namespace puml {

/**
 * A block of nodes in a binary MSH 4.1 file: first the tags, then the coordinates.
 */
struct Msh4NodeBlock {
  std::uint64_t count;
  std::uint64_t tagsOffset;
  std::uint64_t coordinatesOffset;
  // three coordinates and the parametric coordinates, if any
  std::uint64_t valuesPerNode;
};

/**
 * A block of elements of one type in a binary MSH 4.1 file: per element its tag and its nodes.
 */
struct Msh4ElementBlock {
  std::uint64_t count;
  std::uint64_t offset;
  std::int64_t physical;
  /** The MSH element type */
  std::int64_t type;
  std::uint64_t nodesPerElement;
};

/**
 * Where the data of a binary MSH 4.1 file lies.
 */
struct Msh4Index {
  std::uint64_t dataSize = 8;
  bool swapBytes = false;
  // the node tags form the contiguous range [firstNodeTag, firstNodeTag + numNodes)
  std::uint64_t firstNodeTag = 0;
  std::uint64_t numNodes = 0;
  std::vector<Msh4NodeBlock> nodeBlocks;
  std::vector<Msh4ElementBlock> cellBlocks;
  std::vector<Msh4ElementBlock> facetBlocks;
  // whether a cell is of an order higher than one
  bool highOrder = false;
  // pairs of node tags (node, master) from $Periodic
  std::vector<std::pair<std::uint64_t, std::uint64_t>> periodic;
};

/**
 * Reads the structure of a binary MSH 4.1 file (entities, block headers, periodic nodes) and
 * skips the node and element data, recording where they lie.
 */
class Msh4Indexer : public GMSH4Parser {
  public:
  Msh4Indexer() : GMSH4Parser(nullptr) {}

  [[nodiscard]] const Msh4Index& getIndex() const { return index; }

  private:
  void parseNodes() override;
  void parseElements() override;
  void parsePeriodic() override;

  void skipData(std::size_t bytes);

  Msh4Index index;
};

/**
 * Whether the file is a binary MSH 4.1 file; only reads its header.
 */
bool isBinaryMsh4(const std::string& fileName);

/**
 * Whether a binary MSH 4.1 file holds cells of an order higher than one; reads the structure of the
 * file only. False for a file which cannot be indexed, whose errors the reader then reports.
 */
bool hasHighOrderCells(const std::string& fileName);

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_MSH4INDEX_H_
