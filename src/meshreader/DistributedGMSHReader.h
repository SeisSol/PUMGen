// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_DISTRIBUTEDGMSHREADER_H_
#define PUMGEN_SRC_MESHREADER_DISTRIBUTEDGMSHREADER_H_

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Msh4Index.h"
#include "helper/Distributor.h"
#include "helper/MPIConvenience.h"
#include "input/MeshData.h"
#include "mesh/CellType.h"
#include "third_party/MPITraits.h"
#include "utils/logger.h"

namespace puml {

namespace distributed {

/**
 * Reads values of a binary MSH 4.1 file at given offsets.
 */
class RawMshFile {
  public:
  RawMshFile(const std::string& fileName, const Msh4Index& index);
  ~RawMshFile();
  RawMshFile(const RawMshFile&) = delete;
  RawMshFile& operator=(const RawMshFile&) = delete;

  void readSizes(std::uint64_t offset, std::size_t count, std::uint64_t* values);
  void readDoubles(std::uint64_t offset, std::size_t count, double* values);

  private:
  void readBytes(std::uint64_t offset, std::size_t bytes, void* data);

  std::FILE* file = nullptr;
  std::uint64_t dataSize;
  bool swapBytes;
  std::vector<std::uint32_t> narrowSizes;
};

/**
 * Sends outgoing[r] to rank r; returns what each rank has sent to this one.
 */
template <typename T>
std::vector<std::vector<T>> exchange(const std::vector<std::vector<T>>& outgoing, MPI_Comm comm) {
  static_assert(std::is_trivially_copyable_v<T>);
  int size = 1;
  MPI_Comm_size(comm, &size);
  std::vector<std::size_t> sendCounts(size);
  std::vector<std::size_t> recvCounts(size);
  for (int r = 0; r < size; ++r) {
    sendCounts[r] = outgoing[r].size();
  }
  MPI_Alltoall(sendCounts.data(), 1, tndm::mpi_type_t<std::size_t>(), recvCounts.data(), 1,
               tndm::mpi_type_t<std::size_t>(), comm);
  std::vector<std::size_t> sendDisp(size + 1, 0);
  std::vector<std::size_t> recvDisp(size + 1, 0);
  for (int r = 0; r < size; ++r) {
    sendDisp[r + 1] = sendDisp[r] + sendCounts[r];
    recvDisp[r + 1] = recvDisp[r] + recvCounts[r];
  }
  std::vector<T> sendBuffer(sendDisp[size]);
  for (int r = 0; r < size; ++r) {
    std::copy(outgoing[r].begin(), outgoing[r].end(), sendBuffer.begin() + sendDisp[r]);
  }
  std::vector<T> recvBuffer(recvDisp[size]);

  MPI_Datatype type;
  MPI_Type_contiguous(sizeof(T), MPI_BYTE, &type);
  MPI_Type_commit(&type);
  sparseAlltoallv(sendBuffer.data(), sendCounts.data(), sendDisp.data(), type, recvBuffer.data(),
                  recvCounts.data(), recvDisp.data(), type, comm);
  MPI_Type_free(&type);

  std::vector<std::vector<T>> incoming(size);
  for (int r = 0; r < size; ++r) {
    incoming[r].assign(recvBuffer.begin() + recvDisp[r], recvBuffer.begin() + recvDisp[r + 1]);
  }
  return incoming;
}

/**
 * Calls f(block, begin, count) for the parts of the items [first, last) of the concatenated
 * blocks, where begin is relative to the block.
 */
template <typename Block, typename F>
void forBlockRanges(const std::vector<Block>& blocks, std::uint64_t first, std::uint64_t last,
                    F&& f) {
  std::uint64_t blockStart = 0;
  for (const auto& block : blocks) {
    const std::uint64_t blockEnd = blockStart + block.count;
    const std::uint64_t begin = std::max(first, blockStart);
    const std::uint64_t end = std::min(last, blockEnd);
    if (begin < end) {
      f(block, begin - blockStart, end - begin);
    }
    blockStart = blockEnd;
  }
}

std::uint64_t hashFace(const std::array<std::uint64_t, 4>& face);

// items are read in pieces of at most this many
constexpr std::uint64_t Piece = std::uint64_t{1} << 20;

} // namespace distributed

/**
 * Reads a binary MSH 4.1 file with all ranks: rank 0 only locates the data, every rank reads a
 * share of the nodes, cells and boundary faces, and the parts are sent to the ranks which own them
 * in the chunk distribution of getChunksize. The boundary conditions are matched in a distributed
 * table of the boundary faces; only faces whose vertices all lie on the surface mesh are looked up.
 * In a mesh with cells of higher order, the vertices are the nodes which are a vertex of a cell,
 * numbered in the order of the nodes, and every rank fetches the other nodes of its cells from the
 * ranks which read them.
 */
class DistributedGMSHReader {
  public:
  constexpr static bool ProvidesLocalMesh = true;

  explicit DistributedGMSHReader(MPI_Comm comm = MPI_COMM_WORLD) : comm(comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  }

  void open(const char* meshFile) {
    readIndex(meshFile);
    distributed::RawMshFile file(meshFile, index);
    readNodes(file);
    readCells(file);
    matchBoundaryConditions(file);
    if (highOrder) {
      findVertices();
      gatherHighOrderGeometry();
      distributeVertices();
      for (auto& vertex : local.connectivity) {
        vertex = vertexOf(vertex);
      }
    } else {
      numVertices = index.numNodes;
      local.geometry = std::move(nodeCoordinates);
    }
    std::vector<double>().swap(nodeCoordinates);
    identifyPeriodicNodes();
  }

  /**
   * The part of the mesh of this rank; possible only once.
   */
  LocalMesh read() { return std::move(local); }

  private:
  using Face = std::array<std::uint64_t, MaxFaceVertices>;
  static constexpr std::uint64_t NoVertex = std::numeric_limits<std::uint64_t>::max();

  struct NodeRecord {
    std::uint64_t index;
    double x[3];
  };
  struct FaceRecord {
    Face vertices;
    std::int64_t bc;
  };
  struct FaceQuery {
    Face vertices;
    std::uint64_t localFace;
  };
  struct FaceReply {
    std::uint64_t localFace;
    std::int64_t bc;
  };
  struct Identification {
    std::uint64_t vertex;
    std::uint64_t representative;
  };
  constexpr static std::int64_t NoMatch = -1;

  [[nodiscard]] std::uint64_t localFirstNode() const {
    return getChunksum(index.numNodes, rank, size);
  }
  [[nodiscard]] std::uint64_t localNodeCount() const {
    return getChunksize(index.numNodes, rank, size);
  }
  [[nodiscard]] int nodeOwner(std::uint64_t node) const {
    return getChunkOwner(index.numNodes, node, size);
  }
  [[nodiscard]] std::uint64_t localFirstVertex() const {
    return getChunksum(numVertices, rank, size);
  }
  [[nodiscard]] std::uint64_t localVertexCount() const {
    return getChunksize(numVertices, rank, size);
  }
  [[nodiscard]] int vertexOwner(std::uint64_t vertex) const {
    return getChunkOwner(numVertices, vertex, size);
  }

  /** Whether a node is a vertex; every node is one in a mesh of linear cells */
  [[nodiscard]] bool isVertex(std::uint64_t node) const {
    return !highOrder || ((vertexBits[node / 64] >> (node % 64)) & 1) != 0;
  }
  [[nodiscard]] std::uint64_t vertexOf(std::uint64_t node) const {
    if (!highOrder) {
      return node;
    }
    const auto below = vertexBits[node / 64] & ((std::uint64_t{1} << (node % 64)) - 1);
    return verticesBefore[node / 64] + static_cast<std::uint64_t>(__builtin_popcountll(below));
  }

  [[nodiscard]] std::uint64_t nodeIndex(std::uint64_t tag) const {
    if (tag < index.firstNodeTag || tag - index.firstNodeTag >= index.numNodes) {
      logError() << "Unknown node tag" << tag;
    }
    return tag - index.firstNodeTag;
  }

  template <typename T> void broadcastVector(std::vector<T>& values) {
    std::uint64_t count = values.size();
    MPI_Bcast(&count, 1, MPI_UINT64_T, 0, comm);
    values.resize(count);
    MPI_Bcast(values.data(), static_cast<int>(count * sizeof(T)), MPI_BYTE, 0, comm);
  }

  void readIndex(const char* meshFile) {
    if (rank == 0) {
      Msh4Indexer indexer;
      if (!indexer.parseFile(meshFile)) {
        logError() << meshFile << std::endl << indexer.getErrorMessage();
      }
      index = indexer.getIndex();
      local.physicalNames = indexer.getPhysicalNames();
    }
    std::array<std::uint64_t, 6> scalars{
        index.dataSize, index.swapBytes ? 1U : 0U,        index.firstNodeTag,
        index.numNodes, index.periodic.empty() ? 0U : 1U, index.highOrder ? 1U : 0U};
    MPI_Bcast(scalars.data(), scalars.size(), MPI_UINT64_T, 0, comm);
    index.dataSize = scalars[0];
    index.swapBytes = scalars[1] != 0;
    index.firstNodeTag = scalars[2];
    index.numNodes = scalars[3];
    hasIdentify = scalars[4] != 0;
    highOrder = scalars[5] != 0;
    broadcastPhysicalNames(local.physicalNames, 0, comm);
    broadcastVector(index.nodeBlocks);
    broadcastVector(index.cellBlocks);
    broadcastVector(index.facetBlocks);
    numCells = 0;
    for (const auto& block : index.cellBlocks) {
      numCells += block.count;
    }
  }

  void readNodes(distributed::RawMshFile& file) {
    const std::uint64_t first = getChunksum(index.numNodes, rank, size);
    const std::uint64_t last = first + getChunksize(index.numNodes, rank, size);

    std::vector<std::vector<NodeRecord>> outgoing(size);
    std::vector<std::uint64_t> tags;
    std::vector<double> coordinates;
    distributed::forBlockRanges(
        index.nodeBlocks, first, last,
        [&](const Msh4NodeBlock& block, std::uint64_t begin, std::uint64_t count) {
          for (std::uint64_t done = 0; done < count; done += distributed::Piece) {
            const auto n = std::min(distributed::Piece, count - done);
            tags.resize(n);
            coordinates.resize(n * block.valuesPerNode);
            file.readSizes(block.tagsOffset + (begin + done) * index.dataSize, n, tags.data());
            file.readDoubles(block.coordinatesOffset +
                                 (begin + done) * block.valuesPerNode * sizeof(double),
                             coordinates.size(), coordinates.data());
            for (std::uint64_t i = 0; i < n; ++i) {
              const auto vertex = nodeIndex(tags[i]);
              const double* x = &coordinates[i * block.valuesPerNode];
              outgoing[nodeOwner(vertex)].push_back({vertex, {x[0], x[1], x[2]}});
            }
          }
        });

    const auto incoming = distributed::exchange(outgoing, comm);
    std::vector<std::vector<NodeRecord>>().swap(outgoing);

    const auto localFirst = localFirstNode();
    nodeCoordinates.assign(localNodeCount() * 3, 0.0);
    std::vector<bool> defined(localNodeCount(), false);
    std::uint64_t numDefined = 0;
    for (const auto& records : incoming) {
      for (const auto& record : records) {
        const auto localVertex = record.index - localFirst;
        if (defined[localVertex]) {
          logError() << "Duplicate node tag" << record.index + index.firstNodeTag;
        }
        defined[localVertex] = true;
        ++numDefined;
        std::copy_n(record.x, 3, &nodeCoordinates[localVertex * 3]);
      }
    }
    if (numDefined != localNodeCount()) {
      logError() << "Missing nodes: the node tags are not unique";
    }
  }

  void readCells(distributed::RawMshFile& file) {
    const std::uint64_t first = getChunksum(numCells, rank, size);
    const std::uint64_t count = getChunksize(numCells, rank, size);
    local.cellTypes.reserve(count);
    local.groups.reserve(count);
    if (highOrder) {
      local.orders.reserve(count);
    }
    std::uint64_t vertices = 0;
    std::uint64_t otherNodes = 0;
    distributed::forBlockRanges(
        index.cellBlocks, first, first + count,
        [&](const Msh4ElementBlock& block, std::uint64_t /*begin*/, std::uint64_t blockCount) {
          const auto vertexCount = shapeOf(gmshElementType(block.type).cellType).vertexCount;
          vertices += blockCount * vertexCount;
          otherNodes += blockCount * (block.nodesPerElement - vertexCount);
        });
    local.connectivity.reserve(vertices);
    highOrderNodes.reserve(otherNodes);

    std::vector<std::uint64_t> values;
    distributed::forBlockRanges(
        index.cellBlocks, first, first + count,
        [&](const Msh4ElementBlock& block, std::uint64_t begin, std::uint64_t blockCount) {
          const auto element = gmshElementType(block.type);
          const auto vertexCount = shapeOf(element.cellType).vertexCount;
          const std::uint64_t valuesPerCell = 1 + block.nodesPerElement;
          for (std::uint64_t done = 0; done < blockCount; done += distributed::Piece) {
            const auto n = std::min(distributed::Piece, blockCount - done);
            values.resize(n * valuesPerCell);
            file.readSizes(block.offset + (begin + done) * valuesPerCell * index.dataSize,
                           values.size(), values.data());
            for (std::uint64_t e = 0; e < n; ++e) {
              // the vertices of a cell come first, the other nodes follow
              for (std::size_t k = 0; k < block.nodesPerElement; ++k) {
                const auto node = nodeIndex(values[e * valuesPerCell + 1 + k]);
                (k < vertexCount ? local.connectivity : highOrderNodes).push_back(node);
              }
              local.cellTypes.push_back(element.cellType);
              local.groups.push_back(static_cast<int>(block.physical));
              if (highOrder) {
                local.orders.push_back(static_cast<std::uint8_t>(element.order));
              }
            }
          }
        });
  }

  /**
   * Marks the nodes which are a vertex of a cell, in a bit set of all nodes, which all ranks share.
   */
  void findVertices() {
    vertexBits.assign((index.numNodes + 63) / 64, 0);
    for (const auto node : local.connectivity) {
      vertexBits[node / 64] |= std::uint64_t{1} << (node % 64);
    }
    MPI_Allreduce(MPI_IN_PLACE, vertexBits.data(), static_cast<int>(vertexBits.size()),
                  MPI_UINT64_T, MPI_BOR, comm);
    verticesBefore.resize(vertexBits.size());
    numVertices = 0;
    for (std::size_t word = 0; word < vertexBits.size(); ++word) {
      verticesBefore[word] = numVertices;
      numVertices += static_cast<std::uint64_t>(__builtin_popcountll(vertexBits[word]));
    }
  }

  /**
   * Fetches the coordinates of the other nodes of the local cells from the ranks which read them.
   */
  void gatherHighOrderGeometry() {
    std::vector<std::uint64_t> wanted(highOrderNodes);
    std::sort(wanted.begin(), wanted.end());
    wanted.erase(std::unique(wanted.begin(), wanted.end()), wanted.end());
    std::vector<std::vector<std::uint64_t>> requests(size);
    for (const auto node : wanted) {
      requests[nodeOwner(node)].push_back(node);
    }

    const auto firstNode = localFirstNode();
    std::vector<std::vector<NodeRecord>> replies(size);
    const auto received = distributed::exchange(requests, comm);
    for (int r = 0; r < size; ++r) {
      for (const auto node : received[r]) {
        const double* x = &nodeCoordinates[(node - firstNode) * 3];
        replies[r].push_back({node, {x[0], x[1], x[2]}});
      }
    }

    // the replies come in the order of the requests, which is the order of the wanted nodes
    std::vector<double> coordinates(wanted.size() * 3);
    std::size_t position = 0;
    for (const auto& records : distributed::exchange(replies, comm)) {
      for (const auto& record : records) {
        if (record.index != wanted[position]) {
          logError() << "Received the coordinates of node" << record.index << "instead of"
                     << wanted[position];
        }
        std::copy_n(record.x, 3, &coordinates[position * 3]);
        ++position;
      }
    }

    local.highOrderGeometry.reserve(highOrderNodes.size() * 3);
    for (const auto node : highOrderNodes) {
      const auto i = static_cast<std::size_t>(std::lower_bound(wanted.begin(), wanted.end(), node) -
                                              wanted.begin());
      local.highOrderGeometry.insert(local.highOrderGeometry.end(), &coordinates[i * 3],
                                     &coordinates[i * 3] + 3);
    }
    std::vector<std::uint64_t>().swap(highOrderNodes);
  }

  /**
   * Sends the coordinates of the vertices from the ranks which read their nodes to the ranks which
   * own them.
   */
  void distributeVertices() {
    std::vector<std::vector<NodeRecord>> outgoing(size);
    const auto firstNode = localFirstNode();
    for (std::uint64_t i = 0; i < localNodeCount(); ++i) {
      const auto node = firstNode + i;
      if (isVertex(node)) {
        const auto vertex = vertexOf(node);
        const double* x = &nodeCoordinates[i * 3];
        outgoing[vertexOwner(vertex)].push_back({vertex, {x[0], x[1], x[2]}});
      }
    }
    const auto firstVertex = localFirstVertex();
    local.geometry.assign(localVertexCount() * 3, 0.0);
    for (const auto& records : distributed::exchange(outgoing, comm)) {
      for (const auto& record : records) {
        std::copy_n(record.x, 3, &local.geometry[(record.index - firstVertex) * 3]);
      }
    }
  }

  void matchBoundaryConditions(distributed::RawMshFile& file) {
    std::uint64_t numFacets = 0;
    for (const auto& block : index.facetBlocks) {
      numFacets += block.count;
    }
    const std::uint64_t first = getChunksum(numFacets, rank, size);
    const std::uint64_t last = first + getChunksize(numFacets, rank, size);

    // the boundary faces go to the ranks given by their hash; their vertices are marked in a bit
    // set of all vertices, which all ranks share
    std::vector<std::vector<FaceRecord>> outgoingFaces(size);
    std::vector<std::uint64_t> onBoundary((index.numNodes + 63) / 64, 0);
    std::vector<std::uint64_t> values;
    distributed::forBlockRanges(
        index.facetBlocks, first, last,
        [&](const Msh4ElementBlock& block, std::uint64_t begin, std::uint64_t count) {
          const auto faceVertices = gmshElementType(block.type).faceVertices;
          const std::uint64_t valuesPerFacet = 1 + block.nodesPerElement;
          for (std::uint64_t done = 0; done < count; done += distributed::Piece) {
            const auto n = std::min(distributed::Piece, count - done);
            values.resize(n * valuesPerFacet);
            file.readSizes(block.offset + (begin + done) * valuesPerFacet * index.dataSize,
                           values.size(), values.data());
            for (std::uint64_t e = 0; e < n; ++e) {
              Face face;
              face.fill(NoVertex);
              for (std::size_t k = 0; k < faceVertices; ++k) {
                face[k] = nodeIndex(values[e * valuesPerFacet + 1 + k]);
                onBoundary[face[k] / 64] |= std::uint64_t{1} << (face[k] % 64);
              }
              std::sort(face.begin(), face.end());
              outgoingFaces[distributed::hashFace(face) % size].push_back({face, block.physical});
            }
          }
        });

    std::vector<FaceRecord> faceTable;
    for (const auto& faces : distributed::exchange(outgoingFaces, comm)) {
      faceTable.insert(faceTable.end(), faces.begin(), faces.end());
    }
    std::vector<std::vector<FaceRecord>>().swap(outgoingFaces);
    const auto byVertices = [](const FaceRecord& a, const FaceRecord& b) {
      return a.vertices < b.vertices;
    };
    std::sort(faceTable.begin(), faceTable.end(), byVertices);

    MPI_Allreduce(MPI_IN_PLACE, onBoundary.data(), static_cast<int>(onBoundary.size()),
                  MPI_UINT64_T, MPI_BOR, comm);
    const auto vertexOnBoundary = [&](std::uint64_t vertex) {
      return ((onBoundary[vertex / 64] >> (vertex % 64)) & 1) != 0;
    };

    const std::uint64_t numLocalCells = local.cellTypes.size();
    // only faces with all vertices on the surface mesh can be part of it
    std::vector<std::vector<FaceQuery>> queries(size);
    std::uint64_t firstVertex = 0;
    std::uint64_t firstFace = 0;
    for (std::uint64_t cell = 0; cell < numLocalCells; ++cell) {
      const auto& shape = shapeOf(local.cellTypes[cell]);
      const auto* vertices = &local.connectivity[firstVertex];
      const auto cellFaces = firstFace;
      firstVertex += shape.vertexCount;
      firstFace += shape.faceCount;
      std::array<bool, MaxCellVertices> flagged{};
      std::size_t numFlagged = 0;
      for (std::size_t k = 0; k < shape.vertexCount; ++k) {
        flagged[k] = vertexOnBoundary(vertices[k]);
        numFlagged += flagged[k] ? 1 : 0;
      }
      if (numFlagged < 3) {
        continue;
      }
      for (std::size_t f = 0; f < shape.faceCount; ++f) {
        Face face;
        face.fill(NoVertex);
        bool candidate = true;
        for (std::size_t k = 0; k < shape.faceVertexCount[f]; ++k) {
          candidate = candidate && flagged[shape.faceVertices[f][k]];
          face[k] = vertices[shape.faceVertices[f][k]];
        }
        if (candidate) {
          std::sort(face.begin(), face.end());
          queries[distributed::hashFace(face) % size].push_back({face, cellFaces + f});
        }
      }
    }

    const auto received = distributed::exchange(queries, comm);
    std::vector<std::vector<FaceQuery>>().swap(queries);
    std::vector<std::vector<FaceReply>> replies(size);
    for (int r = 0; r < size; ++r) {
      for (const auto& query : received[r]) {
        const FaceRecord key{query.vertices, 0};
        const auto [lower, upper] =
            std::equal_range(faceTable.begin(), faceTable.end(), key, byVertices);
        if (upper - lower > 1) {
          logError() << "A face of an element exists multiple times in the surface mesh.";
        }
        replies[r].push_back({query.localFace, lower != upper ? lower->bc : NoMatch});
      }
    }

    local.boundaries.assign(firstFace, 0);
    for (const auto& answersFrom : distributed::exchange(replies, comm)) {
      for (const auto& reply : answersFrom) {
        if (reply.bc != NoMatch) {
          local.boundaries[reply.localFace] = static_cast<int>(boundaryConditionOf(reply.bc));
        }
      }
    }
  }

  void identifyPeriodicNodes() {
    if (!hasIdentify) {
      return;
    }
    // rank 0 joins the periodic nodes into classes, represented by their smallest vertex; nodes
    // which are no vertices drop out
    std::vector<std::vector<Identification>> outgoing(size);
    if (rank == 0) {
      std::unordered_map<std::uint64_t, std::uint64_t> parent;
      const auto find = [&](std::uint64_t node) {
        auto it = parent.try_emplace(node, node).first;
        while (it->second != it->first) {
          const auto next = parent.find(it->second);
          it->second = next->second;
          it = next;
        }
        return it->first;
      };
      for (const auto& [nodeTag, masterTag] : index.periodic) {
        const auto a = find(nodeIndex(nodeTag));
        const auto b = find(nodeIndex(masterTag));
        parent[std::max(a, b)] = std::min(a, b);
      }
      std::vector<std::uint64_t> nodes;
      nodes.reserve(parent.size());
      for (const auto& entry : parent) {
        if (isVertex(entry.first)) {
          nodes.push_back(entry.first);
        }
      }
      std::unordered_map<std::uint64_t, std::uint64_t> smallest;
      for (const auto node : nodes) {
        const auto [it, inserted] = smallest.try_emplace(find(node), node);
        if (!inserted) {
          it->second = std::min(it->second, node);
        }
      }
      for (const auto node : nodes) {
        const auto vertex = vertexOf(node);
        outgoing[vertexOwner(vertex)].push_back({vertex, vertexOf(smallest[find(node)])});
      }
    }

    const auto localFirst = localFirstVertex();
    local.identify.resize(localVertexCount());
    std::iota(local.identify.begin(), local.identify.end(), localFirst);
    for (const auto& identifications : distributed::exchange(outgoing, comm)) {
      for (const auto& identification : identifications) {
        local.identify[identification.vertex - localFirst] = identification.representative;
      }
    }
  }

  MPI_Comm comm;
  int rank = 0;
  int size = 1;
  Msh4Index index;
  std::uint64_t numCells = 0;
  bool hasIdentify = false;
  bool highOrder = false;
  std::uint64_t numVertices = 0;
  /** The coordinates of the chunk of the nodes this rank read */
  std::vector<double> nodeCoordinates;
  /** The nodes of the local cells but their vertices, one cell after the other */
  std::vector<std::uint64_t> highOrderNodes;
  /** Which nodes are vertices, and how many vertices come before the nodes of each word */
  std::vector<std::uint64_t> vertexBits;
  std::vector<std::uint64_t> verticesBefore;
  LocalMesh local;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_DISTRIBUTEDGMSHREADER_H_
