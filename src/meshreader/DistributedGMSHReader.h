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
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "GMSHBuilder.h"
#include "Msh4Index.h"
#include "helper/Distributor.h"
#include "helper/MPIConvenience.h"
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

std::uint64_t hashFace(const std::array<std::uint64_t, 3>& face);

// items are read in pieces of at most this many
constexpr std::uint64_t Piece = std::uint64_t{1} << 20;

} // namespace distributed

/**
 * Reads a binary MSH 4.1 file with all ranks: rank 0 only locates the data, every rank reads a
 * share of the nodes, cells and boundary faces, and the parts are sent to the ranks which own
 * them in the chunk distribution of getChunksize. The boundary conditions are matched in a
 * distributed table of the boundary faces; only faces whose vertices all lie on the surface mesh
 * are looked up.
 */
template <std::size_t OrderP> class DistributedGMSHReader {
  public:
  constexpr static std::size_t Order = OrderP;
  constexpr static std::size_t Dim = 3;
  constexpr static bool SupportsIdentify = true;
  constexpr static std::size_t NodesPerCell = nodeCount(Dim, Order);
  constexpr static std::size_t NodesPerFacet = nodeCount(Dim - 1, Order);
  constexpr static std::size_t Facet2Nodes[Dim + 1][Dim] = {
      {1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
  constexpr static long BoundaryConditionOffset = 100;

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
    identifyPeriodicNodes();
  }

  [[nodiscard]] std::size_t nVertices() const { return index.numNodes; }
  [[nodiscard]] std::size_t nElements() const { return numCells; }

  // Each read function hands over one part of the local mesh and releases it, so it can be called
  // only once.

  void readVertices(double* vertices) { handOver(geometry, vertices); }
  void readElements(std::size_t* elements) { handOver(connectivity, elements); }
  void readGroups(int* groupsOut) { handOver(groups, groupsOut); }
  void readBoundaries(int* boundariesOut) { handOver(boundaries, boundariesOut); }
  [[nodiscard]] bool hasIdentify() const { return hasIdentify_; }
  void readIdentify(std::size_t* vertices) { handOver(identify, vertices); }

  private:
  struct NodeRecord {
    std::uint64_t index;
    double x[3];
  };
  struct FaceRecord {
    std::array<std::uint64_t, 3> vertices;
    std::int64_t bc;
  };
  struct FaceQuery {
    std::array<std::uint64_t, 3> vertices;
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

  template <typename T, typename U> static void handOver(std::vector<T>& values, U* out) {
    std::copy(values.begin(), values.end(), out);
    std::vector<T>().swap(values);
  }

  [[nodiscard]] std::uint64_t localFirstVertex() const {
    return getChunksum(index.numNodes, rank, size);
  }
  [[nodiscard]] std::uint64_t localVertexCount() const {
    return getChunksize(index.numNodes, rank, size);
  }
  [[nodiscard]] int vertexOwner(std::uint64_t vertex) const {
    return getChunkOwner(index.numNodes, vertex, size);
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
      Msh4Indexer indexer(GMSHSimplexType<Dim, Order>::type, GMSHSimplexType<Dim - 1, Order>::type);
      if (!indexer.parseFile(meshFile)) {
        logError() << meshFile << std::endl << indexer.getErrorMessage();
      }
      index = indexer.getIndex();
    }
    std::array<std::uint64_t, 5> scalars{index.dataSize, index.swapBytes ? 1U : 0U,
                                         index.firstNodeTag, index.numNodes,
                                         index.periodic.empty() ? 0U : 1U};
    MPI_Bcast(scalars.data(), scalars.size(), MPI_UINT64_T, 0, comm);
    index.dataSize = scalars[0];
    index.swapBytes = scalars[1] != 0;
    index.firstNodeTag = scalars[2];
    index.numNodes = scalars[3];
    hasIdentify_ = scalars[4] != 0;
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
              outgoing[vertexOwner(vertex)].push_back({vertex, {x[0], x[1], x[2]}});
            }
          }
        });

    const auto incoming = distributed::exchange(outgoing, comm);
    std::vector<std::vector<NodeRecord>>().swap(outgoing);

    const auto localFirst = localFirstVertex();
    geometry.assign(localVertexCount() * Dim, 0.0);
    std::vector<bool> defined(localVertexCount(), false);
    std::uint64_t numDefined = 0;
    for (const auto& records : incoming) {
      for (const auto& record : records) {
        const auto local = record.index - localFirst;
        if (defined[local]) {
          logError() << "Duplicate node tag" << record.index + index.firstNodeTag;
        }
        defined[local] = true;
        ++numDefined;
        std::copy_n(record.x, Dim, &geometry[local * Dim]);
      }
    }
    if (numDefined != localVertexCount()) {
      logError() << "Missing nodes: the node tags are not unique";
    }
  }

  void readCells(distributed::RawMshFile& file) {
    const std::uint64_t first = getChunksum(numCells, rank, size);
    const std::uint64_t count = getChunksize(numCells, rank, size);
    connectivity.resize(count * NodesPerCell);
    groups.resize(count);

    constexpr std::uint64_t ValuesPerCell = 1 + NodesPerCell;
    std::vector<std::uint64_t> values;
    std::uint64_t cell = 0;
    distributed::forBlockRanges(
        index.cellBlocks, first, first + count,
        [&](const Msh4ElementBlock& block, std::uint64_t begin, std::uint64_t blockCount) {
          for (std::uint64_t done = 0; done < blockCount; done += distributed::Piece) {
            const auto n = std::min(distributed::Piece, blockCount - done);
            values.resize(n * ValuesPerCell);
            file.readSizes(block.offset + (begin + done) * ValuesPerCell * index.dataSize,
                           values.size(), values.data());
            for (std::uint64_t e = 0; e < n; ++e, ++cell) {
              for (std::size_t k = 0; k < NodesPerCell; ++k) {
                connectivity[cell * NodesPerCell + k] =
                    nodeIndex(values[e * ValuesPerCell + 1 + k]);
              }
              groups[cell] = static_cast<int>(block.physical);
            }
          }
        });
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
    constexpr std::uint64_t ValuesPerFacet = 1 + NodesPerFacet;
    std::vector<std::uint64_t> values;
    distributed::forBlockRanges(
        index.facetBlocks, first, last,
        [&](const Msh4ElementBlock& block, std::uint64_t begin, std::uint64_t count) {
          for (std::uint64_t done = 0; done < count; done += distributed::Piece) {
            const auto n = std::min(distributed::Piece, count - done);
            values.resize(n * ValuesPerFacet);
            file.readSizes(block.offset + (begin + done) * ValuesPerFacet * index.dataSize,
                           values.size(), values.data());
            for (std::uint64_t e = 0; e < n; ++e) {
              std::array<std::uint64_t, 3> face{};
              for (std::size_t k = 0; k < Dim; ++k) {
                face[k] = nodeIndex(values[e * ValuesPerFacet + 1 + k]);
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

    const std::uint64_t numLocalCells = groups.size();
    // only faces with all vertices on the surface mesh can be part of it
    std::vector<std::vector<FaceQuery>> queries(size);
    for (std::uint64_t cell = 0; cell < numLocalCells; ++cell) {
      std::array<bool, Dim + 1> flagged{};
      std::size_t numFlagged = 0;
      for (std::size_t k = 0; k < Dim + 1; ++k) {
        flagged[k] = vertexOnBoundary(connectivity[cell * NodesPerCell + k]);
        numFlagged += flagged[k] ? 1 : 0;
      }
      if (numFlagged < Dim) {
        continue;
      }
      for (std::size_t localFace = 0; localFace < Dim + 1; ++localFace) {
        std::array<std::uint64_t, 3> face{};
        bool candidate = true;
        for (std::size_t k = 0; k < Dim; ++k) {
          candidate = candidate && flagged[Facet2Nodes[localFace][k]];
          face[k] = connectivity[cell * NodesPerCell + Facet2Nodes[localFace][k]];
        }
        if (candidate) {
          std::sort(face.begin(), face.end());
          queries[distributed::hashFace(face) % size].push_back(
              {face, cell * (Dim + 1) + localFace});
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

    boundaries.assign(numLocalCells * (Dim + 1), 0);
    for (const auto& answersFrom : distributed::exchange(replies, comm)) {
      for (const auto& reply : answersFrom) {
        if (reply.bc != NoMatch) {
          const auto bc =
              reply.bc >= BoundaryConditionOffset ? reply.bc - BoundaryConditionOffset : reply.bc;
          boundaries[reply.localFace] = static_cast<int>(bc);
        }
      }
    }
  }

  void identifyPeriodicNodes() {
    if (!hasIdentify_) {
      return;
    }
    // rank 0 joins the periodic nodes into classes, represented by their smallest vertex
    std::vector<std::vector<Identification>> outgoing(size);
    if (rank == 0) {
      std::unordered_map<std::uint64_t, std::uint64_t> parent;
      const auto find = [&](std::uint64_t vertex) {
        auto it = parent.try_emplace(vertex, vertex).first;
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
      std::vector<std::uint64_t> vertices;
      vertices.reserve(parent.size());
      for (const auto& entry : parent) {
        vertices.push_back(entry.first);
      }
      for (const auto vertex : vertices) {
        outgoing[vertexOwner(vertex)].push_back({vertex, find(vertex)});
      }
    }

    const auto localFirst = localFirstVertex();
    identify.resize(localVertexCount());
    std::iota(identify.begin(), identify.end(), localFirst);
    for (const auto& identifications : distributed::exchange(outgoing, comm)) {
      for (const auto& identification : identifications) {
        identify[identification.vertex - localFirst] = identification.representative;
      }
    }
  }

  MPI_Comm comm;
  int rank = 0;
  int size = 1;
  Msh4Index index;
  std::uint64_t numCells = 0;
  bool hasIdentify_ = false;
  std::vector<double> geometry;
  std::vector<std::size_t> connectivity;
  std::vector<int> groups;
  std::vector<int> boundaries;
  std::vector<std::size_t> identify;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_DISTRIBUTEDGMSHREADER_H_
