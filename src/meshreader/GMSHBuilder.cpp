// SPDX-FileCopyrightText: 2020 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSHBuilder.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>

#include "helper/Distributor.h"
#include "helper/MPIConvenience.h"
#include "utils/logger.h"

namespace puml {

namespace {
/**
 * Boundary conditions in SeisSol used to start at 100, e.g. 101 = free surface. In the hdf5 format
 * one starts counting from 0, e.g. 1 = free surface. In order to be compatible with legacy gmsh
 * scripts, 100 is subtracted from a boundary condition larger than or equal to 100.
 */
int adjustBoundaryCondition(int bc) {
  constexpr int BoundaryConditionOffset = 100;
  return bc >= BoundaryConditionOffset ? bc - BoundaryConditionOffset : bc;
}

template <typename T> void release(std::vector<T>& values) { std::vector<T>().swap(values); }

std::size_t nodesOf(CellType type, std::uint8_t order) { return lagrangeNodeCount(type, order); }
} // namespace

void GMSHBuilder::resizeIdentifyIfNeeded(std::size_t newSize) {
  const auto oldSize = identify.size();
  if (newSize > oldSize) {
    identify.resize(newSize);
    std::iota(identify.begin() + static_cast<std::ptrdiff_t>(oldSize), identify.end(), oldSize);
  }
}

void GMSHBuilder::setNumVertices(std::size_t numVertices) { vertices.resize(numVertices); }

void GMSHBuilder::setVertex(long id, const std::array<double, 3>& x) {
  vertices[static_cast<std::size_t>(id)] = x;
}

void GMSHBuilder::setNumElements(std::size_t numElements) {
  cellTypes.clear();
  cellOrders.clear();
  cellNodes.clear();
  groups.clear();
  facets.clear();
  bcs.clear();

  // all elements may be cells; as tetrahedra, until more nodes arrive
  constexpr std::size_t TetrahedronNodes = 4;
  cellTypes.reserve(numElements);
  cellOrders.reserve(numElements);
  groups.reserve(numElements);
  cellNodes.reserve(TetrahedronNodes * numElements);
}

void GMSHBuilder::addElement(long type, long tag, long* node, std::size_t numNodes) {
  const auto element = gmshElementType(type);
  if (element.kind == GmshElementType::Kind::Cell) {
    if (!element.complete) {
      throw std::runtime_error("Element type " + std::to_string(type) + " is a " +
                               shapeOf(element.cellType).name +
                               " without all nodes of a Lagrange cell, which is not supported");
    }
    cellTypes.push_back(element.cellType);
    cellOrders.push_back(static_cast<std::uint8_t>(element.order));
    cellNodes.insert(cellNodes.end(), node, node + numNodes);
    groups.push_back(static_cast<int>(tag));
  } else if (element.kind == GmshElementType::Kind::Face) {
    Face face;
    face.fill(NoVertex);
    std::copy_n(node, element.faceVertices, face.begin());
    std::sort(face.begin(), face.end());
    facets.push_back(face);
    bcs.push_back(static_cast<int>(tag));
  }
}

void GMSHBuilder::addVertexLink(std::size_t vertex, std::size_t linkVertex) {
  // needed in case we parse periodic before we parse the vertex/node section
  resizeIdentifyIfNeeded(vertex + 1);
  identify[vertex] = linkVertex;
}

void GMSHBuilder::postprocess() {
  // GMSH does not pick the same node to identify for all
  // thus, stratify (i.e. "path compress" in "union find" terms)

  if (!identify.empty()) {
    // only resize if we have an identification array
    resizeIdentifyIfNeeded(vertices.size());
    for (std::size_t i = 0; i < identify.size(); ++i) {
      if (identify[i] != i) {

        // find until we loop
        std::size_t curr = i;
        std::set<std::size_t> equivalent;
        equivalent.insert(curr);
        while (equivalent.find(identify[curr]) == equivalent.end()) {
          curr = identify[curr];
          equivalent.insert(curr);
        }

        // then set everything to one value
        const auto identifyEquiv = *equivalent.begin();
        for (const auto& index : equivalent) {
          identify[index] = identifyEquiv;
        }
      }
    }
  }
}

bool GMSHBuilder::hasHighOrder() const {
  return std::any_of(cellOrders.begin(), cellOrders.end(),
                     [](std::uint8_t order) { return order > 1; });
}

GlobalMesh prepareMesh(GMSHBuilder& builder) {
  GlobalMesh mesh;
  const std::size_t numNodes = builder.vertices.size();
  const std::size_t numCells = builder.cellTypes.size();
  const bool highOrder = builder.hasHighOrder();
  constexpr std::uint64_t NoId = std::numeric_limits<std::uint64_t>::max();

  // the vertex of every node; with cells of higher order, only the vertices of the cells are
  std::vector<std::uint64_t> vertexOf(numNodes);
  if (highOrder) {
    std::fill(vertexOf.begin(), vertexOf.end(), NoId);
    std::size_t first = 0;
    for (std::size_t cell = 0; cell < numCells; ++cell) {
      const auto type = builder.cellTypes[cell];
      for (std::size_t k = 0; k < shapeOf(type).vertexCount; ++k) {
        vertexOf[builder.cellNodes[first + k]] = 0;
      }
      first += nodesOf(type, builder.cellOrders[cell]);
    }
    std::uint64_t next = 0;
    for (auto& vertex : vertexOf) {
      if (vertex != NoId) {
        vertex = next++;
      }
    }
  } else {
    std::iota(vertexOf.begin(), vertexOf.end(), 0);
  }

  mesh.geometry.reserve(3 * numNodes);
  for (std::size_t node = 0; node < numNodes; ++node) {
    if (vertexOf[node] != NoId) {
      mesh.geometry.insert(mesh.geometry.end(), builder.vertices[node].begin(),
                           builder.vertices[node].end());
    }
  }
  if (!highOrder) {
    release(builder.vertices);
  }

  if (!builder.identify.empty()) {
    // a class of identified nodes is represented by its first vertex
    std::vector<std::uint64_t> representative(numNodes, NoId);
    mesh.identify.resize(mesh.geometry.size() / 3);
    for (std::size_t node = 0; node < numNodes; ++node) {
      if (vertexOf[node] != NoId) {
        auto& first = representative[builder.identify[node]];
        if (first == NoId) {
          first = vertexOf[node];
        }
        mesh.identify[vertexOf[node]] = first;
      }
    }
    release(builder.identify);
  }

  // the boundary conditions of the faces of the cells, from the surface mesh
  {
    std::size_t numFaces = 0;
    for (const auto type : builder.cellTypes) {
      numFaces += shapeOf(type).faceCount;
    }
    mesh.boundaries.assign(numFaces, 0);
    std::vector<std::uint8_t> onBoundary(numNodes, 0);
    std::vector<std::size_t> byVertices(builder.facets.size());
    std::iota(byVertices.begin(), byVertices.end(), 0);
    for (const auto& face : builder.facets) {
      for (const auto vertex : face) {
        if (vertex != GMSHBuilder::NoVertex) {
          onBoundary[vertex] = 1;
        }
      }
    }
    std::sort(byVertices.begin(), byVertices.end(),
              [&](std::size_t a, std::size_t b) { return builder.facets[a] < builder.facets[b]; });
    const auto less = [&](std::size_t facet, const GMSHBuilder::Face& face) {
      return builder.facets[facet] < face;
    };

    std::size_t first = 0;
    std::size_t firstFace = 0;
    for (std::size_t cell = 0; cell < numCells; ++cell) {
      const auto& shape = shapeOf(builder.cellTypes[cell]);
      const auto* nodes = &builder.cellNodes[first];
      first += nodesOf(shape.type, builder.cellOrders[cell]);
      const auto cellFaces = firstFace;
      firstFace += shape.faceCount;

      // only faces with all vertices on the surface mesh can be part of it
      std::size_t verticesOnBoundary = 0;
      for (std::size_t k = 0; k < shape.vertexCount; ++k) {
        verticesOnBoundary += onBoundary[nodes[k]];
      }
      if (verticesOnBoundary < 3) {
        continue;
      }
      for (std::size_t f = 0; f < shape.faceCount; ++f) {
        GMSHBuilder::Face face;
        face.fill(GMSHBuilder::NoVertex);
        bool candidate = true;
        for (std::size_t k = 0; k < shape.faceVertexCount[f]; ++k) {
          face[k] = nodes[shape.faceVertices[f][k]];
          candidate = candidate && onBoundary[face[k]] != 0;
        }
        if (!candidate) {
          continue;
        }
        std::sort(face.begin(), face.end());
        const auto match = std::lower_bound(byVertices.begin(), byVertices.end(), face, less);
        if (match == byVertices.end() || builder.facets[*match] != face) {
          continue;
        }
        if (match + 1 != byVertices.end() && builder.facets[*(match + 1)] == face) {
          logError() << "A face of an element exists multiple times in the surface mesh.";
        }
        mesh.boundaries[cellFaces + f] = adjustBoundaryCondition(builder.bcs[*match]);
      }
    }
    release(builder.facets);
    release(builder.bcs);
  }

  // the cells, with their vertices renumbered and their other nodes as higher-order geometry
  if (highOrder) {
    std::size_t first = 0;
    for (std::size_t cell = 0; cell < numCells; ++cell) {
      const auto type = builder.cellTypes[cell];
      const auto count = nodesOf(type, builder.cellOrders[cell]);
      const auto vertexCount = shapeOf(type).vertexCount;
      for (std::size_t k = 0; k < count; ++k) {
        const auto node = builder.cellNodes[first + k];
        if (k < vertexCount) {
          mesh.connectivity.push_back(vertexOf[node]);
        } else {
          const auto& x = builder.vertices[node];
          mesh.highOrderGeometry.insert(mesh.highOrderGeometry.end(), x.begin(), x.end());
        }
      }
      first += count;
    }
  }
  release(builder.vertices);
  if (highOrder) {
    release(builder.cellNodes);
    mesh.orders = std::move(builder.cellOrders);
  } else {
    mesh.connectivity = std::move(builder.cellNodes);
  }
  mesh.cellTypes = std::move(builder.cellTypes);
  mesh.groups = std::move(builder.groups);
  return mesh;
}

namespace {
template <typename T> MPI_Datatype bytesOf() {
  MPI_Datatype type;
  MPI_Type_contiguous(static_cast<int>(sizeof(T)), MPI_BYTE, &type);
  MPI_Type_commit(&type);
  return type;
}

/**
 * Scatters the rows of an array of width values each, as getChunksize splits them, and releases
 * the array on rank 0.
 */
template <typename T>
void scatterRows(std::vector<T>& global, std::vector<T>& local, std::uint64_t rows,
                 std::size_t width, MPI_Comm comm) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  std::vector<std::size_t> counts(size);
  std::vector<std::size_t> displs(size + 1, 0);
  for (int r = 0; r < size; ++r) {
    counts[r] = width * getChunksize(rows, r, size);
    displs[r + 1] = displs[r] + counts[r];
  }
  local.resize(counts[rank]);
  auto type = bytesOf<T>();
  largeScatterv(global.data(), counts.data(), displs.data(), type, local.data(), local.size(), type,
                0, comm);
  MPI_Type_free(&type);
  release(global);
}

/**
 * Scatters an array holding a different number of values per row, whose count per row rank 0
 * knows, and releases it on rank 0.
 */
template <typename T, typename CountF>
void scatterRagged(std::vector<T>& global, std::vector<T>& local, std::uint64_t rows,
                   CountF&& countOf, MPI_Comm comm) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  std::vector<std::size_t> counts(size, 0);
  std::vector<std::size_t> displs(size + 1, 0);
  if (rank == 0) {
    std::size_t row = 0;
    for (int r = 0; r < size; ++r) {
      for (std::size_t i = 0; i < getChunksize(rows, r, size); ++i, ++row) {
        counts[r] += countOf(row);
      }
      displs[r + 1] = displs[r] + counts[r];
    }
  }
  std::uint64_t mine = 0;
  static_assert(sizeof(std::size_t) == sizeof(std::uint64_t));
  MPI_Scatter(counts.data(), 1, MPI_UINT64_T, &mine, 1, MPI_UINT64_T, 0, comm);
  local.resize(mine);
  auto type = bytesOf<T>();
  largeScatterv(global.data(), counts.data(), displs.data(), type, local.data(), local.size(), type,
                0, comm);
  MPI_Type_free(&type);
  release(global);
}
} // namespace

LocalMesh distributeMesh(GlobalMesh& mesh, MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  std::array<std::uint64_t, 4> sizes{mesh.cellTypes.size(), mesh.geometry.size() / 3,
                                     mesh.identify.empty() ? 0U : 1U,
                                     mesh.orders.empty() ? 0U : 1U};
  MPI_Bcast(sizes.data(), static_cast<int>(sizes.size()), MPI_UINT64_T, 0, comm);
  const auto cells = sizes[0];
  const auto vertices = sizes[1];

  LocalMesh local;
  // the counts of the ragged arrays follow from the kinds (and orders) of the cells, which are
  // distributed first
  int size = 1;
  MPI_Comm_size(comm, &size);
  if (size == 1) {
    local.cellTypes = std::move(mesh.cellTypes);
    local.connectivity = std::move(mesh.connectivity);
    local.geometry = std::move(mesh.geometry);
    local.groups = std::move(mesh.groups);
    local.boundaries = std::move(mesh.boundaries);
    local.identify = std::move(mesh.identify);
    local.orders = std::move(mesh.orders);
    local.highOrderGeometry = std::move(mesh.highOrderGeometry);
    return local;
  }

  std::vector<std::size_t> vertexCounts;
  std::vector<std::size_t> faceCounts;
  std::vector<std::size_t> highOrderCounts;
  if (rank == 0) {
    vertexCounts.resize(cells);
    faceCounts.resize(cells);
    highOrderCounts.resize(sizes[3] != 0 ? cells : 0);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto type = mesh.cellTypes[cell];
      vertexCounts[cell] = shapeOf(type).vertexCount;
      faceCounts[cell] = shapeOf(type).faceCount;
      if (sizes[3] != 0) {
        highOrderCounts[cell] = 3 * (nodesOf(type, mesh.orders[cell]) - vertexCounts[cell]);
      }
    }
  }

  scatterRows(mesh.cellTypes, local.cellTypes, cells, 1, comm);
  scatterRagged(
      mesh.connectivity, local.connectivity, cells,
      [&](std::size_t cell) { return vertexCounts[cell]; }, comm);
  scatterRows(mesh.groups, local.groups, cells, 1, comm);
  scatterRagged(
      mesh.boundaries, local.boundaries, cells, [&](std::size_t cell) { return faceCounts[cell]; },
      comm);
  if (sizes[3] != 0) {
    scatterRows(mesh.orders, local.orders, cells, 1, comm);
    scatterRagged(
        mesh.highOrderGeometry, local.highOrderGeometry, cells,
        [&](std::size_t cell) { return highOrderCounts[cell]; }, comm);
  }
  scatterRows(mesh.geometry, local.geometry, vertices, 3, comm);
  if (sizes[2] != 0) {
    scatterRows(mesh.identify, local.identify, vertices, 1, comm);
  }
  return local;
}

} // namespace puml
