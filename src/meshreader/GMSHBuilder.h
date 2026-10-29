// SPDX-FileCopyrightText: 2020 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_
#define PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "GMSHMeshBuilder.h"
#include "input/MeshData.h"
#include "mesh/CellType.h"

namespace puml {

/**
 * Collects what the parsers read: all nodes, the cells of any kind and order, and the faces of the
 * surface mesh with their boundary conditions.
 */
class GMSHBuilder : public GMSHMeshBuilder {
  public:
  /** Fills up the faces of three vertices */
  static constexpr std::uint64_t NoVertex = std::numeric_limits<std::uint64_t>::max();
  using Face = std::array<std::uint64_t, MaxFaceVertices>;

  std::vector<std::array<double, 3>> vertices;
  std::vector<std::size_t> identify;

  std::vector<CellType> cellTypes;
  std::vector<std::uint8_t> cellOrders;
  /** The nodes of all cells, one cell after the other; a cell has as many as its kind and order */
  std::vector<std::uint64_t> cellNodes;
  std::vector<int> groups;

  /** The vertices of the faces of the surface mesh, sorted */
  std::vector<Face> facets;
  std::vector<int> bcs;

  std::vector<PhysicalName> physicalNames;

  void setNumVertices(std::size_t numVertices) override;
  void setVertex(long id, const std::array<double, 3>& x) override;
  void setNumElements(std::size_t numElements) override;
  void addElement(long type, long tag, long* node, std::size_t numNodes) override;
  void addVertexLink(std::size_t vertex, std::size_t linkVertex) override;
  void addPhysicalName(int dimension, long tag, const std::string& name) override;
  void postprocess() override;

  /** Whether a cell is of an order higher than one */
  [[nodiscard]] bool hasHighOrder() const;

  private:
  void resizeIdentifyIfNeeded(std::size_t newSize);
};

/**
 * A whole mesh, as rank 0 holds it before distributing it.
 */
struct GlobalMesh {
  std::vector<CellType> cellTypes;
  /** The vertices of all cells, one cell after the other */
  std::vector<std::uint64_t> connectivity;
  std::vector<double> geometry;
  std::vector<int> groups;
  /** As many per cell as its kind has faces */
  std::vector<int> boundaries;
  std::vector<std::uint64_t> identify;
  /** Empty if all cells are linear */
  std::vector<std::uint8_t> orders;
  std::vector<double> highOrderGeometry;
  std::vector<PhysicalName> physicalNames;
};

/**
 * Turns what the parser read into the mesh to write, consuming the builder. In a mesh with cells of
 * higher order, the vertices are the nodes that are a vertex of a cell, numbered in the order of
 * the nodes; the other nodes of a cell go to its higher-order geometry. The boundary conditions of
 * the surface mesh are assigned to the faces of the cells.
 */
GlobalMesh prepareMesh(GMSHBuilder& builder);

/**
 * Sends every rank its part of the mesh rank 0 holds: the chunks of getChunksize of the cells and
 * of the vertices. Collective.
 */
LocalMesh distributeMesh(GlobalMesh& mesh, MPI_Comm comm);

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_
