// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PUMGEN_SRC_INPUT_MESHDATA_H_
#define PUMGEN_SRC_INPUT_MESHDATA_H_

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "mesh/CellType.h"

/**
 * The part of a mesh one rank holds: a contiguous range of the cells and a contiguous range of the
 * vertices of the global numbering.
 */
struct LocalMesh {
  /** The kind of every cell */
  std::vector<puml::CellType> cellTypes;
  /** The (global) vertices of all cells, one cell after the other */
  std::vector<std::uint64_t> connectivity;
  /** Three coordinates per vertex */
  std::vector<double> geometry;
  std::vector<int> groups;
  /** The boundary condition of every face, one cell after the other */
  std::vector<int> boundaries;
  /** The vertex every vertex is identified with, for periodic meshes; empty otherwise */
  std::vector<std::uint64_t> identify;
  /** The order of every cell; empty if the mesh has only linear cells */
  std::vector<std::uint8_t> orders;
  /**
   * The coordinates of the nodes of every cell except its vertices, in the order of gmsh, one cell
   * after the other; how many a cell has follows from its kind and order
   */
  std::vector<double> highOrderGeometry;
};

/**
 * Interface for mesh input
 */
class MeshData {
  public:
  virtual ~MeshData() = default;

  [[nodiscard]] virtual std::size_t cellCount() const = 0;
  [[nodiscard]] virtual std::size_t vertexCount() const = 0;

  [[nodiscard]] virtual const std::vector<puml::CellType>& cellTypes() const = 0;
  [[nodiscard]] virtual const std::vector<std::uint64_t>& connectivity() const = 0;
  [[nodiscard]] virtual const std::vector<double>& geometry() const = 0;
  [[nodiscard]] virtual const std::vector<int>& group() const = 0;
  [[nodiscard]] virtual const std::vector<int>& boundary() const = 0;
  [[nodiscard]] virtual const std::vector<std::uint64_t>& identify() const = 0;
  [[nodiscard]] virtual bool hasIdentify() const = 0;
  [[nodiscard]] virtual const std::vector<std::uint8_t>& orders() const = 0;
  [[nodiscard]] virtual const std::vector<double>& highOrderGeometry() const = 0;
};

class FullStorageMeshData : public MeshData {
  public:
  [[nodiscard]] std::size_t cellCount() const override { return cellCountValue; }
  [[nodiscard]] std::size_t vertexCount() const override { return vertexCountValue; }

  [[nodiscard]] const std::vector<puml::CellType>& cellTypes() const override {
    return cellTypeData;
  }
  [[nodiscard]] const std::vector<std::uint64_t>& connectivity() const override {
    return connectivityData;
  }
  [[nodiscard]] const std::vector<double>& geometry() const override { return geometryData; }
  [[nodiscard]] const std::vector<int>& group() const override { return groupData; }
  [[nodiscard]] const std::vector<int>& boundary() const override { return boundaryData; }
  [[nodiscard]] const std::vector<std::uint64_t>& identify() const override { return identifyData; }
  [[nodiscard]] bool hasIdentify() const override { return !identifyData.empty(); }
  [[nodiscard]] const std::vector<std::uint8_t>& orders() const override { return orderData; }
  [[nodiscard]] const std::vector<double>& highOrderGeometry() const override {
    return highOrderGeometryData;
  }

  protected:
  static constexpr std::size_t TetrahedronFaces = 4;

  std::size_t cellCountValue = 0;
  std::size_t vertexCountValue = 0;

  std::vector<puml::CellType> cellTypeData;
  std::vector<std::uint64_t> connectivityData;
  std::vector<double> geometryData;
  std::vector<int> groupData;
  std::vector<int> boundaryData;
  std::vector<std::uint64_t> identifyData;
  std::vector<std::uint8_t> orderData;
  std::vector<double> highOrderGeometryData;

  /**
   * Sets the boundary condition of a face of a mesh of tetrahedra.
   */
  void setBoundary(std::size_t cell, int face, int value) {
    boundaryData[cell * TetrahedronFaces + static_cast<std::size_t>(face)] = value;
  }

  /**
   * Makes room for a mesh of tetrahedra.
   */
  void setup(std::size_t cellCount, std::size_t vertexCount, bool identify = false) {
    cellCountValue = cellCount;
    vertexCountValue = vertexCount;

    cellTypeData.assign(cellCount, puml::CellType::Tetrahedron);
    connectivityData.resize(cellCount * puml::shapeOf(puml::CellType::Tetrahedron).vertexCount);
    geometryData.resize(vertexCount * 3);
    groupData.resize(cellCount);
    boundaryData.assign(cellCount * TetrahedronFaces, 0);

    if (identify) {
      identifyData.resize(vertexCount);
    }
  }

  /**
   * Takes over the part of a mesh a reader distributed to this rank.
   */
  void take(LocalMesh&& mesh) {
    cellCountValue = mesh.cellTypes.size();
    vertexCountValue = mesh.geometry.size() / 3;
    cellTypeData = std::move(mesh.cellTypes);
    connectivityData = std::move(mesh.connectivity);
    geometryData = std::move(mesh.geometry);
    groupData = std::move(mesh.groups);
    boundaryData = std::move(mesh.boundaries);
    identifyData = std::move(mesh.identify);
    orderData = std::move(mesh.orders);
    highOrderGeometryData = std::move(mesh.highOrderGeometry);
  }
};

#endif // PUMGEN_SRC_INPUT_MESHDATA_H_
