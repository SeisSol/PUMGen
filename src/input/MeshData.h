
#ifndef MESH_DATA_H
#define MESH_DATA_H

#include <cstdint>
#include <limits>
#include <vector>

#include "utils/logger.h"

// TODO: to save space, the connectivity, geometry, group, and boundary arrays can also be
// constructed on the fly (for that, change the MeshData data structure to work on an iterator
// instead)

/**
 * Interface for mesh input
 */
class MeshData {
  public:
  virtual ~MeshData() = default;

  virtual std::size_t cellCount() const = 0;
  virtual std::size_t vertexCount() const = 0;

  virtual const std::vector<uint64_t>& connectivity() const = 0;
  virtual const std::vector<double>& geometry() const = 0;
  virtual const std::vector<int>& group() const = 0;
  virtual const std::vector<int64_t>& boundary() const = 0;
  virtual const std::vector<uint64_t>& identify() const = 0;
  virtual bool hasIdentify() const = 0;

  // TODO: generalize?
  virtual std::size_t vertexSize() const { return 3; }
  virtual std::size_t cellSize() const { return 4; }
};

class FullStorageMeshData : public MeshData {
  public:
  ~FullStorageMeshData() override = default;

  std::size_t cellCount() const override { return cellCountValue; }
  std::size_t vertexCount() const override { return vertexCountValue; }

  const std::vector<uint64_t>& connectivity() const override { return connectivityData; }
  const std::vector<double>& geometry() const override { return geometryData; }
  const std::vector<int>& group() const override { return groupData; }
  const std::vector<int64_t>& boundary() const override { return boundaryData; }
  const std::vector<uint64_t>& identify() const override { return identifyData; }

  bool hasIdentify() const override { return !identifyData.empty(); }

  protected:
  FullStorageMeshData(int boundarySize) : bndShift(boundarySize) {}

  std::size_t cellCountValue = 0;
  std::size_t vertexCountValue = 0;

  std::vector<uint64_t> connectivityData;
  std::vector<double> geometryData;
  std::vector<int> groupData;
  std::vector<int64_t> boundaryData;
  std::vector<uint64_t> identifyData;

  int bndShift = 8;

  void setBoundary(std::size_t cell, int face, int64_t value) {
    if (bndShift < 0) {
      boundaryData[cell * (vertexSize() + 1) + face] = value;
    } else {
      constexpr int64_t i32limit = std::numeric_limits<uint8_t>::max();
      constexpr int64_t i64limit = std::numeric_limits<uint16_t>::max();
      if (value < 0 || (value > i32limit && bndShift == 8) ||
          (value > i64limit && bndShift == 16)) {
        logError() << "Cannot handle boundary condition" << value;
      }

      boundaryData[cell] |= static_cast<int64_t>(value) << (face * bndShift);
    }
  }

  void setup(std::size_t cellCount, std::size_t vertexCount, bool identify = false) {
    cellCountValue = cellCount;
    vertexCountValue = vertexCount;

    connectivityData.resize(cellCount * cellSize());
    geometryData.resize(vertexCount * vertexSize());
    groupData.resize(cellCount);
    if (bndShift < 0) {
      boundaryData.resize(cellCount * (vertexSize() + 1));
    } else {
      boundaryData.resize(cellCount);
    }

    // TODO: reconsider the times 4 or 3

    if (identify) {
      identifyData.resize(vertexCount);
    }
  }
};

#endif // MESH_INPUT_H
