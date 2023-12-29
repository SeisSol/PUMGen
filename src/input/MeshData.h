
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

  virtual std::size_t cellCount() = 0;
  virtual std::size_t vertexCount() = 0;

  virtual const std::vector<uint64_t>& connectivity() = 0;
  virtual const std::vector<double>& geometry() = 0;
  virtual const std::vector<int>& group() = 0;
  virtual const std::vector<int64_t>& boundary() = 0;
};

class FullStorageMeshData : public MeshData {
  public:
  virtual ~FullStorageMeshData() = default;

  virtual std::size_t cellCount() { return cellCountValue; }
  virtual std::size_t vertexCount() { return vertexCountValue; }

  virtual const std::vector<uint64_t>& connectivity() { return connectivityData; }
  virtual const std::vector<double>& geometry() { return geometryData; }
  virtual const std::vector<int>& group() { return groupData; }
  virtual const std::vector<int64_t>& boundary() { return boundaryData; }

  protected:
  std::size_t cellCountValue = 0;
  std::size_t vertexCountValue = 0;

  std::vector<uint64_t> connectivityData;
  std::vector<double> geometryData;
  std::vector<int> groupData;
  std::vector<int64_t> boundaryData;

  int bndShift = 8;

  void setBoundary(std::size_t cell, int face, int value) {
    constexpr auto i32limit = std::numeric_limits<unsigned char>::max();
    constexpr auto i64limit = std::numeric_limits<unsigned short>::max();
    if (value < 0 || (value > i32limit && bndShift == 8) || (value > i64limit && bndShift == 16)) {
      logError() << "Cannot handle boundary condition" << value;
    }

    boundaryData[cell] |= static_cast<int64_t>(value) << (face * bndShift);
  }

  void setup(std::size_t cellCount, std::size_t vertexCount) {
    cellCountValue = cellCount;
    vertexCountValue = vertexCount;

    connectivityData.resize(cellCount * 4);
    geometryData.resize(vertexCount * 3);
    groupData.resize(cellCount);
    boundaryData.resize(cellCount);

    // TODO: reconsider the times 4 or 3
  }
};

#endif // MESH_INPUT_H
