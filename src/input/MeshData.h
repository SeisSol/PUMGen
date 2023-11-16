
#ifndef MESH_DATA_H
#define MESH_DATA_H

#include <cstdint>
#include <vector>

/**
 * Interface for mesh input
 */
class MeshData {
  public:
  virtual std::size_t cellCount() = 0;
  virtual std::size_t vertexCount() = 0;

  virtual const std::vector<uint64_t>& connectivity() = 0;
  virtual const std::vector<double>& geometry() = 0;
  virtual const std::vector<int>& group() = 0;
  virtual const std::vector<int64_t>& boundary() = 0;
};

class FullStorageMeshData : public MeshData {
  public:
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
