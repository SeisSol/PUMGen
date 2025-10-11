#ifndef GMSHMESHBUILDER_20200901_H
#define GMSHMESHBUILDER_20200901_H

#include <array>

namespace tndm {

class GMSHMeshBuilder {
  public:
  virtual ~GMSHMeshBuilder() = default;
  virtual void setNumVertices(std::size_t numVertices) = 0;
  virtual void setVertex(long id, std::array<double, 3> const& x) = 0;
  virtual void setNumElements(std::size_t numElements) = 0;
  virtual void addElement(long type, long tag, long* node, std::size_t numNodes) = 0;
  virtual void addVertexLink(std::size_t vertex, std::size_t linkVertex) = 0;
};

} // namespace tndm

#endif // GMSHMESHBUILDER_20200901_H
