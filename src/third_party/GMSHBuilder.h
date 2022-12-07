#ifndef PUMGEN_GMSHBUILDER_H
#define PUMGEN_GMSHBUILDER_H

#include <array>
namespace tndm {

class GMSHMeshBuilder {
  public:
  virtual ~GMSHMeshBuilder() {}
  virtual void setNumVertices(std::size_t numVertices) = 0;
  virtual void setVertex(long id, std::array<double, 3> const& x) = 0;
  virtual void setNumElements(std::size_t numElements) = 0;
  virtual void addElement(long type, long tag, long* node, std::size_t numNodes) = 0;
};

}


#endif // PUMGEN_GMSHBUILDER_H
