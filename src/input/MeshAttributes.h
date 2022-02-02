#ifndef MESHATTRIBUTES_H
#define MESHATTRIBUTES_H

#include <cassert>
#include <list>
#include <string>

#include "split.h"
#include "tinyxml2/tinyxml2.h"
#include "utils/logger.h"

// for PMU_rank()
#include <SimPartitionedMesh.h>

struct Cube {
  double CubeMSize = 0, CubeCenter[3], CubeWidth[3], CubeHeight[3], CubeDepth[3];
};

enum class ElementType { vertex = 0, edge = 1, face = 2, region = 3 };

using namespace tinyxml2;
class MeshAttributes {
  public:
  enum class SmoothingType { Laplacian = 0, Gradient = 1 };
  int surfaceSmoothingLevel = 2;
  SmoothingType surfaceSmoothingType = SmoothingType::Gradient;
  double surfaceFaceRotationLimit = 5.0;
  int surfaceSnap = 0;
  int volumeSmoothingLevel = 1;
  SmoothingType volumeSmoothingType = SmoothingType::Gradient;
  int VolumeMesherOptimization = 1;
  double globalMSize = 0., faultMSize = 0., gradation = 0., vol_AspectRatio = 0.,
         area_AspectRatio = 0.;
  bool useVelocityAwareMeshing = false;
  double targetedFrequency = 0., elementsPerWaveLength = 0;
  std::string easiFileName;
  const std::list<std::pair<std::list<int>, double>>& getMSizeList(ElementType type);

  std::list<std::pair<std::list<int>, double>> lpair_lVertexId_MSize;
  std::list<std::pair<std::list<int>, double>> lpair_lEdgeId_MSize;
  std::list<std::pair<std::list<int>, double>> lpair_lFaceId_MSize;
  std::list<std::pair<std::list<int>, double>> lpair_lRegionId_MSize;

  std::list<Cube> lCube;
  std::list<int> lFaceIdMeshSizePropagation;
  std::list<int> lFaceIdUseDiscreteMesh;
  std::list<int> lFaceIdNoMesh;
  std::list<int> lRegionIdNoMesh;
  int UseDiscreteMesh_noModification;
  double MeshSizePropagationScalingFactor, MeshSizePropagationDistance = 0.0;
  MeshAttributes();
  MeshAttributes(const char*);
  void init(const char*);

  private:
  XMLDocument doc;
  void readXmlFile(const char* xmlFilename);
  void set_MeshRefinementZoneCube();
  void set_SurfaceMeshingAttributes();
  void set_VolumeMeshingAttributes();
  void set_UseDiscreteMesh();
  void set_NoMesh();
  void set_vertexMSize();
  void set_edgeMSize();
  void set_surfaceMSize();
  void set_MeshSizePropagation();
  void set_regionMSize();
  void set_global_mesh_attributes();
  void set_velocity_aware_meshing();
  void update_attribute_from_xml(XMLNode& element, const char* sElementName,
                                 const char* sAttributeName, int& MeshAttributesMember);
  void update_attribute_from_xml(XMLNode& element, const char* sElementName,
                                 const char* sAttributeName, double& MeshAttributesMember);
  std::string retrieve_text_from_xml(XMLNode* element, const char* sElementName);
  std::list<int> fill_list_using_parsed_string(std::string line);
};

#endif
