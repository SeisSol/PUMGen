/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2018 LMU-TUM 
 * @author Thomas Ulrich
 */
#ifndef MESHATTRIBUTES_H
#define MESHATTRIBUTES_H

#include <cassert>
#include <string>
#include <list>

#include "tinyxml2/tinyxml2.h"
#include "utils/logger.h"
#include "split.h"

//for PMU_rank()
#include <SimPartitionedMesh.h>

struct Cube {
  double CubeMSize=0,CubeCenter[3],CubeWidth[3],CubeHeight[3],CubeDepth[3];
};


using namespace tinyxml2;
class MeshAttributes {
  public:
    enum class SmoothingType {
      Laplacian = 0,
      Gradient = 1
    };
    int surfaceSmoothingLevel=2;
    SmoothingType surfaceSmoothingType = SmoothingType::Gradient;
    double surfaceFaceRotationLimit=5.0;
    int surfaceSnap=0;
    int volumeSmoothingLevel=1;
    SmoothingType volumeSmoothingType = SmoothingType::Gradient;
    int VolumeMesherOptimization=1;
    double globalMSize=0., faultMSize=0., gradation=0., vol_AspectRatio=0., area_AspectRatio=0.;
    std::list<double> lsurfaceMSize;
    std::list <std::list<int>> llsurfaceMSizeFaceId;
    std::list<double> lregionMSize;
    std::list <std::list<int>> llregionMSizeRegionId;
    std::list<Cube> lCube;
    std::list<int> lFaceIdMeshSizePropagation;
    std::list<int> lFaceIdUseDiscreteMesh;
    int UseDiscreteMesh_noModification;
    double MeshSizePropagationScalingFactor,MeshSizePropagationDistance=0.0;
    MeshAttributes() {};
    MeshAttributes(const char*);
    void init(const char*);

  private:
    XMLDocument doc;
    void readXmlFile(const char* xmlFilename);
    void set_MeshRefinementZoneCube ();
    void set_SurfaceMeshingAttributes();
    void set_VolumeMeshingAttributes();
    void set_UseDiscreteMesh();
    void set_surfaceMSize();
    void set_MeshSizePropagation();
    void set_regionMSize();
    void set_global_mesh_attributes();
    void update_attribute_from_xml(XMLNode& element, const char* sElementName, const char* sAttributeName, int& MeshAttributesMember);
    void update_attribute_from_xml(XMLNode& element, const char* sElementName, const char* sAttributeName, double& MeshAttributesMember);
    std::string retrieve_text_from_xml(XMLNode* element, const char* sElementName);
    std::list<int>  fill_list_using_parsed_string(std::string line);
};



#endif
