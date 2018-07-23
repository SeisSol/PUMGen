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


using namespace tinyxml2;
class MeshAttributes {
  public:
    int surfaceSmoothingLevel=2;
    int surfaceSmoothingType=1;
    double surfaceFaceRotationLimit=5.0;
    int surfaceSnap=0;
    int volumeSmoothingLevel=1;
    int volumeSmoothingType=1;
    int VolumeMesherOptimization=1;
    double globalMSize=0., faultMSize=0., gradation=0., vol_AspectRatio=0., area_AspectRatio=0.;
    std::list<double> lsurfaceMSize;
    std::list <std::list<int>> llsurfaceMSizeFaceId;
    std::list<double> lregionMSize;
    std::list <std::list<int>> llregionMSizeRegionId;
    double CubeMSize=0,CubeCenter[3],CubeWidth[3],CubeHeight[3],CubeDepth[3];
    std::list<int> lFaceIdMeshSizePropagation;
    std::list<int> lFaceIdUseDiscreteMesh;
    int UseDiscreteMesh_noModification;
    double MeshSizePropagationScalingFactor,MeshSizePropagationDistance=0.0;
    MeshAttributes() {};
    MeshAttributes(const char*, int , int);
    void init(const char*, int , int);

  private:
    XMLDocument doc;
    void readXmlFile(const char* xmlFilename);
    void set_MeshRefinementZoneCube ();
    void set_SurfaceMeshingAttributes();
    void set_VolumeMeshingAttributes();
    void set_UseDiscreteMesh(int numFaces);
    void set_surfaceMSize(int numFaces);
    void set_MeshSizePropagation(int numFaces);
    void set_regionMSize(int numRegions);
    void set_gradation();
    void set_globalMSize();
    void set_area_AspectRatio();
    void set_vol_AspectRatio();
};


#endif
