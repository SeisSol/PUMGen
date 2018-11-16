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

#include "MeshAttributes.h"

MeshAttributes::MeshAttributes(const char* xmlFilename, int numFaces, int numRegions) {
    init(xmlFilename, numFaces, numRegions);
}

void MeshAttributes::init(const char* xmlFilename, int numFaces, int numRegions) {
    readXmlFile(xmlFilename);
    set_MeshRefinementZoneCube ();
    set_SurfaceMeshingAttributes();
    set_VolumeMeshingAttributes();
    set_UseDiscreteMesh(numFaces);
    set_surfaceMSize(numFaces);
    set_MeshSizePropagation(numFaces);
    set_regionMSize(numRegions);
    set_gradation();
    set_globalMSize();
    set_area_AspectRatio();
    set_vol_AspectRatio();
}

void MeshAttributes::readXmlFile(const char* xmlFilename) {
    if(doc.LoadFile(xmlFilename) !=XML_SUCCESS) {
       logError() << "Unable to open or to parse file" << xmlFilename;
    }
}

void MeshAttributes::set_MeshRefinementZoneCube() {
    std::string sval;
    XMLElement* child;
    XMLElement* child2;
    child = doc.FirstChildElement("MeshRefinementZoneCube");
    Cube mycube;
    for (child; child; child = child->NextSiblingElement("MeshRefinementZoneCube"))
       {
       sval = child->Attribute("value");
       mycube.CubeMSize =  atof(sval.c_str());

       child2 = child->FirstChildElement("Center");
       sval = child2->Attribute("x");
       mycube.CubeCenter[0] =  atof(sval.c_str());
       sval = child2->Attribute("y");
       mycube.CubeCenter[1] =  atof(sval.c_str());
       sval = child2->Attribute("z");
       mycube.CubeCenter[2] =  atof(sval.c_str());

       child2 = child->FirstChildElement("Width");
       sval = child2->Attribute("x");
       mycube.CubeWidth[0] =  atof(sval.c_str());
       sval = child2->Attribute("y");
       mycube.CubeWidth[1] =  atof(sval.c_str());
       sval = child2->Attribute("z");
       mycube.CubeWidth[2] =  atof(sval.c_str());

       child2 = child->FirstChildElement("Height");
       sval = child2->Attribute("x");
       mycube.CubeHeight[0] =  atof(sval.c_str());
       sval = child2->Attribute("y");
       mycube.CubeHeight[1] =  atof(sval.c_str());
       sval = child2->Attribute("z");
       mycube.CubeHeight[2] =  atof(sval.c_str());

       child2 = child->FirstChildElement("Depth");
       sval = child2->Attribute("x");
       mycube.CubeDepth[0] =  atof(sval.c_str());
       sval = child2->Attribute("y");
       mycube.CubeDepth[1] =  atof(sval.c_str());
       sval = child2->Attribute("z");
       mycube.CubeDepth[2] =  atof(sval.c_str());
       lCube.push_back(mycube);
    }
}


void MeshAttributes::set_SurfaceMeshingAttributes() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("SurfaceMeshing");
    if(pRoot) {
       sval = pRoot->Attribute("SmoothingLevel");
       surfaceSmoothingLevel  =  atoi(sval.c_str());
       sval = pRoot->Attribute("Snap");
       surfaceSnap  =  atoi(sval.c_str());
       sval = pRoot->Attribute("SmoothingType");
       if (sval.compare("Laplacian")==0) {
          surfaceSmoothingType=0;
       } else if (sval.compare("Gradient")==0) {
          surfaceSmoothingType=1;
       } else {
          logError() << "Unrecognised surfaceSmoothingType (Laplacian or Gradient)" << sval;
       }
       sval = pRoot->Attribute("DiscreteAngle");
       surfaceFaceRotationLimit  =  atof(sval.c_str());
       logInfo(PMU_rank()) << "surface smoothing option: surfaceSmoothingLevel surfaceSmoothingType surfaceFaceRotationLimit Snap"<<
         surfaceSmoothingLevel<<" "<< surfaceSmoothingType<<" "<< surfaceFaceRotationLimit <<" "<<surfaceSnap;
    }
}

void MeshAttributes::set_VolumeMeshingAttributes() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("VolumeMeshing");
    if(pRoot) {
       sval = pRoot->Attribute("SmoothingLevel");
       volumeSmoothingLevel  =  atoi(sval.c_str());
       sval = pRoot->Attribute("SetOptimisation");
       VolumeMesherOptimization =  atoi(sval.c_str());
       sval = pRoot->Attribute("SmoothingType");
       if (sval.compare("Laplacian")==0) {
          volumeSmoothingType=0;
       } else if (sval.compare("Gradient")==0) {
          volumeSmoothingType=1;
       } else {
          logError() << "Unrecognised volumeSmoothingType (Laplacian or Gradient)" << sval;
       }
       logInfo(PMU_rank()) << "volume smoothing option: volumeSmoothingLevel volumrSmoothingType"<< 
         volumeSmoothingLevel<<" "<< volumeSmoothingType;
    }
}

void MeshAttributes::set_UseDiscreteMesh(int numFaces) {
    std::string line,sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("UseDiscreteMesh");
    if(pRoot) {
       sval = pRoot->Attribute("noModification");
       UseDiscreteMesh_noModification =  atoi(sval.c_str());
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           int faceid = std::atoi(tokens[i].c_str());
           assert (faceid<=numFaces);
           lFaceIdUseDiscreteMesh.push_back(faceid);
       }
    }
}

void MeshAttributes::set_surfaceMSize(int numFaces) {
    std::string line,sval;
    XMLElement* pRoot;
    XMLElement* child;
    child = doc.FirstChildElement("surfaceMSize");
    for (child; child; child = child->NextSiblingElement("surfaceMSize"))
    {
       sval = child->Attribute("value");
       double surfaceMSize =  atof(sval.c_str());
       lsurfaceMSize.push_back(surfaceMSize);
       std::vector<std::string> tokens;
       line =  child-> GetText();
       split(tokens, line, ',');
       std::list<int> lFaceId;
       for(int i = 0; i < tokens.size(); i++) {
           int faceid = std::atoi(tokens[i].c_str());
           assert (faceid<=numFaces);
           lFaceId.push_back(faceid);
       }
       llsurfaceMSizeFaceId.push_back(lFaceId);
    }
}

void MeshAttributes::set_MeshSizePropagation(int numFaces) {
    std::string line,sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("MeshSizePropagation");
    if(pRoot) {
       sval = pRoot->Attribute("ScalingFactor");
       MeshSizePropagationScalingFactor =  atof(sval.c_str());
       sval = pRoot->Attribute("Distance");
       MeshSizePropagationDistance =  atof(sval.c_str());
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           int faceid = std::atoi(tokens[i].c_str());
           assert (faceid<=numFaces);
           lFaceIdMeshSizePropagation.push_back(faceid);
       }
    }
}


void MeshAttributes::set_regionMSize(int numRegions) {
    std::string line,sval;
    XMLElement* child;
    XMLElement* pRoot;
    child = doc.FirstChildElement("regionMSize");
    for (child; child; child = child->NextSiblingElement("regionMSize"))
    {
       sval = child->Attribute("value");
       double regionMSize =  atof(sval.c_str());
       lregionMSize.push_back(regionMSize);
       std::vector<std::string> tokens;
       line =  child-> GetText();
       split(tokens, line, ',');
       std::list<int> lregionId;
       for(int i = 0; i < tokens.size(); i++) {
           int regionid = std::atoi(tokens[i].c_str());
           assert (regionid<=numRegions);
           lregionId.push_back(regionid);
       }
       llregionMSizeRegionId.push_back(lregionId);
    }
}


void MeshAttributes::set_gradation() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("gradation");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       gradation =  atof(sval.c_str());
    }
}

void MeshAttributes::set_globalMSize() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("globalMSize");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       globalMSize =  atof(sval.c_str());
    }
}


void MeshAttributes::set_area_AspectRatio() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("area_AspectRatio");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       area_AspectRatio =  atof(sval.c_str());
    }
}

void MeshAttributes::set_vol_AspectRatio() {
    std::string sval;
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("vol_AspectRatio");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       vol_AspectRatio =  atof(sval.c_str());
    }
}


/*
    //removed from previous version:
    //- detecSmallFeatures, which could be useful but does not fall into the scope of this file
    //- write smd which never worked properly

    pRoot = doc.FirstChildElement("detectSmallFeatures");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       detectSmallFeaturesThreshold =  atof(sval.c_str());
    }

    pRoot = doc.FirstChildElement("writeSmd");
    if(pRoot) {
       writeSmd = 1; 
    }

    if (detectSmallFeaturesThreshold>0) {
       pSmallFeatureInfo smallFeats = GM_detectSmallFeatures(model,1,detectSmallFeaturesThreshold,0,0,0);
       pPList lsmallFeats=GM_getSmallFeatures(smallFeats); 
       logInfo(PMU_rank()) << "Number of small features returned: " <<PList_size(lsmallFeats);
    }

    if (writeSmd>0) {
	    logInfo(PMU_rank()) << "writing the smd";
	    GM_setAttManager(model, attMngr);
	    GM_write(model,"model.smd",0,0); // write out the model before the mesh!
	    logInfo(PMU_rank()) << "done writing smd";
    }


*/
