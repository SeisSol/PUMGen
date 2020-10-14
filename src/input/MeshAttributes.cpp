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

MeshAttributes::MeshAttributes(const char* xmlFilename) {
    init(xmlFilename);
}

void MeshAttributes::init(const char* xmlFilename) {
    readXmlFile(xmlFilename);
    set_MeshRefinementZoneCube ();
    set_SurfaceMeshingAttributes();
    set_VolumeMeshingAttributes();
    set_UseDiscreteMesh();
    set_surfaceMSize();
    set_MeshSizePropagation();
    set_regionMSize();
    set_global_mesh_attributes();
}

void MeshAttributes::readXmlFile(const char* xmlFilename) {
    if(doc.LoadFile(xmlFilename) !=XML_SUCCESS) {
       logError() << "Unable to open or to parse file" << xmlFilename;
    }
}

void MeshAttributes::set_MeshRefinementZoneCube() {
    std::string xyz = "xyz";
    for (auto child = doc.FirstChildElement("MeshRefinementZoneCube"); child; child = child->NextSiblingElement("MeshRefinementZoneCube"))
       {
       Cube mycube;
       mycube.CubeMSize =  std::atof(child->Attribute("value"));
       for (std::string::size_type i = 0; i < xyz.size(); i++) {
           update_attribute_from_xml(*child, "Center", &xyz[i], mycube.CubeCenter[i]);
           update_attribute_from_xml(*child, "Width", &xyz[i], mycube.CubeWidth[i]);
           update_attribute_from_xml(*child, "Height", &xyz[i], mycube.CubeHeight[i]);
           update_attribute_from_xml(*child, "Depth", &xyz[i], mycube.CubeDepth[i]);
       }
       lCube.push_back(mycube);
    }
}


void MeshAttributes::set_SurfaceMeshingAttributes() {
    update_attribute_from_xml(doc, "SurfaceMeshing", "DiscreteAngle", surfaceFaceRotationLimit);
    update_attribute_from_xml(doc, "SurfaceMeshing", "SmoothingLevel", surfaceSmoothingLevel);
    update_attribute_from_xml(doc, "SurfaceMeshing", "Snap", surfaceSnap);
    if (auto pRoot = doc.FirstChildElement("SurfaceMeshing")) {
        std::string sSmoothingType = pRoot->Attribute("SmoothingType");
        if (sSmoothingType.compare("Laplacian")==0) {
           surfaceSmoothingType=0;
        } else if (sSmoothingType.compare("Gradient")==0) {
           surfaceSmoothingType=1;
        } else {
           logError() << "Unrecognised surfaceSmoothingType (Laplacian or Gradient)" << sSmoothingType;
        } 
        logInfo(PMU_rank()) << "surface smoothing option: surfaceSmoothingLevel surfaceSmoothingType surfaceFaceRotationLimit Snap"<<
            surfaceSmoothingLevel<<" "<< surfaceSmoothingType<<" "<< surfaceFaceRotationLimit <<" "<<surfaceSnap;
    }
}

void MeshAttributes::set_VolumeMeshingAttributes() {
    update_attribute_from_xml(doc, "VolumeMeshing", "SmoothingLevel", volumeSmoothingLevel);
    update_attribute_from_xml(doc, "VolumeMeshing", "SetOptimisation", VolumeMesherOptimization);
    if (auto pRoot = doc.FirstChildElement("VolumeMeshing")) {
        std::string sSmoothingType = pRoot->Attribute("SmoothingType");
        if (sSmoothingType.compare("Laplacian")==0) {
           volumeSmoothingType = 0;
        } else if (sSmoothingType.compare("Gradient")==0) {
           volumeSmoothingType = 1;
        } else {
           logError() << "Unrecognised volumeSmoothingType (Laplacian or Gradient)" << sSmoothingType;
        }
        logInfo(PMU_rank()) << "volume smoothing option: volumeSmoothingLevel volumeSmoothingType"<< 
            volumeSmoothingLevel<<" "<< volumeSmoothingType;
    }
}

void MeshAttributes::set_UseDiscreteMesh() {
    update_attribute_from_xml(doc, "UseDiscreteMesh", "noModification", UseDiscreteMesh_noModification);
    if (auto child = doc.FirstChildElement("UseDiscreteMesh")) {
       lFaceIdUseDiscreteMesh = fill_list_using_parsed_string(child-> GetText());
    }
}


void MeshAttributes::set_surfaceMSize() {
    for (auto child = doc.FirstChildElement("surfaceMSize"); child; child = child->NextSiblingElement("surfaceMSize"))
    {
       lsurfaceMSize.push_back( std::atof(child->Attribute("value")));
       llsurfaceMSizeFaceId.push_back(fill_list_using_parsed_string(child-> GetText()));
    }
}


void MeshAttributes::set_MeshSizePropagation() {
    update_attribute_from_xml(doc, "MeshSizePropagation", "ScalingFactor", MeshSizePropagationScalingFactor);
    update_attribute_from_xml(doc, "MeshSizePropagation", "Distance", MeshSizePropagationDistance);
    if (auto child = doc.FirstChildElement("MeshSizePropagation")) {
       lFaceIdMeshSizePropagation = fill_list_using_parsed_string(child-> GetText());
    }
}


void MeshAttributes::set_regionMSize() {
    for (auto child = doc.FirstChildElement("regionMSize"); child; child = child->NextSiblingElement("regionMSize"))
    {
       lregionMSize.push_back( std::atof(child->Attribute("value")));
       llregionMSizeRegionId.push_back(fill_list_using_parsed_string(child-> GetText()));
    }
}


void MeshAttributes::set_global_mesh_attributes() {
    update_attribute_from_xml(doc, "gradation", "value", gradation);
    update_attribute_from_xml(doc, "globalMSize", "value", globalMSize);
    update_attribute_from_xml(doc, "area_AspectRatio", "value", area_AspectRatio);
    update_attribute_from_xml(doc, "vol_AspectRatio", "value", vol_AspectRatio);

}


void MeshAttributes::update_attribute_from_xml(XMLNode& element, const char* sElementName, const char* sAttributeName, double& MeshAttributesMember) {
    if (auto child = element.FirstChildElement(sElementName)) {
        const char* attr = child->Attribute(sAttributeName);
        if (attr != NULL) {
            MeshAttributesMember = std::atof(attr);
        } else {
            logError() << "XML Tag " << sElementName << " has no attribute " << sAttributeName;
        }
    }
}

void MeshAttributes::update_attribute_from_xml(XMLNode& element, const char* sElementName, const char* sAttributeName, int& MeshAttributesMember) {
    if (auto child = element.FirstChildElement(sElementName)) {
        const char* attr = child->Attribute(sAttributeName);
        if (attr != NULL) {
            MeshAttributesMember = std::atoi(attr);
        } else {
            logError() << "XML Tag " << sElementName << " has no attribute " << sAttributeName;
        }
    }
}



std::string MeshAttributes::retrieve_text_from_xml(XMLNode* element, const char* sElementName) {
    return element->FirstChildElement(sElementName)-> GetText();
}    

std::list<int>  MeshAttributes::fill_list_using_parsed_string(std::string line) {
    std::vector<std::string> tokens;
    split(tokens, line, ',');
    std::list<int> ltemp;
    for(int i = 0; i < tokens.size(); i++) {
        ltemp.push_back(std::atoi(tokens[i].c_str()));
    }
    return ltemp;
}
