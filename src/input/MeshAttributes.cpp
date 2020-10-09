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
           retrieve_attribute_and_update_class_member(child, "Center", &xyz[i], &mycube.CubeCenter[i]);
           retrieve_attribute_and_update_class_member(child, "Width", &xyz[i], &mycube.CubeWidth[i]);
           retrieve_attribute_and_update_class_member(child, "Height", &xyz[i], &mycube.CubeHeight[i]);
           retrieve_attribute_and_update_class_member(child, "Depth", &xyz[i], &mycube.CubeDepth[i]);
       }
       lCube.push_back(mycube);
    }
}


void MeshAttributes::set_SurfaceMeshingAttributes() {
    retrieve_attribute_and_update_class_member(&doc, "SurfaceMeshing", "DiscreteAngle", &surfaceFaceRotationLimit);
    retrieve_attribute_and_update_class_member(&doc, "SurfaceMeshing", "SmoothingLevel", &surfaceSmoothingLevel);
    retrieve_attribute_and_update_class_member(&doc, "SurfaceMeshing", "Snap", &surfaceSnap);
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
    retrieve_attribute_and_update_class_member(&doc, "VolumeMeshing", "SmoothingLevel", &volumeSmoothingLevel);
    retrieve_attribute_and_update_class_member(&doc, "VolumeMeshing", "SetOptimisation", &VolumeMesherOptimization);
    if (auto pRoot = doc.FirstChildElement("VolumeMeshing")) {
        std::string sSmoothingType = pRoot->Attribute("SmoothingType");
        if (sSmoothingType.compare("Laplacian")==0) {
           volumeSmoothingType=0;
        } else if (sSmoothingType.compare("Gradient")==0) {
           volumeSmoothingType=1;
        } else {
           logError() << "Unrecognised volumeSmoothingType (Laplacian or Gradient)" << sSmoothingType;
        }
        logInfo(PMU_rank()) << "volume smoothing option: volumeSmoothingLevel volumrSmoothingType"<< 
            volumeSmoothingLevel<<" "<< volumeSmoothingType;
    }
}

void MeshAttributes::set_UseDiscreteMesh() {
    retrieve_attribute_and_update_class_member(&doc, "UseDiscreteMesh", "noModification", &UseDiscreteMesh_noModification);
    if (auto child = doc.FirstChildElement("UseDiscreteMesh")) {
       lFaceIdUseDiscreteMesh = parse_string_and_fill_list(child-> GetText());
    }
}


void MeshAttributes::set_surfaceMSize() {
    for (auto child = doc.FirstChildElement("surfaceMSize"); child; child = child->NextSiblingElement("surfaceMSize"))
    {
       lsurfaceMSize.push_back( std::atof(child->Attribute("value")));
       llsurfaceMSizeFaceId.push_back(parse_string_and_fill_list(child-> GetText()));
    }
}


void MeshAttributes::set_MeshSizePropagation() {
    retrieve_attribute_and_update_class_member(&doc, "MeshSizePropagation", "ScalingFactor", &MeshSizePropagationScalingFactor);
    retrieve_attribute_and_update_class_member(&doc, "MeshSizePropagation", "Distance", &MeshSizePropagationDistance);
    if (auto child = doc.FirstChildElement("MeshSizePropagation")) {
       lFaceIdMeshSizePropagation = parse_string_and_fill_list(child-> GetText());
    }
}


void MeshAttributes::set_regionMSize() {
    for (auto child = doc.FirstChildElement("regionMSize"); child; child = child->NextSiblingElement("regionMSize"))
    {
       lregionMSize.push_back( std::atof(child->Attribute("value")));
       llregionMSizeRegionId.push_back(parse_string_and_fill_list(child-> GetText()));
    }
}


void MeshAttributes::set_global_mesh_attributes() {
    retrieve_attribute_and_update_class_member(&doc, "gradation", "value", &gradation);
    retrieve_attribute_and_update_class_member(&doc, "globalMSize", "value", &globalMSize);
    retrieve_attribute_and_update_class_member(&doc, "area_AspectRatio", "value", &area_AspectRatio);
    retrieve_attribute_and_update_class_member(&doc, "vol_AspectRatio", "value", &vol_AspectRatio);

}


void MeshAttributes::retrieve_attribute_and_update_class_member(XMLNode* element, const char* sElementName, const char* sAttributeName, double* MeshAttributesMember) {
    if (auto child = element->FirstChildElement(sElementName)) {
        *MeshAttributesMember = std::atof(child->Attribute(sAttributeName));
    }
}
void MeshAttributes::retrieve_attribute_and_update_class_member(XMLNode* element, const char* sElementName, const char* sAttributeName, int* MeshAttributesMember) {
    if (auto child = element->FirstChildElement(sElementName)) {
        *MeshAttributesMember = std::atoi(child->Attribute(sAttributeName));
    }
}


std::string MeshAttributes::retrieve_text_from_xml(XMLNode* element, const char* sElementName) {
    return element->FirstChildElement(sElementName)-> GetText();
}    

std::list<int>  MeshAttributes::parse_string_and_fill_list(std::string line) {
    std::vector<std::string> tokens;
    split(tokens, line, ',');
    std::list<int> ltemp;
    for(int i = 0; i < tokens.size(); i++) {
        ltemp.push_back(std::atoi(tokens[i].c_str()));
    }
    return ltemp;
}
