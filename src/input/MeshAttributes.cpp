#include "MeshAttributes.h"

MeshAttributes::MeshAttributes()
    : lpair_lVertexId_MSize(), lpair_lEdgeId_MSize(), lpair_lFaceId_MSize(),
      lpair_lRegionId_MSize() {}

MeshAttributes::MeshAttributes(const char* xmlFilename) { init(xmlFilename); }

void MeshAttributes::init(const char* xmlFilename) {
    readXmlFile(xmlFilename);
    set_MeshRefinementZoneCube();
    set_SurfaceMeshingAttributes();
    set_VolumeMeshingAttributes();
    set_UseDiscreteMesh();
    set_vertexMSize();
    set_edgeMSize();
    set_surfaceMSize();
    set_MeshSizePropagation();
    set_regionMSize();
    set_global_mesh_attributes();
}

void MeshAttributes::readXmlFile(const char* xmlFilename) {
    if (doc.LoadFile(xmlFilename) != XML_SUCCESS) {
        logError() << "Unable to open or to parse file" << xmlFilename;
    }
}

void MeshAttributes::set_MeshRefinementZoneCube() {
    for (auto child = doc.FirstChildElement("MeshRefinementZoneCube"); child;
         child = child->NextSiblingElement("MeshRefinementZoneCube")) {
        Cube mycube;
        mycube.CubeMSize = std::atof(child->Attribute("value"));
        for (char i = 0; i < 3; i++) {
            char xyz[] = {static_cast<char>('x' + i), '\0'};
            update_attribute_from_xml(*child, "Center", xyz, mycube.CubeCenter[i]);
            update_attribute_from_xml(*child, "Width", xyz, mycube.CubeWidth[i]);
            update_attribute_from_xml(*child, "Height", xyz, mycube.CubeHeight[i]);
            update_attribute_from_xml(*child, "Depth", xyz, mycube.CubeDepth[i]);
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
        if (sSmoothingType.compare("Laplacian") == 0) {
            surfaceSmoothingType = SmoothingType::Laplacian;
        } else if (sSmoothingType.compare("Gradient") == 0) {
            surfaceSmoothingType = SmoothingType::Gradient;
        } else {
            logError() << "Unrecognised surfaceSmoothingType (Laplacian or Gradient)"
                       << sSmoothingType;
        }
        logInfo(PMU_rank()) << "surface smoothing option: surfaceSmoothingLevel "
                               "surfaceSmoothingType surfaceFaceRotationLimit Snap"
                            << surfaceSmoothingLevel << " "
                            << static_cast<int>(surfaceSmoothingType) << " "
                            << surfaceFaceRotationLimit << " " << surfaceSnap;
    }
}

void MeshAttributes::set_VolumeMeshingAttributes() {
    update_attribute_from_xml(doc, "VolumeMeshing", "SmoothingLevel", volumeSmoothingLevel);
    update_attribute_from_xml(doc, "VolumeMeshing", "SetOptimisation", VolumeMesherOptimization);
    if (auto pRoot = doc.FirstChildElement("VolumeMeshing")) {
        std::string sSmoothingType = pRoot->Attribute("SmoothingType");
        if (sSmoothingType.compare("Laplacian") == 0) {
            volumeSmoothingType = SmoothingType::Laplacian;
        } else if (sSmoothingType.compare("Gradient") == 0) {
            volumeSmoothingType = SmoothingType::Gradient;
        } else {
            logError() << "Unrecognised volumeSmoothingType (Laplacian or Gradient)"
                       << sSmoothingType;
        }
        logInfo(PMU_rank()) << "volume smoothing option: volumeSmoothingLevel volumeSmoothingType"
                            << volumeSmoothingLevel << " " << static_cast<int>(volumeSmoothingType);
    }
}

void MeshAttributes::set_UseDiscreteMesh() {
    update_attribute_from_xml(doc, "UseDiscreteMesh", "noModification",
                              UseDiscreteMesh_noModification);
    if (auto child = doc.FirstChildElement("UseDiscreteMesh")) {
        lFaceIdUseDiscreteMesh = fill_list_using_parsed_string(child->GetText());
    }
}

void MeshAttributes::set_vertexMSize() {
    for (auto child = doc.FirstChildElement("vertexMSize"); child;
         child = child->NextSiblingElement("vertexMSize")) {
        lpair_lVertexId_MSize.push_back(std::make_pair(
            fill_list_using_parsed_string(child->GetText()), std::atof(child->Attribute("value"))));
    }
}

void MeshAttributes::set_edgeMSize() {
    for (auto child = doc.FirstChildElement("edgeMSize"); child;
         child = child->NextSiblingElement("edgeMSize")) {
        lpair_lEdgeId_MSize.push_back(std::make_pair(
            fill_list_using_parsed_string(child->GetText()), std::atof(child->Attribute("value"))));
    }
}

void MeshAttributes::set_surfaceMSize() {
    for (auto child = doc.FirstChildElement("surfaceMSize"); child;
         child = child->NextSiblingElement("surfaceMSize")) {
        lpair_lFaceId_MSize.push_back(std::make_pair(
            fill_list_using_parsed_string(child->GetText()), std::atof(child->Attribute("value"))));
    }
}

void MeshAttributes::set_MeshSizePropagation() {
    update_attribute_from_xml(doc, "MeshSizePropagation", "ScalingFactor",
                              MeshSizePropagationScalingFactor);
    update_attribute_from_xml(doc, "MeshSizePropagation", "Distance", MeshSizePropagationDistance);
    if (auto child = doc.FirstChildElement("MeshSizePropagation")) {
        lFaceIdMeshSizePropagation = fill_list_using_parsed_string(child->GetText());
    }
}

void MeshAttributes::set_regionMSize() {
    for (auto child = doc.FirstChildElement("regionMSize"); child;
         child = child->NextSiblingElement("regionMSize")) {
        lpair_lRegionId_MSize.push_back(std::make_pair(
            fill_list_using_parsed_string(child->GetText()), std::atof(child->Attribute("value"))));
    }
}

void MeshAttributes::set_global_mesh_attributes() {
    update_attribute_from_xml(doc, "gradation", "value", gradation);
    update_attribute_from_xml(doc, "globalMSize", "value", globalMSize);
    update_attribute_from_xml(doc, "area_AspectRatio", "value", area_AspectRatio);
    update_attribute_from_xml(doc, "vol_AspectRatio", "value", vol_AspectRatio);
}

const std::list<std::pair<std::list<int>, double>>& MeshAttributes::getMSizeList(ElementType type) {
    switch (type) {
    case ElementType::vertex:
        return lpair_lVertexId_MSize;
        break;
    case ElementType::edge:
        return lpair_lEdgeId_MSize;
        break;
    case ElementType::face:
        return lpair_lFaceId_MSize;
        break;
    case ElementType::region:
        return lpair_lRegionId_MSize;
        break;
    default:
        __builtin_unreachable();
        break;
    }
}

void MeshAttributes::update_attribute_from_xml(XMLNode& element, const char* sElementName,
                                               const char* sAttributeName,
                                               double& MeshAttributesMember) {
    if (auto child = element.FirstChildElement(sElementName)) {
        const char* attr = child->Attribute(sAttributeName);
        if (attr != NULL) {
            MeshAttributesMember = std::atof(attr);
        } else {
            logError() << "XML Tag " << sElementName << " has no attribute " << sAttributeName;
        }
    }
}

void MeshAttributes::update_attribute_from_xml(XMLNode& element, const char* sElementName,
                                               const char* sAttributeName,
                                               int& MeshAttributesMember) {
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
    return element->FirstChildElement(sElementName)->GetText();
}

std::list<int> MeshAttributes::fill_list_using_parsed_string(std::string line) {
    std::vector<std::string> tokens;
    split(tokens, line, ',');
    std::list<int> ltemp;
    for (int i = 0; i < tokens.size(); i++) {
        ltemp.push_back(std::atoi(tokens[i].c_str()));
    }
    return ltemp;
}
