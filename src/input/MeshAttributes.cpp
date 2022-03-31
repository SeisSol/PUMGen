#include "MeshAttributes.h"

VelocityAwareRefinementSettings::VelocityAwareRefinementSettings(double elementsPerWaveLength,
                                                                 std::string easiFileName)
    : elementsPerWaveLength(elementsPerWaveLength), easiFileName(std::move(easiFileName)) {}

void VelocityAwareRefinementSettings::addRefinementRegion(SimpleCuboid cuboid,
                                                          double targetedFrequency,
                                                          int bypassFindRegionAndUseGroup) {
  refinementRegions.emplace_back(cuboid, targetedFrequency, bypassFindRegionAndUseGroup);
}

bool VelocityAwareRefinementSettings::isVelocityAwareRefinementOn() const {
  return !refinementRegions.empty();
}

const std::string& VelocityAwareRefinementSettings::getEasiFileName() const { return easiFileName; }

double VelocityAwareRefinementSettings::getElementsPerWaveLength() const {
  return elementsPerWaveLength;
}

const std::vector<VelocityRefinementCube>&
VelocityAwareRefinementSettings::getRefinementRegions() const {
  return refinementRegions;
}

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
  set_NoMesh();
  set_vertexMSize();
  set_edgeMSize();
  set_surfaceMSize();
  set_MeshSizePropagation();
  set_regionMSize();
  set_global_mesh_attributes();
  set_velocity_aware_meshing();
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
      update_attribute_from_xml(*child, "HalfWidth", xyz, mycube.CubeWidth[i]);
      update_attribute_from_xml(*child, "HalfHeight", xyz, mycube.CubeHeight[i]);
      update_attribute_from_xml(*child, "HalfDepth", xyz, mycube.CubeDepth[i]);
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
      logError() << "Unrecognised surfaceSmoothingType (Laplacian or Gradient)" << sSmoothingType;
    }
    logInfo(PMU_rank()) << "surface smoothing option: surfaceSmoothingLevel "
                           "surfaceSmoothingType surfaceFaceRotationLimit Snap"
                        << surfaceSmoothingLevel << " " << static_cast<int>(surfaceSmoothingType)
                        << " " << surfaceFaceRotationLimit << " " << surfaceSnap;
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
      logError() << "Unrecognised volumeSmoothingType (Laplacian or Gradient)" << sSmoothingType;
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

void MeshAttributes::set_NoMesh() {
  if (auto child = doc.FirstChildElement("surfaceNoMesh")) {
    lFaceIdNoMesh = fill_list_using_parsed_string(child->GetText());
  }
  if (auto child = doc.FirstChildElement("regionNoMesh")) {
    lRegionIdNoMesh = fill_list_using_parsed_string(child->GetText());
  }
}

void MeshAttributes::set_vertexMSize() {
  for (auto child = doc.FirstChildElement("vertexMSize"); child;
       child = child->NextSiblingElement("vertexMSize")) {
    lpair_lVertexId_MSize.push_back(std::make_pair(fill_list_using_parsed_string(child->GetText()),
                                                   std::atof(child->Attribute("value"))));
  }
}

void MeshAttributes::set_edgeMSize() {
  for (auto child = doc.FirstChildElement("edgeMSize"); child;
       child = child->NextSiblingElement("edgeMSize")) {
    lpair_lEdgeId_MSize.push_back(std::make_pair(fill_list_using_parsed_string(child->GetText()),
                                                 std::atof(child->Attribute("value"))));
  }
}

void MeshAttributes::set_surfaceMSize() {
  for (auto child = doc.FirstChildElement("surfaceMSize"); child;
       child = child->NextSiblingElement("surfaceMSize")) {
    lpair_lFaceId_MSize.push_back(std::make_pair(fill_list_using_parsed_string(child->GetText()),
                                                 std::atof(child->Attribute("value"))));
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
    lpair_lRegionId_MSize.push_back(std::make_pair(fill_list_using_parsed_string(child->GetText()),
                                                   std::atof(child->Attribute("value"))));
  }
}

void MeshAttributes::set_velocity_aware_meshing() {
  int numChilds = 0;
  const auto name = "VelocityAwareMeshing";
  for (auto velocityAwareMeshingElement = doc.FirstChildElement(name); velocityAwareMeshingElement;
       velocityAwareMeshingElement = velocityAwareMeshingElement->NextSiblingElement(name)) {
    const auto easiFileName = velocityAwareMeshingElement->Attribute("easiFile");
    const auto elementsPerWaveLength =
        std::stof(velocityAwareMeshingElement->Attribute("elementsPerWaveLength"));
    velocityAwareRefinementSettings =
        VelocityAwareRefinementSettings(elementsPerWaveLength, easiFileName);
    constexpr auto cuboidName = "VelocityRefinementCuboid";
    logInfo(PMU_rank()) << "Activating velocity aware meshing, using" << elementsPerWaveLength
                        << "elements per wavelength and easi file" << easiFileName;
    for (auto child = velocityAwareMeshingElement->FirstChildElement(cuboidName); child;
         child = child->NextSiblingElement(cuboidName)) {

      const char* attr = child->Attribute("rotationZAnticlockwiseFromX");
      double rotationZAnticlockwiseFromX = 0.0;
      if (attr != nullptr) {
        rotationZAnticlockwiseFromX = std::stof(attr);
      }
      attr = child->Attribute("bypassFindRegionAndUseGroup");
      int bypassFindRegionAndUseGroup = 0;
      if (attr != nullptr) {
        bypassFindRegionAndUseGroup = std::stoi(attr);
      }
      auto cuboid = SimpleCuboid{{
                                     std::stof(child->Attribute("centerX")),
                                     std::stof(child->Attribute("centerY")),
                                     std::stof(child->Attribute("centerZ")),
                                 },
                                 {
                                     std::stof(child->Attribute("halfSizeX")),
                                     std::stof(child->Attribute("halfSizeY")),
                                     std::stof(child->Attribute("halfSizeZ")),
                                 },
                                 {std::cos(rotationZAnticlockwiseFromX * toRadians),
                                  std::sin(rotationZAnticlockwiseFromX * toRadians)},
                                 rotationZAnticlockwiseFromX};
      const auto targetedFrequency = std::stof(child->Attribute("frequency"));

      velocityAwareRefinementSettings.addRefinementRegion(cuboid, targetedFrequency,
                                                          bypassFindRegionAndUseGroup);
      logInfo(PMU_rank()) << "Adding velocity aware refinement region targeting"
                          << targetedFrequency << "Hz, centered at x =" << cuboid.center[0]
                          << "y=" << cuboid.center[1] << "z=" << cuboid.center[2]
                          << "with half sizes"
                          << "x =" << cuboid.halfSize[0] << "y =" << cuboid.halfSize[1]
                          << "z =" << cuboid.halfSize[2];
      if (std::abs(cuboid.rotationZ) > 0.0) {
        logInfo(PMU_rank()) << "rotated around z axis by " << cuboid.rotationZ
                            << "degree(s) counterclockwise from x axis.";
      }
      if (bypassFindRegionAndUseGroup) {
        logInfo(PMU_rank()) << "bypass findRegion and use group =" << bypassFindRegionAndUseGroup;
      }
    }
    if (!velocityAwareRefinementSettings.isVelocityAwareRefinementOn()) {
      logWarning(PMU_rank())
          << "Activated velocity aware meshing but did not specify any refinement region!";
    }
    ++numChilds;
  }
  if (numChilds > 1) {
    logError() << "Multiple definitions of velocityAwareMeshing";
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
    if (attr != nullptr) {
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
    if (attr != nullptr) {
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
