// SPDX-FileCopyrightText: 2023 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_SIMMODSUITEAPF_H_
#define PUMGEN_SRC_INPUT_SIMMODSUITEAPF_H_

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <unordered_map>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <gmi_sim.h>

#include <MeshSim.h>
#include <SimDiscrete.h>
#ifdef BEFORE_SIM_18
#include <SimError.h>
#include <SimErrorCodes.h>
#include <SimMeshingErrorCodes.h>
#else
#include <SimInfo.h>
#include <SimInfoCodes.h>
#include <SimMeshingInfoCodes.h>
#endif
#include <SimModel.h>
#include <SimModelerUtil.h>
#ifdef PARASOLID
#include <SimParasolidKrnl.h>
#endif
#include <SimPartitionedMesh.h>

#include "utils/logger.h"
#include "utils/path.h"
#include "utils/progress.h"

#include "MeshInput.h"

#include "AnalysisAttributes.h"
#include "EasiMeshSize.h"
#include "MeshAttributes.h"

#include <SimDisplay.h>
#include <SimMeshTools.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

/**
 * @todo Currently it is not supported to create more than one instance
 *  of this class
 * @todo Maybe add MS_setMaxEntities to limit the number of elements
 */
class SimModSuiteApf : public ApfMeshInput {
  private:
  EasiMeshSize easiMeshSize;
  pGModel m_model;

  pParMesh m_simMesh;

  /** Enable Simmetrix logging file */
  bool m_log;

  public:
  SimModSuiteApf(const char* modFile, int boundarySize, const char* cadFile = nullptr,
                 const char* licenseFile = nullptr, const char* meshCaseName = "mesh",
                 const char* analysisCaseName = "analysis", int enforceSize = 0,
                 const char* xmlFile = nullptr, const bool analyseAR = false,
                 const char* logFile = nullptr)
      : ApfMeshInput(boundarySize) {
    // Init SimModSuite
    SimModel_start();
    SimPartitionedMesh_start(nullptr, nullptr);
    if (logFile) {
      m_log = true;
      Sim_logOn(logFile);
    } else
      m_log = false;
    Sim_readLicenseFile(licenseFile);
    MS_init();
    SimDiscrete_start(0);
#ifdef PARASOLID
    SimParasolid_start(1);
#endif
    Sim_setMessageHandler(messageHandler);

    // Load CAD
    logInfo() << "Loading model";

    std::string smodFile = modFile;
    if (cadFile != nullptr) {
      loadCAD(modFile, cadFile);
    } else if (smodFile.substr(smodFile.find_last_of(".") + 1) == "smd") {
      loadCAD(modFile, cadFile);
    } else {
      logError() << "unsupported model file";
    }

    // Create group lookup table
    std::unordered_map<pGRegion, int> groupMap;
    GRIter regionIt = GM_regionIter(m_model);
    while (pGRegion region = GRIter_next(regionIt)) {
      groupMap[region] = static_cast<int>(groupMap.size() + 1);
    }

    // Extract cases
    logInfo() << "Extracting cases";
    pACase meshCase, analysisCase;
    MeshAttributes MeshAtt;

    if (xmlFile != nullptr) {

      // Read mesh Attributes from xml file
      int numFaces = GM_numFaces(m_model);
      MeshAtt.init(xmlFile);
      AnalysisAttributes AnalysisAtt(xmlFile, numFaces);
      setCases(m_model, meshCase, analysisCase, MeshAtt, AnalysisAtt, groupMap);
    } else {
      extractCases(m_model, meshCase, meshCaseName, analysisCase, analysisCaseName);
    }

    m_simMesh = PM_new(0, m_model, PMU_size());

    pProgress prog = Progress_new();
    Progress_setCallback(prog, progressHandler);

    // create the mesh
    logInfo() << "Starting the surface mesher";
    pSurfaceMesher surfaceMesher = SurfaceMesher_new(meshCase, m_simMesh);
    if (xmlFile != nullptr) {
      SurfaceMesher_setSmoothing(surfaceMesher, MeshAtt.surfaceSmoothingLevel);
      SurfaceMesher_setSmoothType(surfaceMesher, static_cast<int>(MeshAtt.surfaceSmoothingType));
      SurfaceMesher_setFaceRotationLimit(surfaceMesher, MeshAtt.surfaceFaceRotationLimit);
      SurfaceMesher_setSnapForDiscrete(surfaceMesher, MeshAtt.surfaceSnap);
    }
    SurfaceMesher_execute(surfaceMesher, prog);
    SurfaceMesher_delete(surfaceMesher);

    PM_setTotalNumParts(m_simMesh, PMU_size());

    logInfo() << "Starting the volume mesher";
    pVolumeMesher volumeMesher = VolumeMesher_new(meshCase, m_simMesh);
    if (xmlFile != nullptr) {
      VolumeMesher_setSmoothing(volumeMesher, MeshAtt.volumeSmoothingLevel);
      VolumeMesher_setSmoothType(volumeMesher, static_cast<int>(MeshAtt.volumeSmoothingType));
      VolumeMesher_setOptimization(volumeMesher, MeshAtt.VolumeMesherOptimization);
    }
    VolumeMesher_setEnforceSize(volumeMesher, enforceSize);
    VolumeMesher_execute(volumeMesher, prog);
    VolumeMesher_delete(volumeMesher);

    Progress_delete(prog);

    if (analyseAR) {
      analyse_mesh();
    }

    if (xmlFile == nullptr) {
      // unassociate case for attributes below
      AttCase_unassociate(meshCase);
    }

    // Convert to APF mesh
    apf::Mesh* tmpMesh = apf::createMesh(m_simMesh);
    gmi_register_sim();
    gmi_model* model = gmi_import_sim(m_model);

    logInfo() << "Converting mesh to APF";
    m_mesh = apf::createMdsMesh(model, tmpMesh);
    apf::destroyMesh(tmpMesh);

    // Set the boundary conditions from the geometric model
    AttCase_associate(analysisCase, nullptr);
    apf::MeshTag* boundaryTag = m_mesh->createIntTag("boundary condition", 1);
    apf::MeshIterator* it = m_mesh->begin(2);
    while (apf::MeshEntity* face = m_mesh->iterate(it)) {
      apf::ModelEntity* modelFace = m_mesh->toModel(face);
      if (m_mesh->getModelType(modelFace) != 2)
        continue;

      pGEntity simFace = reinterpret_cast<pGEntity>(modelFace);

      pAttribute attr = GEN_attrib(simFace, "boundaryCondition");
      if (attr) {
        char* image = Attribute_imageClass(attr);
        int boundary = parseBoundary(image);
        Sim_deleteString(image);

        m_mesh->setIntTag(face, boundaryTag, &boundary);
      }
    }
    m_mesh->end(it);

    // Set groups
    apf::MeshTag* groupTag = m_mesh->createIntTag("group", 1);
    it = m_mesh->begin(3);
    while (apf::MeshEntity* element = m_mesh->iterate(it)) {
      apf::ModelEntity* modelRegion = m_mesh->toModel(element);

      pGRegion simRegion = reinterpret_cast<pGRegion>(modelRegion);
      std::unordered_map<pGRegion, int>::const_iterator i = groupMap.find(simRegion);
      if (i == groupMap.end())
        logError() << "Mesh element with unknown region found.";

      m_mesh->setIntTag(element, groupTag, &i->second);
    }
    m_mesh->end(it);

    AttCase_unassociate(analysisCase);

    // Delete cases
    MS_deleteMeshCase(meshCase);
    MS_deleteMeshCase(analysisCase);
  }

  virtual ~SimModSuiteApf() {
    M_release(m_simMesh);
    // We cannot delete the model here because it is still
    // connected to the mesh
    // GM_release(m_model);

    // Finalize SimModSuite
#ifdef PARASOLID
    SimParasolid_stop(1);
#endif
    SimDiscrete_stop(0);
    MS_exit();
    Sim_unregisterAllKeys();
    if (m_log)
      Sim_logOff();
    SimPartitionedMesh_stop();
    SimModel_stop();
  }

  private:
  pACase extractCase(pAManager attMngr, const char* name) {
    pACase acase = AMAN_findCase(attMngr, name);
    if (!acase)
      logError() << "Case" << std::string(name) << "not found.";

    AttCase_setModel(acase, m_model);

    return acase;
  }

  private:
  /**
   * @todo Make this private as soon as APF does no longer need it
   */
  static unsigned int parseBoundary(const char* boundaryCondition) {
    if (strcmp(boundaryCondition, "freeSurface") == 0)
      return 1;
    if (strcmp(boundaryCondition, "dynamicRupture") == 0)
      return 3;
    if (strcmp(boundaryCondition, "absorbing") == 0)
      return 5;

    // check if boundaryCondition starts with pattern
    std::string pattern = "BC";
    std::string sboundaryCondition(boundaryCondition, 11);
    if (sboundaryCondition.find(pattern) == 0) {
      std::string sNumber = sboundaryCondition.substr(pattern.length());
      return stoi(sNumber);
    }
    logError() << "Unknown boundary condition" << boundaryCondition;
    return -1;
  }

  static void messageHandler(int type, const char* msg) {
    switch (type) {
    case Sim_InfoMsg:
      // Show sim info messages as debug messages
      logDebug() << "SimModeler:" << msg;
      break;
    case Sim_DebugMsg:
      // Ignore sim debug messages
      break;
    case Sim_WarningMsg:
      logWarning() << "SimModeler:" << msg;
      break;
    case Sim_ErrorMsg:
      // Use warning because error will abort the program
      logWarning() << "SimModeler:" << msg;
      break;
    }
  }

  static void progressHandler(const char* what, int level, int startVal, int endVal, int currentVal,
                              void* ignore) {
    if (PMU_rank() != 0)
      return;
    if (endVal != 0) {
      switch (currentVal) {
      case -2:
        // task is started, do nothing
        logInfo() << "Progress:" << what << ", 0"
                  << "/" << endVal;
        break;
      case -1:
        // end of the task
        logInfo() << "Progress:" << what << ", done";
        break;
      default:
        logInfo() << "Progress:" << what << "," << currentVal << "/" << endVal;
        break;
      }
    } else {
      if (currentVal == -2)
        logInfo() << "Progress:" << what;
    }
    logDebug() << what << level << startVal << endVal << currentVal;
  }

  private:
  void extractCases(pGModel m_model, pACase& meshCase, const char* meshCaseName,
                    pACase& analysisCase, const char* analysisCaseName) {
    logInfo() << "Extracting cases";

#ifdef BEFORE_SIM_2024
    pAManager attMngr = GM_attManager(m_model);
#else
    pAManager attMngr = GM_attManager(m_model, false);
#endif

    MeshingOptions meshingOptions;
    meshCase = MS_newMeshCase(m_model);

    pACase meshCaseFile = extractCase(attMngr, meshCaseName);
    AttCase_associate(meshCaseFile, nullptr);
    MS_processSimModelerMeshingAtts(meshCaseFile, meshCase, &meshingOptions);
    AttCase_setModel(meshCase, m_model);

    analysisCase = extractCase(attMngr, analysisCaseName);
    pPList children = AttNode_children(analysisCase);
    void* iter = nullptr;
    while (pANode child = static_cast<pANode>(PList_next(children, &iter)))
      AttCase_setModel(static_cast<pACase>(child), m_model);
    PList_delete(children);
  }

  private:
  void setMeshSize(pGModel model, pACase& meshCase, MeshAttributes& MeshAtt) {
    // Set mesh size on vertices, edges, surfaces and regions
    const char* element_name[4] = {"vertex", "edge", "face", "region"};
    const ElementType element_type[4] = {ElementType::vertex, ElementType::edge, ElementType::face,
                                         ElementType::region};

    for (int element_type_id = 0; element_type_id < 4; element_type_id++) {
      std::list<std::pair<std::list<int>, double>> lpair_lId_MSize =
          MeshAtt.getMSizeList(element_type[element_type_id]);
      for (const auto& pair_lId_MSize : lpair_lId_MSize) {
        double MSize = pair_lId_MSize.second;
        std::list<int> lElements = pair_lId_MSize.first;
        for (const auto& element_id : lElements) {
          pGEntity entity = GM_entityByTag(model, element_type_id, element_id);
          if (entity == nullptr) {
            logError() << element_name[element_type_id] << "id:" << element_id
                       << "not found in model.";
          } else {
            MS_setMeshSize(meshCase, entity, 1, MSize, nullptr);
            logInfo() << element_name[element_type_id] << "id:" << element_id
                      << ", MSize =" << MSize;
          }
        }
      }
    }
  }

  void setCases(pGModel model, pACase& meshCase, pACase& analysisCase, MeshAttributes& MeshAtt,
                AnalysisAttributes& AnalysisAtt, std::unordered_map<pGRegion, int> groupMap) {

    logInfo() << "Setting cases";
    // ------------------------------ Set boundary conditions
    // ------------------------------

    // Create a new manager. The manager is responsible for creating
    // new attributes, saving/retrieving to/from file, and to look up
    // cases, AttModels etc.
    pAManager attMngr = AMAN_new();

    // Create a case. A case serves as a grouping mechanism. It will (should)
    // contain all the attributes that make up a complete analysis of
    // some kind, although the system has no way of verifying this.
    analysisCase = AMAN_newCase(attMngr, "analysis", "", (pModel)model);

    // Now we create an attribute information node. This can be seen
    // as an "attribute generator object", as we will later on create an
    // actual attribute for each geometric face from this attribute
    // information node. We name the attribute T1, and give it the
    // information type "boundaryCondition".

    // With this code we can tag any BC, not only 1,3,5
    int numFaces = GM_numFaces(model);
    std::set<int> UniquefaceBound;
    for (auto const& f : AnalysisAtt.faceBound) {
      UniquefaceBound.insert(f.bcType);
    }
    int numBC = UniquefaceBound.size();
    std::map<int, pModelAssoc> aBC;

    for (auto const& ufb : UniquefaceBound) {
      std::string strboundaryCondition = "BC";
      std::ostringstream stringstream;
      stringstream << std::setw(9) << std::setfill('0') << ufb;
      std::string sNumber = stringstream.str();

      pAttInfoVoid iBC = AMAN_newAttInfoVoid(attMngr, "BC", "boundaryCondition");

      strboundaryCondition.append(sNumber);
      AttNode_setImageClass((pANode)iBC, strboundaryCondition.c_str());
      AttCase_addNode(analysisCase, (pANode)iBC);
      aBC[ufb] = AttCase_newModelAssoc(analysisCase, (pANode)iBC);
    }

    for (auto const& fb : AnalysisAtt.faceBound) {
      // Get the face
      pGEntity face = GM_entityByTag(model, 2, fb.faceID + 1);
      if (face == nullptr) {
        logError() << "faceid:" << fb.faceID + 1 << "not found in model.";
      } else {
        // Add the face to the model association. Note that we passed
        // the Attribute Information Node into the Model Association
        // at the time when the Model Association was created. That prepares
        // the creation of the AttributeVoid on the face as soon as the
        // association process is started
        if (fb.bcType != 0) {
          logInfo() << "faceBound[" << fb.faceID + 1 << "] =" << fb.bcType;
          AMA_addGEntity(aBC[fb.bcType], face);
        }
      }
    }

    // ------------------------------ Set meshing parameters
    // ------------------------------

    meshCase = MS_newMeshCase(model);

    // Set global mesh size
    pModelItem modelDomain = GM_domain(model);
    if (MeshAtt.globalMSize > 0) {
      logInfo() << "globalMSize =" << MeshAtt.globalMSize;
      // ( <meshing case>, <entity>, <1=absolute, 2=relative>, <size>, <size
      // expression> )
      MS_setMeshSize(meshCase, modelDomain, 1, MeshAtt.globalMSize, nullptr);
    }

    // Set mesh size on vertices, edges, surfaces and regions
    setMeshSize(model, meshCase, MeshAtt);

    if (MeshAtt.velocityAwareRefinementSettings.isVelocityAwareRefinementOn()) {
      logInfo() << "Enabling velocity aware meshing";
      easiMeshSize = EasiMeshSize(MeshAtt.velocityAwareRefinementSettings, model, groupMap);
      auto easiMeshSizeFunc = [](pSizeAttData sadata, void* userdata) {
        auto* easiMeshSize = static_cast<EasiMeshSize*>(userdata);
        std::array<double, 3> pt{};
        int haspt = SizeAttData_point(sadata, pt.data());
        if (!haspt) {
          logError() << "!haspt";
        }
        return easiMeshSize->getMeshSize(pt);
      };
      // set the user-defined function for isotropic size
      MS_setSizeAttFunc(meshCase, "setCustomMeshSize", easiMeshSizeFunc, &easiMeshSize);
      // Relative anisotropic size is set for the entire model through the
      // function anisoSize
      MS_setMeshSize(meshCase, GM_domain(model), MS_userDefinedType | 1, 0, "setCustomMeshSize");
    }

    if (MeshAtt.gradation > 0) {
      // Set gradation relative
      logInfo() << "Gradation rate =" << MeshAtt.gradation;
      MS_setGlobalSizeGradationRate(meshCase, MeshAtt.gradation);
    }
    if (MeshAtt.vol_AspectRatio > 0) {
      // Set target equivolume AspectRatio
      logInfo() << "Target equivolume AspectRatio =" << MeshAtt.vol_AspectRatio;
      MS_setVolumeShapeMetric(meshCase, modelDomain, ShapeMetricType_AspectRatio,
                              MeshAtt.vol_AspectRatio);
    }
    if (MeshAtt.area_AspectRatio > 0) {
      // Set target equiarea AspectRatio
      logInfo() << "Target equiarea AspectRatio =" << MeshAtt.area_AspectRatio;
#ifdef BEFORE_SIM_11
      MS_setSurfaceShapeMetric(meshCase, modelDomain, ShapeMetricType_AspectRatio,
                               MeshAtt.area_AspectRatio);
#else
      MS_setSurfaceShapeMetric(meshCase, modelDomain, ShapeMetricType_AspectRatio,
                               ElementType_Triangle, MeshAtt.area_AspectRatio);
#endif
    }
    for (auto& mycube : MeshAtt.lCube) {
      logInfo() << "Cube mesh refinement: " << mycube.CubeMSize;
      logInfo() << "Center" << mycube.CubeCenter[0] << " " << mycube.CubeCenter[1] << " "
                << mycube.CubeCenter[2];
      logInfo() << "Width" << mycube.CubeWidth[0] << " " << mycube.CubeWidth[1] << " "
                << mycube.CubeWidth[2];
      logInfo() << "Height" << mycube.CubeHeight[0] << " " << mycube.CubeHeight[1] << " "
                << mycube.CubeHeight[2];
      logInfo() << "Depth" << mycube.CubeDepth[0] << " " << mycube.CubeDepth[1] << " "
                << mycube.CubeDepth[2];
      MS_addCubeRefinement(meshCase, mycube.CubeMSize, &mycube.CubeCenter[0], &mycube.CubeWidth[0],
                           &mycube.CubeHeight[0], &mycube.CubeDepth[0]);
    }

    if (MeshAtt.MeshSizePropagationDistance > 0.0) {
      for (auto& iElem : MeshAtt.lFaceIdMeshSizePropagation) {
        pGEntity face = GM_entityByTag(model, 2, iElem);
        if (face == nullptr) {
          logError() << "MeshSizeProp faceid:" << iElem << "not found in model.";
        } else {
          logInfo() << "MeshSizeProp faceid:" << iElem
                    << ", distance =" << MeshAtt.MeshSizePropagationDistance
                    << ", scaling factor =" << MeshAtt.MeshSizePropagationScalingFactor;
#ifdef BEFORE_SIM_15
          MS_setMeshSizePropagation(meshCase, face, 1, MeshAtt.MeshSizePropagationDistance,
                                    MeshAtt.MeshSizePropagationScalingFactor);
#else
          MS_setMeshSizePropagation(meshCase, face, 2, 1, MeshAtt.MeshSizePropagationDistance,
                                    MeshAtt.MeshSizePropagationScalingFactor);
#endif
        }
      }
    }
    for (auto& iElem : MeshAtt.lFaceIdUseDiscreteMesh) {
      pGEntity face = GM_entityByTag(model, 2, iElem);
      if (face == nullptr) {
        logError() << "UseDiscreteMesh; faceid:" << iElem << "not found in model.";
      } else {
        logInfo() << "UseDiscreteMesh; faceid, noModification:" << iElem
                  << MeshAtt.UseDiscreteMesh_noModification;
        MS_useDiscreteGeometryMesh(meshCase, face, MeshAtt.UseDiscreteMesh_noModification);
        // MS_limitSurfaceMeshModification(meshCase,face,UseDiscreteMesh_noModification);
      }
    }
    for (auto& iElem : MeshAtt.lFaceIdNoMesh) {
      pGEntity face = GM_entityByTag(model, 2, iElem);
      if (face == nullptr) {
        logError() << "No Mesh; faceid:" << iElem << "not found in model.";
      } else {
        logInfo() << "No Mesh; faceid:" << iElem;
        MS_setNoMesh(meshCase, face, 1);
      }
    }
    for (auto& iElem : MeshAtt.lRegionIdNoMesh) {
      pGEntity region = GM_entityByTag(model, 3, iElem);
      if (region == nullptr) {
        logError() << "No Mesh; regionid:" << iElem << "not found in model.";
      } else {
        logInfo() << "No Mesh; regionid:" << iElem;
        MS_setNoMesh(meshCase, region, 1);
      }
    }
  }

  private:
  void loadCAD(const char* modFile, const char* cadFile) {
    // Load CAD
    std::string sCadFile;
    if (cadFile)
      sCadFile = cadFile;
    else {
      sCadFile = modFile;
      utils::StringUtils::replaceLast(sCadFile, ".smd", "_nat.x_t");
    }
    pNativeModel nativeModel = nullptr;
#ifdef PARASOLID
    if (utils::Path(sCadFile).exists())
      nativeModel = ParasolidNM_createFromFile(sCadFile.c_str(), 0);
#endif
    m_model = GM_load(modFile, nativeModel, nullptr);
    nativeModel = GM_nativeModel(m_model);

    if (nativeModel)
      NM_release(nativeModel);

#ifdef BEFORE_SIM_18
    // check for model errors
    pPList modelErrors = PList_new();
    if (!GM_isValid(m_model, 1, modelErrors)) {
      void* iter = nullptr;
      while (pSimError err = static_cast<pSimError>(PList_next(modelErrors, &iter))) {
        logInfo() << "  Error code: " << SimError_code(err) << std::endl;
        logInfo() << "  Error string: " << SimError_toString(err) << std::endl;
      }
      logError() << "Input model is not valid";
    }
    PList_delete(modelErrors);
#else
    // check for model infos
    pPList modelInfos = PList_new();
    if (!GM_isValid(m_model, 1, modelInfos)) {
      void* iter = nullptr;
      while (pSimInfo info = static_cast<pSimInfo>(PList_next(modelInfos, &iter))) {
        logInfo() << "  Info code: " << SimInfo_code(info) << std::endl;
        logInfo() << "  Info string: " << SimInfo_toString(info) << std::endl;
      }
      logError() << "Input model is not valid";
    }
    PList_delete(modelInfos);
#endif

    // check for self-intersecting CAD mesh
    pProgress prog = Progress_new();
    Progress_setCallback(prog, progressHandler);
    pMesh mesh = DM_getMesh(static_cast<pDiscreteModel>(m_model));
    pPList entities = PList_new();
    if (MS_checkMeshIntersections(mesh, entities, prog) > 0) {

      const char* element_name[4] = {"vertex", "edge", "face", "region"};

      for (int i = 0; i < PList_size(entities); i += 2) {
        pEntity firstEnt = (pEntity)PList_item(entities, i);
        pEntity secondEnt = (pEntity)PList_item(entities, i + 1);
        logInfo() << "Self-intersection between" << element_name[EN_whatInType(firstEnt)]
                  << GEN_tag(EN_whatIn(firstEnt)) << "and" << element_name[EN_whatInType(secondEnt)]
                  << GEN_tag(EN_whatIn(secondEnt));
        double centroid[3], centroid2[3];
        EN_centroid(firstEnt, &centroid[0]);
        EN_centroid(secondEnt, &centroid2[0]);
        logInfo() << "centroids" << centroid[0] << centroid[1] << centroid[2] << "and "
                  << centroid2[0] << centroid2[1] << centroid2[2] << std::flush;
      }
      logError() << PList_size(entities) / 2 << "Self-intersection(s) detected in CAD mesh";
    }
    PList_delete(entities);
    Progress_delete(prog);
  }

  private:
  void analyse_mesh() {
    int num_bins = 8;
    double AR[8] = {0, 2, 4, 6, 10, 20, 40, 100};
    long int AR_vol_bins[num_bins];
    long int skew_area_bins[num_bins];
    memset(AR_vol_bins, 0, num_bins * sizeof(long int));
    memset(skew_area_bins, 0, num_bins * sizeof(long int));

    // Accumulate equivolume AspectRatio data
    RIter reg_it;
    int num_partMeshes = PM_numParts(m_simMesh);
    pRegion reg;
    double maxAR = 0, AR_global = 0, elAR;
    int kAR;
    for (int i = 0; i < num_partMeshes; i++) {
      reg_it = M_regionIter(PM_mesh(m_simMesh, i));
      while ((reg = RIter_next(reg_it))) {
        kAR = num_bins;
        for (int k = 1; k < num_bins; k++) {
          elAR = R_aspectRatio(reg);
          if (elAR < AR[k]) {
            kAR = k;
            break;
          }
        }
        AR_vol_bins[kAR - 1]++;
        maxAR = std::max(maxAR, elAR);
      }
    }
    RIter_delete(reg_it);

    // Print the statistics
    logInfo() << "AR statistics:";
    MPI_Allreduce(&maxAR, &AR_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    logInfo() << "AR max:" << AR_global;
    logInfo() << "AR (target: < ~10):";
    long int bin_global;
    for (int i = 0; i < num_bins - 1; i++) {
      MPI_Allreduce(&AR_vol_bins[i], &bin_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      logInfo() << std::fixed << std::setprecision(2) << "[" << AR[i] << "," << AR[i + 1]
                << "):" << bin_global;
    }
    MPI_Allreduce(&AR_vol_bins[num_bins - 1], &bin_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    logInfo() << std::fixed << std::setprecision(2) << "[" << AR[num_bins - 1]
              << ",inf):" << bin_global;
  }
};

#endif // PUMGEN_SRC_INPUT_SIMMODSUITEAPF_H_
