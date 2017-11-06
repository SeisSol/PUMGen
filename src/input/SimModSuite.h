/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2017 Technical University of Munich
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 */

#ifndef SIM_MOD_SUITE_H
#define SIM_MOD_SUITE_H

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

#include <SimModel.h>
#include <SimDiscrete.h>
#include <SimParasolidKrnl.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimError.h>
#include <SimErrorCodes.h>
#include <SimMeshingErrorCodes.h>
#include <SimModelerUtil.h>

#include "utils/logger.h"
#include "utils/path.h"
#include "utils/progress.h"

#include "MeshInput.h"



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <SimMeshTools.h>
#include <SimDisplay.h>
#include "tinyxml2/tinyxml2.h"
#include <list>
#include <SimExport.h>

//forward declare
pAManager SModel_attManager(pModel model);

/**
 * @todo Currently it is not supported to create more than one instance
 *  of this class
 * @todo Maybe add MS_setMaxEntities to limit the number of elements
 */
class SimModSuite : public MeshInput
{
private:
	pGModel m_model;

	pParMesh m_simMesh;

	/** Enable Simmetrix logging file */
	bool m_log;
        int surfaceSmoothingLevel=2;
        int surfaceSmoothingType=1;
        double surfaceFaceRotationLimit=5.0;
        int surfaceSnap=0;
        int volumeSmoothingLevel=1;
        int volumeSmoothingType=1;
        int VolumeMesherOptimization=1;

public:
	SimModSuite(const char* modFile, const char* cadFile = 0L,
			const char* licenseFile = 0L,
			const char* meshCaseName = "mesh",
			const char* analysisCaseName = "analysis",
			int enforceSize = 0,
         const char* stl_ParFile=0L,
         const char* export_sxp_file=0L,
         const bool probe_faces=false,
         const bool analyseAR=false,
			const char* logFile = 0L)
	{
		// Init SimModSuite
		SimPartitionedMesh_start(0L, 0L);
		if (logFile) {
			m_log = true;
			Sim_logOn(logFile);
		} else
			m_log = false;
		Sim_readLicenseFile(licenseFile);
		MS_init();
		SimDiscrete_start(0);
		SimParasolid_start(1);
		Sim_setMessageHandler(messageHandler);

		// Load CAD
		logInfo(PMU_rank()) << "Loading model";

      std::string smodFile = modFile;
      if (cadFile!=0L) {
         loadCAD(modFile, cadFile);
      } else if(smodFile.substr(smodFile.find_last_of(".") + 1) == "smd") {
         loadCAD(modFile, cadFile);
      } else {
         loadSTL(modFile);
      }


      // Probe faces
      if(probe_faces) {
         probeFaceCoords(m_model);
      }

      // Extract cases
      logInfo(PMU_rank()) << "Extracting cases";
      pACase meshCase, analysisCase;
      if(stl_ParFile != NULL) {
        setCases(m_model, meshCase, analysisCase, stl_ParFile);
      } else {
        extractCases(m_model, meshCase, meshCaseName, analysisCase, analysisCaseName);
      }

		//if (nativeModel)
			m_simMesh = PM_new(0, m_model, PMU_size());
		//else
			// Discrete model
			//m_simMesh = PM_new(0, m_model, 1);

		pProgress prog = Progress_new();
		Progress_setCallback(prog, progressHandler);

		// create the mesh
		logInfo(PMU_rank()) << "Starting the surface mesher";
		pSurfaceMesher surfaceMesher = SurfaceMesher_new(meshCase, m_simMesh);
                if(stl_ParFile != NULL) {
                   SurfaceMesher_setSmoothing(surfaceMesher, surfaceSmoothingLevel);
                   SurfaceMesher_setSmoothType(surfaceMesher, surfaceSmoothingType);
                   SurfaceMesher_setFaceRotationLimit(surfaceMesher, surfaceFaceRotationLimit);
                   SurfaceMesher_setSnapForDiscrete(surfaceMesher,surfaceSnap);
                }
		progressBar.setTotal(26);
		SurfaceMesher_execute(surfaceMesher, prog);
		SurfaceMesher_delete(surfaceMesher);
      if(export_sxp_file != NULL) {
          SimExport_start();

          // Instantiate Exporter
          pExporter exporter = Exporter_new();

          // Load Nastran pattern file
          pCompiledPattern pattern = CompiledPattern_createFromFile(export_sxp_file);

          // Set inputs:
          //  1. mesh
          //  2. AttCase (set to null)
          pMesh mesh = M_createFromParMesh(m_simMesh,2,prog);
          Exporter_setInputs(exporter, mesh, 0);

          // Set outputs:
          //  1. output file
          //  2. output directory (set to null - will use local directory
          Exporter_setOutputs(exporter, "out.ts", "");

          // Run Exporter
          Exporter_executeCompiledPattern(exporter, pattern, prog);
          logInfo(PMU_rank()) <<"Export complete\n";

          // Release resources
          CompiledPattern_release(pattern);
          Exporter_delete(exporter);
          M_release(mesh);

         exit(0);
      }

		//if (!nativeModel)
			// Discrete model
			PM_setTotalNumParts(m_simMesh, PMU_size());

		logInfo(PMU_rank()) << "Starting the volume mesher";
		pVolumeMesher volumeMesher = VolumeMesher_new(meshCase, m_simMesh);
                if(stl_ParFile != NULL) {
                   VolumeMesher_setSmoothing(volumeMesher, volumeSmoothingLevel);
                   VolumeMesher_setSmoothType(volumeMesher, volumeSmoothingType);
                   VolumeMesher_setOptimization(volumeMesher,VolumeMesherOptimization);
                }
		VolumeMesher_setEnforceSize(volumeMesher, enforceSize);
		progressBar.setTotal(6);
		VolumeMesher_execute(volumeMesher, prog);
		VolumeMesher_delete(volumeMesher);

		Progress_delete(prog);

      if (analyseAR) {
        analyse_mesh();
      }
      
		// Convert to APF mesh
		apf::Mesh* tmpMesh = apf::createMesh(m_simMesh);
		gmi_register_sim();
		gmi_model* model = gmi_import_sim(m_model);

		logInfo(PMU_rank()) << "Converting mesh to APF";
		m_mesh = apf::createMdsMesh(model, tmpMesh);
		apf::destroyMesh(tmpMesh);

		// Set the boundary conditions from the geometric model
		AttCase_associate(analysisCase, 0L);
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

		// Create group lookup table
		std::unordered_map<pGRegion, int> groupTable;
		GRIter regionIt = GM_regionIter(m_model);
		while (pGRegion region = GRIter_next(regionIt)) {
			int id = groupTable.size() + 1;
			groupTable[region] = id;
		}

		// Set groups
		apf::MeshTag* groupTag = m_mesh->createIntTag("group", 1);
		it = m_mesh->begin(3);
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			apf::ModelEntity* modelRegion = m_mesh->toModel(element);

			pGRegion simRegion = reinterpret_cast<pGRegion>(modelRegion);
			std::unordered_map<pGRegion, int>::const_iterator i = groupTable.find(simRegion);
			if (i == groupTable.end())
				logError() << "Mesh element with unknown region found.";

			m_mesh->setIntTag(element, groupTag, &i->second);
		}
		m_mesh->end(it);

		AttCase_unassociate(analysisCase);

		// Delete cases
		MS_deleteMeshCase(meshCase);
		MS_deleteMeshCase(analysisCase);
	}

	virtual ~SimModSuite()
	{
		M_release(m_simMesh);
		// We cannot delete the model here because it is still
		// connected to the mesh
		//GM_release(m_model);

		// Finalize SimModSuite
		SimParasolid_stop(1);
		SimDiscrete_stop(0);
		MS_exit();
		Sim_unregisterAllKeys();
		if (m_log)
			Sim_logOff();
		SimPartitionedMesh_stop();
	}

private:
	pACase extractCase(pAManager attMngr, const char* name)
	{
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
	static unsigned int parseBoundary(const char* boundaryCondition)
	{
		if (strcmp(boundaryCondition, "freeSurface") == 0)
			return 1;
		if (strcmp(boundaryCondition, "dynamicRupture") == 0)
			return 3;
		if (strcmp(boundaryCondition, "absorbing") == 0)
			return 5;

		logError() << "Unknown boundary condition" << boundaryCondition;
		return -1;
	}

	static void messageHandler(int type, const char* msg)
	{
		switch (type) {
		case Sim_InfoMsg:
			// Show sim info messages as debug messages
			logDebug(PMU_rank()) << "SimModeler:" << msg;
			break;
		case Sim_DebugMsg:
			// Ignore sim debug messages
			break;
		case Sim_WarningMsg:
			logWarning(PMU_rank()) << "SimModeler:" << msg;
			break;
		case Sim_ErrorMsg:
			// Use warning because error will abort the program
			logWarning() << "SimModeler:" << msg;
			break;
		}
	}

	static void progressHandler(const char* what, int level, int startVal, int endVal, int currentVal, void *ignore)
	{
		if (PMU_rank() != 0)
			return;

		switch (level) {
		case 0:
			if (currentVal == -2)
				progressBar.update(0);
			else
				progressBar.clear();
			break;
		case 1:
			if (currentVal == 0)
				progressBar.update();
			else
				progressBar.increment();
			break;
		}

		logDebug() << what << level << startVal << endVal << currentVal;
	}



std::vector<std::string> &split(std::vector<std::string> &elems, const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
private:
void extractCases(pGModel m_model, pACase &meshCase, const char *meshCaseName, pACase &analysisCase, const char *analysisCaseName) {
    logInfo(PMU_rank()) << "Extracting cases";
    pAManager attMngr = SModel_attManager(m_model);

    MeshingOptions meshingOptions;
    meshCase = MS_newMeshCase(m_model);
    MS_processSimModelerMeshingAtts(extractCase(attMngr, meshCaseName),
meshCase, &meshingOptions); 

    analysisCase = extractCase(attMngr, analysisCaseName);
    pPList children = AttNode_children(analysisCase);
    void* iter = 0L;
    while (pANode child = static_cast<pANode>(PList_next(children, &iter)))
        AttCase_setModel(static_cast<pACase>(child), m_model);
    PList_delete(children);
}
private:
void setCases(pGModel model, pACase &meshCase, pACase &analysisCase, const char* stl_ParFile) {

    //Read settings from xml file

    std::ifstream casefile(stl_ParFile);
    std::string line,sval;
    int numFaces = GM_numFaces(model);
    int * faceBound = new int[numFaces];
    std::fill_n(faceBound, numFaces, 0);
    double globalMSize=0., faultMSize=0., gradation=0., vol_AspectRatio=0., area_AspectRatio=0.;
    std::list<double> lsurfaceMSize;
    std::list <std::list<int>> llsurfaceMSizeFaceId;
    int numRegions = GM_numRegions(model);
    std::list<double> lregionMSize;
    std::list <std::list<int>> llregionMSizeRegionId;
    double CubeMSize,CubeCenter[3],CubeWidth[3],CubeHeight[3],CubeDepth[3];
    CubeMSize=0;
    std::list<int> lFaceIdMeshSizePropagation;
    std::list<int> lFaceIdUseDiscreteMesh;
    int UseDiscreteMesh_noModification;
    double MeshSizePropagationScalingFactor,MeshSizePropagationDistance=0.0;
    int writeSmd;
    double detectSmallFeaturesThreshold=0;

    //TU 2.05.2007 use xml file for stl_ParFile
    using namespace tinyxml2;
    XMLDocument doc;
    XMLElement* child;
    if(doc.LoadFile(stl_ParFile) ==XML_SUCCESS)
    {
    XMLElement* pRoot;
    pRoot = doc.FirstChildElement("gradation");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       gradation =  atof(sval.c_str());
    }
    pRoot = doc.FirstChildElement("writeSmd");
    if(pRoot) {
       writeSmd = 1; 
    }
    pRoot = doc.FirstChildElement("globalMSize");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       globalMSize =  atof(sval.c_str());
    }
    pRoot = doc.FirstChildElement("area_AspectRatio");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       area_AspectRatio =  atof(sval.c_str());
    }
    pRoot = doc.FirstChildElement("vol_AspectRatio");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       vol_AspectRatio =  atof(sval.c_str());
    }

    pRoot = doc.FirstChildElement("detectSmallFeatures");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       detectSmallFeaturesThreshold =  atof(sval.c_str());
    }

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

    pRoot = doc.FirstChildElement("MeshRefinementZoneCube");
    if(pRoot) {
       sval = pRoot->Attribute("value");
       CubeMSize =  atof(sval.c_str());

       child = pRoot->FirstChildElement("Center");
       sval = child->Attribute("x");
       CubeCenter[0] =  atof(sval.c_str());
       sval = child->Attribute("y");
       CubeCenter[1] =  atof(sval.c_str());
       sval = child->Attribute("z");
       CubeCenter[2] =  atof(sval.c_str());

       child = pRoot->FirstChildElement("Width");
       sval = child->Attribute("x");
       CubeWidth[0] =  atof(sval.c_str());
       sval = child->Attribute("y");
       CubeWidth[1] =  atof(sval.c_str());
       sval = child->Attribute("z");
       CubeWidth[2] =  atof(sval.c_str());

       child = pRoot->FirstChildElement("Height");
       sval = child->Attribute("x");
       CubeHeight[0] =  atof(sval.c_str());
       sval = child->Attribute("y");
       CubeHeight[1] =  atof(sval.c_str());
       sval = child->Attribute("z");
       CubeHeight[2] =  atof(sval.c_str());

       child = pRoot->FirstChildElement("Depth");
       sval = child->Attribute("x");
       CubeDepth[0] =  atof(sval.c_str());
       sval = child->Attribute("y");
       CubeDepth[1] =  atof(sval.c_str());
       sval = child->Attribute("z");
       CubeDepth[2] =  atof(sval.c_str());
    }

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
   
    pRoot = doc.FirstChildElement("freeSurface");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound[std::atoi(tokens[i].c_str())-1] = 1;
       }
    }
    pRoot = doc.FirstChildElement("dynamicRupture");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound[std::atoi(tokens[i].c_str())-1] = 3;
       }
    }
    pRoot = doc.FirstChildElement("absorbing");
    if (pRoot) {
       std::vector<std::string> tokens;
       line =  pRoot-> GetText();
       split(tokens, line, ',');
       for(int i = 0; i < tokens.size(); i++) {
           faceBound[std::atoi(tokens[i].c_str())-1] = 5;
       }
    }

    } else {
        logError() << "Unable to open or to parse file" << stl_ParFile;
        return;
    }
    if (detectSmallFeaturesThreshold>0) {
       pSmallFeatureInfo smallFeats = GM_detectSmallFeatures(model,1,detectSmallFeaturesThreshold,0,0,0);
       pPList lsmallFeats=GM_getSmallFeatures(smallFeats); 
       logInfo(PMU_rank()) << "Number of small features returned: " <<PList_size(lsmallFeats);
    }
    logInfo(PMU_rank()) << "Setting cases";
    // ------------------------------ Set boundary conditions ------------------------------

    // Create a new manager. The manager is responsible for creating
    // new attributes, saving/retrieving to/from file, and to look up
    // cases, AttModels etc.
    pAManager attMngr = AMAN_new();

    // Create a case. A case serves as a grouping mechanism. It will (should)
    // contain all the attributes that make up a complete analysis of
    // some kind, although the system has no way of verifying this.
    analysisCase = AMAN_newCase(attMngr,"analysis","",(pModel)model);

    // Now we create an attribute information node. This can be seen
    // as an "attribute generator object", as we will later on create an
    // actual attribute for each geometric face from this attribute
    // information node. We name the attribute T1, and give it the
    // information type "boundaryCondition".
    pAttInfoVoid iSurf = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iSurf,"freeSurface");
    pAttInfoVoid iDynRup = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iDynRup,"dynamicRupture");
    pAttInfoVoid iAbsorb = AMAN_newAttInfoVoid(attMngr,"BC","boundaryCondition");
    AttNode_setImageClass((pANode)iAbsorb,"absorbing");

    // We need to add the Attribute Information Nodes to the case
    AttCase_addNode(analysisCase,(pANode)iSurf);
    AttCase_addNode(analysisCase,(pANode)iDynRup);
    AttCase_addNode(analysisCase,(pANode)iAbsorb);

    // We need to associate the Attribute Information Nodes with the
    // model. To do that we have to create a Model Association for the case
    pModelAssoc aSurf = AttCase_newModelAssoc(analysisCase,(pANode)iSurf);
    pModelAssoc aDynRup = AttCase_newModelAssoc(analysisCase,(pANode)iDynRup);
    pModelAssoc aAbsorb = AttCase_newModelAssoc(analysisCase,(pANode)iAbsorb);

    pGEntity face;
    for (int i = 0; i < numFaces; i++) {
        // Get the face
        face = GM_entityByTag(model, 2, i + 1);

        // Add the face to the model association. Note that we passed
        // the Attribute Information Node into the Model Association
        // at the time when the Model Association was created. That prepares
        // the creation of the AttributeVoid on the face as soon as the
        // association process is started
        if (faceBound[i]!=0)
           logInfo(PMU_rank()) << "faceBound[" << i + 1 <<"] =" << faceBound[i];
        switch (faceBound[i]) {
            case 0:
                break;
            case 1:
                AMA_addGEntity(aSurf,face);
                break;
            case 3:
                AMA_addGEntity(aDynRup,face);
                break;
            case 5:
                AMA_addGEntity(aAbsorb,face);
                break;
            default:
                logError() << "Unrecognised boundary type integer";
        }
    }

    // ------------------------------ Set meshing parameters ------------------------------

    meshCase = MS_newMeshCase(model);

    // Set global mesh size
    pModelItem modelDomain = GM_domain(model);
    if (globalMSize>0) {
       logInfo(PMU_rank()) << "globalMSize =" << globalMSize;
       // ( <meshing case>, <entity>, <1=absolute, 2=relative>, <size>, <size expression> )
       MS_setMeshSize(meshCase, modelDomain, 1, globalMSize, NULL);
    }
    
    // Set mesh size on surfaces
    std::list< std::list<int>>::iterator itr;
    std::list<double>::iterator itMeshSize;
    itMeshSize = lsurfaceMSize.begin();
    for (itr=llsurfaceMSizeFaceId.begin(); itr != llsurfaceMSizeFaceId.end(); itr++)
    {
       double surfaceMSize = *itMeshSize;
       std::list<int>tl=*itr;
       std::list<int>::iterator it;
       for (it=tl.begin(); it != tl.end(); it++)
       {
           face = GM_entityByTag(model, 2, *it);
           MS_setMeshSize(meshCase, face, 1, surfaceMSize, NULL);
           logInfo(PMU_rank()) << "faceid:"<<*it <<", surfaceMSize =" << surfaceMSize;
       }
       itMeshSize++;
    }

    // Set mesh size on region
    itMeshSize = lregionMSize.begin();
    for (itr=llregionMSizeRegionId.begin(); itr != llregionMSizeRegionId.end(); itr++)
    {
       double regionMSize = *itMeshSize;
       std::list<int>tl=*itr;
       std::list<int>::iterator it;
       for (it=tl.begin(); it != tl.end(); it++)
       {
           pGRegion region = (pGRegion) GM_entityByTag(model, 3, *it);
           MS_setMeshSize(meshCase, region, 1, regionMSize, NULL);
           logInfo(PMU_rank()) << "regionid:"<<*it <<", regionMSize =" << regionMSize;
       }
       itMeshSize++;
    }
    if (gradation>0) {
       // Set gradation relative
       logInfo(PMU_rank()) << "Gradation rate =" << gradation;
       MS_setGlobalSizeGradationRate(meshCase, gradation);
    }
    if (vol_AspectRatio>0) {
       // Set target equivolume AspectRatio
       logInfo(PMU_rank()) << "Target equivolume AspectRatio =" << vol_AspectRatio;
       MS_setVolumeShapeMetric(meshCase, modelDomain, ShapeMetricType_AspectRatio, vol_AspectRatio);
    }
    if (area_AspectRatio>0) {
       // Set target equiarea AspectRatio
       logInfo(PMU_rank()) << "Target equiarea AspectRatio =" << area_AspectRatio;
       for (int i = 0; i < numFaces; i++) {
           face = GM_entityByTag(model, 2, i + 1);
           MS_setSurfaceShapeMetric(meshCase, face, ShapeMetricType_AspectRatio, area_AspectRatio);
       }
    }
    if (CubeMSize>0.0) {
       logInfo(PMU_rank()) << "Cube mesh refinement: " <<CubeMSize;
       logInfo(PMU_rank()) << "Center"<< CubeCenter[0]<<" "<<CubeCenter[1]<<" "<<CubeCenter[2];
       logInfo(PMU_rank()) << "Width"<< CubeWidth[0]<<" "<<CubeWidth[1]<<" "<<CubeWidth[2];
       logInfo(PMU_rank()) << "Height"<< CubeHeight[0]<<" "<<CubeHeight[1]<<" "<<CubeHeight[2];
       logInfo(PMU_rank()) << "Depth"<< CubeDepth[0]<<" "<<CubeDepth[1]<<" "<<CubeDepth[2];
       MS_addCubeRefinement (meshCase, CubeMSize, &CubeCenter[0], &CubeWidth[0], &CubeHeight[0], &CubeDepth[0]);
    }
    if (MeshSizePropagationDistance>0.0) {
       std::list<int>::iterator it;
       for (it=lFaceIdMeshSizePropagation.begin(); it != lFaceIdMeshSizePropagation.end(); it++)
       {
           face = GM_entityByTag(model, 2, *it);
           logInfo(PMU_rank()) << "MeshSizeProp faceid:"<<*it <<", distance =" << MeshSizePropagationDistance << ", scaling factor =" << MeshSizePropagationScalingFactor;
           MS_setMeshSizePropagation(meshCase,face,1,MeshSizePropagationDistance,MeshSizePropagationScalingFactor);
       }
    }
     if (lFaceIdUseDiscreteMesh.size()>0) {
       std::list<int>::iterator it;
       for (it=lFaceIdUseDiscreteMesh.begin(); it != lFaceIdUseDiscreteMesh.end(); it++)
       {
           face = GM_entityByTag(model, 2, *it);
           logInfo(PMU_rank()) << "UseDiscreteMesh; faceid, noModification:"<<*it, UseDiscreteMesh_noModification;
           MS_useDiscreteGeometryMesh(meshCase,face,UseDiscreteMesh_noModification);
           //MS_limitSurfaceMeshModification(meshCase,face,UseDiscreteMesh_noModification);
       }
    } 
    if (writeSmd>0) {
	    logInfo(PMU_rank()) << "writing the smd";
	    GM_setAttManager(model, attMngr);
	    GM_write(model,"model.smd",0,0); // write out the model before the mesh!
	    logInfo(PMU_rank()) << "done writing smd";
    }
}

private:
void loadSTL(const char *filename){
    pMesh mesh = M_new(0,0);
    pDiscreteModel d_model = 0;
    if(M_importFromSTLFile(mesh, filename, 0L)) { //check for error
        logError() << "Error importing file";
        M_release(mesh);
        return;
    }
    logInfo(PMU_rank()) <<"done importing stl";
    // check the input mesh for intersections
    // this call must occur before the discrete model is created
    if(MS_checkMeshIntersections(mesh, 0, 0L)) {
        logError() << "There are intersections in the input mesh";
        M_release(mesh);
        return;
    }
    logInfo(PMU_rank()) <<"done checking for intersections";

    // create the discrete model
    d_model = DM_createFromMesh(mesh, 1, 0L);
    if(!d_model) { //check for error
        logError() << "Error creating Discrete model from mesh";
        M_release(mesh);
        return;
    }
    logInfo(PMU_rank()) <<"done creating the discrete mesh";


    // define the discrete model
    //DM_findEdgesByFaceNormalsDegrees(d_model, 70, 0L);
    DM_eliminateDanglingEdges(d_model, 0L);
    logInfo(PMU_rank()) <<"done eliminating Dangling edges";

    if(DM_completeTopology(d_model, 0L)) { //check for error
        logError() << "Error completing Discrete model topology";
        M_release(mesh);
        GM_release(d_model);
        return;
    }
    logInfo(PMU_rank()) <<"done checking topology";

    // Print out information about the model
    logInfo(PMU_rank()) << "Number of model vertices: " << GM_numVertices(d_model);
    logInfo(PMU_rank()) << "Number of model edges: " << GM_numEdges(d_model);
    logInfo(PMU_rank()) << "Number of model faces: " << GM_numFaces(d_model);
    logInfo(PMU_rank()) << "Number of model regions: " << GM_numRegions(d_model);

    m_model = d_model;

    // Since we told the Discrete model to use the input mesh, we release our
    // pointer to it.  It will be fully released when the Discrete model is released.
    M_release(mesh);
}

private:
void loadCAD(const char* modFile, const char* cadFile){
    // Load CAD
    std::string sCadFile;
    if (cadFile)
        sCadFile = cadFile;
    else {
        sCadFile = modFile;
        utils::StringUtils::replaceLast(sCadFile, ".smd", "_nat.x_t");
    }
    pNativeModel nativeModel = 0L;
    if (utils::Path(sCadFile).exists())
        nativeModel = ParasolidNM_createFromFile(sCadFile.c_str(), 0);

    m_model = GM_load(modFile, nativeModel, 0L);

    if (nativeModel)
        NM_release(nativeModel);

    // check for model errors
    pPList modelErrors = PList_new();
    if (!GM_isValid(m_model, 0, modelErrors))
            // TODO print more detail about errors
            logError() << "Input model is not valid";
    PList_delete(modelErrors);
}
private:
void analyse_mesh() {
    int num_bins = 8;
    double AR[8]={0,2,4,6,10,20,40,100};
    long int AR_vol_bins[num_bins];
    long int skew_area_bins[num_bins];
    memset(AR_vol_bins, 0, num_bins*sizeof(long int));
    memset(skew_area_bins, 0, num_bins*sizeof(long int));

    // Accumulate equivolume AspectRatio data
    RIter reg_it;
    int num_partMeshes = PM_numParts(m_simMesh);
    pRegion reg;
    double maxAR=0, AR_global=0,elAR;
    int kAR;
    for(int i = 0; i < num_partMeshes; i++) {
        reg_it = M_regionIter(PM_mesh(m_simMesh, i));
        while (reg = RIter_next(reg_it)) {
            kAR=num_bins;
            for(int k = 1; k < num_bins; k++) {
               elAR=R_aspectRatio(reg);
               if (elAR<AR[k]) {
               kAR=k;
               break;
               }
            }
            AR_vol_bins[kAR-1]++;
            maxAR=std::max(maxAR,elAR);
        }
    }
    RIter_delete(reg_it);
    /*
    // Accumulate equiarea skewness data
    FIter face_it;
    pFace face;
    for(int i = 0; i < num_partMeshes; i++) {
        face_it = M_faceIter(PM_mesh(m_simMesh, i));
        while (face = FIter_next(face_it)) {
            skew_area_bins[(int)(F_equiareaSkewness(face) * num_bins)]++; // because skewness lies in [0,1]
        }
    }
    FIter_delete(face_it);
    */
    // Print the statistics
    logInfo(PMU_rank()) << "AR statistics:";
    MPI_Allreduce(&maxAR, &AR_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    logInfo(PMU_rank()) << "AR max:"<< AR_global;
    logInfo(PMU_rank()) << "AR (target: < ~10):";
    long int bin_global;
    for(int i = 0; i < num_bins-1; i++) {
        MPI_Allreduce(&AR_vol_bins[i], &bin_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        logInfo(PMU_rank()) << std::fixed << std::setprecision(2) << "[" << AR[i] << "," << AR[i+1] << "):" << bin_global;
    }
    MPI_Allreduce(&AR_vol_bins[num_bins-1], &bin_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    logInfo(PMU_rank()) << std::fixed << std::setprecision(2) << "[" << AR[num_bins-1] << ",inf):" << bin_global;
    /*
    logInfo(PMU_rank()) << "Equiarea skewness (target: < 0.8):";
    for(int i = 0; i < num_bins; i++) {
        logInfo(PMU_rank()) << std::fixed << std::setprecision(2) << "[" << i * 1.0 / num_bins << "," << (i + 1) * 1.0 / num_bins << "):" << skew_area_bins[i];
    }*/
}

private:
// Method for probing the locations of the model faces to facilitate parameter setup
void probeFaceCoords(pGModel model) {
    GFIter modelFaces;
    pGFace modelFace;
    int ID;
    pPList edgeList;  // Edges bounding a face
    pGEdge thisEdge;
    pSimPolygons poly; // tessellation of a face
    const int maxPolyPoints = 100;
    int polypoint[maxPolyPoints];   // ID of the points of a polygon
    double pntlocation[3];
    double pntnormal[3];
    modelFaces = GM_faceIter(model);
    logInfo(PMU_rank()) << "Face information:";
    while(modelFace=GFIter_next(modelFaces)) { // get the next model face

        ID = GEN_tag(modelFace);
        poly = GF_displayRep(modelFace);
        int npolys = SimPolygons_numPolys(poly);
        int npolypnts = SimPolygons_numPoints(poly);
        logInfo(PMU_rank()) << "There are" << npolys << "polygons and" << npolypnts << "points on model face" << ID << ", e.g.:";

        int j;
        for (j=0; j<1; j++) { // loop over the polygons

            int myPoints = SimPolygons_polySize(poly, j);
            SimPolygons_poly(poly, j, polypoint);

            std::stringstream polygon_str;
            polygon_str << "Polygon " << j << " has the following points:";
            int k;
            for (k=0; k<myPoints; k++)
                polygon_str << " " << polypoint[k];
            logInfo(PMU_rank()) << polygon_str.str().c_str();

            for (k=0; k<myPoints; k++) {
                int hasnorm = SimPolygons_pointData(poly, polypoint[k], pntlocation, pntnormal);
                logInfo(PMU_rank()) << " Point" << polypoint[k] << ": (" << pntlocation[0] << "," << pntlocation[1] << "," << pntlocation[2] << ")";
            }
        }
        SimPolygons_delete(poly); // cleanup
    }
    GFIter_delete(modelFaces); // cleanup
}


private:
	static utils::Progress progressBar;

};

#endif // SIM_MOD_SUITE_H
