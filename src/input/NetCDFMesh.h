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

#ifndef NETCDF_MESH_H
#define NETCDF_MESH_H

#include <mpi.h>

#include <netcdf.h>
#include <netcdf_par.h>

#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <gmi_null.h>

#include "utils/logger.h"

#include "MeshInput.h"
#include "NetCDFPartition.h"
#include "ParallelVertexFilter.h"

/**
 * Read PUMGen generated mesh files
 */
class NetCDFMesh : public MeshInput
{
public:
	NetCDFMesh(const char* meshFile, MPI_Comm comm = MPI_COMM_WORLD)
	{
		int rank = 0;
		int nProcs = 1;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &nProcs);

		gmi_register_null();
		gmi_model* model = gmi_load(".null");
		m_mesh = apf::makeEmptyMdsMesh(model, 3, false);

		int ncFile;
		checkNcError(nc_open_par(meshFile, NC_MPIIO,
				comm, MPI_INFO_NULL, &ncFile));

		// Get number of partitions
		int ncDimPart;
		checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
		size_t nPartitions;
		checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &nPartitions));

		// Local partitions
		unsigned int nMaxLocalPart = (nPartitions + nProcs - 1) / nProcs;
		unsigned int nLocalPart = nMaxLocalPart;
		if (nPartitions < (rank+1) * nMaxLocalPart)
			nLocalPart = std::max(0, static_cast<int>(nPartitions - rank * nMaxLocalPart));

		MPI_Comm commIO;
		MPI_Comm_split(MPI_COMM_WORLD, (nLocalPart > 0 ? 0 : MPI_UNDEFINED), 0, &commIO);

		// Reopen netCDF file with correct communicator
		checkNcError(nc_close(ncFile));

		if (nLocalPart > 0)
			checkNcError(nc_open_par(meshFile, NC_MPIIO,
					commIO, MPI_INFO_NULL, &ncFile));

		PCU_Switch_Comm(commIO);

		unsigned int nElements = 0;
		unsigned int nVertices = 0;
		int* elements = 0L;
		double* vertices = 0L;
		int* boundaries = 0L;
		int* groups = 0L;

		if (nLocalPart > 0) {
			// Create netCDF variables
			int ncVarElemSize;
			checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));
			collectiveAccess(ncFile, ncVarElemSize);

			int ncVarElemVertices;
			checkNcError(nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices));
			collectiveAccess(ncFile, ncVarElemVertices);

			int ncVarElemBoundaries;
			checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
			collectiveAccess(ncFile, ncVarElemBoundaries);

			int ncVarElemGroup;
			bool useGroups = true;
			if (nc_inq_varid(ncFile, "element_group", &ncVarElemGroup) != NC_NOERR) {
				useGroups = false;
				logWarning() << "No group found, using group 0 for all elements";
			}else{
				collectiveAccess(ncFile, ncVarElemGroup);
			}

			int ncVarVrtxSize;
			checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
			collectiveAccess(ncFile, ncVarVrtxSize);

			int ncVarVrtxCoords;
			checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
			collectiveAccess(ncFile, ncVarVrtxCoords);

			Partition* partitions = new Partition[nLocalPart];

			// Read elements
			logInfo(rank) << "Reading netCDF file";
			for (unsigned int i = 0; i < nMaxLocalPart; i++) {
				unsigned int j = i % nLocalPart;

				size_t start[3] = {j + rank*nMaxLocalPart, 0, 0};

				// Element size
				unsigned int size;
				checkNcError(nc_get_var1_uint(ncFile, ncVarElemSize, start, &size));
				partitions[j].setElemSize(size);

				size_t count[3] = {1, size, 4};

				// Elements
				checkNcError(nc_get_vara_int(ncFile, ncVarElemVertices, start, count,
						partitions[j].elements()));

				// Boundaries and group
				checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count,
						partitions[j].boundaries()));
				if (useGroups)
					checkNcError(nc_get_vara_int(ncFile, ncVarElemGroup, start, count,
						partitions[j].groups()));

				// Vertex size
				checkNcError(nc_get_var1_uint(ncFile, ncVarVrtxSize, start, &size));
				partitions[j].setVrtxSize(size);

				// Vertices
				count[1] = size;
				count[2] = 3;

				checkNcError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count,
						partitions[j].vertices()));
			}

			checkNcError(nc_close(ncFile));

			for (unsigned int i = 0; i < nLocalPart; i++) {
				nElements += partitions[i].nElements();
				nVertices += partitions[i].nVertices();
			}

			// Copy to the buffer
			unsigned int* elementsLocal = new unsigned int[nElements*4];
			elements = new int[nElements*4];
			vertices = new double[nVertices*3];

			boundaries = new int[nElements*4];
			groups = new int[nElements];

			unsigned int elementOffset = 0;
			unsigned int vertexOffset = 0;
			for (unsigned int i = 0; i < nLocalPart; i++) {
#ifdef _OPENMP
				#pragma omp parallel
#endif
				for (unsigned int j = 0; j < partitions[i].nElements()*4; j++)
					elementsLocal[elementOffset*4 + j] = partitions[i].elements()[j] + vertexOffset;

				memcpy(&vertices[vertexOffset*3], partitions[i].vertices(),
						partitions[i].nVertices()*3*sizeof(double));

				partitions[i].convertBoundary();
				memcpy(&boundaries[elementOffset*4], partitions[i].boundaries(),
						partitions[i].nElements()*4*sizeof(int));
				memcpy(&groups[elementOffset], partitions[i].groups(),
						partitions[i].nElements()*sizeof(int));

				elementOffset += partitions[i].nElements();
				vertexOffset += partitions[i].nVertices();
			}

			logInfo(rank) << "Running vertex filter";
			ParallelVertexFilter filter(commIO);
			filter.filter(nVertices, vertices);

			// Create filtered vertex list
			delete [] vertices;

			nVertices = filter.numLocalVertices();
			vertices = new double[nVertices*3];
			memcpy(vertices, filter.localVertices(), nVertices*3*sizeof(double));

			logInfo(rank) << "Converting local to global vertex identifier";
#ifdef _OPENMP
			#pragma omp parallel
#endif
			for (unsigned int i = 0; i < nElements*4; i++)
				elements[i] = filter.globalIds()[elementsLocal[i]];

			delete [] partitions;
		}

		logInfo(rank) << "Constructing the mesh";
		apf::GlobalToVert vertMap;
		apf::construct(m_mesh, elements, nElements, apf::Mesh::TET, vertMap);
		delete [] elements;

		apf::alignMdsRemotes(m_mesh);
		apf::deriveMdsModel(m_mesh);

		logInfo(rank) << "Set coordinates in APF";
		apf::setCoords(m_mesh, vertices, nVertices, vertMap);
		delete [] vertices;

		// Set boundaries
		apf::MeshTag* boundaryTag = m_mesh->createIntTag("boundary condition", 1);
		apf::MeshIterator* it = m_mesh->begin(3);
		unsigned int i = 0;
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			apf::Adjacent adjacent;
			m_mesh->getAdjacent(element, 2, adjacent);

			for (unsigned int j = 0; j < 4; j++) {
				if (!boundaries[i*4 + j])
					continue;

				m_mesh->setIntTag(adjacent[j], boundaryTag, &boundaries[i*4 + j]);
			}

			i++;
		}
		m_mesh->end(it);
		delete [] boundaries;

		// Set groups
		apf::MeshTag* groupTag = m_mesh->createIntTag("group", 1);
		it = m_mesh->begin(3);
		i = 0;
		while (apf::MeshEntity* element = m_mesh->iterate(it)) {
			m_mesh->setIntTag(element, groupTag, &groups[i]);
			i++;
		}
		m_mesh->end(it);
		delete [] groups;

		PCU_Switch_Comm(MPI_COMM_WORLD);
	}

private:
	/**
	 * Switch to collective access for a netCDf variable
	 */
	static void collectiveAccess(int ncFile, int ncVar)
	{
		checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
	}

	static void checkNcError(int error)
	{
		if (error != NC_NOERR)
			logError() << "Error while reading netCDF file:" << nc_strerror(error);
	}
};

#endif // NETCDF_MESH_H
