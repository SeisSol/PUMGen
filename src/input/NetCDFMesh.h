// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#ifndef PUMGEN_SRC_INPUT_NETCDFMESH_H_
#define PUMGEN_SRC_INPUT_NETCDFMESH_H_

#include <mpi.h>

#include <netcdf.h>
#include <netcdf_par.h>

#include "MeshData.h"
#include "utils/logger.h"

#include "MeshData.h"
#include "NetCDFPartition.h"
#include "ParallelVertexFilter.h"

/**
 * Read PUMGen generated mesh files
 */
class NetCDFMesh : public FullStorageMeshData {
  public:
  virtual ~NetCDFMesh() = default;

  NetCDFMesh(const char* meshFile, int boundarySize, MPI_Comm comm = MPI_COMM_WORLD)
      : FullStorageMeshData(boundarySize) {
    int rank = 0;
    int nProcs = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nProcs);

    int ncFile;
    checkNcError(nc_open_par(meshFile, NC_MPIIO, comm, MPI_INFO_NULL, &ncFile));

    // Get number of partitions
    int ncDimPart;
    checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
    size_t nPartitions;
    checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &nPartitions));

    // Local partitions
    const std::size_t nMaxLocalPart = (nPartitions + nProcs - 1) / nProcs;
    std::size_t nLocalPart = nMaxLocalPart;
    if (nPartitions < (rank + 1) * nMaxLocalPart && nPartitions >= rank * nMaxLocalPart) {
      nLocalPart = static_cast<std::size_t>(nPartitions - rank * nMaxLocalPart);
    }

    MPI_Comm commIO;
    MPI_Comm_split(MPI_COMM_WORLD, (nLocalPart > 0 ? 0 : MPI_UNDEFINED), 0, &commIO);

    // Reopen netCDF file with correct communicator
    checkNcError(nc_close(ncFile));

    if (nLocalPart > 0) {
      checkNcError(nc_open_par(meshFile, NC_MPIIO, commIO, MPI_INFO_NULL, &ncFile));
    }

    std::size_t nElements = 0;
    std::size_t nVertices = 0;

    if (nLocalPart > 0) {
      std::vector<Partition> partitions(nLocalPart);

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
      } else {
        collectiveAccess(ncFile, ncVarElemGroup);
      }

      int ncVarVrtxSize;
      checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
      collectiveAccess(ncFile, ncVarVrtxSize);

      int ncVarVrtxCoords;
      checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
      collectiveAccess(ncFile, ncVarVrtxCoords);

      // Read elements
      logInfo() << "Reading netCDF file";
      for (std::size_t i = 0; i < nMaxLocalPart; i++) {
        std::size_t j = i % nLocalPart;

        // for now, each partition stays limited to about 2^31 maximum elements

        size_t start[3] = {j + rank * nMaxLocalPart, 0, 0};

        // Element size
        unsigned int size;
        checkNcError(nc_get_var1_uint(ncFile, ncVarElemSize, start, &size));
        partitions[j].setElemSize(size);

        size_t count[3] = {1, size, 4};

        // Elements
        checkNcError(
            nc_get_vara_int(ncFile, ncVarElemVertices, start, count, partitions[j].elements()));

        // Boundaries and group
        checkNcError(
            nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, partitions[j].boundaries()));
        if (useGroups)
          checkNcError(
              nc_get_vara_int(ncFile, ncVarElemGroup, start, count, partitions[j].groups()));

        // Vertex size
        checkNcError(nc_get_var1_uint(ncFile, ncVarVrtxSize, start, &size));
        partitions[j].setVrtxSize(size);

        // Vertices
        count[1] = size;
        count[2] = 3;

        checkNcError(
            nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count, partitions[j].vertices()));
      }

      checkNcError(nc_close(ncFile));

      for (std::size_t i = 0; i < nLocalPart; i++) {
        nElements += partitions[i].nElements();
        nVertices += partitions[i].nVertices();
      }

      // Copy to the buffer
      std::vector<std::size_t> elementsLocal(nElements * 4);
      std::vector<double> verticesLocal(nVertices * 3);

      std::size_t vertexOffset = 0;
      for (std::size_t i = 0; i < nLocalPart; i++) {
        std::copy_n(partitions[i].vertices(), partitions[i].nVertices() * 3,
                    verticesLocal.begin() + vertexOffset * 3);
        vertexOffset += partitions[i].nVertices();
      }

      logInfo() << "Running vertex filter";
      ParallelVertexFilter filter(commIO);
      filter.filter(nVertices, verticesLocal);

      nVertices = filter.numLocalVertices();

      setup(nElements, nVertices);

      std::copy_n(filter.localVertices().begin(), nVertices * 3, geometryData.begin());

      std::size_t elementOffset = 0;
      vertexOffset = 0;
      for (std::size_t i = 0; i < nLocalPart; i++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t j = 0; j < partitions[i].nElements() * 4; j++) {
          elementsLocal[elementOffset * 4 + j] = partitions[i].elements()[j] + vertexOffset;
        }

        partitions[i].convertBoundary();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t j = 0; j < partitions[i].nElements(); ++j) {
          for (int k = 0; k < 4; ++k) {
            setBoundary(elementOffset + j, k, partitions[i].boundaries()[j * 4 + k]);
          }
        }

        std::copy_n(partitions[i].groups(), partitions[i].nElements(),
                    groupData.begin() + elementOffset);

        elementOffset += partitions[i].nElements();
        vertexOffset += partitions[i].nVertices();
      }

      logInfo() << "Converting local to global vertex identifier";
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t i = 0; i < nElements * 4; i++) {
        connectivityData[i] = filter.globalIds()[elementsLocal[i]];
      }
    }
  }

  private:
  /**
   * Switch to collective access for a netCDf variable
   */
  static void collectiveAccess(int ncFile, int ncVar) {
    checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
  }

  static void checkNcError(int error) {
    if (error != NC_NOERR)
      logError() << "Error while reading netCDF file:" << nc_strerror(error);
  }
};

#endif // PUMGEN_SRC_INPUT_NETCDFMESH_H_
