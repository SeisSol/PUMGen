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

#include <H5Ppublic.h>
#include <H5Zpublic.h>
#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <string>

#include <hdf5.h>

#include <apfMesh2.h>
#include <apfNumbering.h>
// #include <apfZoltan.h>
#include <maMesh.h>

#include "third_party/MPITraits.h"
#include "utils/args.h"
#include "utils/logger.h"

#include "input/SerialMeshFile.h"
#ifdef USE_NETCDF
#include "input/NetCDFMesh.h"
#endif // USE_NETCDF
#include "input/ApfNative.h"
#ifdef USE_SIMMOD
#include "input/SimModSuite.h"
#endif // USE_SIMMOD
#include "meshreader/GMSH4Parser.h"
#include "meshreader/ParallelFidapReader.h"
#include "meshreader/ParallelGMSHReader.h"
#include "meshreader/ParallelGambitReader.h"
#include "third_party/GMSH2Parser.h"

template <typename TT> static TT _checkH5Err(TT&& status, const char* file, int line) {
  if (status < 0) {
    logError() << utils::nospace << "An HDF5 error occurred (" << file << ": " << line << ")";
  }
  return std::forward<TT>(status);
}

#define checkH5Err(...) _checkH5Err(__VA_ARGS__, __FILE__, __LINE__)

static int ilog(std::size_t value, int expbase = 1) {
  int count = 0;
  while (value > 0) {
    value >>= expbase;
    ++count;
  }
  return count;
}

template <typename F> static void iterate(apf::Mesh* mesh, int dim, F&& function) {
  apf::MeshIterator* it = mesh->begin(dim);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    std::invoke(function, element);
  }
  mesh->end(it);
}

int main(int argc, char* argv[]) {
  int rank = 0;
  int processes = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  PCU_Comm_Init();

  // Parse command line arguments
  utils::Args args;
  const char* source[] = {"gambit", "fidap", "msh2", "msh4", "netcdf", "apf", "simmodsuite"};
  args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);
  args.addOption("dump", 'd', "Dump APF mesh before partitioning it", utils::Args::Required, false);
  args.addOption("model", 0, "Dump/Load a specific model file", utils::Args::Required, false);
  args.addOption("vtk", 0, "Dump mesh to VTK files", utils::Args::Required, false);

  const char* filters[] = {"none",     "scaleoffset", "deflate1", "deflate2",
                           "deflate3", "deflate4",    "deflate5", "deflate6",
                           "deflate7", "deflate8",    "deflate9"};
  args.addOption("compactify-datatypes", 0,
                 "Compress index and group data types to minimum byte size (no HDF5 filters)",
                 utils::Args::Required, false);
  args.addEnumOption("filter-enable", filters, 0,
                     "Apply HDF5 filters (i.e. compression). Disabled by default.", false);
  args.addOption("filter-chunksize", 0, "Chunksize for filters (default=4096).",
                 utils::Args::Required, false);

  args.addOption("license", 'l', "License file (only used by SimModSuite)", utils::Args::Required,
                 false);
#ifdef PARASOLID
  args.addOption("cad", 'c', "CAD file (only used by SimModSuite)", utils::Args::Required, false);
#endif
  args.addOption("mesh", 0, "Mesh attributes name (only used by SimModSuite, default: \"mesh\")",
                 utils::Args::Required, false);
  args.addOption("analysis", 0,
                 "Analysis attributes name (only used by SimModSuite, default: "
                 "\"analysis\")",
                 utils::Args::Required, false);
  args.addOption("xml", 0,
                 "Use mesh and attributes parameters from a xml file (only "
                 "used by SimModSuite)",
                 utils::Args::Required, false);
  args.addOption("analyseAR", 0, "produce an histogram of AR", utils::Args::No, false);
  const char* forces[] = {"0", "1", "2"};
  args.addEnumOption("enforce-size", forces, 0,
                     "Enforce mesh size (only used by SimModSuite, default: 0)", false);
  args.addOption("sim_log", 0, "Create SimModSuite log", utils::Args::Required, false);
  args.addAdditionalOption("input", "Input file (mesh or model)");
  args.addAdditionalOption("output", "Output parallel unstructured mesh file", false);

  if (args.parse(argc, argv, rank == 0) != utils::Args::Success)
    return 1;

  const char* inputFile = args.getAdditionalArgument<const char*>("input");

  std::string outputFile;
  if (args.isSetAdditional("output")) {
    outputFile = args.getAdditionalArgument<std::string>("output");
  } else {
    // Compute default output filename
    outputFile = inputFile;
    size_t dotPos = outputFile.find_last_of('.');
    if (dotPos != std::string::npos)
      outputFile.erase(dotPos);
    outputFile.append(".puml.h5");
  }

  bool reduceInts = args.isSet("compactify-datatypes");
  int filterEnable = args.getArgument("filter-enable", 0);
  hsize_t filterChunksize = args.getArgument<hsize_t>("filter-chunksize", 4096);
  bool applyFilters = filterEnable > 0;
  if (reduceInts) {
    logInfo(rank) << "Using compact integer types.";
  }
  if (filterEnable == 0) {
    logInfo(rank) << "No filtering enabled (contiguous storage)";
  } else {
    logInfo(rank) << "Using filtering. Chunksize:" << filterChunksize;
    if (filterEnable == 1) {
      logInfo(rank) << "Compression: scale-offset compression for integers (disabled for floats)";
    } else if (filterEnable < 11) {
      logInfo(rank) << "Compression: deflate, strength" << filterEnable - 1
                    << "(note: this requires HDF5 to be compiled with GZIP support; this applies "
                       "to SeisSol as well)";
    }
  }

  std::string xdmfFile = outputFile;
  if (utils::StringUtils::endsWith(outputFile, ".puml.h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".puml.h5", ".xdmf");
  }
  if (utils::StringUtils::endsWith(outputFile, ".h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".h5", ".xdmf");
  } else {
    xdmfFile.append(".xdmf");
  }

  // Create/read the mesh
  MeshInput* meshInput = 0L;
  apf::Mesh2* mesh = 0L;
  switch (args.getArgument<int>("source", 0)) {
  case 0:
    logInfo(rank) << "Using Gambit mesh";
    meshInput = new SerialMeshFile<puml::ParallelGambitReader>(inputFile);
    break;
  case 1:
    logInfo(rank) << "Using Fidap mesh";
    meshInput = new SerialMeshFile<ParallelFidapReader>(inputFile);
    break;
  case 2:
    logInfo(rank) << "Using GMSH mesh format 2 (msh2) mesh";
    meshInput = new SerialMeshFile<puml::ParallelGMSHReader<tndm::GMSH2Parser>>(inputFile);
    break;
  case 3:
    logInfo(rank) << "Using GMSH mesh format 4 (msh4) mesh";
    meshInput = new SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH4Parser>>(inputFile);
    break;
  case 4:
#ifdef USE_NETCDF
    logInfo(rank) << "Using netCDF mesh";
    meshInput = new NetCDFMesh(inputFile);
#else  // USE_NETCDF
    logError() << "netCDF is not supported in this version";
#endif // USE_NETCDF
    break;
  case 5:
    logInfo(rank) << "Using APF native format";
    meshInput = new ApfNative(inputFile, args.getArgument<const char*>("model", 0L));
    break;
  case 6:
#ifdef USE_SIMMOD
    logInfo(rank) << "Using SimModSuite";

    meshInput = new SimModSuite(
        inputFile, args.getArgument<const char*>("cad", 0L),
        args.getArgument<const char*>("license", 0L), args.getArgument<const char*>("mesh", "mesh"),
        args.getArgument<const char*>("analysis", "analysis"),
        args.getArgument<int>("enforce-size", 0), args.getArgument<const char*>("xml", 0L),
        args.isSet("analyseAR"), args.getArgument<const char*>("sim_log", 0L));
#else  // USE_SIMMOD
    logError() << "SimModSuite is not supported in this version";
#endif // USE_SIMMOD
    break;
  default:
    logError() << "Unknown source";
  }

  mesh = meshInput->getMesh();

  // Check mesh
  if (alignMdsMatches(mesh))
    logWarning() << "Fixed misaligned matches";
  mesh->verify();

  // Dump mesh for later usage?
  const char* dumpFile = args.getArgument<const char*>("dump", 0L);
  if (dumpFile) {
    logInfo(PCU_Comm_Self()) << "Writing native APF mesh";
    mesh->writeNative(dumpFile);

    const char* modelFile = args.getArgument<const char*>("model", 0L);
    if (modelFile)
      gmi_write_dmg(mesh->getModel(), modelFile);
  }

  // Dump VTK mesh
  const char* vtkPrefix = args.getArgument<const char*>("vtk", 0L);
  if (vtkPrefix) {
    logInfo(PCU_Comm_Self()) << "Writing VTK mesh";
    apf::writeVtkFiles(vtkPrefix, mesh);
  }

  apf::Sharing* sharing = apf::getSharing(mesh);

  // oriented at the apf::countOwned method ... But with size_t
  auto countOwnedLong = [sharing, mesh](int dim) {
    std::size_t counter = 0;

    apf::MeshIterator* it = mesh->begin(dim);
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      if (sharing->isOwned(element)) {
        ++counter;
      }
    }
    mesh->end(it);

    return counter;
  };

  // TODO(David): replace by apf::countOwned again once extended

  // Get local/global size
  std::size_t localSize[2] = {countOwnedLong(3), countOwnedLong(0)};
  std::size_t globalSize[2] = {localSize[0], localSize[1]};
  MPI_Allreduce(MPI_IN_PLACE, globalSize, 2, tndm::mpi_type_t<std::size_t>(), MPI_SUM,
                MPI_COMM_WORLD);

  logInfo(rank) << "Mesh size:" << globalSize[0];

  // Compute min insphere radius
  double min = std::numeric_limits<double>::max();
  apf::MeshIterator* it = mesh->begin(3);
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    min = std::min(min, ma::getInsphere(mesh, element));
  }
  mesh->end(it);
  MPI_Reduce((rank == 0 ? MPI_IN_PLACE : &min), &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  logInfo(rank) << "Minimum insphere found:" << min;

  // Get offsets
  std::size_t offsets[2] = {localSize[0], localSize[1]};
  MPI_Scan(MPI_IN_PLACE, offsets, 2, tndm::mpi_type_t<std::size_t>(), MPI_SUM, MPI_COMM_WORLD);
  offsets[0] -= localSize[0];
  offsets[1] -= localSize[1];

  // Create numbering for the vertices/elements
  apf::GlobalNumbering* vertexNum = apf::makeGlobal(apf::numberOwnedNodes(mesh, "vertices"));
  apf::synchronize(vertexNum);

  // Create the H5 file
  hid_t h5falist = H5Pcreate(H5P_FILE_ACCESS);
  checkH5Err(h5falist);
#ifdef H5F_LIBVER_V18
  checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
  checkH5Err(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
  checkH5Err(H5Pset_fapl_mpio(h5falist, MPI_COMM_WORLD, MPI_INFO_NULL));
  hid_t h5file = H5Fcreate(outputFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5falist);
  checkH5Err(h5file);
  checkH5Err(H5Pclose(h5falist));

  hid_t h5dxlist = H5Pcreate(H5P_DATASET_XFER);
  checkH5Err(h5dxlist);
  checkH5Err(H5Pset_dxpl_mpio(h5dxlist, H5FD_MPIO_COLLECTIVE));

  // Write cells
  std::size_t connectSize = 8;
  logInfo(rank) << "Writing cells";
  {
    hsize_t sizes[2] = {globalSize[0], 4};
    hid_t h5space = H5Screate_simple(2, sizes, 0L);
    checkH5Err(h5space);

    hid_t connectType = H5T_STD_U64LE;
    if (reduceInts) {
      connectType = H5Tcopy(H5T_STD_U64LE);
      std::size_t bits = ilog(globalSize[0]);
      std::size_t unsignedSize = (bits + 7) / 8;
      checkH5Err(connectType);
      checkH5Err(H5Tset_size(connectType, unsignedSize));
      checkH5Err(
          H5Tcommit(h5file, "/connectType", connectType, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    hid_t connectFilter = H5P_DEFAULT;
    if (applyFilters) {
      connectFilter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
      hsize_t chunk[2] = {std::min(filterChunksize, sizes[0]), 4};
      checkH5Err(H5Pset_chunk(connectFilter, 2, chunk));
      if (filterEnable == 1) {
        checkH5Err(H5Pset_scaleoffset(connectFilter, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT));
      } else if (filterEnable < 11) {
        int deflateStrength = filterEnable - 1;
        checkH5Err(H5Pset_deflate(connectFilter, deflateStrength));
      }
    }

    hid_t h5connect = H5Dcreate(h5file, "/connect", connectType, h5space, H5P_DEFAULT,
                                connectFilter, H5P_DEFAULT);
    checkH5Err(h5connect);

    hsize_t start[2] = {offsets[0], 0};
    hsize_t count[2] = {localSize[0], 4};
    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

    sizes[0] = localSize[0];
    hid_t h5memspace = H5Screate_simple(2, sizes, 0L);
    checkH5Err(h5memspace);

    std::vector<unsigned long> connect(localSize[0] * 4);
    it = mesh->begin(3);
    std::size_t index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      apf::NewArray<long> vn;
      apf::getElementNumbers(vertexNum, element, vn);

      for (int i = 0; i < 4; i++) {
        connect[index * 4 + i] = vn[i];
      }

      index++;
    }
    mesh->end(it);

    checkH5Err(
        H5Dwrite(h5connect, H5T_NATIVE_ULONG, h5memspace, h5space, h5dxlist, connect.data()));

    if (applyFilters) {
      checkH5Err(H5Pclose(connectFilter));
    }
    if (reduceInts) {
      checkH5Err(H5Tclose(connectType));
    }
    checkH5Err(H5Sclose(h5space));
    checkH5Err(H5Sclose(h5memspace));
    checkH5Err(H5Dclose(h5connect));
  }

  // Vertices
  logInfo(rank) << "Writing vertices";
  {
    hsize_t sizes[2] = {0, 0};
    hsize_t start[2] = {0, 0};
    hsize_t count[2] = {0, 0};
    sizes[0] = globalSize[1];
    sizes[1] = 3;
    hid_t h5space = H5Screate_simple(2, sizes, 0L);
    checkH5Err(h5space);

    hid_t geometryFilter = H5P_DEFAULT;
    if (applyFilters && filterEnable > 1) {
      geometryFilter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
      hsize_t chunk[2] = {std::min(filterChunksize, sizes[0]), 3};
      checkH5Err(H5Pset_chunk(geometryFilter, 2, chunk));
      if (filterEnable == 1) {
        // float compression disabled at the moment (would be lossy)
        // checkH5Err(H5Pset_scaleoffset(geometryFilter, H5Z_SO_FLOAT_DSCALE,
        // H5Z_SO_INT_MINBITS_DEFAULT));
      } else if (filterEnable < 11) {
        int deflateStrength = filterEnable - 1;
        checkH5Err(H5Pset_deflate(geometryFilter, deflateStrength));
      }
    }

    hid_t h5geometry = H5Dcreate(h5file, "/geometry", H5T_IEEE_F64LE, h5space, H5P_DEFAULT,
                                 geometryFilter, H5P_DEFAULT);
    checkH5Err(h5geometry);

    start[0] = offsets[1];
    count[0] = localSize[1];
    count[1] = 3;
    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

    sizes[0] = localSize[1];
    hid_t h5memspace = H5Screate_simple(2, sizes, 0L);
    checkH5Err(h5memspace);

    std::vector<double> geometry(localSize[1] * 3);
    it = mesh->begin(0);
    std::size_t index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      if (!sharing->isOwned(element)) {
        continue;
      }

      long gid = apf::getNumber(vertexNum, apf::Node(element, 0));

      if (gid != static_cast<long>(offsets[1] + index)) {
        logError() << "Global vertex numbering is incorrect";
      }

      apf::Vector3 point;
      mesh->getPoint(element, 0, point);
      point.toArray(&geometry[index * 3]);

      index++;
    }
    mesh->end(it);

    checkH5Err(
        H5Dwrite(h5geometry, H5T_NATIVE_DOUBLE, h5memspace, h5space, h5dxlist, geometry.data()));

    if (applyFilters && filterEnable > 1) {
      checkH5Err(H5Pclose(geometryFilter));
    }

    checkH5Err(H5Sclose(h5space));
    checkH5Err(H5Sclose(h5memspace));
    checkH5Err(H5Dclose(h5geometry));
  }

  // Group information
  apf::MeshTag* groupTag = mesh->findTag("group");
  std::size_t groupSize = 4;
  if (groupTag) {
    logInfo(rank) << "Writing group information";

    hsize_t sizes[1] = {0};
    hsize_t start[1] = {0};
    hsize_t count[1] = {0};

    sizes[0] = globalSize[0];
    hid_t h5space = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5space);

    std::vector<int> group(localSize[0]);
    it = mesh->begin(3);
    std::size_t index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      assert(mesh->hasTag(element, groupTag));

      mesh->getIntTag(element, groupTag, &group[index]);
      index++;
    }
    mesh->end(it);

    hid_t groupType = H5T_STD_I32LE;
    if (reduceInts) {
      groupType = H5Tcopy(H5T_STD_I32LE);
      uint32_t minvalue = std::max(-(*std::min_element(group.begin(), group.end()) + 1), 0);
      uint32_t maxvalue = *std::max_element(group.begin(), group.end());
      uint32_t totalmax = std::max(minvalue, maxvalue);
      MPI_Allreduce(MPI_IN_PLACE, &totalmax, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
      std::size_t bits = ilog(totalmax);
      std::size_t signedSize = (bits + 7 + 1) / 8;
      checkH5Err(groupType);
      checkH5Err(H5Tset_size(groupType, signedSize));
      checkH5Err(H5Tcommit(h5file, "/groupType", groupType, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    hid_t groupFilter = H5P_DEFAULT;
    if (applyFilters) {
      groupFilter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
      hsize_t chunk[1] = {std::min(filterChunksize, sizes[0])};
      checkH5Err(H5Pset_chunk(groupFilter, 1, chunk));
      if (filterEnable == 1) {
        checkH5Err(H5Pset_scaleoffset(groupFilter, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT));
      } else if (filterEnable < 11) {
        int deflateStrength = filterEnable - 1;
        checkH5Err(H5Pset_deflate(groupFilter, deflateStrength));
      }
    }

    hid_t h5group =
        H5Dcreate(h5file, "/group", H5T_STD_I32LE, h5space, H5P_DEFAULT, groupFilter, H5P_DEFAULT);
    checkH5Err(h5group);

    start[0] = offsets[0];
    count[0] = localSize[0];
    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

    sizes[0] = localSize[0];
    hid_t h5memspace = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5memspace);

    checkH5Err(H5Dwrite(h5group, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, group.data()));

    if (applyFilters) {
      checkH5Err(H5Pclose(groupFilter));
    }
    if (reduceInts) {
      checkH5Err(H5Tclose(groupType));
    }
    checkH5Err(H5Sclose(h5space));
    checkH5Err(H5Sclose(h5memspace));
    checkH5Err(H5Dclose(h5group));
  } else {
    logInfo() << "No group information found in mesh";
  }

  // Write boundary condition
  logInfo(rank) << "Writing boundary condition";
  {
    apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
    assert(boundaryTag);

    hsize_t sizes[1] = {0};
    hsize_t start[1] = {0};
    hsize_t count[1] = {0};
    sizes[0] = globalSize[0];
    hid_t h5space = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5space);

    hid_t boundaryFilter = H5P_DEFAULT;
    if (applyFilters) {
      boundaryFilter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
      hsize_t chunk[1] = {std::min(filterChunksize, sizes[0])};
      checkH5Err(H5Pset_chunk(boundaryFilter, 1, chunk));
      if (filterEnable == 1) {
        checkH5Err(H5Pset_scaleoffset(boundaryFilter, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT));
      } else if (filterEnable < 11) {
        int deflateStrength = filterEnable - 1;
        checkH5Err(H5Pset_deflate(boundaryFilter, deflateStrength));
      }
    }

    hid_t h5boundary = H5Dcreate(h5file, "/boundary", H5T_STD_I32LE, h5space, H5P_DEFAULT,
                                 boundaryFilter, H5P_DEFAULT);
    checkH5Err(h5boundary);

    start[0] = offsets[0];
    count[0] = localSize[0];
    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

    sizes[0] = localSize[0];
    hid_t h5memspace = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5memspace);

    std::vector<int> boundary(localSize[0]);

    it = mesh->begin(3);
    std::size_t index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      apf::Downward faces;
      mesh->getDownward(element, 2, faces);

      for (int i = 0; i < 4; i++) {
        if (mesh->hasTag(faces[i], boundaryTag)) {
          int b;
          mesh->getIntTag(faces[i], boundaryTag, &b);

          if (b <= 0 || b > std::numeric_limits<unsigned char>::max())
            logError() << "Cannot handle boundary condition" << b;

          boundary[index] += b << (i * 8);
        }
      }

      index++;
    }
    mesh->end(it);

    checkH5Err(
        H5Dwrite(h5boundary, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, boundary.data()));

    if (boundaryFilter) {
      checkH5Err(H5Pclose(boundaryFilter));
    }
    checkH5Err(H5Sclose(h5space));
    checkH5Err(H5Sclose(h5memspace));
    checkH5Err(H5Dclose(h5boundary));
  }

  checkH5Err(H5Pclose(h5dxlist));

  // Writing XDMF file
  if (rank == 0) {
    logInfo() << "Writing XDMF file";

    // Strip all pathes from the filename
    std::string basename = outputFile;
    size_t pathPos = basename.find_last_of('/');
    if (pathPos != std::string::npos)
      basename = basename.substr(pathPos + 1);

    std::ofstream xdmf(xdmfFile.c_str());

    xdmf << "<?xml version=\"1.0\" ?>" << std::endl
         << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
         << "<Xdmf Version=\"2.0\">" << std::endl
         << " <Domain>" << std::endl
         << "  <Grid Name=\"puml mesh\" GridType=\"Uniform\">" << std::endl
         << "   <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"" << globalSize[0]
         << "\">"
         << std::endl
         // This should be UInt but for some reason this does not work with
         // binary data
         << "    <DataItem NumberType=\"Int\" Precision=\"" << connectSize
         << "\" Format=\"HDF\" "
            "Dimensions=\""
         << globalSize[0] << " 4\">" << basename << ":/connect</DataItem>" << std::endl
         << "   </Topology>" << std::endl
         << "   <Geometry name=\"geo\" GeometryType=\"XYZ\" NumberOfElements=\"" << globalSize[1]
         << "\">" << std::endl
         << "    <DataItem NumberType=\"Float\" Precision=\"" << sizeof(double)
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[1] << " 3\">" << basename
         << ":/geometry</DataItem>" << std::endl
         << "   </Geometry>" << std::endl
         << "   <Attribute Name=\"group\" Center=\"Cell\">" << std::endl
         << "    <DataItem  NumberType=\"Int\" Precision=\"" << groupSize
         << "\" Format=\"HDF\" "
            "Dimensions=\""
         << globalSize[0] << "\">" << basename << ":/group</DataItem>" << std::endl
         << "   </Attribute>" << std::endl
         << "   <Attribute Name=\"boundary\" Center=\"Cell\">" << std::endl
         << "    <DataItem NumberType=\"Int\" Precision=\"4\" Format=\"HDF\" "
            "Dimensions=\""
         << globalSize[0] << "\">" << basename << ":/boundary</DataItem>" << std::endl
         << "   </Attribute>" << std::endl
         << "  </Grid>" << std::endl
         << " </Domain>" << std::endl
         << "</Xdmf>" << std::endl;
  }

  checkH5Err(H5Fclose(h5file));

  delete sharing;
  delete mesh;
  delete meshInput;

  logInfo(rank) << "Finished successfully";

  PCU_Comm_Free();

  MPI_Finalize();
  return 0;
}
