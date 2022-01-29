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

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>

#include <hdf5.h>

#include <apfMesh2.h>
#include <apfNumbering.h>
// #include <apfZoltan.h>
#include <maMesh.h>

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
#include "meshreader/ParallelFidapReader.h"
#include "meshreader/ParallelGMSHReader.h"
#include "meshreader/ParallelGambitReader.h"

template <typename TT> static void _checkH5Err(TT status, const char* file, int line) {
  if (status < 0)
    logError() << utils::nospace << "An HDF5 error occurred (" << file << ": " << line << ")";
}

#define checkH5Err(...) _checkH5Err(__VA_ARGS__, __FILE__, __LINE__)

int main(int argc, char* argv[]) {
  int rank = 0;
  int processes = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  PCU_Comm_Init();

  // Parse command line arguments
  utils::Args args;
  const char* source[] = {"gambit", "fidap", "msh2", "netcdf", "apf", "simmodsuite"};
  args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);
  args.addOption("dump", 'd', "Dump APF mesh before partitioning it", utils::Args::Required, false);
  args.addOption("model", 0, "Dump/Load a specific model file", utils::Args::Required, false);
  args.addOption("vtk", 0, "Dump mesh to VTK files", utils::Args::Required, false);
  args.addOption("license", 'l', "License file (only used by SimModSuite)", utils::Args::Required,
                 false);
  args.addOption("features", 'f', "SimModSuite features (available in License file,only used by SimModSuite)", utils::Args::Required,
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

  std::string xdmfFile = outputFile;
  if (utils::StringUtils::endsWith(outputFile, ".puml.h5"))
    utils::StringUtils::replaceLast(xdmfFile, ".puml.h5", ".xdmf");
  if (utils::StringUtils::endsWith(outputFile, ".h5"))
    utils::StringUtils::replaceLast(xdmfFile, ".h5", ".xdmf");
  else
    xdmfFile.append(".xdmf");

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
    meshInput = new SerialMeshFile<puml::ParallelGMSHReader>(inputFile);
    break;
  case 3:
#ifdef USE_NETCDF
    logInfo(rank) << "Using netCDF mesh";
    meshInput = new NetCDFMesh(inputFile);
#else  // USE_NETCDF
    logError() << "netCDF is not supported in this version";
#endif // USE_NETCDF
    break;
  case 4:
    logInfo(rank) << "Using APF native format";
    meshInput = new ApfNative(inputFile, args.getArgument<const char*>("model", 0L));
    break;
  case 5:
#ifdef USE_SIMMOD
    logInfo(rank) << "Using SimModSuite";

    meshInput = new SimModSuite(
        inputFile, args.getArgument<const char*>("cad", 0L),
        args.getArgument<const char*>("license", 0L), args.getArgument<const char*>("features", 0L), 
        args.getArgument<const char*>("mesh", "mesh"),
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

  // Get local/global size
  unsigned int localSize[2] = {static_cast<unsigned int>(apf::countOwned(mesh, 3)),
                               static_cast<unsigned int>(apf::countOwned(mesh, 0))};
  unsigned long globalSize[2] = {localSize[0], localSize[1]};
  MPI_Allreduce(MPI_IN_PLACE, globalSize, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

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
  unsigned long offsets[2] = {localSize[0], localSize[1]};
  MPI_Scan(MPI_IN_PLACE, offsets, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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

  // Write cells
  logInfo(rank) << "Writing cells";

  hsize_t sizes[2] = {globalSize[0], 4};
  hid_t h5space = H5Screate_simple(2, sizes, 0L);
  checkH5Err(h5space);

  hid_t h5connect =
      H5Dcreate(h5file, "/connect", H5T_STD_U64LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  checkH5Err(h5connect);

  hsize_t start[2] = {offsets[0], 0};
  hsize_t count[2] = {localSize[0], 4};
  checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

  sizes[0] = localSize[0];
  hid_t h5memspace = H5Screate_simple(2, sizes, 0L);
  checkH5Err(h5memspace);

  hid_t h5dxlist = H5Pcreate(H5P_DATASET_XFER);
  checkH5Err(h5dxlist);
  checkH5Err(H5Pset_dxpl_mpio(h5dxlist, H5FD_MPIO_COLLECTIVE));

  unsigned long* connect = new unsigned long[localSize[0] * 4];
  it = mesh->begin(3);
  unsigned int index = 0;
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    apf::NewArray<long> vn;
    apf::getElementNumbers(vertexNum, element, vn);

    for (unsigned int i = 0; i < 4; i++)
      connect[index * 4 + i] = vn[i];

    index++;
  }
  mesh->end(it);

  checkH5Err(H5Dwrite(h5connect, H5T_NATIVE_ULONG, h5memspace, h5space, h5dxlist, connect));

  checkH5Err(H5Sclose(h5space));
  checkH5Err(H5Sclose(h5memspace));
  checkH5Err(H5Dclose(h5connect));

  delete[] connect;

  // Vertices
  logInfo(rank) << "Writing vertices";

  sizes[0] = globalSize[1];
  sizes[1] = 3;
  h5space = H5Screate_simple(2, sizes, 0L);
  checkH5Err(h5space);

  hid_t h5geometry = H5Dcreate(h5file, "/geometry", H5T_IEEE_F64LE, h5space, H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT);
  checkH5Err(h5geometry);

  start[0] = offsets[1];
  count[0] = localSize[1];
  count[1] = 3;
  checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

  sizes[0] = localSize[1];
  h5memspace = H5Screate_simple(2, sizes, 0L);
  checkH5Err(h5memspace);

  apf::Sharing* sharing = apf::getSharing(mesh);

  double* geometry = new double[localSize[1] * 3];
  it = mesh->begin(0);
  index = 0;
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    if (!sharing->isOwned(element))
      continue;

    long gid = apf::getNumber(vertexNum, apf::Node(element, 0));

    if (gid != static_cast<long>(offsets[1] + index))
      logError() << "Global vertex numbering is incorrect";

    apf::Vector3 point;
    mesh->getPoint(element, 0, point);
    point.toArray(&geometry[index * 3]);

    index++;
  }
  mesh->end(it);

  checkH5Err(H5Dwrite(h5geometry, H5T_NATIVE_DOUBLE, h5memspace, h5space, h5dxlist, geometry));

  checkH5Err(H5Sclose(h5space));
  checkH5Err(H5Sclose(h5memspace));
  checkH5Err(H5Dclose(h5geometry));

  delete[] geometry;

  // Group information
  apf::MeshTag* groupTag = mesh->findTag("group");
  if (groupTag) {
    logInfo(rank) << "Writing group information";

    sizes[0] = globalSize[0];
    h5space = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5space);

    hid_t h5group =
        H5Dcreate(h5file, "/group", H5T_STD_I32LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    checkH5Err(h5group);

    start[0] = offsets[0];
    count[0] = localSize[0];
    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

    sizes[0] = localSize[0];
    h5memspace = H5Screate_simple(1, sizes, 0L);
    checkH5Err(h5memspace);

    int* group = new int[localSize[0]];
    it = mesh->begin(3);
    index = 0;
    while (apf::MeshEntity* element = mesh->iterate(it)) {
      assert(mesh->hasTag(element, groupTag));

      mesh->getIntTag(element, groupTag, &group[index]);
      index++;
    }
    mesh->end(it);

    checkH5Err(H5Dwrite(h5group, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, group));

    checkH5Err(H5Sclose(h5space));
    checkH5Err(H5Sclose(h5memspace));
    checkH5Err(H5Dclose(h5group));

    delete[] group;
  } else {
    logInfo() << "No group information found in mesh";
  }

  // Write boundary condition
  logInfo(rank) << "Writing boundary condition";
  apf::MeshTag* boundaryTag = mesh->findTag("boundary condition");
  assert(boundaryTag);

  int* boundary = new int[localSize[0]];
  memset(boundary, 0, localSize[0] * sizeof(int));

  sizes[0] = globalSize[0];
  h5space = H5Screate_simple(1, sizes, 0L);
  checkH5Err(h5space);

  hid_t h5boundary =
      H5Dcreate(h5file, "/boundary", H5T_STD_I32LE, h5space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  checkH5Err(h5boundary);

  start[0] = offsets[0];
  count[0] = localSize[0];
  checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, 0L, count, 0L));

  sizes[0] = localSize[0];
  h5memspace = H5Screate_simple(1, sizes, 0L);
  checkH5Err(h5memspace);

  it = mesh->begin(3);
  index = 0;
  while (apf::MeshEntity* element = mesh->iterate(it)) {
    apf::Downward faces;
    mesh->getDownward(element, 2, faces);

    for (unsigned int i = 0; i < 4; i++) {
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

  checkH5Err(H5Dwrite(h5boundary, H5T_NATIVE_INT, h5memspace, h5space, h5dxlist, boundary));

  checkH5Err(H5Sclose(h5space));
  checkH5Err(H5Sclose(h5memspace));
  checkH5Err(H5Dclose(h5boundary));

  delete[] boundary;

  checkH5Err(H5Pclose(h5dxlist));

  delete meshInput;

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
         << "    <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" "
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
         << "    <DataItem  NumberType=\"Int\" Precision=\"4\" Format=\"HDF\" "
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

  logInfo(rank) << "Finished successfully";

  PCU_Comm_Free();

  MPI_Finalize();
  return 0;
}
