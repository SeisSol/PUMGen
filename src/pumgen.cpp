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

#include "meshreader/GMSHBuilder.h"
#include <H5Tpublic.h>
#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <string>

#include <hdf5.h>

#ifdef USE_SCOREC
#include <PCU.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
// #include <apfZoltan.h>
#include "input/ApfNative.h"
#include <maMesh.h>
#ifdef USE_SIMMOD
#include "input/SimModSuiteApf.h"
#endif
#endif
#include <type_traits>

#include "third_party/MPITraits.h"
#include "utils/args.h"
#include "utils/logger.h"

#include "input/SerialMeshFile.h"
#ifdef USE_NETCDF
#include "input/NetCDFMesh.h"
#endif // USE_NETCDF
#ifdef USE_SIMMOD
#include "input/SimModSuite.h"
#endif // USE_SIMMOD
#include "meshreader/GMSH4Parser.h"
#include "meshreader/ParallelGMSHReader.h"
#include "meshreader/ParallelGambitReader.h"
#include "third_party/GMSH2Parser.h"

#include "helper/InsphereCalculator.h"

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

constexpr std::size_t NoSecondDim = 0;

template <typename T, typename F>
static void writeH5Data(const F& handler, hid_t h5file, const std::string& name, void* mesh,
                        int meshdim, hid_t h5memtype, hid_t h5outtype, std::size_t chunk,
                        std::size_t localSize, std::size_t globalSize, bool reduceInts,
                        int filterEnable, std::size_t filterChunksize, std::size_t secondDim) {
  const std::size_t secondSize = std::max(secondDim, static_cast<std::size_t>(1));
  const std::size_t dimensions = secondDim == 0 ? 1 : 2;
  const std::size_t chunkSize = chunk / secondSize / sizeof(T);
  const std::size_t bufferSize = std::min(localSize, chunkSize);
  std::vector<T> data(secondSize * bufferSize);

  std::size_t rounds = (localSize + chunkSize - 1) / chunkSize;

  MPI_Allreduce(MPI_IN_PLACE, &rounds, 1, tndm::mpi_type_t<decltype(rounds)>(), MPI_MAX,
                MPI_COMM_WORLD);

  hsize_t globalSizes[2] = {globalSize, secondDim};
  hid_t h5space = H5Screate_simple(dimensions, globalSizes, nullptr);

  hsize_t sizes[2] = {bufferSize, secondDim};
  hid_t h5memspace = H5Screate_simple(dimensions, sizes, nullptr);

  std::size_t offset = localSize;

  MPI_Scan(MPI_IN_PLACE, &offset, 1, tndm::mpi_type_t<decltype(offset)>(), MPI_SUM, MPI_COMM_WORLD);

  offset -= localSize;

  hsize_t start[2] = {offset, 0};
  hsize_t count[2] = {bufferSize, secondDim};

  hid_t h5dxlist = H5Pcreate(H5P_DATASET_XFER);
  checkH5Err(h5dxlist);
  checkH5Err(H5Pset_dxpl_mpio(h5dxlist, H5FD_MPIO_COLLECTIVE));

  hid_t h5type = h5outtype;
  if (reduceInts && std::is_integral_v<T>) {
    h5type = H5Tcopy(h5outtype);
    std::size_t bits = ilog(globalSize);
    std::size_t unsignedSize = (bits + 7) / 8;
    checkH5Err(h5type);
    checkH5Err(H5Tset_size(h5type, unsignedSize));
    checkH5Err(H5Tcommit(h5file, (std::string("/") + name + std::string("Type")).c_str(), h5type,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  }

  hid_t h5filter = H5P_DEFAULT;
  if (filterEnable > 0) {
    h5filter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
    hsize_t chunk[2] = {std::min(filterChunksize, bufferSize), secondDim};
    checkH5Err(H5Pset_chunk(h5filter, 2, chunk));
    if (filterEnable == 1 && std::is_integral_v<T>) {
      checkH5Err(H5Pset_scaleoffset(h5filter, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT));
    } else if (filterEnable < 11) {
      int deflateStrength = filterEnable - 1;
      checkH5Err(H5Pset_deflate(h5filter, deflateStrength));
    }
  }

  hid_t h5data = H5Dcreate(h5file, (std::string("/") + name).c_str(), h5type, h5space, H5P_DEFAULT,
                           h5filter, H5P_DEFAULT);

  std::size_t written = 0;

  hsize_t nullstart[2] = {0, 0};

  for (std::size_t i = 0; i < rounds; ++i) {
    start[0] = offset + written;
    count[0] = std::min(localSize - written, bufferSize);

    std::copy_n(handler.begin() + written * secondSize, count[0] * secondSize, data.begin());

    checkH5Err(H5Sselect_hyperslab(h5memspace, H5S_SELECT_SET, nullstart, nullptr, count, nullptr));

    checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, nullptr, count, nullptr));

    checkH5Err(H5Dwrite(h5data, h5memtype, h5memspace, h5space, h5dxlist, data.data()));

    written += count[0];
  }

  if (filterEnable > 0) {
    checkH5Err(H5Pclose(h5filter));
  }
  if (reduceInts) {
    checkH5Err(H5Tclose(h5type));
  }
  checkH5Err(H5Sclose(h5space));
  checkH5Err(H5Sclose(h5memspace));
  checkH5Err(H5Dclose(h5data));
  checkH5Err(H5Pclose(h5dxlist));
}

void addAttribute(hid_t h5file, const std::string& name, const std::string& value) {
  hid_t attrSpace = checkH5Err(H5Screate(H5S_SCALAR));
  hid_t attrType = checkH5Err(H5Tcopy(H5T_C_S1));
  checkH5Err(H5Tset_size(attrType, H5T_VARIABLE));
  hid_t attrBoundary =
      checkH5Err(H5Acreate(h5file, name.data(), attrType, attrSpace, H5P_DEFAULT, H5P_DEFAULT));
  const void* stringData = value.data();
  checkH5Err(H5Awrite(attrBoundary, attrType, &stringData));
  checkH5Err(H5Aclose(attrBoundary));
  checkH5Err(H5Sclose(attrSpace));
  checkH5Err(H5Tclose(attrType));
}

template <std::size_t Order>
using SMF2 = SerialMeshFile<puml::ParallelGMSHReader<tndm::GMSH2Parser, Order>>;
template <std::size_t Order>
using SMF4 = SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH4Parser, Order>>;

int main(int argc, char* argv[]) {
  int rank = 0;
  int processes = 1;
  int mpithreadstate = MPI_THREAD_SINGLE;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpithreadstate);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

#ifdef USE_SCOREC
  PCU_Comm_Init();
#endif

  // Parse command line arguments
  utils::Args args;
  const char* source[] = {"gambit", "msh2",        "msh4",           "netcdf",
                          "apf",    "simmodsuite", "simmodsuite-apf"};
  args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);

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
  args.addOption("chunksize", 0, "Chunksize for writing (default=1 GiB).", utils::Args::Required,
                 false);
  const char* boundarytypes[] = {"int32", "int64", "i32", "i64", "i32x4"};
  args.addEnumOption("boundarytype", boundarytypes, 0,
                     "Type for writing out boundary data (default: i32).", false);

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
  args.addOption("order", 'o', "Mesh order", utils::Args::Optional, false);

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
    if (dotPos != std::string::npos) {
      outputFile.erase(dotPos);
    }
    outputFile.append(".puml.h5");
  }

  hsize_t chunksize = args.getArgument<hsize_t>("filter-chunksize", static_cast<hsize_t>(1) << 30);

  bool reduceInts = args.isSet("compactify-datatypes");
  int filterEnable = args.getArgument("filter-enable", 0);
  hsize_t filterChunksize = args.getArgument<hsize_t>("filter-chunksize", 4096);
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

  int boundaryType = args.getArgument("boundarytype", 0);
  hid_t boundaryDatatype;
  int faceOffset;
  int secondShape = NoSecondDim;
  std::string boundaryFormatAttr = "";
  int boundaryPrecision = 0;
  std::string secondDimBoundary = "";
  if (boundaryType == 0 || boundaryType == 2) {
    boundaryDatatype = H5T_STD_I32LE;
    faceOffset = 8;
    secondShape = NoSecondDim;
    boundaryFormatAttr = "i32";
    boundaryPrecision = 4;
    logInfo(rank) << "Using 32-bit integer boundary type conditions, or 8 bit per face (i32).";
  } else if (boundaryType == 1 || boundaryType == 3) {
    boundaryDatatype = H5T_STD_I64LE;
    faceOffset = 16;
    secondShape = NoSecondDim;
    boundaryFormatAttr = "i64";
    boundaryPrecision = 8;
    logInfo(rank) << "Using 64-bit integer boundary type conditions, or 16 bit per face (i64).";
  } else if (boundaryType == 4) {
    boundaryDatatype = H5T_STD_I32LE;
    secondShape = 4;
    faceOffset = -1;
    boundaryFormatAttr = "i32x4";
    boundaryPrecision = 4;
    secondDimBoundary = " 4";
    logInfo(rank) << "Using 32-bit integer per boundary face (i32x4).";
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

  const auto meshOrder = args.getArgument<int>("order", 1);

  // Create/read the mesh
  MeshData* meshInput = nullptr;
  switch (args.getArgument<int>("source", 0)) {
  case 0:
    logInfo(rank) << "Using Gambit mesh";
    meshInput = new SerialMeshFile<puml::ParallelGambitReader>(inputFile, faceOffset);
    break;
  case 1:
    logInfo(rank) << "Using GMSH mesh format 2 (msh2) mesh";
    meshInput = puml::makePointer<MeshData, SMF2>(meshOrder, inputFile, faceOffset);
    break;
  case 2:
    logInfo(rank) << "Using GMSH mesh format 4 (msh4) mesh";
    meshInput = puml::makePointer<MeshData, SMF4>(meshOrder, inputFile, faceOffset);
    break;
  case 3:
#ifdef USE_NETCDF
    logInfo(rank) << "Using netCDF mesh";
    meshInput = new NetCDFMesh(inputFile, faceOffset);
#else  // USE_NETCDF
    logError() << "netCDF is not supported in this version";
#endif // USE_NETCDF
    break;
  case 4:
    logInfo(rank) << "Using APF native format";
#ifdef USE_SCOREC
    meshInput = new ApfNative(inputFile, faceOffset, args.getArgument<const char*>("input", 0L));
    (dynamic_cast<ApfMeshInput*>(meshInput))->generate();
#else
    logError() << "This version of PUMgen has been compiled without SCOREC. Hence, the APF format "
                  "is not available.";
#endif
    break;
  case 5:
#ifdef USE_SIMMOD
    logInfo(rank) << "Using SimModSuite";

    meshInput = new SimModSuite(
        inputFile, faceOffset, args.getArgument<const char*>("cad", 0L),
        args.getArgument<const char*>("license", 0L), args.getArgument<const char*>("mesh", "mesh"),
        args.getArgument<const char*>("analysis", "analysis"),
        args.getArgument<int>("enforce-size", 0), args.getArgument<const char*>("xml", 0L),
        args.isSet("analyseAR"), args.getArgument<const char*>("sim_log", 0L));
#else  // USE_SIMMOD
    logError() << "SimModSuite is not supported in this version.";
#endif // USE_SIMMOD
    break;
  case 6:
#ifdef USE_SIMMOD
#ifdef USE_SCOREC
    logInfo(rank) << "Using SimModSuite with APF (deprecated)";

    meshInput = new SimModSuiteApf(
        inputFile, faceOffset, args.getArgument<const char*>("cad", 0L),
        args.getArgument<const char*>("license", 0L), args.getArgument<const char*>("mesh", "mesh"),
        args.getArgument<const char*>("analysis", "analysis"),
        args.getArgument<int>("enforce-size", 0), args.getArgument<const char*>("xml", 0L),
        args.isSet("analyseAR"), args.getArgument<const char*>("sim_log", 0L));
    (dynamic_cast<ApfMeshInput*>(meshInput))->generate();
#else
    logError() << "This version of PUMgen has been compiled without SCOREC. Hence, this reader for "
                  "the SimModSuite is not available here.";
#endif
#else
    logError() << "SimModSuite is not supported in this version.";
#endif
    break;
  default:
    logError() << "Unknown source.";
  }

  logInfo(rank) << "Parsed mesh successfully, writing output...";

  void* mesh = nullptr;

  // Get local/global size
  std::size_t localSize[2] = {meshInput->cellCount(), meshInput->vertexCount()};
  std::size_t globalSize[2] = {localSize[0], localSize[1]};
  MPI_Allreduce(MPI_IN_PLACE, globalSize, 2, tndm::mpi_type_t<std::size_t>(), MPI_SUM,
                MPI_COMM_WORLD);

  logInfo(rank) << "Coordinates in vertex:" << meshInput->vertexSize();
  logInfo(rank) << "Vertices in cell:" << meshInput->cellSize();
  logInfo(rank) << "Total cell count:" << globalSize[0];
  logInfo(rank) << "Total vertex count:" << globalSize[1];

  auto inspheres =
      calculateInsphere(meshInput->connectivity(), meshInput->geometry(), MPI_COMM_WORLD);
  double min = *std::min_element(inspheres.begin(), inspheres.end());
  MPI_Reduce((rank == 0 ? MPI_IN_PLACE : &min), &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  logInfo(rank) << "Minimum insphere found:" << min;

  // Get offsets
  std::size_t offsets[2] = {localSize[0], localSize[1]};
  MPI_Scan(MPI_IN_PLACE, offsets, 2, tndm::mpi_type_t<std::size_t>(), MPI_SUM, MPI_COMM_WORLD);
  offsets[0] -= localSize[0];
  offsets[1] -= localSize[1];

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
  std::size_t connectBytesPerData = 8;
  logInfo(rank) << "Writing cells";
  writeH5Data<uint64_t>(meshInput->connectivity(), h5file, "connect", mesh, 3, H5T_NATIVE_UINT64,
                        H5T_STD_U64LE, chunksize, localSize[0], globalSize[0], reduceInts,
                        filterEnable, filterChunksize, meshInput->cellSize());

  // Vertices
  logInfo(rank) << "Writing vertices";
  writeH5Data<double>(meshInput->geometry(), h5file, "geometry", mesh, 0, H5T_IEEE_F64LE,
                      H5T_IEEE_F64LE, chunksize, localSize[1], globalSize[1], reduceInts,
                      filterEnable, filterChunksize, meshInput->vertexSize());

  // Group information

  std::size_t groupBytesPerData = 4;
  logInfo(rank) << "Writing group information";
  writeH5Data<int32_t>(meshInput->group(), h5file, "group", mesh, 3, H5T_NATIVE_INT32,
                       H5T_STD_I32LE, chunksize, localSize[0], globalSize[0], reduceInts,
                       filterEnable, filterChunksize, NoSecondDim);

  // Write boundary condition
  logInfo(rank) << "Writing boundary condition";
  if (secondShape != NoSecondDim) {
    // TODO: a bit ugly, but it works
    secondShape = meshInput->vertexSize() + 1;
  }
  if (boundaryFormatAttr == "i32") {
    writeH5Data<int32_t>(meshInput->boundary(), h5file, "boundary", mesh, 3, H5T_NATIVE_INT32,
                         boundaryDatatype, chunksize, localSize[0], globalSize[0], reduceInts,
                         filterEnable, filterChunksize, secondShape);
  } else {
    writeH5Data<int64_t>(meshInput->boundary(), h5file, "boundary", mesh, 3, H5T_NATIVE_INT64,
                         boundaryDatatype, chunksize, localSize[0], globalSize[0], reduceInts,
                         filterEnable, filterChunksize, secondShape);
  }

  addAttribute(h5file, "boundary-format", boundaryFormatAttr);

  if (meshInput->hasIdentify()) {
    logInfo(rank) << "Writing vertex topology identification";
    writeH5Data<uint64_t>(meshInput->identify(), h5file, "identify", mesh, 0, H5T_NATIVE_UINT64,
                          H5T_STD_U64LE, chunksize, localSize[1], globalSize[1], reduceInts,
                          filterEnable, filterChunksize, NoSecondDim);
    addAttribute(h5file, "topology-format", "identify-vertex");
  } else {
    addAttribute(h5file, "topology-format", "geometric");
  }

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
         << "    <DataItem NumberType=\"Int\" Precision=\"" << connectBytesPerData
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
         << "    <DataItem  NumberType=\"Int\" Precision=\"" << groupBytesPerData
         << "\" Format=\"HDF\" "
            "Dimensions=\""
         << globalSize[0] << "\">" << basename << ":/group</DataItem>" << std::endl
         << "   </Attribute>" << std::endl
         << "   <Attribute Name=\"boundary\" Center=\"Cell\">" << std::endl
         << "    <DataItem NumberType=\"Int\" Precision=\"" << boundaryPrecision
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0] << secondDimBoundary << "\">"
         << basename << ":/boundary</DataItem>" << std::endl
         << "   </Attribute>" << std::endl;
    if (meshInput->hasIdentify()) {
      xdmf << "   <Attribute Name=\"boundary\" Center=\"Node\">" << std::endl
           << "    <DataItem NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Dimensions=\""
           << globalSize[1] << "\">" << basename << ":/identify</DataItem>" << std::endl
           << "   </Attribute>" << std::endl;
    }
    xdmf << "  </Grid>" << std::endl << " </Domain>" << std::endl << "</Xdmf>" << std::endl;
  }

  checkH5Err(H5Fclose(h5file));

  delete meshInput;

  logInfo(rank) << "Finished successfully";

#ifdef USE_SCOREC
  PCU_Comm_Free();
#endif

  MPI_Finalize();
  return 0;
}
