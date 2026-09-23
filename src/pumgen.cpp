// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

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
#include <sstream>
#include <string>
#include <type_traits>

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
#include "meshreader/DistributedGMSHReader.h"
#include "meshreader/GMSH2Parser.h"
#include "meshreader/GMSH4Parser.h"
#include "meshreader/ParallelGMSHReader.h"
#include "meshreader/ParallelGambitReader.h"

#include "helper/InsphereCalculator.h"
#include "helper/IntegerWidth.h"

#ifdef USE_EASI
#include "sizing/VelocityAwareMeshSize.h"
#include "sizing/VelocityAwareSettings.h"
#include "sizing/VelocityCheck.h"
#include "tinyxml2/tinyxml2.h"
#endif

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

/**
 * The precision to declare in XDMF for integers of the given size: XDMF knows 1, 2, 4 and 8 bytes,
 * and HDF5 widens the sizes in between which the compaction may choose.
 */
static std::size_t xdmfPrecision(std::size_t bytes) {
  std::size_t precision = 1;
  while (precision < bytes) {
    precision *= 2;
  }
  return precision;
}

/**
 * The smallest size in bytes of an integer type with the given signedness which holds all values
 * of the (distributed) data, among the sizes of 1, 2, 4 and 8 bytes. Data with negative values is
 * not reduced.
 */
template <typename F>
static std::size_t requiredIntegerBytes(const F& data, std::size_t count, bool isSigned) {
  using ValueT = typename F::value_type;
  unsigned long long largest = 0;
  int negative = 0;
  for (std::size_t i = 0; i < count; ++i) {
    if constexpr (std::is_signed_v<ValueT>) {
      if (data[i] < 0) {
        negative = 1;
        continue;
      }
    }
    largest = std::max(largest, static_cast<unsigned long long>(data[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE, &largest, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &negative, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (negative != 0) {
    return sizeof(ValueT);
  }
  return compactIntegerBytes(ilog(largest) + (isSigned ? 1 : 0));
}

/**
 * Writes a (distributed) dataset; returns the size of its values in the file.
 */
template <typename T, typename F>
static std::size_t
writeH5Data(const F& handler, hid_t h5file, const std::string& name, hid_t h5memtype,
            hid_t h5outtype, std::size_t chunk, std::size_t localSize, std::size_t globalSize,
            bool reduceInts, int filterEnable, std::size_t filterChunksize, std::size_t secondDim) {
  const std::size_t secondSize = std::max(secondDim, static_cast<std::size_t>(1));
  const std::size_t dimensions = secondDim == 0 ? 1 : 2;
  // at least one row per round, however small the chunk is
  const std::size_t chunkSize = std::max<std::size_t>(1, chunk / secondSize / sizeof(T));
  const std::size_t bufferSize = std::min(localSize, chunkSize);
  // data of the memory type is written from where it is, other data is converted round by round
  constexpr bool Direct = std::is_same_v<typename F::value_type, T>;
  std::vector<T> data(Direct ? 0 : secondSize * bufferSize);

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
  if constexpr (std::is_integral_v<T>) {
    if (reduceInts) {
      const std::size_t bytes = requiredIntegerBytes(handler, localSize * secondSize,
                                                     H5Tget_sign(h5outtype) == H5T_SGN_2);
      if (bytes < H5Tget_size(h5outtype)) {
        h5type = checkH5Err(H5Tcopy(h5outtype));
        checkH5Err(H5Tset_size(h5type, bytes));
      }
    }
  }

  hid_t h5filter = H5P_DEFAULT;
  if (filterEnable > 0) {
    h5filter = checkH5Err(H5Pcreate(H5P_DATASET_CREATE));
    // the dataset creation is collective, so the chunk shape must not depend on the local size
    const hsize_t chunkRows =
        std::max(static_cast<hsize_t>(1), std::min<hsize_t>(filterChunksize, globalSize));
    hsize_t chunk[2] = {chunkRows, secondDim};
    checkH5Err(H5Pset_chunk(h5filter, dimensions, chunk));
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

    const T* source = data.data();
    if constexpr (Direct) {
      source = handler.data() + written * secondSize;
    } else {
      std::copy_n(handler.begin() + written * secondSize, count[0] * secondSize, data.begin());
    }

    if (count[0] > 0) {
      checkH5Err(
          H5Sselect_hyperslab(h5memspace, H5S_SELECT_SET, nullstart, nullptr, count, nullptr));
      checkH5Err(H5Sselect_hyperslab(h5space, H5S_SELECT_SET, start, nullptr, count, nullptr));
    } else {
      // ranks without (further) data still take part in the collective write
      checkH5Err(H5Sselect_none(h5memspace));
      checkH5Err(H5Sselect_none(h5space));
    }

    checkH5Err(H5Dwrite(h5data, h5memtype, h5memspace, h5space, h5dxlist, source));

    written += count[0];
  }

  const std::size_t valueBytes = H5Tget_size(h5type);
  if (filterEnable > 0) {
    checkH5Err(H5Pclose(h5filter));
  }
  if (h5type != h5outtype) {
    checkH5Err(H5Tclose(h5type));
  }
  checkH5Err(H5Sclose(h5space));
  checkH5Err(H5Sclose(h5memspace));
  checkH5Err(H5Dclose(h5data));
  checkH5Err(H5Pclose(h5dxlist));
  return valueBytes;
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
using SMF2 = SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH2Parser, Order>>;
template <std::size_t Order>
using SMF4 = SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH4Parser, Order>>;
template <std::size_t Order>
using SMF4Distributed = SerialMeshFile<puml::DistributedGMSHReader<Order>>;

int main(int argc, char* argv[]) {
  int rank = 0;
  int processes = 1;
  int mpithreadstate = MPI_THREAD_SINGLE;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpithreadstate);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  utils::Logger::setRank(rank);

#ifdef USE_SCOREC
  PCU_Comm_Init();
#endif

  // Parse command line arguments
  utils::Args args("A mesh converter to the PUML format, as used in SeisSol.");
  const auto source = std::vector<std::string>{"gambit", "msh2",        "msh4",           "netcdf",
                                               "apf",    "simmodsuite", "simmodsuite-apf"};
  args.addEnumOption("source", source, 's', "Mesh source (default: gambit)", false);

  const auto filters = std::vector<std::string>{"none",     "scaleoffset", "deflate1", "deflate2",
                                                "deflate3", "deflate4",    "deflate5", "deflate6",
                                                "deflate7", "deflate8",    "deflate9"};
  args.addOption("compactify-datatypes", 0,
                 "Compress index and group data types to minimum byte size (no HDF5 filters)",
                 utils::Args::No, false);
  args.addEnumOption("filter-enable", filters, 0,
                     "Apply HDF5 filters (i.e. compression). Disabled by default.", false);
  args.addOption("filter-chunksize", 0, "Chunksize for filters (default=4096).",
                 utils::Args::Required, false);
  args.addOption("chunksize", 0, "Chunksize for writing (default=1 GiB).", utils::Args::Required,
                 false);
  const auto boundarytypes = std::vector<std::string>{"int32", "int64", "i32", "i64", "i32x4"};
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
  const auto forces = std::vector<std::string>{"0", "1", "2"};
  args.addEnumOption("enforce-size", forces, 0,
                     "Enforce mesh size (only used by SimModSuite, default: 0)", false);
  args.addOption("sim_log", 0, "Create SimModSuite log", utils::Args::Required, false);
  args.addAdditionalOption("input", "Input file (mesh or model)");
  args.addAdditionalOption("output", "Output parallel unstructured mesh file", false);
  args.addOption("order", 'o', "Mesh order (default: 1)", utils::Args::Required, false);
  const auto gmshReaders = std::vector<std::string>{"auto", "serial"};
  args.addEnumOption("gmsh-reader", gmshReaders, 0,
                     "Reader for msh4 meshes (default: auto, which reads binary MSH 4.1 files "
                     "with all ranks; serial: rank 0 reads the file)",
                     false);
  args.addOption("velocity-check", 0,
                 "Compare the cell sizes with the velocity-aware mesh size of the "
                 "VelocityAwareMeshing element of a mesh attributes file (XML)",
                 utils::Args::Required, false);

  const auto parsed = args.parse(argc, argv, rank == 0);
  if (parsed != utils::Args::Success) {
    MPI_Finalize();
    return parsed == utils::Args::Help ? 0 : 1;
  }

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

  hsize_t chunksize = args.getArgument<hsize_t>("chunksize", static_cast<hsize_t>(1) << 30);

  bool reduceInts = args.isSet("compactify-datatypes");
  int filterEnable = args.getArgument("filter-enable", 0);
  hsize_t filterChunksize = args.getArgument<hsize_t>("filter-chunksize", 4096);
  if (reduceInts) {
    logInfo() << "Using compact integer types.";
  }
  if (filterEnable == 0) {
    logInfo() << "No filtering enabled (contiguous storage)";
  } else {
    logInfo() << "Using filtering. Chunksize:" << filterChunksize;
    if (filterEnable == 1) {
      logInfo() << "Compression: scale-offset compression for integers (disabled for floats)";
    } else if (filterEnable < 11) {
      logInfo() << "Compression: deflate, strength" << filterEnable - 1
                << "(note: this requires HDF5 to be compiled with GZIP support; this applies "
                   "to SeisSol as well)";
    }
  }

  int boundaryType = args.getArgument("boundarytype", 0);
  hid_t boundaryDatatype;
  int faceOffset;
  int secondShape = NoSecondDim;
  std::string boundaryFormatAttr = "";
  std::string secondDimBoundary = "";
  if (boundaryType == 0 || boundaryType == 2) {
    boundaryDatatype = H5T_STD_I32LE;
    faceOffset = 8;
    secondShape = NoSecondDim;
    boundaryFormatAttr = "i32";
    logInfo() << "Using 32-bit integer boundary type conditions, or 8 bit per face (i32).";
  } else if (boundaryType == 1 || boundaryType == 3) {
    boundaryDatatype = H5T_STD_I64LE;
    faceOffset = 16;
    secondShape = NoSecondDim;
    boundaryFormatAttr = "i64";
    logInfo() << "Using 64-bit integer boundary type conditions, or 16 bit per face (i64).";
  } else if (boundaryType == 4) {
    boundaryDatatype = H5T_STD_I32LE;
    secondShape = 4;
    faceOffset = -1;
    boundaryFormatAttr = "i32x4";
    secondDimBoundary = " 4";
    logInfo() << "Using 32-bit integer per boundary face (i32x4).";
  }

  std::string xdmfFile = outputFile;
  if (utils::StringUtils::endsWith(outputFile, ".puml.h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".puml.h5", ".xdmf");
  } else if (utils::StringUtils::endsWith(outputFile, ".h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".h5", ".xdmf");
  } else {
    xdmfFile.append(".xdmf");
  }

  const auto meshOrder = args.getArgument<int>("order", 1);
  if (meshOrder < 1) {
    logError() << "The mesh order has to be at least 1, got" << meshOrder;
  }

  // Create/read the mesh
  std::unique_ptr<MeshData> meshInput;
  switch (args.getArgument<int>("source", 0)) {
  case 0:
    logInfo() << "Using Gambit mesh";
    meshInput = std::make_unique<SerialMeshFile<puml::ParallelGambitReader>>(inputFile, faceOffset);
    break;
  case 1:
    logInfo() << "Using GMSH mesh format 2 (msh2) mesh";
    meshInput.reset(puml::makePointer<MeshData, SMF2>(meshOrder, inputFile, faceOffset));
    break;
  case 2: {
    int distributed = 0;
    if (rank == 0) {
      distributed = args.getArgument<int>("gmsh-reader", 0) == 0 && puml::isBinaryMsh4(inputFile);
    }
    MPI_Bcast(&distributed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (distributed != 0) {
      logInfo() << "Using GMSH mesh format 4 (msh4) mesh, binary, read by all ranks";
      meshInput.reset(
          puml::makePointer<MeshData, SMF4Distributed>(meshOrder, inputFile, faceOffset));
    } else {
      logInfo() << "Using GMSH mesh format 4 (msh4) mesh";
      meshInput.reset(puml::makePointer<MeshData, SMF4>(meshOrder, inputFile, faceOffset));
    }
    break;
  }
  case 3:
#ifdef USE_NETCDF
    logInfo() << "Using netCDF mesh";
    meshInput = std::make_unique<NetCDFMesh>(inputFile, faceOffset);
#else  // USE_NETCDF
    logError() << "netCDF is not supported in this version";
#endif // USE_NETCDF
    break;
  case 4:
    logInfo() << "Using APF native format";
#ifdef USE_SCOREC
    meshInput = std::make_unique<ApfNative>(inputFile, faceOffset,
                                            args.getArgument<const char*>("input", 0L));
    (dynamic_cast<ApfMeshInput*>(meshInput.get()))->generate();
#else
    logError() << "This version of PUMgen has been compiled without SCOREC. Hence, the APF format "
                  "is not available.";
#endif
    break;
  case 5:
#ifdef USE_SIMMOD
    logInfo() << "Using SimModSuite";

    meshInput = std::make_unique<SimModSuite>(
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
    logInfo() << "Using SimModSuite with APF (deprecated)";

    meshInput = std::make_unique<SimModSuiteApf>(
        inputFile, faceOffset, args.getArgument<const char*>("cad", 0L),
        args.getArgument<const char*>("license", 0L), args.getArgument<const char*>("mesh", "mesh"),
        args.getArgument<const char*>("analysis", "analysis"),
        args.getArgument<int>("enforce-size", 0), args.getArgument<const char*>("xml", 0L),
        args.isSet("analyseAR"), args.getArgument<const char*>("sim_log", 0L));
    (dynamic_cast<ApfMeshInput*>(meshInput.get()))->generate();
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

  logInfo() << "Parsed mesh successfully, writing output...";

  // Get local/global size
  std::size_t localSize[2] = {meshInput->cellCount(), meshInput->vertexCount()};
  std::size_t globalSize[2] = {localSize[0], localSize[1]};
  MPI_Allreduce(MPI_IN_PLACE, globalSize, 2, tndm::mpi_type_t<std::size_t>(), MPI_SUM,
                MPI_COMM_WORLD);

  logInfo() << "Coordinates in vertex:" << meshInput->vertexSize();
  logInfo() << "Vertices in cell:" << meshInput->cellSize();
  logInfo() << "Total cell count:" << globalSize[0];
  logInfo() << "Total vertex count:" << globalSize[1];
  logInfo() << "Vertex identification:" << meshInput->hasIdentify();

  const CellVertices cellVertices(meshInput->connectivity(), meshInput->geometry(),
                                  meshInput->cellSize(), MPI_COMM_WORLD);
  double min = std::numeric_limits<double>::infinity();
  {
    const auto inspheres = calculateInsphere(cellVertices);
    if (!inspheres.empty()) {
      min = *std::min_element(inspheres.begin(), inspheres.end());
    }
  }
  MPI_Reduce((rank == 0 ? MPI_IN_PLACE : &min), &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  logInfo() << "Minimum insphere found:" << min;

  if (args.isSet("velocity-check")) {
#ifdef USE_EASI
    const auto* settingsFile = args.getArgument<const char*>("velocity-check");
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(settingsFile) != tinyxml2::XML_SUCCESS) {
      logError() << "Could not read" << settingsFile;
    }
    const auto settings = readVelocityAwareSettings(doc);
    if (!settings.isVelocityAwareRefinementOn()) {
      logError() << "No refinement cuboid in" << settingsFile;
    }
    const VelocityAwareMeshSize meshSize(settings);
    checkVelocityAwareMeshSize(meshSize, cellVertices, meshInput->group(), MPI_COMM_WORLD);
#else
    logError() << "This version of PUMGen has been compiled without easi; hence, --velocity-check "
                  "is not available.";
#endif
  }

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
  logInfo() << "Writing cells";
  const auto connectBytes =
      writeH5Data<uint64_t>(meshInput->connectivity(), h5file, "connect", H5T_NATIVE_UINT64,
                            H5T_STD_U64LE, chunksize, localSize[0], globalSize[0], reduceInts,
                            filterEnable, filterChunksize, meshInput->cellSize());

  // Vertices
  logInfo() << "Writing vertices";
  writeH5Data<double>(meshInput->geometry(), h5file, "geometry", H5T_IEEE_F64LE, H5T_IEEE_F64LE,
                      chunksize, localSize[1], globalSize[1], reduceInts, filterEnable,
                      filterChunksize, meshInput->vertexSize());

  // Group information

  logInfo() << "Writing group information";
  const auto groupBytes = writeH5Data<int32_t>(
      meshInput->group(), h5file, "group", H5T_NATIVE_INT32, H5T_STD_I32LE, chunksize, localSize[0],
      globalSize[0], reduceInts, filterEnable, filterChunksize, NoSecondDim);

  // Write boundary condition
  logInfo() << "Writing boundary condition";
  if (secondShape != NoSecondDim) {
    // TODO: a bit ugly, but it works
    secondShape = meshInput->vertexSize() + 1;
  }
  std::size_t boundaryBytes = 0;
  if (boundaryFormatAttr == "i32") {
    boundaryBytes = writeH5Data<int32_t>(
        meshInput->boundary(), h5file, "boundary", H5T_NATIVE_INT32, boundaryDatatype, chunksize,
        localSize[0], globalSize[0], reduceInts, filterEnable, filterChunksize, secondShape);
  } else {
    boundaryBytes = writeH5Data<int64_t>(
        meshInput->boundary(), h5file, "boundary", H5T_NATIVE_INT64, boundaryDatatype, chunksize,
        localSize[0], globalSize[0], reduceInts, filterEnable, filterChunksize, secondShape);
  }

  addAttribute(h5file, "boundary-format", boundaryFormatAttr);

  std::size_t identifyBytes = 0;
  if (meshInput->hasIdentify()) {
    logInfo() << "Writing vertex topology identification";
    identifyBytes = writeH5Data<uint64_t>(
        meshInput->identify(), h5file, "identify", H5T_NATIVE_UINT64, H5T_STD_U64LE, chunksize,
        localSize[1], globalSize[1], reduceInts, filterEnable, filterChunksize, NoSecondDim);
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

    // the cells are shown with their vertices; for higher-order meshes, the other nodes follow
    // the vertices in the connectivity and are left out
    const auto cellSize = meshInput->cellSize();
    std::ostringstream connect;
    connect << "<DataItem NumberType=\"Int\" Precision=\"" << xdmfPrecision(connectBytes)
            << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0] << " " << cellSize << "\">"
            << basename << ":/connect</DataItem>";

    xdmf << "<?xml version=\"1.0\" ?>" << std::endl
         << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
         << "<Xdmf Version=\"2.0\">" << std::endl
         << " <Domain>" << std::endl
         << "  <Grid Name=\"puml mesh\" GridType=\"Uniform\">" << std::endl
         << "   <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"" << globalSize[0]
         << "\">" << std::endl;
    if (cellSize == CellVertices::VerticesPerCell) {
      xdmf << "    " << connect.str() << std::endl;
    } else {
      xdmf << "    <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << globalSize[0] << " "
           << CellVertices::VerticesPerCell << "\" Type=\"HyperSlab\">" << std::endl
           << "     <DataItem Dimensions=\"3 2\" Format=\"XML\">0 0 1 1 " << globalSize[0] << " "
           << CellVertices::VerticesPerCell << "</DataItem>" << std::endl
           << "     " << connect.str() << std::endl
           << "    </DataItem>" << std::endl;
    }
    xdmf << "   </Topology>" << std::endl
         << "   <Geometry name=\"geo\" GeometryType=\"XYZ\" NumberOfElements=\"" << globalSize[1]
         << "\">" << std::endl
         << "    <DataItem NumberType=\"Float\" Precision=\"" << sizeof(double)
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[1] << " 3\">" << basename
         << ":/geometry</DataItem>" << std::endl
         << "   </Geometry>" << std::endl
         << "   <Attribute Name=\"group\" Center=\"Cell\">" << std::endl
         << "    <DataItem  NumberType=\"Int\" Precision=\"" << xdmfPrecision(groupBytes)
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0] << "\">" << basename
         << ":/group</DataItem>" << std::endl
         << "   </Attribute>" << std::endl
         << "   <Attribute Name=\"boundary\" Center=\"Cell\">" << std::endl
         << "    <DataItem NumberType=\"Int\" Precision=\"" << xdmfPrecision(boundaryBytes)
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0] << secondDimBoundary << "\">"
         << basename << ":/boundary</DataItem>" << std::endl
         << "   </Attribute>" << std::endl;
    if (meshInput->hasIdentify()) {
      xdmf << "   <Attribute Name=\"identify\" Center=\"Node\">" << std::endl
           << "    <DataItem NumberType=\"Int\" Precision=\"" << xdmfPrecision(identifyBytes)
           << "\" Format=\"HDF\" Dimensions=\"" << globalSize[1] << "\">" << basename
           << ":/identify</DataItem>" << std::endl
           << "   </Attribute>" << std::endl;
    }
    xdmf << "  </Grid>" << std::endl << " </Domain>" << std::endl << "</Xdmf>" << std::endl;
  }

  checkH5Err(H5Fclose(h5file));

  meshInput.reset();

  logInfo() << "Finished successfully";

#ifdef USE_SCOREC
  PCU_Comm_Free();
#endif

  MPI_Finalize();
  return 0;
}
