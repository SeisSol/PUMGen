// SPDX-FileCopyrightText: 2017 SeisSol Group
// SPDX-FileCopyrightText: 2017 Technical University of Munich
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileContributor: Sebastian Rettenberger <sebastian.rettenberger@tum.de>

#include <H5Tpublic.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
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

#include "helper/BoundaryFormat.h"
#include "helper/InsphereCalculator.h"
#include "helper/IntegerWidth.h"
#include "mesh/CellType.h"

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

#ifndef PUMGEN_VERSION
#define PUMGEN_VERSION "unknown"
#endif

static void addIntegerAttribute(hid_t object, const std::string& name,
                                const std::vector<std::int64_t>& values) {
  const hsize_t count = values.size();
  hid_t space = checkH5Err(H5Screate_simple(1, &count, nullptr));
  hid_t attribute =
      checkH5Err(H5Acreate(object, name.c_str(), H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT));
  checkH5Err(H5Awrite(attribute, H5T_NATIVE_INT64, values.data()));
  checkH5Err(H5Aclose(attribute));
  checkH5Err(H5Sclose(space));
}

/**
 * A string of fixed length, as VTK expects the Type of a VTKHDF file.
 */
static void addFixedStringAttribute(hid_t object, const std::string& name,
                                    const std::string& value) {
  hid_t space = checkH5Err(H5Screate(H5S_SCALAR));
  hid_t type = checkH5Err(H5Tcopy(H5T_C_S1));
  checkH5Err(H5Tset_size(type, value.size()));
  checkH5Err(H5Tset_strpad(type, H5T_STR_NULLPAD));
  hid_t attribute =
      checkH5Err(H5Acreate(object, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT));
  checkH5Err(H5Awrite(attribute, type, value.data()));
  checkH5Err(H5Aclose(attribute));
  checkH5Err(H5Tclose(type));
  checkH5Err(H5Sclose(space));
}

/**
 * Names the values of a dataset by the physical groups of a dimension: the attribute ids holds the
 * values as the dataset holds them, the attribute names their names in the same order.
 */
static void addPhysicalNames(hid_t h5file, const char* dataset, int dimension,
                             const std::vector<puml::PhysicalName>& names, bool boundary) {
  std::vector<std::int32_t> ids;
  std::vector<const char*> strings;
  for (const auto& name : names) {
    if (name.dimension == dimension) {
      ids.push_back(
          static_cast<std::int32_t>(boundary ? puml::boundaryConditionOf(name.tag) : name.tag));
      strings.push_back(name.name.c_str());
    }
  }
  if (ids.empty()) {
    return;
  }
  hid_t data = checkH5Err(H5Dopen(h5file, dataset, H5P_DEFAULT));
  const hsize_t count = ids.size();
  hid_t space = checkH5Err(H5Screate_simple(1, &count, nullptr));
  hid_t idAttribute =
      checkH5Err(H5Acreate(data, "ids", H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT));
  checkH5Err(H5Awrite(idAttribute, H5T_NATIVE_INT32, ids.data()));
  checkH5Err(H5Aclose(idAttribute));
  hid_t type = checkH5Err(H5Tcopy(H5T_C_S1));
  checkH5Err(H5Tset_size(type, H5T_VARIABLE));
  checkH5Err(H5Tset_cset(type, H5T_CSET_UTF8));
  hid_t nameAttribute = checkH5Err(H5Acreate(data, "names", type, space, H5P_DEFAULT, H5P_DEFAULT));
  checkH5Err(H5Awrite(nameAttribute, type, strings.data()));
  checkH5Err(H5Aclose(nameAttribute));
  checkH5Err(H5Tclose(type));
  checkH5Err(H5Sclose(space));
  checkH5Err(H5Dclose(data));
}

static std::size_t shapeIndex(puml::CellType type) {
  return static_cast<std::size_t>(&puml::shapeOf(type) - puml::CellShapes.data());
}

/**
 * Makes the file a VTKHDF file as well: the group /VTKHDF describes the mesh as an unstructured
 * grid of a single partition through links to the datasets, so that VTK and ParaView read the file
 * directly.
 */
static void writeVtkHdfView(hid_t h5file, std::uint64_t points, std::uint64_t cells,
                            std::uint64_t connectivityIds, bool identify, bool highOrder, int rank,
                            std::size_t chunksize) {
  hid_t group = checkH5Err(H5Gcreate(h5file, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  addIntegerAttribute(group, "Version", {1, 0});
  addFixedStringAttribute(group, "Type", "UnstructuredGrid");
  checkH5Err(H5Gclose(group));

  const std::array<std::pair<const char*, std::uint64_t>, 3> counts = {
      {{"VTKHDF/NumberOfPoints", points},
       {"VTKHDF/NumberOfCells", cells},
       {"VTKHDF/NumberOfConnectivityIds", connectivityIds}}};
  for (const auto& [name, count] : counts) {
    const std::vector<std::int64_t> value{static_cast<std::int64_t>(count)};
    writeH5Data<std::int64_t>(value, h5file, name, H5T_NATIVE_INT64, H5T_STD_I64LE, chunksize,
                              rank == 0 ? 1 : 0, 1, false, 0, 0, NoSecondDim);
  }

  for (const char* name : {"/VTKHDF/CellData", "/VTKHDF/PointData"}) {
    checkH5Err(
        H5Gclose(checkH5Err(H5Gcreate(h5file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT))));
  }
  const auto link = [&](const char* target, const char* name) {
    checkH5Err(H5Lcreate_hard(h5file, target, h5file, name, H5P_DEFAULT, H5P_DEFAULT));
  };
  link("/geometry", "/VTKHDF/Points");
  link("/connect", "/VTKHDF/Connectivity");
  link("/connect_offsets", "/VTKHDF/Offsets");
  link("/cell_type", "/VTKHDF/Types");
  link("/group", "/VTKHDF/CellData/group");
  link("/boundary", "/VTKHDF/CellData/boundary");
  if (highOrder) {
    link("/order", "/VTKHDF/CellData/order");
  }
  if (identify) {
    link("/identify", "/VTKHDF/PointData/identify");
  }
}

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
                     "Type for writing out boundary data (default: i32 for cells of at most four "
                     "faces, one integer per face otherwise; i32x4 is one 32-bit integer per face)",
                     false);

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

  std::string xdmfFile = outputFile;
  if (utils::StringUtils::endsWith(outputFile, ".puml.h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".puml.h5", ".xdmf");
  } else if (utils::StringUtils::endsWith(outputFile, ".h5")) {
    utils::StringUtils::replaceLast(xdmfFile, ".h5", ".xdmf");
  } else {
    xdmfFile.append(".xdmf");
  }

  // Create/read the mesh
  const auto sourceIndex = args.getArgument<int>("source", 0);
  std::unique_ptr<MeshData> meshInput;
  switch (sourceIndex) {
  case 0:
    logInfo() << "Using Gambit mesh";
    meshInput = std::make_unique<SerialMeshFile<puml::ParallelGambitReader>>(inputFile);
    break;
  case 1:
    logInfo() << "Using GMSH mesh format 2 (msh2) mesh";
    meshInput =
        std::make_unique<SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH2Parser>>>(inputFile);
    break;
  case 2: {
    // binary MSH 4.1 files of linear cells are read by all ranks
    int distributed = 0;
    if (rank == 0 && args.getArgument<int>("gmsh-reader", 0) == 0 &&
        puml::isBinaryMsh4(inputFile)) {
      distributed = puml::hasHighOrderCells(inputFile) ? 2 : 1;
    }
    MPI_Bcast(&distributed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (distributed == 1) {
      logInfo() << "Using GMSH mesh format 4 (msh4) mesh, binary, read by all ranks";
      meshInput = std::make_unique<SerialMeshFile<puml::DistributedGMSHReader>>(inputFile);
    } else {
      if (distributed == 2) {
        logInfo() << "The binary msh4 mesh has cells of higher order, which rank 0 reads";
      }
      logInfo() << "Using GMSH mesh format 4 (msh4) mesh";
      meshInput =
          std::make_unique<SerialMeshFile<puml::ParallelGMSHReader<puml::GMSH4Parser>>>(inputFile);
    }
    break;
  }
  case 3:
#ifdef USE_NETCDF
    logInfo() << "Using netCDF mesh";
    meshInput = std::make_unique<NetCDFMesh>(inputFile);
#else  // USE_NETCDF
    logError() << "netCDF is not supported in this version";
#endif // USE_NETCDF
    break;
  case 4:
    logInfo() << "Using APF native format";
#ifdef USE_SCOREC
    meshInput = std::make_unique<ApfNative>(inputFile, args.getArgument<const char*>("input", 0L));
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
        inputFile, args.getArgument<const char*>("cad", 0L),
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
        inputFile, args.getArgument<const char*>("cad", 0L),
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

  const std::uint64_t localCells = meshInput->cellCount();
  const std::uint64_t localVertices = meshInput->vertexCount();
  std::array<std::uint64_t, 2> globalSize{localCells, localVertices};
  MPI_Allreduce(MPI_IN_PLACE, globalSize.data(), 2, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

  // what the whole mesh holds: the kinds of its cells, cells of higher order, identification
  unsigned kinds = 0;
  for (const auto type : meshInput->cellTypes()) {
    kinds |= 1U << shapeIndex(type);
  }
  MPI_Allreduce(MPI_IN_PLACE, &kinds, 1, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
  const puml::CellShape* shape = &puml::CellShapes[0];
  std::size_t kindCount = 0;
  std::size_t maxFaces = 0;
  for (std::size_t k = 0; k < puml::CellShapes.size(); ++k) {
    if ((kinds & (1U << k)) != 0) {
      shape = &puml::CellShapes[k];
      ++kindCount;
      maxFaces = std::max(maxFaces, puml::CellShapes[k].faceCount);
    }
  }
  const bool mixed = kindCount > 1;
  int highOrder = meshInput->orders().empty() ? 0 : 1;
  MPI_Allreduce(MPI_IN_PLACE, &highOrder, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  int identify = meshInput->hasIdentify() ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &identify, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  logInfo() << "Cells:" << (mixed ? "mixed" : shape->name)
            << (highOrder != 0 ? "(some of higher order)" : "");
  logInfo() << "Total cell count:" << globalSize[0];
  logInfo() << "Total vertex count:" << globalSize[1];
  logInfo() << "Vertex identification:" << identify;

  const CellVertices cellVertices(meshInput->cellTypes(), meshInput->connectivity(),
                                  meshInput->geometry(), MPI_COMM_WORLD);
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

  // the boundary format: i32 unless given, and one integer per face for cells of more faces than
  // the packed formats hold
  int boundaryType = args.getArgument("boundarytype", 0);
  if (!args.isSet("boundarytype") && maxFaces > 4) {
    boundaryType = 4;
  }
  const std::size_t facesPerCell = std::max<std::size_t>(maxFaces, 4);
  int bitsPerFace = 0;
  std::string boundaryFormat;
  if (boundaryType == 0 || boundaryType == 2) {
    bitsPerFace = 8;
    boundaryFormat = "i32";
    logInfo() << "Using 32-bit integer boundary type conditions, or 8 bit per face (i32).";
  } else if (boundaryType == 1 || boundaryType == 3) {
    bitsPerFace = 16;
    boundaryFormat = "i64";
    logInfo() << "Using 64-bit integer boundary type conditions, or 16 bit per face (i64).";
  } else {
    boundaryFormat = "i32x" + std::to_string(facesPerCell);
    logInfo() << "Using 32-bit integer per boundary face (" << boundaryFormat << ").";
  }
  if (bitsPerFace > 0 && maxFaces > 4) {
    logError() << "The boundary format" << boundaryFormat
               << "packs four faces per cell, but the mesh has cells with" << maxFaces
               << "faces; --boundarytype=i32x4 stores one integer per face.";
  }

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

  // Write cells: a rectangle for cells of one kind; for several kinds, the layout of a VTKHDF
  // unstructured grid, a flat connectivity with the offsets of the cells in it and their kinds
  logInfo() << "Writing cells";
  const auto& cellTypes = meshInput->cellTypes();
  std::size_t connectBytes = 0;
  std::uint64_t connectivityLength = meshInput->connectivity().size();
  if (mixed) {
    std::uint64_t first = 0;
    MPI_Exscan(&connectivityLength, &first, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
      first = 0;
    }
    std::uint64_t totalLength = connectivityLength;
    MPI_Allreduce(MPI_IN_PLACE, &totalLength, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    connectBytes = writeH5Data<uint64_t>(
        meshInput->connectivity(), h5file, "connect", H5T_NATIVE_UINT64, H5T_STD_U64LE, chunksize,
        connectivityLength, totalLength, reduceInts, filterEnable, filterChunksize, NoSecondDim);

    // the last rank adds the end of the last cell
    std::vector<std::uint64_t> offsets;
    offsets.reserve(localCells + 1);
    for (const auto type : cellTypes) {
      offsets.push_back(first);
      first += puml::shapeOf(type).vertexCount;
    }
    if (rank == processes - 1) {
      offsets.push_back(first);
    }
    writeH5Data<uint64_t>(offsets, h5file, "connect_offsets", H5T_NATIVE_UINT64, H5T_STD_U64LE,
                          chunksize, offsets.size(), globalSize[0] + 1, reduceInts, filterEnable,
                          filterChunksize, NoSecondDim);

    std::vector<std::uint8_t> kindValues(localCells);
    for (std::size_t cell = 0; cell < localCells; ++cell) {
      kindValues[cell] = static_cast<std::uint8_t>(cellTypes[cell]);
    }
    writeH5Data<std::uint8_t>(kindValues, h5file, "cell_type", H5T_NATIVE_UINT8, H5T_STD_U8LE,
                              chunksize, localCells, globalSize[0], false, filterEnable,
                              filterChunksize, NoSecondDim);
    connectivityLength = totalLength;
  } else {
    connectBytes = writeH5Data<uint64_t>(
        meshInput->connectivity(), h5file, "connect", H5T_NATIVE_UINT64, H5T_STD_U64LE, chunksize,
        localCells, globalSize[0], reduceInts, filterEnable, filterChunksize, shape->vertexCount);
  }

  // Vertices
  logInfo() << "Writing vertices";
  writeH5Data<double>(meshInput->geometry(), h5file, "geometry", H5T_IEEE_F64LE, H5T_IEEE_F64LE,
                      chunksize, localVertices, globalSize[1], reduceInts, filterEnable,
                      filterChunksize, 3);

  // Group information
  logInfo() << "Writing group information";
  const auto groupBytes = writeH5Data<int32_t>(
      meshInput->group(), h5file, "group", H5T_NATIVE_INT32, H5T_STD_I32LE, chunksize, localCells,
      globalSize[0], reduceInts, filterEnable, filterChunksize, NoSecondDim);

  // Write boundary condition
  logInfo() << "Writing boundary condition";
  const auto& faces = meshInput->boundary();
  std::size_t boundaryBytes = 0;
  if (bitsPerFace > 0) {
    // the packed bits are written as they are, as signed integers
    const auto packed = packBoundaries(faces, cellTypes, bitsPerFace);
    if (bitsPerFace == 8) {
      boundaryBytes = writeH5Data<int32_t>(packed, h5file, "boundary", H5T_NATIVE_INT32,
                                           H5T_STD_I32LE, chunksize, localCells, globalSize[0],
                                           reduceInts, filterEnable, filterChunksize, NoSecondDim);
    } else {
      boundaryBytes = writeH5Data<int64_t>(packed, h5file, "boundary", H5T_NATIVE_INT64,
                                           H5T_STD_I64LE, chunksize, localCells, globalSize[0],
                                           reduceInts, filterEnable, filterChunksize, NoSecondDim);
    }
  } else {
    const auto values = boundariesPerFace(faces, cellTypes, facesPerCell);
    boundaryBytes = writeH5Data<int32_t>(values, h5file, "boundary", H5T_NATIVE_INT32,
                                         H5T_STD_I32LE, chunksize, localCells, globalSize[0],
                                         reduceInts, filterEnable, filterChunksize, facesPerCell);
  }

  // the names of the groups and boundary conditions, where the source has them
  addPhysicalNames(h5file, "/group", 3, meshInput->physicalNames(), false);
  addPhysicalNames(h5file, "/boundary", 2, meshInput->physicalNames(), true);

  std::size_t identifyBytes = 0;
  if (identify != 0) {
    logInfo() << "Writing vertex topology identification";
    identifyBytes = writeH5Data<uint64_t>(
        meshInput->identify(), h5file, "identify", H5T_NATIVE_UINT64, H5T_STD_U64LE, chunksize,
        localVertices, globalSize[1], reduceInts, filterEnable, filterChunksize, NoSecondDim);
  }

  // the nodes of the cells of higher order but their vertices, in the order of gmsh, with the
  // offsets of the cells in them (in values) and the order of every cell
  if (highOrder != 0) {
    logInfo() << "Writing the geometry of higher order";
    const auto& orders = meshInput->orders();
    std::uint64_t length = meshInput->highOrderGeometry().size();
    std::uint64_t first = 0;
    MPI_Exscan(&length, &first, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
      first = 0;
    }
    std::uint64_t totalLength = length;
    MPI_Allreduce(MPI_IN_PLACE, &totalLength, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    writeH5Data<double>(meshInput->highOrderGeometry(), h5file, "geometry_ho", H5T_IEEE_F64LE,
                        H5T_IEEE_F64LE, chunksize, length, totalLength, false, filterEnable,
                        filterChunksize, NoSecondDim);

    std::vector<std::uint64_t> offsets;
    offsets.reserve(localCells + 1);
    for (std::size_t cell = 0; cell < localCells; ++cell) {
      offsets.push_back(first);
      const auto type = cellTypes[cell];
      first += 3 * (puml::lagrangeNodeCount(type, orders[cell]) - puml::shapeOf(type).vertexCount);
    }
    if (rank == processes - 1) {
      offsets.push_back(first);
    }
    writeH5Data<uint64_t>(offsets, h5file, "geometry_ho_offsets", H5T_NATIVE_UINT64, H5T_STD_U64LE,
                          chunksize, offsets.size(), globalSize[0] + 1, reduceInts, filterEnable,
                          filterChunksize, NoSecondDim);
    writeH5Data<std::uint8_t>(orders, h5file, "order", H5T_NATIVE_UINT8, H5T_STD_U8LE, chunksize,
                              localCells, globalSize[0], false, filterEnable, filterChunksize,
                              NoSecondDim);

    hid_t dataset = checkH5Err(H5Dopen(h5file, "/geometry_ho", H5P_DEFAULT));
    addAttribute(dataset, "node-ordering", "gmsh");
    checkH5Err(H5Dclose(dataset));
  }

  // what the file holds, how it is laid out, and where it comes from
  addAttribute(h5file, "boundary-format", boundaryFormat);
  addAttribute(h5file, "topology-format", identify != 0 ? "identify-vertex" : "geometric");
  addAttribute(h5file, "cell-type", mixed ? "mixed" : shape->name);
  addIntegerAttribute(h5file, "format-version", {2, 0});
  addAttribute(h5file, "generator", std::string("PUMGen ") + PUMGEN_VERSION);
  {
    std::string command;
    for (int i = 0; i < argc; ++i) {
      command += (i > 0 ? " " : "") + std::string(argv[i]);
    }
    addAttribute(h5file, "command", command);
  }
  addAttribute(h5file, "source-format", source[static_cast<std::size_t>(sourceIndex)]);
  {
    std::array<char, 32> created{};
    if (rank == 0) {
      const std::time_t now = std::time(nullptr);
      std::strftime(created.data(), created.size(), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    }
    MPI_Bcast(created.data(), static_cast<int>(created.size()), MPI_CHAR, 0, MPI_COMM_WORLD);
    addAttribute(h5file, "created", created.data());
  }

  if (mixed) {
    writeVtkHdfView(h5file, globalSize[1], globalSize[0], connectivityLength, identify != 0,
                    highOrder != 0, rank, chunksize);
    logInfo() << "The cells are of several kinds, which XDMF cannot describe; the mesh file is a "
                 "VTKHDF file as well, which ParaView opens directly";
  } else if (rank == 0) {
    // Writing XDMF file
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
         << "   <Topology TopologyType=\"" << shape->xdmfName << "\" NumberOfElements=\""
         << globalSize[0] << "\">" << std::endl
         << "    <DataItem NumberType=\"Int\" Precision=\"" << xdmfPrecision(connectBytes)
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0] << " " << shape->vertexCount << "\">"
         << basename << ":/connect</DataItem>" << std::endl
         << "   </Topology>" << std::endl
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
         << "\" Format=\"HDF\" Dimensions=\"" << globalSize[0]
         << (bitsPerFace > 0 ? std::string() : " " + std::to_string(facesPerCell)) << "\">"
         << basename << ":/boundary</DataItem>" << std::endl
         << "   </Attribute>" << std::endl;
    if (highOrder != 0) {
      xdmf << "   <Attribute Name=\"order\" Center=\"Cell\">" << std::endl
           << "    <DataItem NumberType=\"UInt\" Precision=\"1\" Format=\"HDF\" Dimensions=\""
           << globalSize[0] << "\">" << basename << ":/order</DataItem>" << std::endl
           << "   </Attribute>" << std::endl;
    }
    if (identify != 0) {
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
