// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

// Writes gmsh Structured fields for velocity-aware meshing (one per refinement cuboid) and a .geo
// file which combines them into the background mesh size.

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "helper/Distributor.h"
#include "sizing/SizeField.h"
#include "sizing/VelocityAwareMeshSize.h"
#include "sizing/VelocityAwareSettings.h"
#include "tinyxml2/tinyxml2.h"
#include "utils/args.h"
#include "utils/logger.h"

namespace {

/**
 * Writes the binary format of gmsh's Structured field: origin and spacing (3 doubles each), node
 * counts (3 4-byte ints), then the values. Each rank writes the nodes with x index from first.
 */
void writeGrid(const std::string& fileName, const StructuredGrid& grid,
               const std::vector<double>& values, std::size_t first, MPI_Comm comm) {
  MPI_File file;
  if (MPI_File_open(comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                    &file) != MPI_SUCCESS) {
    logError() << "Could not open" << fileName;
  }
  MPI_File_set_size(file, 0);

  constexpr std::size_t HeaderSize = 6 * sizeof(double) + 3 * sizeof(std::int32_t);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    char header[HeaderSize];
    const double spacing[3] = {grid.spacing, grid.spacing, grid.spacing};
    std::int32_t count[3];
    for (int d = 0; d < 3; ++d) {
      count[d] = static_cast<std::int32_t>(grid.count[d]);
    }
    std::memcpy(header, grid.origin.data(), 3 * sizeof(double));
    std::memcpy(header + 3 * sizeof(double), spacing, 3 * sizeof(double));
    std::memcpy(header + 6 * sizeof(double), count, 3 * sizeof(std::int32_t));
    MPI_File_write_at(file, 0, header, HeaderSize, MPI_BYTE, MPI_STATUS_IGNORE);
  }

  // in pieces, since MPI counts are ints
  constexpr std::size_t Piece = std::size_t{1} << 26;
  const MPI_Offset offset = HeaderSize + first * grid.count[1] * grid.count[2] * sizeof(double);
  for (std::size_t done = 0; done < values.size(); done += Piece) {
    const auto count = static_cast<int>(std::min(Piece, values.size() - done));
    MPI_File_write_at(file, offset + static_cast<MPI_Offset>(done * sizeof(double)),
                      values.data() + done, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
  MPI_File_close(&file);
}

std::string baseName(const std::string& path) {
  const auto slash = path.find_last_of('/');
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

void writeGeo(const std::string& fileName, const std::string& settingsFile,
              const std::vector<std::string>& gridFiles, const std::vector<int>& groups,
              int firstField) {
  std::ofstream geo(fileName);
  geo << "// Velocity-aware mesh size written by pumgen-sizefield from " << settingsFile << ":\n"
      << "// one Structured field per refinement cuboid, restricted to the physical volume of\n"
      << "// the cuboid's group, and their minimum as background field. The grid files are\n"
      << "// expected next to this file. Recommended options:\n"
      << "// Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0;\n"
      << "// Mesh.MeshSizeFromCurvature = 0;\n";
  std::string restricted;
  for (std::size_t c = 0; c < gridFiles.size(); ++c) {
    const int grid = firstField + 2 * static_cast<int>(c);
    const int restrict = grid + 1;
    geo << "Field[" << grid << "] = Structured;\n"
        << "Field[" << grid << "].FileName = StrCat(CurrentDirectory, \"" << baseName(gridFiles[c])
        << "\");\n"
        << "Field[" << grid << "].TextFormat = 0;\n"
        << "Field[" << grid << "].SetOutsideValue = 1;\n"
        << "Field[" << grid << "].OutsideValue = " << UnconstrainedSize << ";\n"
        << "Field[" << restrict << "] = Restrict;\n"
        << "Field[" << restrict << "].InField = " << grid << ";\n"
        << "Field[" << restrict << "].VolumesList = {Physical Volume{" << groups[c] << "}};\n"
        << "Field[" << restrict << "].IncludeBoundary = 1;\n";
    restricted += (c == 0 ? "" : ", ") + std::to_string(restrict);
  }
  const int minimum = firstField + 2 * static_cast<int>(gridFiles.size());
  geo << "Field[" << minimum << "] = Min;\n"
      << "Field[" << minimum << "].FieldsList = {" << restricted << "};\n"
      << "Background Field = " << minimum << ";\n";
}

} // namespace

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  utils::Logger::setRank(rank);

  utils::Args args("Writes gmsh Structured fields for velocity-aware meshing, one per refinement "
                   "cuboid, and a .geo file which restricts each field to the physical volume of "
                   "the cuboid's group (bypassFindRegionAndUseGroup) and uses their minimum as "
                   "background field.");
  args.addOption("spacing", 0, "Grid spacing; about the smallest mesh size", utils::Args::Required,
                 true);
  args.addOption("first-field", 0, "Number of the first gmsh field (default: 1000)",
                 utils::Args::Required, false);
  args.addAdditionalOption("settings", "Mesh attributes (XML) with a VelocityAwareMeshing element");
  args.addAdditionalOption("output", "Prefix of the output files");
  if (args.parse(argc, argv, rank == 0) != utils::Args::Success) {
    MPI_Finalize();
    return 1;
  }
  const auto settingsFile = args.getAdditionalArgument<std::string>("settings");
  const auto prefix = args.getAdditionalArgument<std::string>("output");
  const auto spacing = args.getArgument<double>("spacing");
  const auto firstField = args.getArgument<int>("first-field", 1000);
  if (!(spacing > 0)) {
    logError() << "The grid spacing has to be positive";
  }

  tinyxml2::XMLDocument doc;
  if (doc.LoadFile(settingsFile.c_str()) != tinyxml2::XML_SUCCESS) {
    logError() << "Could not read" << settingsFile;
  }
  const auto settings = readVelocityAwareSettings(doc);
  const auto& regions = settings.getRefinementRegions();
  if (regions.empty()) {
    logError() << "No refinement cuboid in" << settingsFile;
  }
  for (std::size_t c = 0; c < regions.size(); ++c) {
    if (regions[c].bypassFindRegionAndUseGroup <= 0) {
      logError() << "Refinement cuboid" << c
                 << "has no group; pumgen-sizefield needs bypassFindRegionAndUseGroup";
    }
  }

  const VelocityAwareMeshSize meshSize(settings);
  std::vector<std::string> gridFiles;
  std::vector<int> groups;
  for (std::size_t c = 0; c < regions.size(); ++c) {
    const auto grid = gridForCuboid(regions[c].cuboid, spacing);
    for (const auto count : grid.count) {
      if (count > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        logError() << "Too many grid nodes for refinement cuboid" << c;
      }
    }
    logInfo() << "Refinement cuboid" << c << "(group" << regions[c].bypassFindRegionAndUseGroup
              << ", " << regions[c].targetedFrequency << "Hz): grid of" << grid.count[0] << "x"
              << grid.count[1] << "x" << grid.count[2] << "nodes";

    const std::size_t first = getChunksum(grid.count[0], rank, size);
    const std::size_t last = first + getChunksize(grid.count[0], rank, size);
    const auto values = cuboidMeshSizes(meshSize, regions[c], grid, first, last);

    double smallest =
        values.empty() ? UnconstrainedSize : *std::min_element(values.begin(), values.end());
    MPI_Allreduce(MPI_IN_PLACE, &smallest, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    logInfo() << "Refinement cuboid" << c << ": smallest mesh size" << smallest;

    const std::string gridFile = prefix + "-" + std::to_string(c) + ".bin";
    writeGrid(gridFile, grid, values, first, MPI_COMM_WORLD);
    gridFiles.push_back(gridFile);
    groups.push_back(regions[c].bypassFindRegionAndUseGroup);
  }

  if (rank == 0) {
    writeGeo(prefix + ".geo", settingsFile, gridFiles, groups, firstField);
    logInfo() << "Wrote" << (prefix + ".geo").c_str();
  }

  MPI_Finalize();
  return 0;
}
