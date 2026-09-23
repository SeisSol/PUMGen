// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

// Reads a mesh file of PUMGen with PUML and checks what PUML makes of it against the file: the
// numbers of cells, vertices, faces and faces on the surface, which are counted in the file from
// the faces of the cells, and the Euler characteristic of a mesh of a ball. Every face on the
// surface has to have a boundary condition in the file.

#include <hdf5.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "CellType.h"
#include "Hdf5Reader.h"
#include "PUML.h"
#include "Upward.h"

namespace {

long globalSum(long value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  return value;
}

// the elements a rank owns: those it shares with no rank of a smaller number
template <typename T> long owned(const std::vector<T>& elements, int rank) {
  long count = 0;
  for (const auto& element : elements) {
    if (element.shared().empty() || element.shared()[0] > rank) {
      ++count;
    }
  }
  return count;
}

std::string stringAttribute(hid_t file, const char* name) {
  const hid_t attribute = H5Aopen(file, name, H5P_DEFAULT);
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, H5T_VARIABLE);
  char* text = nullptr;
  H5Aread(attribute, type, static_cast<void*>(&text));
  std::string value = text != nullptr ? text : "";
  H5free_memory(text);
  H5Tclose(type);
  H5Aclose(attribute);
  return value;
}

std::array<hsize_t, 2> extent(hid_t file, const char* name) {
  const hid_t data = H5Dopen(file, name, H5P_DEFAULT);
  const hid_t space = H5Dget_space(data);
  std::array<hsize_t, 2> dims{0, 1};
  H5Sget_simple_extent_dims(space, dims.data(), nullptr);
  H5Sclose(space);
  H5Dclose(data);
  return dims;
}

template <typename T> std::vector<T> readAll(hid_t file, const char* name, hid_t type) {
  const auto dims = extent(file, name);
  std::vector<T> values(dims[0] * dims[1]);
  const hid_t data = H5Dopen(file, name, H5P_DEFAULT);
  H5Dread(data, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
  H5Dclose(data);
  return values;
}

constexpr std::size_t MaxFaces = 6;

// the faces of the cells by their VTK cell types, as PUML and PUMGen number them
const std::map<int, std::vector<std::vector<std::size_t>>> Faces = {
    {10, {{1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}},
    {12, {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}}},
    {13, {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}}},
    {14, {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}}}};
const std::map<int, std::size_t> VertexCounts = {{10, 4}, {12, 8}, {13, 6}, {14, 5}};
const std::map<std::string, int> KindsByName = {
    {"tetrahedron", 10}, {"hexahedron", 12}, {"wedge", 13}, {"pyramid", 14}};

struct FileCounts {
  long cells = 0;
  long vertices = 0;
  long faces = 0;
  long surface = 0;
  long surfaceWithoutCondition = 0;
  bool mixed = false;
  std::string cellType;
};

FileCounts readFile(const std::string& path) {
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  FileCounts counts;
  counts.mixed = H5Lexists(file, "/connect_offsets", H5P_DEFAULT) > 0;
  counts.cells = static_cast<long>(counts.mixed ? extent(file, "/connect_offsets")[0] - 1
                                                : extent(file, "/connect")[0]);
  counts.vertices = static_cast<long>(extent(file, "/geometry")[0]);
  counts.cellType = stringAttribute(file, "cell-type");

  // the boundary conditions, packed into one integer per cell or one integer per face
  const auto format = stringAttribute(file, "boundary-format");
  std::vector<std::int64_t> boundary(static_cast<std::size_t>(counts.cells) * MaxFaces, 0);
  if (format == "i32" || format == "i64") {
    const int bits = format == "i32" ? 8 : 16;
    const std::uint64_t mask = (std::uint64_t{1} << bits) - 1;
    const auto values = readAll<std::int64_t>(file, "/boundary", H5T_NATIVE_INT64);
    for (std::size_t cell = 0; cell < values.size(); ++cell) {
      const auto packed = format == "i32" ? static_cast<std::uint32_t>(values[cell])
                                          : static_cast<std::uint64_t>(values[cell]);
      for (std::size_t face = 0; face < 4; ++face) {
        boundary[cell * MaxFaces + face] =
            static_cast<std::int64_t>((packed >> (face * bits)) & mask);
      }
    }
  } else {
    const auto columns = extent(file, "/boundary")[1];
    const auto values = readAll<std::int32_t>(file, "/boundary", H5T_NATIVE_INT32);
    for (std::size_t i = 0; i < values.size(); ++i) {
      boundary[(i / columns) * MaxFaces + i % columns] = values[i];
    }
  }

  // the cells: their kinds and where their vertices begin in the connectivity
  const auto connect = readAll<std::uint64_t>(file, "/connect", H5T_NATIVE_UINT64);
  std::vector<int> kinds(static_cast<std::size_t>(counts.cells));
  std::vector<std::uint64_t> offsets(kinds.size() + 1);
  if (counts.mixed) {
    const auto types = readAll<std::uint8_t>(file, "/cell_type", H5T_NATIVE_UINT8);
    std::copy(types.begin(), types.end(), kinds.begin());
    offsets = readAll<std::uint64_t>(file, "/connect_offsets", H5T_NATIVE_UINT64);
  } else {
    std::fill(kinds.begin(), kinds.end(), KindsByName.at(counts.cellType));
    for (std::size_t cell = 0; cell <= kinds.size(); ++cell) {
      offsets[cell] = cell * VertexCounts.at(KindsByName.at(counts.cellType));
    }
  }
  H5Fclose(file);

  // every face by its sorted vertices, with the faces of the cells which have it
  std::map<std::array<std::uint64_t, 4>, std::vector<std::size_t>> faces;
  for (std::size_t cell = 0; cell < kinds.size(); ++cell) {
    const auto& cellFaces = Faces.at(kinds[cell]);
    for (std::size_t f = 0; f < cellFaces.size(); ++f) {
      std::array<std::uint64_t, 4> key{};
      key.fill(std::numeric_limits<std::uint64_t>::max());
      for (std::size_t k = 0; k < cellFaces[f].size(); ++k) {
        key[k] = connect[offsets[cell] + cellFaces[f][k]];
      }
      std::sort(key.begin(), key.end());
      faces[key].push_back(cell * MaxFaces + f);
    }
  }
  counts.faces = static_cast<long>(faces.size());
  for (const auto& entry : faces) {
    if (entry.second.size() == 1) {
      ++counts.surface;
      counts.surfaceWithoutCondition += boundary[entry.second[0]] == 0 ? 1 : 0;
    }
  }
  return counts;
}

} // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (argc != 2) {
    if (rank == 0) {
      std::fprintf(stderr, "Usage: %s <mesh file>\n", argv[0]);
    }
    MPI_Finalize();
    return 2;
  }
  const std::string path = argv[1];
  const auto file = readFile(path);

  // a file of one kind of cells holds a rectangle, which PUML reads as the kind it is given
  const std::map<std::string, PUML::CellType> kinds = {{"tetrahedron", PUML::CellType::Tetrahedron},
                                                       {"hexahedron", PUML::CellType::Hexahedron},
                                                       {"wedge", PUML::CellType::Wedge},
                                                       {"pyramid", PUML::CellType::Pyramid},
                                                       {"mixed", PUML::CellType::Tetrahedron}};
  bool ok = false;
  {
    PUML::MIXEDPUML puml;
    puml.setComm(MPI_COMM_WORLD);
    {
      PUML::Hdf5Reader<PUML::MIXED> reader(puml, kinds.at(file.cellType));
      reader.openMixed(path + ":/connect", path + ":/connect_offsets", path + ":/cell_type",
                       path + ":/geometry");
    }
    puml.generateMesh();

    const long cells = globalSum(static_cast<long>(puml.cells().size()));
    const long vertices = globalSum(owned(puml.vertices(), rank));
    const long faces = globalSum(owned(puml.faces(), rank));
    const long edges = globalSum(owned(puml.edges(), rank));

    // the faces of a single cell lie on the surface
    long surface = 0;
    for (const auto& face : puml.faces()) {
      std::array<PUML::LocalId, 2> adjacent{};
      PUML::Upward::cells(puml, face, adjacent.data());
      if (adjacent[0] != PUML::InvalidLocalId && adjacent[1] == PUML::InvalidLocalId &&
          face.shared().empty()) {
        ++surface;
      }
    }
    surface = globalSum(surface);
    const long euler = vertices - edges + faces - cells;

    ok = cells == file.cells && vertices == file.vertices && faces == file.faces &&
         surface == file.surface && file.surfaceWithoutCondition == 0 && euler == 1;
    if (rank == 0) {
      std::printf("%s: cells %ld (file %ld), vertices %ld (file %ld), faces %ld (file %ld), on the "
                  "surface %ld (file %ld, without a boundary condition %ld), edges %ld, Euler "
                  "characteristic %ld: %s\n",
                  path.c_str(), cells, file.cells, vertices, file.vertices, faces, file.faces,
                  surface, file.surface, file.surfaceWithoutCondition, edges, euler,
                  ok ? "ok" : "WRONG");
    }
  }
  MPI_Finalize();
  return ok ? 0 : 1;
}
