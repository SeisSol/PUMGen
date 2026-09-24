// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESH_PHYSICALNAME_H_
#define PUMGEN_SRC_MESH_PHYSICALNAME_H_

#include <mpi.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace puml {

/**
 * The name of a physical group of gmsh: of cells (dimension 3), whose tag is their group, or of
 * faces (dimension 2), whose tag gives their boundary condition.
 */
struct PhysicalName {
  int dimension = 0;
  long tag = 0;
  std::string name;
};

inline bool operator==(const PhysicalName& a, const PhysicalName& b) {
  return a.dimension == b.dimension && a.tag == b.tag && a.name == b.name;
}

/**
 * Boundary conditions in SeisSol used to start at 100, e.g. 101 = free surface. In the hdf5 format
 * one starts counting from 0, e.g. 1 = free surface. In order to be compatible with legacy gmsh
 * scripts, 100 is subtracted from a boundary condition larger than or equal to 100.
 */
constexpr long boundaryConditionOf(long tag) {
  constexpr long Offset = 100;
  return tag >= Offset ? tag - Offset : tag;
}

/**
 * Sends the names which the rank root holds to all ranks. Collective.
 */
inline void broadcastPhysicalNames(std::vector<PhysicalName>& names, int root, MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  // the dimension, tag and length of every name, and then all characters
  std::vector<std::int64_t> header;
  std::string characters;
  if (rank == root) {
    for (const auto& name : names) {
      header.insert(header.end(),
                    {name.dimension, name.tag, static_cast<std::int64_t>(name.name.size())});
      characters += name.name;
    }
  }
  std::array<std::uint64_t, 2> sizes{header.size(), characters.size()};
  MPI_Bcast(sizes.data(), 2, MPI_UINT64_T, root, comm);
  header.resize(sizes[0]);
  characters.resize(sizes[1]);
  MPI_Bcast(header.data(), static_cast<int>(sizes[0]), MPI_INT64_T, root, comm);
  MPI_Bcast(characters.data(), static_cast<int>(sizes[1]), MPI_CHAR, root, comm);
  if (rank != root) {
    names.clear();
    std::size_t position = 0;
    for (std::size_t i = 0; i + 2 < header.size(); i += 3) {
      const auto length = static_cast<std::size_t>(header[i + 2]);
      names.push_back({static_cast<int>(header[i]), static_cast<long>(header[i + 1]),
                       characters.substr(position, length)});
      position += length;
    }
  }
}

} // namespace puml

#endif // PUMGEN_SRC_MESH_PHYSICALNAME_H_
