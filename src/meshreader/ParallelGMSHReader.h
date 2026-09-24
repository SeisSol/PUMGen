// SPDX-FileCopyrightText: 2020 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_
#define PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_

#include <mpi.h>

#include "GMSHBuilder.h"
#include "input/MeshData.h"
#include "utils/logger.h"

namespace puml {

/**
 * Reads an MSH file on rank 0 with the parser P and distributes the mesh.
 */
template <typename P> class ParallelGMSHReader {
  public:
  constexpr static bool ProvidesLocalMesh = true;

  explicit ParallelGMSHReader(MPI_Comm comm = MPI_COMM_WORLD) : comm(comm) {}

  void open(const char* meshFile) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      GMSHBuilder builder;
      P parser(&builder);
      if (!parser.parseFile(meshFile)) {
        logError() << meshFile << std::endl << parser.getErrorMessage();
      }
      builder.postprocess();
      if (builder.hasHighOrder()) {
        logInfo() << "The mesh has cells of higher order; the vertices of the cells become the "
                     "vertices of the mesh";
      }
      mesh = prepareMesh(builder);
    }
  }

  /**
   * The part of the mesh of this rank. Collective, and possible only once.
   */
  LocalMesh read() { return distributeMesh(mesh, comm); }

  private:
  MPI_Comm comm;
  GlobalMesh mesh;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_PARALLELGMSHREADER_H_
