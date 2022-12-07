#include "ParallelGMSHReader.h"
#include "third_party/GMSHParser.h"

#include <array>
#include <cstddef>
#include <iterator>

namespace puml {

template <> void ParallelGMSHReader<gmsh_version::v2>::open(char const* meshFile) {
  int rank;
  MPI_Comm_rank(comm_, &rank);

  if (rank == 0) {
    tndm::GMSHParser parser(&builder_);
    bool ok = parser.parseFile(meshFile);
    if (!ok) {
      logError() << meshFile << std::endl << parser.getErrorMessage();
    }
    convertBoundaryConditions();

    nVertices_ = builder_.vertices.size();
    nElements_ = builder_.elements.size();
  }

  MPI_Bcast(&nVertices_, 1, tndm::mpi_type_t<decltype(nVertices_)>(), 0, comm_);
  MPI_Bcast(&nElements_, 1, tndm::mpi_type_t<decltype(nVertices_)>(), 0, comm_);
}

template <> void ParallelGMSHReader<gmsh_version::v4>::open(char const* meshFile) {
  int rank;
  MPI_Comm_rank(comm_, &rank);

  if (rank == 0) {
    tndm::GMSHParser parser(&builder_);
    bool ok = parser.parseFile(meshFile);
    if (!ok) {
      logError() << meshFile << std::endl << parser.getErrorMessage();
    }
    convertBoundaryConditions();

    nVertices_ = builder_.vertices.size();
    nElements_ = builder_.elements.size();
  }

  MPI_Bcast(&nVertices_, 1, tndm::mpi_type_t<decltype(nVertices_)>(), 0, comm_);
  MPI_Bcast(&nElements_, 1, tndm::mpi_type_t<decltype(nVertices_)>(), 0, comm_);
}

} // namespace puml
