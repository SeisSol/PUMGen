#include "ParallelGMSHReader.h"
#include "third_party/GMSHParser.h"
#include "utils/logger.h"

#include <array>
#include <cstddef>
#include <iterator>

namespace puml {

void ParallelGMSHReader::open(char const *meshFile) {
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

void ParallelGMSHReader::convertBoundaryConditions() {
  const auto nVerts = builder_.vertices.size();
  const auto nFacets = builder_.facets.size();
  const auto nElements = builder_.elements.size();

  bcs_.clear();
  bcs_.resize(nElements);

  constexpr auto nFacetsPerTet = Dim + 1u;
  constexpr auto nVertsPerTri = Dim;

  // vertex to triangle map
  std::vector<std::vector<std::size_t>> vertex2facets(nVerts);
  for (std::size_t fctNo = 0; fctNo < nFacets; ++fctNo) {
    for (auto const &vtxNo : builder_.facets[fctNo]) {
      vertex2facets[vtxNo].push_back(fctNo);
    }
  }
  for (auto &v2f : vertex2facets) {
    std::sort(v2f.begin(), v2f.end());
  }

  /**
   * GMSH stores boundary conditions on a surface mesh whereas SeisSol expects
   * boundary conditions to be stored per element. In the following we convert
   * the boundary condition representation by matching the 4 faces of an element
   * with the surface mesh.
   */
  for (unsigned elNo = 0; elNo < nElements; ++elNo) {
    for (unsigned localFctNo = 0; localFctNo < nFacetsPerTet; ++localFctNo) {
      unsigned nodes[nVertsPerTri];
      for (unsigned localNo = 0; localNo < nVertsPerTri; ++localNo) {
        const auto localNoElement = Facet2Nodes[localFctNo][localNo];
        nodes[localNo] = builder_.elements[elNo][localNoElement];
      }
      std::vector<unsigned> intersect[nVertsPerTri - 1];
      std::set_intersection(
          vertex2facets[nodes[0]].begin(), vertex2facets[nodes[0]].end(),
          vertex2facets[nodes[1]].begin(), vertex2facets[nodes[1]].end(),
          std::back_inserter(intersect[0]));
      for (unsigned node = 2; node < nVertsPerTri; ++node) {
        std::set_intersection(intersect[node - 2].begin(),
                              intersect[node - 2].end(),
                              vertex2facets[nodes[node]].begin(),
                              vertex2facets[nodes[node]].end(),
                              std::back_inserter(intersect[node - 1]));
      }
      if (!intersect[nVertsPerTri - 2].empty()) {
        if (intersect[nVertsPerTri - 2].size() > 1) {
          logError() << "A face of an element exists multiple times in the "
                        "surface mesh.";
        }
        const auto fctNo = intersect[nVertsPerTri - 2][0];
        bcs_[elNo][localFctNo] = adjustBoundaryCondition(builder_.bcs[fctNo]);
      } else {
        bcs_[elNo][localFctNo] = 0;
      }
    }
  }
}

} // namespace puml
