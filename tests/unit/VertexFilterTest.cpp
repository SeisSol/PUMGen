// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "input/ParallelVertexFilter.h"
#include "third_party/MPITraits.h"

MPI_TEST_CASE("The vertex filter merges duplicates, also with a rank without vertices", 3) {
  // vertex k lies at (k, 2k, 3k); rank 1 has none, rank 2 repeats two of rank 0
  const std::vector<std::vector<int>> keys = {{0, 1, 2}, {}, {1, 3, 0, 4}};
  std::vector<double> vertices;
  for (const int key : keys[test_rank]) {
    vertices.insert(vertices.end(), {1.0 * key, 2.0 * key, 3.0 * key});
  }

  ParallelVertexFilter filter(test_comm);
  filter.filter(keys[test_rank].size(), vertices);

  std::size_t unique = filter.numLocalVertices();
  MPI_Allreduce(MPI_IN_PLACE, &unique, 1, tndm::mpi_type_t<std::size_t>(), MPI_SUM, test_comm);
  CHECK(unique == 5);

  // the same vertex has the same global id on all ranks, and different vertices different ones
  std::array<std::size_t, 8> mine{};
  mine.fill(~std::size_t{0});
  for (std::size_t i = 0; i < keys[test_rank].size(); ++i) {
    mine[2 * i] = keys[test_rank][i];
    mine[2 * i + 1] = filter.globalIds()[i];
  }
  std::array<std::size_t, 24> all{};
  MPI_Allgather(mine.data(), 8, tndm::mpi_type_t<std::size_t>(), all.data(), 8,
                tndm::mpi_type_t<std::size_t>(), test_comm);
  std::map<std::size_t, std::size_t> idOf;
  std::set<std::size_t> ids;
  for (std::size_t i = 0; i < all.size(); i += 2) {
    if (all[i] == ~std::size_t{0}) {
      continue;
    }
    const auto [it, inserted] = idOf.emplace(all[i], all[i + 1]);
    CHECK(it->second == all[i + 1]);
    ids.insert(all[i + 1]);
  }
  CHECK(idOf.size() == 5);
  CHECK(ids.size() == 5);
  CHECK(*ids.rbegin() == 4);
}
