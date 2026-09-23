// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "doctest/extensions/doctest_mpi.h"

#include <cstddef>

#include "helper/Distributor.h"

TEST_CASE("The chunks of all ranks cover the items exactly once") {
  for (const std::size_t total : {0, 1, 2, 3, 7, 64, 1000, 1597}) {
    for (const int size : {1, 2, 3, 4, 7, 16}) {
      CAPTURE(total);
      CAPTURE(size);
      std::size_t sum = 0;
      for (int rank = 0; rank < size; ++rank) {
        REQUIRE(getChunksum(total, rank, size) == sum);
        const auto chunk = getChunksize(total, rank, size);
        // the chunks differ by at most one item, the larger ones come first
        CHECK(chunk + 1 >= getChunksize(total, 0, size));
        for (std::size_t item = sum; item < sum + chunk; ++item) {
          REQUIRE(getChunkOwner(total, item, size) == rank);
        }
        sum += chunk;
      }
      CHECK(sum == total);
      CHECK(getChunksum(total, size, size) == total);
    }
  }
}
