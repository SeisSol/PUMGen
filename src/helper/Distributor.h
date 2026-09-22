// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_DISTRIBUTOR_H_
#define PUMGEN_SRC_HELPER_DISTRIBUTOR_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>

// no namespace for now

constexpr std::size_t getChunksize(std::size_t total, int rank, int size) {
  return (total / size) + (total % size > rank);
}

constexpr std::size_t getChunksum(std::size_t total, int until, int size) {
  std::size_t base = (total / size) * until;
  std::size_t rest = total % size;
  std::size_t addon = std::min(static_cast<std::size_t>(until), rest);
  return base + addon;
}

/**
 * The rank whose chunk (as given by getChunksize/getChunksum) contains the item index.
 */
constexpr int getChunkOwner(std::size_t total, std::size_t index, int size) {
  const std::size_t base = total / size;
  const std::size_t rest = total % size;
  const std::size_t largeChunks = rest * (base + 1);
  if (index < largeChunks) {
    return static_cast<int>(index / (base + 1));
  }
  return static_cast<int>(rest + (index - largeChunks) / base);
}

#endif // PUMGEN_SRC_HELPER_DISTRIBUTOR_H_
