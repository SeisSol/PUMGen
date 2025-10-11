#ifndef PUMGEN_AUX_DISTRIBUTOR_H_
#define PUMGEN_AUX_DISTRIBUTOR_H_

#include <cmath>
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

#endif // PUMGEN_AUX_DISTRIBUTOR_H_
