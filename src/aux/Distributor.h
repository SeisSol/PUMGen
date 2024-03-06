#ifndef PUMGEN_AUX_DISTRIBUTOR_H_
#define PUMGEN_AUX_DISTRIBUTOR_H_

#include <cstdlib>

// no namespace for now

constexpr std::size_t getChunksize(std::size_t total, int rank, int size) {
  return (total / size) + (total % size > rank);
}

#endif // PUMGEN_AUX_DISTRIBUTOR_H_
