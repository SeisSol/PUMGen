// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_INTEGERWIDTH_H_
#define PUMGEN_SRC_HELPER_INTEGERWIDTH_H_

#include <cstddef>

/**
 * The size of the smallest integer of 1, 2, 4 or 8 bytes with at least the given number of bits.
 * Readers such as XDMF and VTK know integers of these sizes only.
 */
constexpr std::size_t compactIntegerBytes(std::size_t bits) {
  std::size_t bytes = 1;
  while (bytes * 8 < bits) {
    bytes *= 2;
  }
  return bytes;
}

#endif // PUMGEN_SRC_HELPER_INTEGERWIDTH_H_
