// SPDX-FileCopyrightText: 2025 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_THIRD_PARTY_HASH_H_
#define PUMGEN_SRC_THIRD_PARTY_HASH_H_

#include <cstdint>
#include <string>

namespace tndm {

// https://de.wikipedia.org/wiki/FNV_(Informatik)
constexpr uint64_t fnv1a0() { return 0xcbf29ce484222325; }
constexpr uint64_t fnv1a_step(uint64_t hash, char c) { return (hash ^ c) * 0x00000100000001b3; }
constexpr uint64_t fnv1a(char const* s, std::size_t len) {
  return len > 0 ? fnv1a_step(fnv1a(s, len - 1), s[len - 1]) : fnv1a0();
}
uint64_t fnv1a(std::string const& s) { return fnv1a(s.data(), s.size()); }
constexpr uint64_t operator""_fnv1a(char const* s, std::size_t len) { return tndm::fnv1a(s, len); }

} // namespace tndm

#endif // PUMGEN_SRC_THIRD_PARTY_HASH_H_
