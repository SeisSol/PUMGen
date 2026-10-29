// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "DistributedGMSHReader.h"

#include <algorithm>
#include <array>
#include <cstring>

namespace puml::distributed {

namespace {
template <typename T> T swapBytesOf(T value) {
  std::array<unsigned char, sizeof(T)> bytes{};
  std::memcpy(bytes.data(), &value, sizeof(T));
  std::reverse(bytes.begin(), bytes.end());
  std::memcpy(&value, bytes.data(), sizeof(T));
  return value;
}

std::uint64_t mix(std::uint64_t value) {
  // splitmix64 finalizer
  value ^= value >> 30;
  value *= 0xbf58476d1ce4e5b9ULL;
  value ^= value >> 27;
  value *= 0x94d049bb133111ebULL;
  value ^= value >> 31;
  return value;
}
} // namespace

RawMshFile::RawMshFile(const std::string& fileName, const Msh4Index& index)
    : dataSize(index.dataSize), swapBytes(index.swapBytes) {
  file = std::fopen(fileName.c_str(), "rb");
  if (file == nullptr) {
    logError() << "Unable to open MSH file" << fileName;
  }
}

RawMshFile::~RawMshFile() {
  if (file != nullptr) {
    std::fclose(file);
  }
}

void RawMshFile::readBytes(std::uint64_t offset, std::size_t bytes, void* data) {
  if (std::fseek(file, static_cast<long>(offset), SEEK_SET) != 0 ||
      std::fread(data, 1, bytes, file) != bytes) {
    logError() << "Unexpected end of file in binary data";
  }
}

void RawMshFile::readSizes(std::uint64_t offset, std::size_t count, std::uint64_t* values) {
  if (dataSize == sizeof(std::uint64_t)) {
    readBytes(offset, count * sizeof(std::uint64_t), values);
    if (swapBytes) {
      std::transform(values, values + count, values, swapBytesOf<std::uint64_t>);
    }
  } else {
    narrowSizes.resize(count);
    readBytes(offset, count * sizeof(std::uint32_t), narrowSizes.data());
    for (std::size_t i = 0; i < count; ++i) {
      values[i] = swapBytes ? swapBytesOf(narrowSizes[i]) : narrowSizes[i];
    }
  }
}

void RawMshFile::readDoubles(std::uint64_t offset, std::size_t count, double* values) {
  readBytes(offset, count * sizeof(double), values);
  if (swapBytes) {
    std::transform(values, values + count, values, swapBytesOf<double>);
  }
}

std::uint64_t hashFace(const std::array<std::uint64_t, 4>& face) {
  return mix(mix(mix(mix(face[0]) ^ face[1]) ^ face[2]) ^ face[3]);
}

} // namespace puml::distributed
