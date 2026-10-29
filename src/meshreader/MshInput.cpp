// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "MshInput.h"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace puml {

MshInput::MshInput(const std::string& fileName, std::size_t bufferSize)
    : fileName(fileName), buffer(std::max(bufferSize, 2 * MaxTokenLength) + 1) {
  file = std::fopen(fileName.c_str(), "rb");
  if (file == nullptr) {
    throw std::runtime_error("Unable to open MSH file " + fileName);
  }
  buffer[0] = '\0';
}

MshInput::~MshInput() {
  if (file != nullptr) {
    std::fclose(file);
  }
}

std::string_view MshInput::readToken() {
  if (!skipWhitespace()) {
    return {};
  }
  ensure(MaxTokenLength);
  std::size_t stop = begin;
  while (stop < end && !isSpace(buffer[stop]) && stop - begin < MaxTokenLength) {
    ++stop;
  }
  const std::string_view token(buffer.data() + begin, stop - begin);
  begin = stop;
  return token;
}

std::optional<std::string> MshInput::readQuoted() {
  if (!skipWhitespace() || buffer[begin] != '"') {
    return std::nullopt;
  }
  ++begin;
  std::string text;
  while (true) {
    while (begin < end) {
      const char c = buffer[begin];
      ++begin;
      if (c == '"') {
        return text;
      }
      if (c == '\n') {
        return std::nullopt;
      }
      text.push_back(c);
    }
    if (!refill()) {
      return std::nullopt;
    }
  }
}

bool MshInput::readLineBreak() {
  ensure(2);
  if (begin < end && buffer[begin] == '\r') {
    ++begin;
  }
  if (begin < end && buffer[begin] == '\n') {
    ++begin;
    return true;
  }
  return false;
}

bool MshInput::readRaw(void* data, std::size_t bytes) {
  auto* out = static_cast<char*>(data);
  while (bytes > 0) {
    if (begin == end && bytes >= (buffer.size() - 1) / 2) {
      // large reads bypass the (empty) buffer
      bufferOffset += end;
      begin = end = 0;
      const std::size_t read = std::fread(out, 1, bytes, file);
      bufferOffset += read;
      if (read < bytes) {
        eof = true;
        return false;
      }
      return true;
    }
    if (begin == end && !refill()) {
      return false;
    }
    const std::size_t count = std::min(bytes, end - begin);
    std::memcpy(out, buffer.data() + begin, count);
    begin += count;
    out += count;
    bytes -= count;
  }
  return true;
}

bool MshInput::seek(std::size_t offset) {
  if (offset >= bufferOffset && offset <= bufferOffset + end) {
    begin = offset - bufferOffset;
    return true;
  }
  if (std::fseek(file, 0, SEEK_END) != 0) {
    return false;
  }
  const auto size = static_cast<std::size_t>(std::ftell(file));
  if (offset > size || std::fseek(file, static_cast<long>(offset), SEEK_SET) != 0) {
    return false;
  }
  bufferOffset = offset;
  begin = end = 0;
  buffer[0] = '\0';
  eof = false;
  return true;
}

bool MshInput::skipTo(std::string_view marker) {
  while (true) {
    const std::string_view data(buffer.data() + begin, end - begin);
    const auto found = data.find(marker);
    if (found != std::string_view::npos) {
      begin += found + marker.size();
      return true;
    }
    // the tail may hold the beginning of the marker
    begin = end - std::min(data.size(), marker.size() - 1);
    if (!refill()) {
      begin = end;
      return false;
    }
  }
}

std::pair<std::size_t, std::size_t> MshInput::location(std::size_t offset) const {
  std::size_t line = 1;
  std::size_t column = 1;
  std::FILE* scan = std::fopen(fileName.c_str(), "rb");
  if (scan == nullptr) {
    return {line, column};
  }
  std::vector<char> chunk(std::size_t{1} << 20);
  std::size_t position = 0;
  while (position < offset) {
    const std::size_t read =
        std::fread(chunk.data(), 1, std::min(chunk.size(), offset - position), scan);
    if (read == 0) {
      break;
    }
    for (std::size_t i = 0; i < read; ++i) {
      if (chunk[i] == '\n') {
        ++line;
        column = 1;
      } else {
        ++column;
      }
    }
    position += read;
  }
  std::fclose(scan);
  return {line, column};
}

bool MshInput::refill() {
  if (eof) {
    return false;
  }
  const std::size_t remaining = end - begin;
  if (begin > 0) {
    std::memmove(buffer.data(), buffer.data() + begin, remaining);
    bufferOffset += begin;
    begin = 0;
    end = remaining;
  }
  const std::size_t capacity = buffer.size() - 1;
  if (end == capacity) {
    return false;
  }
  const std::size_t read = std::fread(buffer.data() + end, 1, capacity - end, file);
  end += read;
  buffer[end] = '\0';
  if (read == 0) {
    eof = true;
  }
  return read > 0;
}

} // namespace puml
