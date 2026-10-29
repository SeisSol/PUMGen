// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_MSHINPUT_H_
#define PUMGEN_SRC_MESHREADER_MSHINPUT_H_

#include <charconv>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace puml {

/**
 * Buffered sequential reading of an MSH file: whitespace-separated text tokens for the ASCII
 * parts and raw bytes for the binary parts.
 */
class MshInput {
  public:
  static constexpr std::size_t DefaultBufferSize = std::size_t{16} << 20;
  // numbers and section names are shorter; the buffer always holds this many bytes ahead
  static constexpr std::size_t MaxTokenLength = 128;

  /**
   * @throws std::runtime_error if the file cannot be opened
   */
  explicit MshInput(const std::string& fileName, std::size_t bufferSize = DefaultBufferSize);
  ~MshInput();
  MshInput(const MshInput&) = delete;
  MshInput& operator=(const MshInput&) = delete;

  /**
   * Skips whitespace; returns false at the end of the file.
   */
  bool skipWhitespace() {
    while (true) {
      while (begin < end && isSpace(buffer[begin])) {
        ++begin;
      }
      if (begin < end) {
        return true;
      }
      if (!refill()) {
        return false;
      }
    }
  }

  /**
   * The first character of the next token, or '\0' at the end of the file.
   */
  char peek() { return skipWhitespace() ? buffer[begin] : '\0'; }

  /**
   * The next whitespace-separated token (at most MaxTokenLength characters); empty at the end of
   * the file. The view is valid until the next read.
   */
  std::string_view readToken();

  /**
   * Reads a whitespace-separated integer; std::nullopt (without consuming anything but
   * whitespace) if the next token is not an integer of type T.
   */
  template <typename T> std::optional<T> readInteger() {
    if (!skipWhitespace()) {
      return std::nullopt;
    }
    ensure(MaxTokenLength);
    const char* first = buffer.data() + begin;
    const char* last = buffer.data() + end;
    if (*first == '+') {
      ++first;
    }
    T value{};
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc() || (ptr != last && !isSpace(*ptr))) {
      return std::nullopt;
    }
    begin = static_cast<std::size_t>(ptr - buffer.data());
    return value;
  }

  /**
   * Reads a whitespace-separated floating-point number; std::nullopt (without consuming anything
   * but whitespace) if the next token is not a number.
   */
  std::optional<double> readReal() {
    if (!skipWhitespace()) {
      return std::nullopt;
    }
    ensure(MaxTokenLength);
    const char* first = buffer.data() + begin;
    const char* last = buffer.data() + end;
    if (*first == '+') {
      ++first;
    }
    double value = 0;
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
    const auto [ptr, ec] = std::from_chars(first, last, value);
    const bool parsed = ec == std::errc();
#else
    // the buffer is null-terminated behind the valid data
    char* ptr = nullptr;
    value = std::strtod(first, &ptr);
    const bool parsed = ptr != first;
#endif
    if (!parsed || (ptr != last && !isSpace(*ptr))) {
      return std::nullopt;
    }
    begin = static_cast<std::size_t>(ptr - buffer.data());
    return value;
  }

  /**
   * Reads a text in double quotes, which ends on its line, and returns it without the quotes.
   */
  std::optional<std::string> readQuoted();

  /**
   * Consumes a single line break ("\n" or "\r\n"), as it separates the header line of a binary
   * section from its data.
   */
  bool readLineBreak();

  /**
   * Reads the next bytes as they are.
   */
  bool readRaw(void* data, std::size_t bytes);

  /**
   * Continues reading at the given file offset; returns false if the offset lies behind the end
   * of the file.
   */
  bool seek(std::size_t offset);

  /**
   * Skips everything up to and including the next occurrence of the marker.
   */
  bool skipTo(std::string_view marker);

  /**
   * The offset of the next unread byte in the file.
   */
  [[nodiscard]] std::size_t offset() const { return bufferOffset + begin; }

  /**
   * Line and column (both starting at 1) of a byte offset; reads the file up to that offset.
   */
  [[nodiscard]] std::pair<std::size_t, std::size_t> location(std::size_t offset) const;

  private:
  static bool isSpace(char c) { return c == ' ' || (c >= '\t' && c <= '\r'); }

  /**
   * Moves the unread data to the front of the buffer and appends more data from the file;
   * returns false if nothing could be appended.
   */
  bool refill();

  /**
   * Makes at least count bytes available unless the end of the file comes first.
   */
  void ensure(std::size_t count) {
    while (end - begin < count && !eof) {
      refill();
    }
  }

  std::string fileName;
  std::FILE* file = nullptr;
  // one byte more than the capacity for the terminating null character
  std::vector<char> buffer;
  // the unread data is buffer[begin, end)
  std::size_t begin = 0;
  std::size_t end = 0;
  // the file offset of buffer[0]
  std::size_t bufferOffset = 0;
  bool eof = false;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_MSHINPUT_H_
