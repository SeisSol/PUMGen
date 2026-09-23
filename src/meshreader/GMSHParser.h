// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSHPARSER_H_
#define PUMGEN_SRC_MESHREADER_GMSHPARSER_H_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

#include "GMSHMeshBuilder.h"
#include "MshInput.h"

namespace puml {

namespace detail {
template <std::size_t N> constexpr std::size_t maxEntry(const std::size_t (&values)[N]) {
  std::size_t result = 0;
  for (std::size_t i = 0; i < N; ++i) {
    result = values[i] > result ? values[i] : result;
  }
  return result;
}
} // namespace detail

/**
 * Common parts of the MSH parsers: element types, error handling and reading of text values.
 */
class GMSHParser {
  public:
  // Look-up table from gmsh type to number of nodes; 0 for the type numbers gmsh does not define
  static constexpr std::size_t NumNodes[] = {
      2,    // 1: Line 2
      3,    // 2: Triangle 3
      4,    // 3: Quadrilateral 4
      4,    // 4: Tetrahedron 4
      8,    // 5: Hexahedron 8
      6,    // 6: Prism 6
      5,    // 7: Pyramid 5
      3,    // 8: Line 3
      6,    // 9: Triangle 6
      9,    // 10: Quadrilateral 9
      10,   // 11: Tetrahedron 10
      27,   // 12: Hexahedron 27
      18,   // 13: Prism 18
      14,   // 14: Pyramid 14
      1,    // 15: Point
      8,    // 16: Quadrilateral 8
      20,   // 17: Hexahedron 20
      15,   // 18: Prism 15
      13,   // 19: Pyramid 13
      9,    // 20: Triangle 9
      10,   // 21: Triangle 10
      12,   // 22: Triangle 12
      15,   // 23: Triangle 15
      15,   // 24: Triangle 15I
      21,   // 25: Triangle 21
      4,    // 26: Line 4
      5,    // 27: Line 5
      6,    // 28: Line 6
      20,   // 29: Tetrahedron 20
      35,   // 30: Tetrahedron 35
      56,   // 31: Tetrahedron 56
      22,   // 32: Tetrahedron 22
      28,   // 33: Tetrahedron 28
      0,    // 34: Polygon
      0,    // 35: Polyhedron
      16,   // 36: Quadrilateral 16
      25,   // 37: Quadrilateral 25
      36,   // 38: Quadrilateral 36
      12,   // 39: Quadrilateral 12
      16,   // 40: Quadrilateral 16I
      20,   // 41: Quadrilateral 20
      28,   // 42: Triangle 28
      36,   // 43: Triangle 36
      45,   // 44: Triangle 45
      55,   // 45: Triangle 55
      66,   // 46: Triangle 66
      49,   // 47: Quadrilateral 49
      64,   // 48: Quadrilateral 64
      81,   // 49: Quadrilateral 81
      100,  // 50: Quadrilateral 100
      121,  // 51: Quadrilateral 121
      18,   // 52: Triangle 18
      21,   // 53: Triangle 21I
      24,   // 54: Triangle 24
      27,   // 55: Triangle 27
      30,   // 56: Triangle 30
      24,   // 57: Quadrilateral 24
      28,   // 58: Quadrilateral 28
      32,   // 59: Quadrilateral 32
      36,   // 60: Quadrilateral 36I
      40,   // 61: Quadrilateral 40
      7,    // 62: Line 7
      8,    // 63: Line 8
      9,    // 64: Line 9
      10,   // 65: Line 10
      11,   // 66: Line 11
      0,    // 67: undefined
      0,    // 68: undefined
      0,    // 69: Polygon Border
      0,    // 70: undefined
      84,   // 71: Tetrahedron 84
      120,  // 72: Tetrahedron 120
      165,  // 73: Tetrahedron 165
      220,  // 74: Tetrahedron 220
      286,  // 75: Tetrahedron 286
      0,    // 76: undefined
      0,    // 77: undefined
      0,    // 78: undefined
      34,   // 79: Tetrahedron 34
      40,   // 80: Tetrahedron 40
      46,   // 81: Tetrahedron 46
      52,   // 82: Tetrahedron 52
      58,   // 83: Tetrahedron 58
      1,    // 84: Line 1
      1,    // 85: Triangle 1
      1,    // 86: Quadrilateral 1
      1,    // 87: Tetrahedron 1
      1,    // 88: Hexahedron 1
      1,    // 89: Prism 1
      40,   // 90: Prism 40
      75,   // 91: Prism 75
      64,   // 92: Hexahedron 64
      125,  // 93: Hexahedron 125
      216,  // 94: Hexahedron 216
      343,  // 95: Hexahedron 343
      512,  // 96: Hexahedron 512
      729,  // 97: Hexahedron 729
      1000, // 98: Hexahedron 1000
      32,   // 99: Hexahedron 32
      44,   // 100: Hexahedron 44
      56,   // 101: Hexahedron 56
      68,   // 102: Hexahedron 68
      80,   // 103: Hexahedron 80
      92,   // 104: Hexahedron 92
      104,  // 105: Hexahedron 104
      126,  // 106: Prism 126
      196,  // 107: Prism 196
      288,  // 108: Prism 288
      405,  // 109: Prism 405
      550,  // 110: Prism 550
      24,   // 111: Prism 24
      33,   // 112: Prism 33
      42,   // 113: Prism 42
      51,   // 114: Prism 51
      60,   // 115: Prism 60
      69,   // 116: Prism 69
      78,   // 117: Prism 78
      30,   // 118: Pyramid 30
      55,   // 119: Pyramid 55
      91,   // 120: Pyramid 91
      140,  // 121: Pyramid 140
      204,  // 122: Pyramid 204
      285,  // 123: Pyramid 285
      385,  // 124: Pyramid 385
      21,   // 125: Pyramid 21
      29,   // 126: Pyramid 29
      37,   // 127: Pyramid 37
      45,   // 128: Pyramid 45
      53,   // 129: Pyramid 53
      61,   // 130: Pyramid 61
      69,   // 131: Pyramid 69
      1,    // 132: Pyramid 1
      0,    // 133: Point Xfem
      0,    // 134: Line Xfem
      0,    // 135: Triangle Xfem
      0,    // 136: Tetrahedron Xfem
      16,   // 137: Tetrahedron 16
  };
  static constexpr std::size_t NumTypes = sizeof(NumNodes) / sizeof(NumNodes[0]);
  static constexpr std::size_t MaxNodesPerElement = detail::maxEntry(NumNodes);

  explicit GMSHParser(GMSHMeshBuilder* builder,
                      std::size_t bufferSize = MshInput::DefaultBufferSize)
      : builder(builder), bufferSize(bufferSize) {}
  virtual ~GMSHParser() = default;

  /**
   * Parses the file; stops at the first error, whose description getErrorMessage() returns.
   */
  bool parseFile(const std::string& fileName);

  [[nodiscard]] std::string_view getErrorMessage() const { return errorMsg; }

  protected:
  struct MeshFormat {
    double version;
    long fileType;
    std::size_t dataSize;
  };

  /**
   * Thrown after an error has been recorded in errorMsg.
   */
  struct ParseError {};

  /**
   * Records an error at the given file offset and aborts parsing.
   */
  [[noreturn]] void failAt(std::size_t offset, std::string_view message);

  /**
   * Records an error at the next token and aborts parsing.
   */
  [[noreturn]] void fail(std::string_view message);

  /**
   * Records an error concerning the whole file and aborts parsing.
   */
  [[noreturn]] void failFile(std::string_view message);

  std::size_t expectSize() {
    if (const auto value = input->readInteger<std::size_t>()) {
      return *value;
    }
    fail("Expected non-negative integer");
  }

  long expectInteger() {
    if (const auto value = input->readInteger<long>()) {
      return *value;
    }
    fail("Expected integer");
  }

  double expectNumber() {
    if (const auto value = input->readReal()) {
      return *value;
    }
    fail("Expected number");
  }

  void expectToken(std::string_view token);

  /**
   * Reads the $MeshFormat section up to (and including) the data size.
   */
  MeshFormat parseMeshFormatHeader();

  /**
   * Returns the name of the next section (e.g. "$Nodes"), or an empty view at the end of the file.
   */
  std::string_view nextSection();

  /**
   * Skips the rest of a section which is not needed.
   */
  void skipSection(std::string_view section);

  virtual void parse_() = 0;

  GMSHMeshBuilder* builder;
  std::unique_ptr<MshInput> input;
  // errors in binary data are located by byte offset instead of line and column
  bool binaryData = false;

  private:
  std::size_t bufferSize;
  std::string errorMsg;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSHPARSER_H_
