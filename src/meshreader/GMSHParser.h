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
  // Look-up table from gmsh type to number of nodes
  static constexpr std::size_t NumNodes[] = {
      2,    // MSH_LIN_2
      3,    // MSH_TRI_3
      4,    // MSH_QUA_4
      4,    // MSH_TET_4
      8,    // MSH_HEX_8
      6,    // MSH_PRI_6
      5,    // MSH_PYR_5
      3,    // MSH_LIN_3
      6,    // MSH_TRI_6
      9,    // MSH_QUA_9
      10,   // MSH_TET_10
      27,   // MSH_HEX_27
      18,   // MSH_PRI_18
      14,   // MSH_PYR_14
      1,    // MSH_PNT
      8,    // MSH_QUA_8
      20,   // MSH_HEX_20
      15,   // MSH_PRI_15
      13,   // MSH_PYR_13
      9,    // MSH_TRI_9
      10,   // MSH_TRI_10
      12,   // MSH_TRI_12
      15,   // MSH_TRI_15
      15,   // MSH_TRI_15I
      21,   // MSH_TRI_21
      4,    // MSH_LIN_4
      5,    // MSH_LIN_5
      6,    // MSH_LIN_6
      20,   // MSH_TET_20
      35,   // MSH_TET_35
      56,   // MSH_TET_56
      22,   // MSH_TET_22
      28,   // MSH_TET_28
      0,    // MSH_POLYG_
      0,    // MSH_POLYH_
      16,   // MSH_QUA_16
      25,   // MSH_QUA_25
      36,   // MSH_QUA_36
      12,   // MSH_QUA_12
      16,   // MSH_QUA_16I
      20,   // MSH_QUA_20
      28,   // MSH_TRI_28
      36,   // MSH_TRI_36
      45,   // MSH_TRI_45
      55,   // MSH_TRI_55
      66,   // MSH_TRI_66
      49,   // MSH_QUA_49
      64,   // MSH_QUA_64
      81,   // MSH_QUA_81
      100,  // MSH_QUA_100
      121,  // MSH_QUA_121
      18,   // MSH_TRI_18
      21,   // MSH_TRI_21I
      24,   // MSH_TRI_24
      27,   // MSH_TRI_27
      30,   // MSH_TRI_30
      24,   // MSH_QUA_24
      28,   // MSH_QUA_28
      32,   // MSH_QUA_32
      36,   // MSH_QUA_36I
      40,   // MSH_QUA_40
      7,    // MSH_LIN_7
      8,    // MSH_LIN_8
      9,    // MSH_LIN_9
      10,   // MSH_LIN_10
      11,   // MSH_LIN_11
      0,    // MSH_LIN_B
      0,    // MSH_TRI_B
      0,    // MSH_POLYG_B
      0,    // MSH_LIN_C
      84,   // MSH_TET_84
      120,  // MSH_TET_120
      165,  // MSH_TET_165
      220,  // MSH_TET_220
      286,  // MSH_TET_286
      34,   // MSH_TET_34
      40,   // MSH_TET_40
      46,   // MSH_TET_46
      52,   // MSH_TET_52
      58,   // MSH_TET_58
      1,    // MSH_LIN_1
      1,    // MSH_TRI_1
      1,    // MSH_QUA_1
      1,    // MSH_TET_1
      1,    // MSH_HEX_1
      1,    // MSH_PRI_1
      40,   // MSH_PRI_40
      75,   // MSH_PRI_75
      64,   // MSH_HEX_64
      125,  // MSH_HEX_125
      216,  // MSH_HEX_216
      343,  // MSH_HEX_343
      512,  // MSH_HEX_512
      729,  // MSH_HEX_729
      1000, // MSH_HEX_1000
      32,   // MSH_HEX_32
      44,   // MSH_HEX_44
      56,   // MSH_HEX_56
      68,   // MSH_HEX_68
      80,   // MSH_HEX_80
      92,   // MSH_HEX_92
      104,  // MSH_HEX_104
      126,  // MSH_PRI_126
      196,  // MSH_PRI_196
      288,  // MSH_PRI_288
      405,  // MSH_PRI_405
      550,  // MSH_PRI_550
      24,   // MSH_PRI_24
      33,   // MSH_PRI_33
      42,   // MSH_PRI_42
      51,   // MSH_PRI_51
      60,   // MSH_PRI_60
      69,   // MSH_PRI_69
      78,   // MSH_PRI_78
      30,   // MSH_PYR_30
      55,   // MSH_PYR_55
      91,   // MSH_PYR_91
      140,  // MSH_PYR_140
      204,  // MSH_PYR_204
      285,  // MSH_PYR_285
      385,  // MSH_PYR_385
      21,   // MSH_PYR_21
      29,   // MSH_PYR_29
      37,   // MSH_PYR_37
      45,   // MSH_PYR_45
      53,   // MSH_PYR_53
      61,   // MSH_PYR_61
      69,   // MSH_PYR_69
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

  private:
  std::size_t bufferSize;
  std::string errorMsg;
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSHPARSER_H_
