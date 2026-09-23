// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

// Parses the fixtures with a tiny read buffer, so that tokens and section markers straddle buffer
// boundaries, and compares the result with the one of the default buffer.

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "meshreader/GMSH2Parser.h"
#include "meshreader/GMSH4Parser.h"
#include "meshreader/GMSHBuilder.h"
#include "meshreader/MshInput.h"

namespace {

constexpr std::size_t TinyBuffer = 2 * puml::MshInput::MaxTokenLength;

template <typename Parser, typename Builder>
bool parse(const std::string& file, std::size_t bufferSize, Builder& builder) {
  Parser parser(&builder, bufferSize);
  if (!parser.parseFile(file)) {
    std::fprintf(stderr, "%s (buffer %zu): %s", file.c_str(), bufferSize,
                 std::string(parser.getErrorMessage()).c_str());
    return false;
  }
  builder.postprocess();
  return true;
}

template <typename Parser> bool check(const std::string& file) {
  puml::GMSHBuilder reference;
  puml::GMSHBuilder tiny;
  if (!parse<Parser>(file, puml::MshInput::DefaultBufferSize, reference) ||
      !parse<Parser>(file, TinyBuffer, tiny)) {
    return false;
  }
  const bool same = reference.vertices == tiny.vertices && reference.cellNodes == tiny.cellNodes &&
                    reference.cellTypes == tiny.cellTypes &&
                    reference.cellOrders == tiny.cellOrders && reference.groups == tiny.groups &&
                    reference.facets == tiny.facets && reference.bcs == tiny.bcs &&
                    reference.identify == tiny.identify &&
                    reference.physicalNames == tiny.physicalNames && !reference.cellTypes.empty();
  std::printf("%s: %s\n", file.c_str(), same ? "ok" : "DIFFERENT");
  return same;
}

std::string readFile(const std::string& file) {
  std::ifstream in(file, std::ios::binary);
  std::stringstream content;
  content << in.rdbuf();
  return content.str();
}

// a section to be skipped, longer than the tiny buffer and with partial end markers
bool checkSkippedSection(const std::string& fixtures) {
  std::string content = readFile(fixtures + "/tiny-v41.msh");
  std::string junk = "$Comments\n";
  for (int i = 0; i < 200; ++i) {
    junk += "$EndComment $EndCommentz ";
  }
  junk += "\n$EndComments\n";
  content.insert(content.find("$Nodes"), junk);
  const std::string file = "reader-test-skipped-section.msh";
  std::ofstream(file, std::ios::binary) << content;

  puml::GMSHBuilder reference;
  puml::GMSHBuilder skipped;
  const bool ok = parse<puml::GMSH4Parser>(fixtures + "/tiny-v41.msh",
                                           puml::MshInput::DefaultBufferSize, reference) &&
                  parse<puml::GMSH4Parser>(file, TinyBuffer, skipped) &&
                  reference.vertices == skipped.vertices &&
                  reference.cellNodes == skipped.cellNodes && reference.groups == skipped.groups;
  std::remove(file.c_str());
  std::printf("skipped section: %s\n", ok ? "ok" : "DIFFERENT");
  return ok;
}

// names in double quotes, with spaces, across the ends of the tiny buffer
bool checkPhysicalNames(const std::string& file) {
  puml::GMSHBuilder builder;
  const bool parsed = parse<puml::GMSH4Parser>(file, TinyBuffer, builder);
  const puml::PhysicalName surface{2, 101, "free surface"};
  const bool ok = parsed && builder.physicalNames.size() == 7 &&
                  std::find(builder.physicalNames.begin(), builder.physicalNames.end(), surface) !=
                      builder.physicalNames.end();
  std::printf("physical names: %s\n", ok ? "ok" : "DIFFERENT");
  return ok;
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "Usage: %s <fixture directory>\n", argv[0]);
    return 2;
  }
  const std::string fixtures = argv[1];
  bool ok = true;
  ok &= check<puml::GMSH4Parser>(fixtures + "/layered-v41.msh");
  ok &= check<puml::GMSH2Parser>(fixtures + "/layered-v22.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/periodic-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/periodic-parametric-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/coarse-o2-v41.msh");
  ok &= check<puml::GMSH2Parser>(fixtures + "/coarse-o2-v22.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/layered-binary-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/periodic-binary-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/coarse-binary-bigendian-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/coarse-binary-size4-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/coarse-o2-binary-v41.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/mixed-v41.msh");
  ok &= check<puml::GMSH2Parser>(fixtures + "/mixed-v22.msh");
  ok &= check<puml::GMSH4Parser>(fixtures + "/mixed-binary-v41.msh");
  ok &= checkPhysicalNames(fixtures + "/mixed-v41.msh");
  ok &= checkSkippedSection(fixtures);
  return ok ? 0 : 1;
}
