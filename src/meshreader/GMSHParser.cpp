// SPDX-FileCopyrightText: 2022 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#include "GMSHParser.h"

#include <exception>
#include <sstream>
#include <string>

namespace puml {

bool GMSHParser::parseFile(const std::string& fileName) {
  errorMsg.clear();
  try {
    input = std::make_unique<MshInput>(fileName, bufferSize);
    parse_();
    input.reset();
    return true;
  } catch (const ParseError&) {
  } catch (const std::exception& error) {
    errorMsg += "GMSH parser error:\n\t";
    errorMsg += error.what();
    errorMsg += '\n';
  }
  input.reset();
  return false;
}

void GMSHParser::failAt(std::size_t offset, std::string_view message) {
  const auto [line, column] = input->location(offset);
  std::stringstream stream;
  stream << "GMSH parser error in line " << line << " in column " << column << ":\n";
  stream << '\t' << message << '\n';
  errorMsg += stream.str();
  throw ParseError{};
}

void GMSHParser::fail(std::string_view message) {
  input->skipWhitespace();
  failAt(input->offset(), message);
}

void GMSHParser::failFile(std::string_view message) {
  errorMsg += "GMSH parser error:\n\t";
  errorMsg += message;
  errorMsg += '\n';
  throw ParseError{};
}

void GMSHParser::expectToken(std::string_view token) {
  input->skipWhitespace();
  const auto offset = input->offset();
  if (input->readToken() != token) {
    failAt(offset, "Expected " + std::string(token));
  }
}

GMSHParser::MeshFormat GMSHParser::parseMeshFormatHeader() {
  expectToken("$MeshFormat");
  MeshFormat format{};
  format.version = expectNumber();
  format.fileType = expectInteger();
  format.dataSize = expectSize();
  return format;
}

std::string_view GMSHParser::nextSection() {
  if (!input->skipWhitespace()) {
    return {};
  }
  const auto offset = input->offset();
  const auto section = input->readToken();
  if (section.empty() || section.front() != '$') {
    failAt(offset, "Expected a section");
  }
  return section;
}

void GMSHParser::skipSection(std::string_view section) {
  const std::string endMarker = "$End" + std::string(section.substr(1));
  if (!input->skipTo(endMarker)) {
    failFile("Missing " + endMarker);
  }
}

} // namespace puml
