#include "GMSHLexer.h"
#include "third_party/Hash.h"

namespace puml {
void GMSHLexer::advance() {
  in->get(lastChar);

  if (lastChar == '\n' || lastChar == '\r') {
    ++loc.line;
    loc.col = 1;
  } else {
    ++loc.col;
  }
}

puml::GMSHToken GMSHLexer::getToken() {
  if (in == nullptr) {
    return puml::GMSHToken::eof;
  }
  in->peek();
  if (!in->good()) {
    return puml::GMSHToken::eof;
  }

  while (isspace(lastChar)) {
    advance();
  }

  if (lastChar == '$') {
    advance();
    auto hash = tndm::fnv1a0();
    while (isalpha(lastChar)) {
      hash = tndm::fnv1a_step(hash, lastChar);
      advance();
    }
    puml::GMSHToken token = puml::GMSHToken::unknown_section;
    using tndm::operator""_fnv1a;
    switch (hash) {
    case "MeshFormat"_fnv1a:
      token = puml::GMSHToken::mesh_format;
      break;
    case "EndMeshFormat"_fnv1a:
      token = puml::GMSHToken::end_mesh_format;
      break;
    case "Entities"_fnv1a:
      token = puml::GMSHToken::entities;
      break;
    case "EndEntities"_fnv1a:
      token = puml::GMSHToken::end_entities;
      break;
    case "Nodes"_fnv1a:
      token = puml::GMSHToken::nodes;
      break;
    case "EndNodes"_fnv1a:
      token = puml::GMSHToken::end_nodes;
      break;
    case "Elements"_fnv1a:
      token = puml::GMSHToken::elements;
      break;
    case "EndElements"_fnv1a:
      token = puml::GMSHToken::end_elements;
      break;
    default:
      break;
    }
    return token;
  }

  auto mustbereal = [](char c) { return c == '.' || c == 'e' || c == 'E'; };
  auto isnumber = [&mustbereal](char c) {
    return isdigit(c) || c == '+' || c == '-' || mustbereal(c);
  };

  if (isnumber(lastChar)) {
    char buf[MaxNumberLength + 1];
    int pos = 0;
    bool isreal = false;
    do {
      buf[pos++] = lastChar;
      isreal = isreal || mustbereal(lastChar);
      advance();
    } while (isnumber(lastChar) && pos < MaxNumberLength);
    buf[pos] = 0;
    if (pos == MaxNumberLength) {
      throw std::runtime_error("Too large number encountered in GMSHLexer: " + std::string(buf));
    }
    if (isreal) {
      real = std::strtod(buf, 0);
      return puml::GMSHToken::real;
    }
    integer = std::strtol(buf, 0, 10);
    return puml::GMSHToken::integer;
  }

  if (lastChar == '"') {
    do {
      advance();
    } while (lastChar != '"');
    advance();
    return puml::GMSHToken::string;
  }
  return puml::GMSHToken::unknown_token;
}

} // namespace puml
