#ifndef PUMGEN_GMSH4LEXER_H
#define PUMGEN_GMSH4LEXER_H

#include <fstream>

namespace puml {

enum class GMSHToken {
  eof,
  integer,
  real,
  string,
  mesh_format,
  end_mesh_format,
  nodes,
  end_nodes,
  elements,
  end_elements,
  unknown_section,
  unknown_token
};

struct GMSHSourceLocation {
  std::size_t line;
  std::size_t col;
};

class GMSH4Lexer {
  public:
  void setIStream(std::istream* istream) {
    in = istream;
    loc = {1, 1};
  };

  private:
  std::istream* in = nullptr;
  GMSHSourceLocation loc = {1, 1};
};
} // namespace puml

#endif // PUMGEN_GMSH4LEXER_H
