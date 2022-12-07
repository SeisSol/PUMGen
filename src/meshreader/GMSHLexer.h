#ifndef PUMGEN_GMSHLEXER_H
#define PUMGEN_GMSHLEXER_H

#include <fstream>

namespace puml {
enum class GMSHToken {
  eof,
  integer,
  real,
  string,
  mesh_format,
  end_mesh_format,
  entities,
  end_entities,
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

class GMSHLexer {
  public:
  static constexpr std::size_t MaxNumberLength = 128;

  void setIStream(std::istream* istream) {
    in = istream;
    loc = {1, 1};
  };

  virtual GMSHToken getToken() = 0;
  GMSHSourceLocation getSourceLoc() const { return loc; }
  uint64_t getIdentifier() const { return identifier; }
  long getInteger() const { return integer; }
  double getReal() const { return real; }

  protected:
  void advance();
  std::istream* in = nullptr;
  GMSHSourceLocation loc = {1, 1};
  uint64_t identifier;
  long integer;
  double real;
  char lastChar = ' ';
};
} // namespace puml

#endif // PUMGEN_GMSHLEXER_H
