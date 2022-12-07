#ifndef PUMGEN_GMSH4LEXER_H
#define PUMGEN_GMSH4LEXER_H

#include "GMSHLexer.h"

namespace puml {

class GMSH4Lexer : public GMSHLexer {
  puml::GMSHToken getToken() override{};
};
} // namespace puml

#endif // PUMGEN_GMSH4LEXER_H
