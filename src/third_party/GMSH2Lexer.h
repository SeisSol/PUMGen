#ifndef GMSH2LEXER_20200901_H
#define GMSH2LEXER_20200901_H

#include <cstdint>
#include <istream>

#include "meshreader/GMSHLexer.h"

namespace tndm {

class GMSH2Lexer : public puml::GMSHLexer {
public:
    puml::GMSHToken getToken() override;
};

} // namespace tndm

#endif // GMSH2LEXER_20200901_H
