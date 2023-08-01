#ifndef GMSH2PARSER_20200901_H
#define GMSH2PARSER_20200901_H

#include "GMSHParser.h"
#include "GMSHLexer.h"
#include "GMSHMeshBuilder.h"

#include <array>
#include <optional>
#include <streambuf>
#include <string>
#include <string_view>

namespace tndm {

class GMSH2Parser: public GMSHParser {
public:
    using GMSHParser::GMSHParser;
private:
    bool parseNodes();
    bool parseElements();
    virtual bool parse_() override;
};

}; // namespace tndm

#endif // GMSH2PARSER_20200901_H
