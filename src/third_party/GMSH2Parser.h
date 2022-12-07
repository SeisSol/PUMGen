#ifndef GMSH2PARSER_20200901_H
#define GMSH2PARSER_20200901_H

#include "GMSH2Lexer.h"
#include "GMSHBuilder.h"

#include <array>
#include <optional>
#include <streambuf>
#include <string>
#include <string_view>

namespace tndm {

struct membuf : std::streambuf {
    membuf(char* b, char* e) { this->setg(b, b, e); }
};

class GMSH2Parser {
private:
    GMSHMeshBuilder* builder;
    GMSHToken curTok;
    GMSHSourceLocation curLoc;
    GMSH2Lexer lexer;
    std::string errorMsg;

    GMSHToken getNextToken() {
        curLoc = lexer.getSourceLoc();
        return curTok = lexer.getToken();
    }

    template <typename T> T logError(std::string_view msg);
    template <typename T> T logErrorAnnotated(std::string_view msg);
    std::optional<double> getNumber();
    double parseMeshFormat();
    bool parseNodes();
    bool parseElements();
    bool parse_();

public:
    // Look-up table from gmsh type to number of nodes
    static constexpr std::size_t NumNodes[] = {
        2,  // line
        3,  // triangle,
        4,  // quadrangle,
        4,  // tetrahedron
        8,  // hexahedron
        6,  // prism
        5,  // pyramid
        3,  // P1 line
        6,  // P1 triangle
        9,  // P1 quadrangle
        10, // P1 tetrahedron
        27, // P2 hexahedron
        18, // P2 prism
        14, // P2 pyramid
        1,  // point
    };

    GMSH2Parser(GMSHMeshBuilder* builder) : builder(builder) {}

    bool parse(std::string& msh);
    bool parse(char* msh, std::size_t len);
    bool parseFile(std::string const& fileName);

    std::string_view getErrorMessage() const { return errorMsg; }
};

}; // namespace tndm

#endif // GMSH2PARSER_20200901_H
