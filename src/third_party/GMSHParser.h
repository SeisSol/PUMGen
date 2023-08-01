#ifndef GMSHPARSER_20200901_H
#define GMSHPARSER_20200901_H

#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>

#include "GMSHMeshBuilder.h"
#include "GMSHLexer.h"

namespace tndm {

class GMSHParser {
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
    static constexpr const char* ElementTypes[] = {
        "line",
        "triangle",
        "quadrangle",
        "tetrahedron",
        "hexahedron",
        "prism",
        "pyramid",
        "P1 line",
        "P1 triangle",
        "P1 quadrangle",
        "P1 tetrahedron",
        "P2 hexahedron",
        "P2 prism",
        "P2 pyramid",
        "point",
    };

    explicit GMSHParser(GMSHMeshBuilder* builder) : builder(builder){};
    [[nodiscard]] std::string_view getErrorMessage() const { return errorMsg; }

    bool parseFile(std::string const& fileName) {
        std::ifstream in(fileName);
        if (!in.is_open()) {
            return logError<bool>("Unable to open MSH file");
        }
        lexer.setIStream(&in);
        return parse_();
    }

protected:
    template <typename T> T logErrorAnnotated(std::string_view msg) {
        std::stringstream ss;
        ss << "GMSH parser error in line " << curLoc.line << " in column " << curLoc.col << ":\n";
        ss << '\t' << msg << '\n';
        errorMsg += ss.str();
        return {};
    }

    template <typename T> T logError(std::string_view msg) {
        errorMsg += "GMSH parser error:\n\t";
        errorMsg += msg;
        errorMsg += '\n';
        return {};
    }

    GMSHMeshBuilder* builder;
    GMSHSourceLocation curLoc = {0, 0};
    GMSHToken curTok;
    GMSHLexer lexer;
    std::string errorMsg;

    GMSHToken getNextToken() {
        curLoc = lexer.getSourceLoc();
        return curTok = lexer.getToken();
    }

    std::optional<double> getNumber() {
        if (curTok == GMSHToken::integer) {
            return {lexer.getInteger()};
        } else if (curTok == GMSHToken::real) {
            return {lexer.getReal()};
        }
        return std::nullopt;
    }

    virtual bool parse_() = 0;

    double parseMeshFormat();
};
} // namespace tndm

#endif // GMSHPARSER_20200901_H
