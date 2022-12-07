#include "GMSH2Lexer.h"
#include "Hash.h"

#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace tndm {

puml::GMSHToken GMSH2Lexer::getToken() {
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
        auto hash = fnv1a0();
        while (isalpha(lastChar)) {
            hash = fnv1a_step(hash, lastChar);
            advance();
        }
        puml::GMSHToken token = puml::GMSHToken::unknown_section;
        switch (hash) {
        case "MeshFormat"_fnv1a:
            token = puml::GMSHToken::mesh_format;
            break;
        case "EndMeshFormat"_fnv1a:
            token = puml::GMSHToken::end_mesh_format;
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
            throw std::runtime_error("Too large number encountered in GMSH2Lexer: " +
                                     std::string(buf));
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

} // namespace tndm
