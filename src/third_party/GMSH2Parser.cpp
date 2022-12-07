#include "GMSH2Parser.h"
#include "meshreader/GMSHLexer.h"

#include <cstdio>

namespace tndm {

bool GMSH2Parser::parse_() {
    errorMsg.clear();
    getNextToken();

    double version = parseMeshFormat();
    if (version < 2.0 || version >= 3.0) {
        char buf[128];
        sprintf(buf, "Unsupported MSH version %.1lf", version);
        return logError<bool>(buf);
    }

    bool hasNodes = false;
    bool hasElements = false;

    while (curTok != puml::GMSHToken::eof) {
        switch (curTok) {
        case puml::GMSHToken::nodes:
            hasNodes = parseNodes();
            break;
        case puml::GMSHToken::elements:
            hasElements = parseElements();
            break;
        default:
            getNextToken();
            break;
        }
    }

    return hasNodes && hasElements;
}

bool GMSH2Parser::parseNodes() {
    getNextToken();
    if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 0) {
        return logErrorAnnotated<bool>("Expected non-zero integer");
    }
    std::size_t numVertices = lexer.getInteger();
    builder->setNumVertices(numVertices);

    for (std::size_t i = 0; i < numVertices; ++i) {
        getNextToken();
        if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 1 ||
            lexer.getInteger() > numVertices) {
            char buf[128];
            sprintf(buf, "Expected node-tag with 1 <= node-tag <= %zu", numVertices);
            return logErrorAnnotated<bool>(buf);
        }
        std::size_t id = lexer.getInteger() - 1;

        std::array<double, 3> x;
        for (std::size_t i = 0; i < 3; ++i) {
            getNextToken();
            auto coord = getNumber();
            if (!coord) {
                return logErrorAnnotated<bool>("Expected coordinate");
            }
            x[i] = *coord;
        }
        builder->setVertex(id, x);
    }
    getNextToken();
    if (curTok != puml::GMSHToken::end_nodes) {
        return logErrorAnnotated<bool>("Expected $EndNodes");
    }
    getNextToken();
    return true;
}

bool GMSH2Parser::parseElements() {
    getNextToken();
    if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 0) {
        return logErrorAnnotated<bool>("Expected positive integer");
    }
    const auto numElements = lexer.getInteger();

    constexpr std::size_t MaxNodes = sizeof(NumNodes) / sizeof(std::size_t);
    long tag;
    std::array<long, MaxNodes> nodes;

    builder->setNumElements(numElements);

    for (std::size_t i = 0; i < numElements; ++i) {
        getNextToken();
        if (curTok != puml::GMSHToken::integer) {
            return logErrorAnnotated<bool>("Expected element-tag");
        }

        getNextToken();
        if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 1 ||
            lexer.getInteger() > MaxNodes) {
            char buf[128];
            sprintf(buf, "Expected element-type with 1 <= element-type <= %zu", MaxNodes);
            return logErrorAnnotated<bool>(buf);
        }
        long type = lexer.getInteger();

        getNextToken();
        if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 0) {
            return logErrorAnnotated<bool>("Expected number of tags");
        }
        long numTags = lexer.getInteger();
        for (long i = 0; i < numTags; ++i) {
            getNextToken();
            if (curTok != puml::GMSHToken::integer) {
                return logErrorAnnotated<bool>("Expected tag (integer)");
            }
            if (i == 0) {
                tag = lexer.getInteger();
            }
        }

        for (std::size_t i = 0; i < NumNodes[type - 1]; ++i) {
            getNextToken();
            if (curTok != puml::GMSHToken::integer || lexer.getInteger() < 1) {
                char buf[128];
                sprintf(buf, "Expected node number > 0 (%zu/%zu for type %li)", i + 1,
                        NumNodes[type - 1], type);
                return logErrorAnnotated<bool>(buf);
            }
            nodes[i] = lexer.getInteger() - 1;
        }

        builder->addElement(type, tag, nodes.data(), NumNodes[type - 1]);
    }
    getNextToken();
    if (curTok != puml::GMSHToken::end_elements) {
        return logErrorAnnotated<bool>("Expected $EndElements");
    }
    getNextToken();

    return true;
}

} // namespace tndm
