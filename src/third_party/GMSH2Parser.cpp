#include "GMSH2Parser.h"
#include "GMSHLexer.h"

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

    while (curTok != GMSHToken::eof) {
        switch (curTok) {
        case GMSHToken::nodes:
            hasNodes = parseNodes();
            break;
        case GMSHToken::elements:
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
    if (curTok != GMSHToken::integer || lexer.getInteger() < 0) {
        return logErrorAnnotated<bool>("Expected non-zero integer");
    }
    std::size_t numVertices = lexer.getInteger();
    builder->setNumVertices(numVertices);

    for (std::size_t i = 0; i < numVertices; ++i) {
        getNextToken();
        if (curTok != GMSHToken::integer || lexer.getInteger() < 1 ||
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
    if (curTok != GMSHToken::end_nodes) {
        return logErrorAnnotated<bool>("Expected $EndNodes");
    }
    getNextToken();
    return true;
}

bool GMSH2Parser::parseElements() {
    getNextToken();
    if (curTok != GMSHToken::integer || lexer.getInteger() < 0) {
        return logErrorAnnotated<bool>("Expected positive integer");
    }
    const auto numElements = lexer.getInteger();

    constexpr std::size_t MaxNodes = sizeof(NumNodes) / sizeof(std::size_t);
    long tag{};
    std::array<long, MaxNodes> nodes;

    builder->setNumElements(numElements);

    for (std::size_t elementIdx = 0; elementIdx < numElements; ++elementIdx) {
        getNextToken();
        if (curTok != GMSHToken::integer) {
            return logErrorAnnotated<bool>("Expected element-tag");
        }

        getNextToken();
        if (curTok != GMSHToken::integer || lexer.getInteger() < 1 ||
            lexer.getInteger() > MaxNodes) {
            char buf[128];
            sprintf(buf, "Expected element-type with 1 <= element-type <= %zu", MaxNodes);
            return logErrorAnnotated<bool>(buf);
        }
        long type = lexer.getInteger();

        getNextToken();
        if (curTok != GMSHToken::integer || lexer.getInteger() < 0) {
            return logErrorAnnotated<bool>("Expected number of tags");
        }
        long numTags = lexer.getInteger();
        for (long tagIdx = 0; tagIdx < numTags; ++tagIdx) {
            getNextToken();
            if (curTok != GMSHToken::integer) {
                return logErrorAnnotated<bool>("Expected tag (integer)");
            }
            if (tagIdx == 0) {
                tag = lexer.getInteger();
            }
        }

        for (std::size_t nodeIdx = 0; nodeIdx < NumNodes[type - 1]; ++nodeIdx) {
            getNextToken();
            if (curTok != GMSHToken::integer || lexer.getInteger() < 1) {
                char buf[128];
                sprintf(buf, "Expected node number > 0 (%zu/%zu for type %li)", nodeIdx + 1,
                        NumNodes[type - 1], type);
                return logErrorAnnotated<bool>(buf);
            }
            nodes[nodeIdx] = lexer.getInteger() - 1;
        }

        builder->addElement(type, tag, nodes.data(), NumNodes[type - 1]);
    }
    getNextToken();
    if (curTok != GMSHToken::end_elements) {
        return logErrorAnnotated<bool>("Expected $EndElements");
    }
    getNextToken();

    return true;
}

} // namespace tndm
