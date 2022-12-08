#ifndef GMSHLEXER_20200901_H
#define GMSHLEXER_20200901_H

#include <fstream>

namespace tndm {

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
private:
    uint64_t identifier;
    long integer;
    double real;
    char lastChar = ' ';
    std::istream* in = nullptr;
    GMSHSourceLocation loc = {1, 1};
protected:
    void advance();
public:
    static constexpr std::size_t MaxNumberLength = 128;

    void setIStream(std::istream* istream) {
        in = istream;
        loc = {1, 1};
    }

    GMSHToken getToken();
    uint64_t getIdentifier() const { return identifier; }
    long getInteger() const { return integer; }
    double getReal() const { return real; }
    GMSHSourceLocation getSourceLoc() const { return loc; }

};
} // namespace tndm

#endif // GMSHLEXER_20200901_H
