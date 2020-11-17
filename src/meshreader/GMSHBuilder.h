#ifndef GMSHBUILDER_20201014_H
#define GMSHBUILDER_20201014_H

#include "third_party/GMSHParser.h"

#include <array>
#include <cassert>
#include <vector>

namespace puml {

template <std::size_t D> struct GMSHSimplexType {};
template <> struct GMSHSimplexType<1u> { static constexpr long type = 1; }; // 1 -> line
template <> struct GMSHSimplexType<2u> { static constexpr long type = 2; }; // 2 -> triangle
template <> struct GMSHSimplexType<3u> { static constexpr long type = 4; }; // 4 -> tetrahedron

template <std::size_t D> class GMSHBuilder : public tndm::GMSHMeshBuilder {
public:
    using vertex_t = std::array<double, D>;
    using element_t = std::array<int, D + 1u>;
    using facet_t = std::array<int, D>;

    std::vector<vertex_t> vertices;
    std::vector<element_t> elements;
    std::vector<int> groups;
    std::vector<facet_t> facets;
    std::vector<int> bcs;

    void setNumVertices(std::size_t numVertices) { vertices.resize(numVertices); }
    void setVertex(long id, std::array<double, 3> const& x) {
        for (std::size_t i = 0; i < D; ++i) {
            vertices[id][i] = x[i];
        }
    }
    void setNumElements(std::size_t numElements) {
        elements.clear();
        elements.reserve(numElements);
        groups.clear();
        groups.reserve(numElements);
        facets.clear();
        facets.reserve(numElements);
        bcs.clear();
        bcs.reserve(numElements);
    }
    void addElement(long type, long tag, long* node, std::size_t numNodes) {
        if (type == GMSHSimplexType<D>::type) {
            assert(numNodes == D + 1u);
            element_t elem;
            std::copy(node, node + D + 1u, elem.begin());
            elements.emplace_back(elem);
            groups.emplace_back(static_cast<int>(tag));
        } else if (type == GMSHSimplexType<D - 1u>::type) {
            assert(numNodes == D);
            facet_t elem;
            std::copy(node, node + D, elem.begin());
            facets.emplace_back(elem);
            bcs.emplace_back(static_cast<int>(tag));
        }
    }
};

} // namespace puml

#endif // GMSHBUILDER_20201014_H
