// SPDX-FileCopyrightText: 2020 SeisSol Group
// SPDX-FileCopyrightText: 2020 Ludwig-Maximilians-Universität München
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_
#define PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <set>
#include <vector>

#include "utils/logger.h"

#include "third_party/GMSHMeshBuilder.h"

namespace puml {

constexpr std::size_t nodeCount(std::size_t D, std::size_t Order) {
  std::size_t count = 1;
  for (std::size_t d = 1; d <= D; ++d) {
    count *= (Order + D - d + 1);
    count /= d;
  }
  return count;
}

template <std::size_t P, typename BaseT, template <std::size_t> typename ParamT, typename... Args>
BaseT* makePointerI(std::size_t param, Args&&... args) {
  if constexpr (P < 11) {
    if (param == P) {
      return new ParamT<P>(std::forward<Args>(args)...);
    } else {
      return makePointerI<P + 1, BaseT, ParamT>(param, std::forward<Args>(args)...);
    }
  } else {
    logError() << "Order too high:" << param;
    return nullptr;
  }
}

template <typename BaseT, template <std::size_t> typename ParamT, typename... Args>
BaseT* makePointer(std::size_t param, Args&&... args) {
  return makePointerI<1, BaseT, ParamT>(param, std::forward<Args>(args)...);
}

template <std::size_t D, std::size_t Order> struct GMSHSimplexType {};

// clang format seems to mess up the following lines from time to time, hence we disable it
template <std::size_t Order> struct GMSHSimplexType<0U, Order> {
  static constexpr long type = 15;
};
template <> struct GMSHSimplexType<1U, 1U> {
  static constexpr long type = 1;
};
template <> struct GMSHSimplexType<2U, 1U> {
  static constexpr long type = 2;
};
template <> struct GMSHSimplexType<3U, 1U> {
  static constexpr long type = 4;
};
template <> struct GMSHSimplexType<1U, 2U> {
  static constexpr long type = 8;
};
template <> struct GMSHSimplexType<2U, 2U> {
  static constexpr long type = 9;
};
template <> struct GMSHSimplexType<3U, 2U> {
  static constexpr long type = 11;
};
template <> struct GMSHSimplexType<1U, 3U> {
  static constexpr long type = 26;
};
template <> struct GMSHSimplexType<1U, 4U> {
  static constexpr long type = 27;
};
template <> struct GMSHSimplexType<1U, 5U> {
  static constexpr long type = 28;
};
template <> struct GMSHSimplexType<1U, 6U> {
  static constexpr long type = 62;
};
template <> struct GMSHSimplexType<1U, 7U> {
  static constexpr long type = 63;
};
template <> struct GMSHSimplexType<1U, 8U> {
  static constexpr long type = 64;
};
template <> struct GMSHSimplexType<1U, 9U> {
  static constexpr long type = 65;
};
template <> struct GMSHSimplexType<1U, 10U> {
  static constexpr long type = 66;
};
template <> struct GMSHSimplexType<2U, 3U> {
  static constexpr long type = 21;
};
template <> struct GMSHSimplexType<2U, 4U> {
  static constexpr long type = 23;
};
template <> struct GMSHSimplexType<2U, 5U> {
  static constexpr long type = 25;
};
template <> struct GMSHSimplexType<2U, 6U> {
  static constexpr long type = 42;
};
template <> struct GMSHSimplexType<2U, 7U> {
  static constexpr long type = 43;
};
template <> struct GMSHSimplexType<2U, 8U> {
  static constexpr long type = 44;
};
template <> struct GMSHSimplexType<2U, 9U> {
  static constexpr long type = 45;
};
template <> struct GMSHSimplexType<2U, 10U> {
  static constexpr long type = 46;
};
template <> struct GMSHSimplexType<3U, 3U> {
  static constexpr long type = 29;
};
template <> struct GMSHSimplexType<3U, 4U> {
  static constexpr long type = 30;
};
template <> struct GMSHSimplexType<3U, 5U> {
  static constexpr long type = 31;
};
template <> struct GMSHSimplexType<3U, 6U> {
  static constexpr long type = 71;
};
template <> struct GMSHSimplexType<3U, 7U> {
  static constexpr long type = 72;
};
template <> struct GMSHSimplexType<3U, 8U> {
  static constexpr long type = 73;
};
template <> struct GMSHSimplexType<3U, 9U> {
  static constexpr long type = 74;
};
template <> struct GMSHSimplexType<3U, 10U> {
  static constexpr long type = 75;
};

template <std::size_t D, std::size_t Order> class GMSHBuilder : public tndm::GMSHMeshBuilder {
  public:
  using vertex_t = std::array<double, D>;
  using element_t = std::array<std::size_t, nodeCount(D, Order)>;
  using facet_t = std::array<int, nodeCount(D - 1u, Order)>;

  std::vector<vertex_t> vertices;
  std::vector<std::size_t> identify;
  std::vector<element_t> elements;
  std::vector<int> groups;
  std::vector<facet_t> facets;
  std::vector<int> bcs;

  void resizeIdentifyIfNeeded(std::size_t newSize) {
    const auto oldSize = identify.size();
    if (newSize > oldSize) {
      identify.resize(newSize);
      std::iota(identify.begin() + oldSize, identify.end(), oldSize);
    }
  }

  void setNumVertices(std::size_t numVertices) { vertices.resize(numVertices); }

  void setVertex(long id, std::array<double, D> const& x) {
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
    if (type == GMSHSimplexType<D, Order>::type) {
      element_t elem;
      assert(numNodes == elem.size());
      std::copy_n(node, elem.size(), elem.begin());
      elements.emplace_back(elem);
      groups.emplace_back(static_cast<int>(tag));
    } else if (type == GMSHSimplexType<D - 1u, Order>::type) {
      facet_t elem;
      assert(numNodes == elem.size());
      std::copy_n(node, elem.size(), elem.begin());
      facets.emplace_back(elem);
      bcs.emplace_back(static_cast<int>(tag));
    }
  }
  void addVertexLink(std::size_t vertex, std::size_t linkVertex) {
    // needed in case we parse periodic before we parse the vertex/node section
    resizeIdentifyIfNeeded(vertex + 1);
    identify[vertex] = linkVertex;
  }

  void postprocess() {
    // GMSH does not pick the same node to identify for all
    // thus, stratify (i.e. "path compress" in "union find" terms)

    if (!identify.empty()) {
      // only resize if we have an identification array
      resizeIdentifyIfNeeded(vertices.size());
      for (std::size_t i = 0; i < identify.size(); ++i) {
        if (identify[i] != i) {

          // find until we loop
          std::size_t curr = i;
          std::set<std::size_t> equivalent;
          equivalent.insert(curr);
          while (equivalent.find(identify[curr]) == equivalent.end()) {
            curr = identify[curr];
            equivalent.insert(curr);
          }

          // then set everything to one value
          const auto identifyEquiv = *equivalent.begin();
          for (const auto& index : equivalent) {
            identify[index] = identifyEquiv;
          }
        }
      }
    }
  }
};

} // namespace puml

#endif // PUMGEN_SRC_MESHREADER_GMSHBUILDER_H_
