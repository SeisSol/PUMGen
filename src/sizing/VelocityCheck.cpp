// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#include "VelocityCheck.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include "utils/logger.h"

namespace {

// log-spaced histogram of the ratios between MinRatio and MaxRatio
constexpr double MinRatio = 1e-3;
constexpr double MaxRatio = 1e3;
constexpr std::size_t NumBins = 6000;
constexpr std::size_t Batch = std::size_t{1} << 20;

std::size_t binOf(double ratio) {
  const double position =
      std::log(ratio / MinRatio) / std::log(MaxRatio / MinRatio) * static_cast<double>(NumBins);
  return static_cast<std::size_t>(std::clamp(position, 0.0, static_cast<double>(NumBins - 1)));
}

double upperEdgeOf(std::size_t bin) {
  return MinRatio * std::pow(MaxRatio / MinRatio, static_cast<double>(bin + 1) / NumBins);
}

struct RatioStatistics {
  std::vector<std::uint64_t> histogram = std::vector<std::uint64_t>(NumBins, 0);
  std::uint64_t above = 0;
  double maximum = 0;

  void add(double ratio) {
    ++histogram[binOf(ratio)];
    above += ratio > 1 ? 1 : 0;
    maximum = std::max(maximum, ratio);
  }

  void reduce(MPI_Comm comm) {
    MPI_Allreduce(MPI_IN_PLACE, histogram.data(), NumBins, MPI_UINT64_T, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &above, 1, MPI_UINT64_T, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_DOUBLE, MPI_MAX, comm);
  }

  /**
   * The quantile, rounded up to the next bin edge (0.23 % apart).
   */
  [[nodiscard]] double quantile(double q, std::uint64_t total) const {
    const auto target = static_cast<std::uint64_t>(std::ceil(q * static_cast<double>(total)));
    std::uint64_t cumulative = 0;
    for (std::size_t bin = 0; bin < NumBins; ++bin) {
      cumulative += histogram[bin];
      if (cumulative >= target && cumulative > 0) {
        return std::min(upperEdgeOf(bin), maximum);
      }
    }
    return maximum;
  }

  [[nodiscard]] std::string summary(std::uint64_t total) const {
    std::stringstream stream;
    stream << std::setprecision(3) << "p50 " << quantile(0.5, total) << ", p90 "
           << quantile(0.9, total) << ", p99 " << quantile(0.99, total) << ", max " << maximum
           << "; " << above << " cells above 1";
    return stream.str();
  }
};

double distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}

} // namespace

void checkVelocityAwareMeshSize(const VelocityAwareMeshSize& meshSize, const CellVertices& cells,
                                const std::vector<int>& groups, MPI_Comm comm) {
  RatioStatistics meanEdge;
  RatioStatistics longestEdge;
  std::uint64_t checked = 0;

  std::vector<std::array<double, 3>> barycenters;
  std::vector<int> cellGroups;
  std::vector<double> frequencies;
  std::vector<std::array<double, 2>> edges;
  for (std::size_t first = 0; first < cells.numCells(); first += Batch) {
    const std::size_t count = std::min(Batch, cells.numCells() - first);
    barycenters.resize(count);
    cellGroups.resize(count);
    frequencies.resize(count);
    edges.resize(count);
    for (std::size_t c = 0; c < count; ++c) {
      const auto vertices = cells(first + c);
      std::array<double, 3> barycenter{};
      double sum = 0;
      double longest = 0;
      for (std::size_t a = 0; a < vertices.size(); ++a) {
        for (int d = 0; d < 3; ++d) {
          barycenter[d] += vertices[a][d] / static_cast<double>(vertices.size());
        }
        for (std::size_t b = a + 1; b < vertices.size(); ++b) {
          const double length = distance(vertices[a], vertices[b]);
          sum += length;
          longest = std::max(longest, length);
        }
      }
      barycenters[c] = barycenter;
      cellGroups[c] = groups[first + c];
      frequencies[c] = meshSize.targetedFrequencyAndGroup(barycenter).first;
      edges[c] = {sum / 6, longest};
    }

    const auto sizes = meshSize.meshSizes(barycenters, cellGroups, frequencies);
    for (std::size_t c = 0; c < count; ++c) {
      if (sizes[c] == std::numeric_limits<double>::max()) {
        continue;
      }
      ++checked;
      meanEdge.add(edges[c][0] / sizes[c]);
      longestEdge.add(edges[c][1] / sizes[c]);
    }
  }

  std::uint64_t total = cells.numCells();
  MPI_Allreduce(MPI_IN_PLACE, &checked, 1, MPI_UINT64_T, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_UINT64_T, MPI_SUM, comm);
  meanEdge.reduce(comm);
  longestEdge.reduce(comm);

  logInfo() << "Velocity check:" << checked << "of" << total
            << "cells lie in refinement regions (ratios of edge length and mesh size)";
  if (checked > 0) {
    logInfo() << "Velocity check, mean edge:" << meanEdge.summary(checked).c_str();
    logInfo() << "Velocity check, longest edge:" << longestEdge.summary(checked).c_str();
  }
}
