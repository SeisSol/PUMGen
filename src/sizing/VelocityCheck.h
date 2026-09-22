// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_SIZING_VELOCITYCHECK_H_
#define PUMGEN_SRC_SIZING_VELOCITYCHECK_H_

#include <mpi.h>
#include <vector>

#include "VelocityAwareMeshSize.h"
#include "helper/InsphereCalculator.h"

/**
 * Compares the edge lengths of the cells in the refinement regions with the velocity-aware mesh
 * size at their barycentre (for the material of the cell's group) and logs the distribution of the
 * ratios, as mean and as longest edge length over mesh size.
 */
void checkVelocityAwareMeshSize(const VelocityAwareMeshSize& meshSize, const CellVertices& cells,
                                const std::vector<int>& groups, MPI_Comm comm);

#endif // PUMGEN_SRC_SIZING_VELOCITYCHECK_H_
