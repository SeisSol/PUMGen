// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
#ifndef PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
#define PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_

#include <cstddef>
#include <mpi.h>
#include <vector>

// assumes a contiguous distribution of all vertices over all processes

std::vector<double> calculateInsphere(const std::vector<std::size_t>& connectivity,
                                      const std::vector<double>& geometry, MPI_Comm comm);

#endif // PUMGEN_SRC_HELPER_INSPHERECALCULATOR_H_
