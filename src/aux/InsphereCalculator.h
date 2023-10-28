#ifndef PUMGEN_AUX_INSPHERE_CALCULATOR_H_
#define PUMGEN_AUX_INSPHERE_CALCULATOR_H_

#include <mpi.h>
#include <vector>

// assumes a contiguous distribution of all vertices over all processes

std::vector<double> calculateInsphere(const std::vector<std::size_t>& connectivity,
                                      const std::vector<double>& geometry, MPI_Comm comm);

#endif // PUMGEN_AUX_INSPHERE_CALCULATOR_H_
