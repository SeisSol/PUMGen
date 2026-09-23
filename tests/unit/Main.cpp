// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/extensions/doctest_mpi.h"

#include <mpi.h>

int main(int argc, char** argv) {
  doctest::mpi_init_thread(argc, argv, MPI_THREAD_MULTIPLE);

  doctest::Context context;
  context.setOption("reporters", "MpiConsoleReporter");
  context.applyCommandLine(argc, argv);
  const int result = context.run();

  doctest::mpi_finalize();
  return result;
}
