# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

find_program(H5DIFF_EXECUTABLE NAMES h5diff DOC "h5diff, used to compare the test output")
if(NOT H5DIFF_EXECUTABLE)
  message(WARNING "h5diff not found: the tests only check whether pumgen succeeds")
endif()

set(PUMGEN_TEST_FIXTURES ${CMAKE_CURRENT_SOURCE_DIR}/fixtures)
set(PUMGEN_TEST_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${PUMGEN_TEST_OUTPUT})

# pumgen_add_test(<name> RANKS <n> SOURCE <source> INPUT <fixture>
#                 [REFERENCE <fixture> [TOLERANCE <tolerance>]] [EXPECT_ERROR <regex>]
#                 [ARGS <argument>...])
function(pumgen_add_test name)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "RANKS;SOURCE;INPUT;REFERENCE;TOLERANCE;EXPECT_ERROR"
                        "ARGS")
  set(output ${PUMGEN_TEST_OUTPUT}/${name}.puml.h5)
  set(command
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG}
      ${TEST_RANKS}
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pumgen>
      ${MPIEXEC_POSTFLAGS}
      -s
      ${TEST_SOURCE}
      ${TEST_ARGS}
      ${PUMGEN_TEST_FIXTURES}/${TEST_INPUT}
      ${output})
  string(REPLACE ";" "|" command "${command}")

  set(script_args -DCOMMAND=${command} -DOUTPUT=${output})
  if(DEFINED TEST_EXPECT_ERROR)
    list(APPEND script_args "-DEXPECT_ERROR=${TEST_EXPECT_ERROR}")
  elseif(DEFINED TEST_REFERENCE AND H5DIFF_EXECUTABLE)
    list(APPEND script_args -DREFERENCE=${PUMGEN_TEST_FIXTURES}/${TEST_REFERENCE}
         -DH5DIFF=${H5DIFF_EXECUTABLE})
    if(DEFINED TEST_TOLERANCE)
      list(APPEND script_args -DTOLERANCE=${TEST_TOLERANCE})
    endif()
  endif()

  add_test(NAME pumgen.${name} COMMAND ${CMAKE_COMMAND} ${script_args} -P
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RunPumgenTest.cmake)
  set_tests_properties(
    pumgen.${name} PROPERTIES PROCESSORS ${TEST_RANKS} TIMEOUT 120 ENVIRONMENT
                              "ASAN_OPTIONS=detect_leaks=0")
endfunction()

# The fixtures are created by fixtures/generate.py; the reference files are computed from the
# gmsh model directly. The layered mesh has element and vertex counts which do not divide evenly
# among 3 or 4 ranks.

foreach(ranks 1 2 3 4)
  pumgen_add_test(msh4.layered.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT layered-v41.msh
                  REFERENCE layered.puml.h5)
endforeach()
foreach(ranks 1 2 3)
  pumgen_add_test(msh2.layered.np${ranks} RANKS ${ranks} SOURCE msh2 INPUT layered-v22.msh
                  REFERENCE layered.puml.h5)
endforeach()
foreach(ranks 1 2 3 4)
  pumgen_add_test(gambit.layered.np${ranks} RANKS ${ranks} SOURCE gambit INPUT layered.neu
                  REFERENCE layered.puml.h5 TOLERANCE 1e-9)
endforeach()
pumgen_add_test(msh4.coarse.np2 RANKS 2 SOURCE msh4 INPUT coarse-v41.msh REFERENCE coarse.puml.h5)
pumgen_add_test(msh2.coarse-o2.np1 RANKS 1 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 ARGS --order=2)
pumgen_add_test(msh4.tiny.np1 RANKS 1 SOURCE msh4 INPUT tiny-v41.msh REFERENCE tiny.puml.h5)
pumgen_add_test(msh2.tiny.np1 RANKS 1 SOURCE msh2 INPUT tiny-v22.msh REFERENCE tiny.puml.h5)

# more ranks than cells
foreach(ranks 2 3)
  pumgen_add_test(msh4.tiny.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT tiny-v41.msh
                  REFERENCE tiny.puml.h5)
endforeach()

# HDF5 filters
pumgen_add_test(msh4.layered.deflate.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --filter-enable=deflate1 --filter-chunksize=64)
pumgen_add_test(msh4.tiny.scaleoffset.np3 RANKS 3 SOURCE msh4 INPUT tiny-v41.msh
                REFERENCE tiny.puml.h5 ARGS --filter-enable=scaleoffset)

# command line options
pumgen_add_test(msh2.coarse-o2.shortopt.np1 RANKS 1 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 ARGS -o 2)
pumgen_add_test(msh4.layered.chunksize.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --chunksize=4096)
pumgen_add_test(msh4.tiny.order0 RANKS 1 SOURCE msh4 INPUT tiny-v41.msh ARGS --order=0
                EXPECT_ERROR "mesh order has to be at least 1")
