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
#                 [REFERENCE <fixture> [TOLERANCE <tolerance>] [IDENTIFY]]
#                 [EXPECT_OUTPUT <regex>] [EXPECT_ERROR <regex>] [ARGS <argument>...])
function(pumgen_add_test name)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "IDENTIFY"
                        "RANKS;SOURCE;INPUT;REFERENCE;TOLERANCE;EXPECT_OUTPUT;EXPECT_ERROR" "ARGS")
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
  if(DEFINED TEST_EXPECT_OUTPUT)
    list(APPEND script_args "-DEXPECT_OUTPUT=${TEST_EXPECT_OUTPUT}")
  endif()
  if(DEFINED TEST_EXPECT_ERROR)
    list(APPEND script_args "-DEXPECT_ERROR=${TEST_EXPECT_ERROR}")
  elseif(DEFINED TEST_REFERENCE AND H5DIFF_EXECUTABLE)
    list(APPEND script_args -DREFERENCE=${PUMGEN_TEST_FIXTURES}/${TEST_REFERENCE}
         -DH5DIFF=${H5DIFF_EXECUTABLE})
    if(DEFINED TEST_TOLERANCE)
      list(APPEND script_args -DTOLERANCE=${TEST_TOLERANCE})
    endif()
    if(TEST_IDENTIFY)
      list(APPEND script_args -DIDENTIFY=ON)
    endif()
  endif()

  add_test(NAME pumgen.${name} COMMAND ${CMAKE_COMMAND} ${script_args} -P
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RunPumgenTest.cmake)
  set_tests_properties(
    pumgen.${name} PROPERTIES PROCESSORS ${TEST_RANKS} TIMEOUT 120 ENVIRONMENT
                              "ASAN_OPTIONS=detect_leaks=0")
endfunction()

# the GMSH parsers with a tiny read buffer
add_executable(
  pumgen_reader_test
  ReaderTest.cpp ${PROJECT_SOURCE_DIR}/src/meshreader/MshInput.cpp
  ${PROJECT_SOURCE_DIR}/src/meshreader/GMSHParser.cpp
  ${PROJECT_SOURCE_DIR}/src/meshreader/GMSH2Parser.cpp
  ${PROJECT_SOURCE_DIR}/src/meshreader/GMSH4Parser.cpp)
target_include_directories(pumgen_reader_test PRIVATE ${PROJECT_SOURCE_DIR}/src
                                                      ${PROJECT_SOURCE_DIR}/submodules)
target_compile_definitions(pumgen_reader_test PRIVATE LOG_LEVEL=${LOG_LEVEL})
add_test(NAME reader.buffer-boundaries COMMAND pumgen_reader_test ${PUMGEN_TEST_FIXTURES})
set_tests_properties(reader.buffer-boundaries PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")

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

# integer types reduced to the value range
pumgen_add_test(msh4.layered.compact.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --compactify-datatypes)
pumgen_add_test(msh2.coarse-o2.compact.np2 RANKS 2 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 ARGS --compactify-datatypes --order=2)

# malformed or unsupported input
pumgen_add_test(msh4.tiny-truncated RANKS 1 SOURCE msh4 INPUT tiny-truncated-v41.msh
                EXPECT_ERROR "GMSH parser error")
pumgen_add_test(msh4.tiny-badnumber RANKS 1 SOURCE msh4 INPUT tiny-badnumber-v41.msh
                EXPECT_ERROR "Expected number")
pumgen_add_test(msh4.tiny-unknowntype RANKS 1 SOURCE msh4 INPUT tiny-unknowntype-v41.msh
                EXPECT_ERROR "Unknown element type 999")
pumgen_add_test(msh4.tiny-v40 RANKS 1 SOURCE msh4 INPUT tiny-v40.msh
                EXPECT_ERROR "MSH 4.0 is not supported")
pumgen_add_test(msh4.coarse-binary RANKS 1 SOURCE msh4 INPUT coarse-binary-v41.msh
                EXPECT_ERROR "Binary MSH files are not supported")
pumgen_add_test(msh2.tiny-as-msh4 RANKS 1 SOURCE msh2 INPUT tiny-v41.msh
                EXPECT_ERROR "Unsupported MSH version 4.1")

# periodic meshes
foreach(ranks 1 2 3)
  pumgen_add_test(msh4.periodic.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT periodic-v41.msh
                  REFERENCE periodic.puml.h5 IDENTIFY)
endforeach()

# node tags
pumgen_add_test(msh4.coarse-offset.np2 RANKS 2 SOURCE msh4 INPUT coarse-offset-v41.msh
                REFERENCE coarse.puml.h5)
pumgen_add_test(msh4.periodic-parametric.np2 RANKS 2 SOURCE msh4 INPUT periodic-parametric-v41.msh
                REFERENCE periodic.puml.h5 IDENTIFY)
pumgen_add_test(msh4.tiny-sparse RANKS 1 SOURCE msh4 INPUT tiny-sparse-v41.msh
                EXPECT_ERROR "Non-contiguous node tags are not supported")
pumgen_add_test(msh4.tiny-duplicate RANKS 1 SOURCE msh4 INPUT tiny-duplicate-v41.msh
                EXPECT_ERROR "Duplicate node tag 2")
pumgen_add_test(msh4.tiny-badnode RANKS 1 SOURCE msh4 INPUT tiny-badnode-v41.msh
                EXPECT_ERROR "Unknown node tag 9")
pumgen_add_test(msh2.tiny-duplicate RANKS 1 SOURCE msh2 INPUT tiny-duplicate-v22.msh
                EXPECT_ERROR "Duplicate node tag 2")
pumgen_add_test(msh2.tiny-badnode RANKS 1 SOURCE msh2 INPUT tiny-badnode-v22.msh
                EXPECT_ERROR "Unknown node tag")

# higher-order meshes (with the inspheres of the linear mesh) and partitioned files
foreach(ranks 1 2)
  pumgen_add_test(msh4.coarse-o2.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT coarse-o2-v41.msh
                  REFERENCE coarse-o2.puml.h5 ARGS --order=2
                  EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
endforeach()
pumgen_add_test(msh2.coarse-o2.np2 RANKS 2 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 ARGS --order=2
                EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
pumgen_add_test(msh4.coarse.insphere RANKS 1 SOURCE msh4 INPUT coarse-v41.msh
                EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
pumgen_add_test(msh4.coarse-partitioned RANKS 1 SOURCE msh4 INPUT coarse-partitioned-v41.msh
                EXPECT_ERROR "Partitioned MSH files are not supported")
pumgen_add_test(msh4.tiny-nophysical RANKS 1 SOURCE msh4 INPUT tiny-nophysical-v41.msh
                EXPECT_OUTPUT "1 volume[(]s[)] and 0 surface[(]s[)] without a physical group")
