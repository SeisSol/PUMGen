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
#                 [REFERENCE <fixture> [TOLERANCE <tolerance>] [IDENTIFY] [HIGHORDER] [MIXED]]
#                 [EXPECT_OUTPUT <regex>] [EXPECT_ERROR <regex>]
#                 [OUTPUT <file>] [EXPECT_XDMF <regex>] [ARGS <argument>...])
function(pumgen_add_test name)
  cmake_parse_arguments(PARSE_ARGV 1 TEST "IDENTIFY;HIGHORDER;MIXED"
                        "RANKS;SOURCE;INPUT;OUTPUT;REFERENCE;TOLERANCE;EXPECT_OUTPUT;EXPECT_ERROR;EXPECT_XDMF" "ARGS")
  set(output ${PUMGEN_TEST_OUTPUT}/${name}.puml.h5)
  if(DEFINED TEST_OUTPUT)
    set(output ${TEST_OUTPUT})
  endif()
  if(DEFINED TEST_EXPECT_ERROR AND TEST_RANKS EQUAL 1)
    # started as MPI singleton, since mpiexec may drop the output of a process ending with
    # MPI_Abort
    set(command $<TARGET_FILE:pumgen>)
  else()
    set(command
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG}
        ${TEST_RANKS}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pumgen>
        ${MPIEXEC_POSTFLAGS})
  endif()
  set(input ${TEST_INPUT})
  if(NOT IS_ABSOLUTE ${input})
    set(input ${PUMGEN_TEST_FIXTURES}/${input})
  endif()
  list(APPEND command -s ${TEST_SOURCE} ${TEST_ARGS} ${input} ${output})
  string(REPLACE ";" "|" command "${command}")

  set(script_args -DCOMMAND=${command} -DOUTPUT=${output})
  if(DEFINED TEST_EXPECT_OUTPUT)
    list(APPEND script_args "-DEXPECT_OUTPUT=${TEST_EXPECT_OUTPUT}")
  endif()
  if(DEFINED TEST_EXPECT_XDMF)
    list(APPEND script_args "-DEXPECT_XDMF=${TEST_EXPECT_XDMF}")
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
    if(TEST_HIGHORDER)
      list(APPEND script_args -DHIGHORDER=ON)
    endif()
    if(TEST_MIXED)
      list(APPEND script_args -DMIXED=ON)
    endif()
  endif()

  add_test(NAME pumgen.${name} COMMAND ${CMAKE_COMMAND} ${script_args} -P
                                       ${CMAKE_CURRENT_SOURCE_DIR}/RunPumgenTest.cmake)
  set_tests_properties(
    pumgen.${name} PROPERTIES PROCESSORS ${TEST_RANKS} TIMEOUT 120 ENVIRONMENT
                              "ASAN_OPTIONS=detect_leaks=0")
endfunction()

# the GMSH parsers with a tiny read buffer
add_executable(pumgen_reader_test ReaderTest.cpp)
target_link_libraries(pumgen_reader_test PRIVATE pumgen-core)
add_test(NAME reader.buffer-boundaries COMMAND pumgen_reader_test ${PUMGEN_TEST_FIXTURES})
set_tests_properties(reader.buffer-boundaries PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")

# unit tests; the ones needing more ranks than started with are skipped
add_executable(pumgen_unit_tests unit/Main.cpp unit/DistributorTest.cpp unit/BoundaryTest.cpp
                                 unit/VertexFilterTest.cpp unit/GMSHTypesTest.cpp
                                 unit/CellTypeTest.cpp)
target_include_directories(pumgen_unit_tests SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/submodules/doctest)
target_link_libraries(pumgen_unit_tests PRIVATE pumgen-core)
foreach(ranks 1 3)
  add_test(NAME unit.np${ranks}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:pumgen_unit_tests> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit.np${ranks} PROPERTIES PROCESSORS ${ranks} TIMEOUT 120 ENVIRONMENT
                                                  "ASAN_OPTIONS=detect_leaks=0")
endforeach()

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
                REFERENCE coarse-o2.puml.h5 HIGHORDER)
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
pumgen_add_test(msh4.layered.chunksize.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --chunksize=4096)
# a chunk smaller than a row
pumgen_add_test(msh4.layered.chunksize1.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --chunksize=1)

# integer types reduced to the value range
pumgen_add_test(msh4.layered.compact.np3 RANKS 3 SOURCE msh4 INPUT layered-v41.msh
                REFERENCE layered.puml.h5 ARGS --compactify-datatypes)
pumgen_add_test(msh2.coarse-o2.compact.np2 RANKS 2 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 HIGHORDER ARGS --compactify-datatypes)

# malformed or unsupported input
pumgen_add_test(msh4.tiny-truncated RANKS 1 SOURCE msh4 INPUT tiny-truncated-v41.msh
                EXPECT_ERROR "GMSH parser error")
pumgen_add_test(msh4.tiny-badnumber RANKS 1 SOURCE msh4 INPUT tiny-badnumber-v41.msh
                EXPECT_ERROR "Expected number")
pumgen_add_test(msh4.tiny-unknowntype RANKS 1 SOURCE msh4 INPUT tiny-unknowntype-v41.msh
                EXPECT_ERROR "Unknown element type 999")
pumgen_add_test(msh4.tiny-v40 RANKS 1 SOURCE msh4 INPUT tiny-v40.msh
                EXPECT_ERROR "MSH 4.0 is not supported")
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

# higher-order meshes: the vertices of the cells are the vertices of the linear mesh, the other
# nodes the geometry of higher order; the inspheres are the ones of the linear mesh. The reference
# is written by fixtures/generate.py from the node tags in the file.
foreach(ranks 1 2)
  pumgen_add_test(msh4.coarse-o2.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT coarse-o2-v41.msh
                  REFERENCE coarse-o2.puml.h5 HIGHORDER
                  EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
endforeach()
pumgen_add_test(msh2.coarse-o2.np2 RANKS 2 SOURCE msh2 INPUT coarse-o2-v22.msh
                REFERENCE coarse-o2.puml.h5 HIGHORDER
                EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
pumgen_add_test(msh4.coarse.insphere RANKS 1 SOURCE msh4 INPUT coarse-v41.msh
                EXPECT_OUTPUT "Minimum insphere found: 0[.]03184[67]")
pumgen_add_test(msh4.coarse-partitioned RANKS 1 SOURCE msh4 INPUT coarse-partitioned-v41.msh
                EXPECT_ERROR "Partitioned MSH files are not supported")
pumgen_add_test(msh4.tiny-nophysical RANKS 1 SOURCE msh4 INPUT tiny-nophysical-v41.msh
                EXPECT_OUTPUT "1 volume[(]s[)] and 0 surface[(]s[)] without a physical group")

# binary MSH 4.1; the coordinates of the references come from the ASCII files
foreach(ranks 1 2 3)
  pumgen_add_test(msh4.layered-binary.np${ranks} RANKS ${ranks} SOURCE msh4
                  INPUT layered-binary-v41.msh REFERENCE layered.puml.h5 TOLERANCE 1e-12)
endforeach()
foreach(variant binary binary-bigendian binary-size4)
  pumgen_add_test(msh4.coarse-${variant}.np2 RANKS 2 SOURCE msh4 INPUT coarse-${variant}-v41.msh
                  REFERENCE coarse.puml.h5 TOLERANCE 1e-12)
endforeach()
pumgen_add_test(msh4.coarse-o2-binary.np2 RANKS 2 SOURCE msh4 INPUT coarse-o2-binary-v41.msh
                REFERENCE coarse-o2.puml.h5 TOLERANCE 1e-12 HIGHORDER
                EXPECT_OUTPUT "higher order, which rank 0 reads")
pumgen_add_test(msh4.periodic-binary.np2 RANKS 2 SOURCE msh4 INPUT periodic-binary-v41.msh
                REFERENCE periodic.puml.h5 TOLERANCE 1e-12 IDENTIFY)
pumgen_add_test(msh4.coarse-binary-truncated RANKS 1 SOURCE msh4
                INPUT coarse-binary-truncated-v41.msh EXPECT_ERROR "Unexpected end of file")

# binary MSH 4.1 files are read by all ranks unless --gmsh-reader=serial
foreach(ranks 4 7)
  pumgen_add_test(msh4.layered-binary.np${ranks} RANKS ${ranks} SOURCE msh4
                  INPUT layered-binary-v41.msh REFERENCE layered.puml.h5 TOLERANCE 1e-12)
endforeach()
pumgen_add_test(msh4.layered-binary-serial.np3 RANKS 3 SOURCE msh4 INPUT layered-binary-v41.msh
                REFERENCE layered.puml.h5 TOLERANCE 1e-12 ARGS --gmsh-reader=serial)
pumgen_add_test(msh4.periodic-binary-serial.np2 RANKS 2 SOURCE msh4 INPUT periodic-binary-v41.msh
                REFERENCE periodic.puml.h5 TOLERANCE 1e-12 IDENTIFY ARGS --gmsh-reader=serial)
pumgen_add_test(msh4.periodic-binary.np3 RANKS 3 SOURCE msh4 INPUT periodic-binary-v41.msh
                REFERENCE periodic.puml.h5 TOLERANCE 1e-12 IDENTIFY)
pumgen_add_test(msh4.tiny-binary.np3 RANKS 3 SOURCE msh4 INPUT tiny-binary-v41.msh
                REFERENCE tiny.puml.h5)
pumgen_add_test(msh4.tiny-binary-duplicate RANKS 1 SOURCE msh4 INPUT tiny-binary-duplicate-v41.msh
                EXPECT_ERROR "Duplicate node tag 2")
pumgen_add_test(msh4.tiny-binary-badnode RANKS 1 SOURCE msh4 INPUT tiny-binary-badnode-v41.msh
                EXPECT_ERROR "Unknown node tag 9")
pumgen_add_test(msh4.tiny-binary-sparse RANKS 1 SOURCE msh4 INPUT tiny-binary-sparse-v41.msh
                EXPECT_ERROR "Non-contiguous node tags are not supported")

# a Gambit mesh whose elements are in no group, derived from the layered one
file(READ ${PUMGEN_TEST_FIXTURES}/layered.neu gambit)
string(REPLACE "      1597         2         3" "      1597         0         3" gambit "${gambit}")
while(TRUE)
  string(FIND "${gambit}" "       ELEMENT GROUP" begin)
  if(begin EQUAL -1)
    break()
  endif()
  string(SUBSTRING "${gambit}" ${begin} -1 tail)
  string(FIND "${tail}" "ENDOFSECTION\n" length)
  math(EXPR end "${begin} + ${length} + 13")
  string(SUBSTRING "${gambit}" 0 ${begin} head)
  string(SUBSTRING "${gambit}" ${end} -1 rest)
  set(gambit "${head}${rest}")
endwhile()
file(WRITE ${PUMGEN_TEST_OUTPUT}/layered-nogroups.neu "${gambit}")
pumgen_add_test(gambit.nogroups.np3 RANKS 3 SOURCE gambit
                INPUT ${PUMGEN_TEST_OUTPUT}/layered-nogroups.neu
                EXPECT_OUTPUT "1597 elements belong to no group and get group 0")

# the XDMF file: the cells of a higher-order mesh with their vertices only, the precision of
# compacted integers, and the file name next to an output in a directory whose name ends in .h5
pumgen_add_test(
  msh4.coarse-o2.xdmf RANKS 1 SOURCE msh4 INPUT coarse-o2-v41.msh
  EXPECT_XDMF "Tetrahedron. NumberOfElements=.320.[^#]*Dimensions=.320 4.>[^#]*Name=.order.")
pumgen_add_test(
  msh4.coarse.compact-xdmf RANKS 2 SOURCE msh4 INPUT coarse-v41.msh ARGS --compactify-datatypes
  EXPECT_XDMF "Precision=.1. Format=.HDF. Dimensions=.320 4.>[^#]*Precision=.8. Format=.HDF. Dimensions=.109 3.")
file(MAKE_DIRECTORY ${PUMGEN_TEST_OUTPUT}/directory.h5)
pumgen_add_test(msh4.tiny.xdmf-path RANKS 1 SOURCE msh4 INPUT tiny-v41.msh
                OUTPUT ${PUMGEN_TEST_OUTPUT}/directory.h5/tiny.puml.h5 EXPECT_XDMF "tiny.puml.h5:/connect")

# a face of the surface mesh given twice
pumgen_add_test(msh2.tiny-duplicateface RANKS 1 SOURCE msh2 INPUT tiny-duplicateface-v22.msh
                EXPECT_ERROR "exists multiple times in the surface mesh")

# velocity-aware sizing, with a material of shear wave speed 1 in group 1 and 2 in group 2
if(EASI)
  configure_file(fixtures/velocity-settings.xml.in ${CMAKE_CURRENT_BINARY_DIR}/velocity-settings.xml
                 @ONLY)
  configure_file(fixtures/velocity-sizefield.xml.in
                 ${CMAKE_CURRENT_BINARY_DIR}/velocity-sizefield.xml @ONLY)

  # the counts are computed independently from the reference of the coarse mesh
  pumgen_add_test(
    msh4.coarse.velocity-check RANKS 2 SOURCE msh4 INPUT coarse-v41.msh
    ARGS --velocity-check=${CMAKE_CURRENT_BINARY_DIR}/velocity-settings.xml
    EXPECT_OUTPUT
      "320 of 320 cells lie in refinement regions.*mean edge:[^\n]*; 4 cells above 1.*longest edge:[^\n]*; 143 cells above 1"
  )

  # the expected grids are verified independently; they must not depend on the number of ranks
  foreach(ranks 1 3)
    set(prefix ${PUMGEN_TEST_OUTPUT}/sizefield-np${ranks})
    set(command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:pumgen-sizefield> ${MPIEXEC_POSTFLAGS} --spacing 0.25
                ${CMAKE_CURRENT_BINARY_DIR}/velocity-sizefield.xml ${prefix})
    string(REPLACE ";" "|" command "${command}")
    add_test(NAME sizefield.np${ranks}
             COMMAND ${CMAKE_COMMAND} -DCOMMAND=${command} -DPREFIX=${prefix}
                     -DEXPECTED=${PUMGEN_TEST_FIXTURES}/velocity-sizefield -DGRIDS=2 -P
                     ${CMAKE_CURRENT_SOURCE_DIR}/RunSizefieldTest.cmake)
    set_tests_properties(sizefield.np${ranks} PROPERTIES PROCESSORS ${ranks} TIMEOUT 120)
  endforeach()
endif()
if(NOT EASI)
  pumgen_add_test(msh4.tiny.velocity-check-without-easi RANKS 1 SOURCE msh4 INPUT tiny-v41.msh
                  ARGS --velocity-check=unused.xml EXPECT_ERROR "compiled without easi")
endif()

# a mesh of all four kinds of cells, in the layout of a VTKHDF unstructured grid; the reference is
# written by fixtures/generate.py from the cells it builds
foreach(ranks 1 2 3)
  pumgen_add_test(msh4.mixed.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT mixed-v41.msh
                  REFERENCE mixed.puml.h5 MIXED)
endforeach()
pumgen_add_test(msh4.mixed-binary.np3 RANKS 3 SOURCE msh4 INPUT mixed-binary-v41.msh
                REFERENCE mixed.puml.h5 MIXED EXPECT_OUTPUT "binary, read by all ranks")
pumgen_add_test(msh4.mixed-binary-serial.np2 RANKS 2 SOURCE msh4 INPUT mixed-binary-v41.msh
                REFERENCE mixed.puml.h5 MIXED ARGS --gmsh-reader=serial)
pumgen_add_test(msh2.mixed.np2 RANKS 2 SOURCE msh2 INPUT mixed-v22.msh
                REFERENCE mixed-v22.puml.h5 MIXED)
pumgen_add_test(msh4.mixed.compact.np2 RANKS 2 SOURCE msh4 INPUT mixed-v41.msh
                REFERENCE mixed.puml.h5 MIXED ARGS --compactify-datatypes)
pumgen_add_test(msh4.mixed.packed RANKS 1 SOURCE msh4 INPUT mixed-v41.msh ARGS --boundarytype=i32
                EXPECT_ERROR "packs four faces per cell")

# more ranks than vertices: the ranks without vertices take part in writing the identification
foreach(input periodic-tiny periodic-tiny-binary)
  pumgen_add_test(msh4.${input}.np9 RANKS 9 SOURCE msh4 INPUT ${input}-v41.msh
                  REFERENCE periodic-tiny.puml.h5 IDENTIFY)
endforeach()

# cells of several kinds and orders: the hexahedron and the wedges of order 2, the others linear
foreach(ranks 1 3)
  pumgen_add_test(msh4.mixed-ho.np${ranks} RANKS ${ranks} SOURCE msh4 INPUT mixed-ho-v41.msh
                  REFERENCE mixed-ho.puml.h5 MIXED HIGHORDER)
endforeach()
pumgen_add_test(msh4.mixed-ho-binary.np3 RANKS 3 SOURCE msh4 INPUT mixed-ho-binary-v41.msh
                REFERENCE mixed-ho.puml.h5 MIXED HIGHORDER
                EXPECT_OUTPUT "higher order, which rank 0 reads")
pumgen_add_test(msh4.mixed-ho-binary-serial.np2 RANKS 2 SOURCE msh4 INPUT mixed-ho-binary-v41.msh
                REFERENCE mixed-ho.puml.h5 MIXED HIGHORDER ARGS --gmsh-reader=serial)
