# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

# Runs pumgen for one test case and checks the result.
#
# Variables:
#   COMMAND       pumgen command line, with '|' as the argument separator
#   OUTPUT        file written by pumgen
#   REFERENCE     (optional) reference PUML file the output has to match
#   TOLERANCE     (optional) absolute tolerance for the vertex coordinates
#   IDENTIFY      (optional) also compare the vertex identification of periodic meshes
#   HIGHORDER     (optional) also compare the geometry of higher order and the orders
#   MIXED         (optional) also compare the offsets and kinds of the cells, and the VTKHDF view
#   VTKHDF        (optional) also compare the VTKHDF view of a mesh of cells of one kind
#   H5DIFF        h5diff executable, required together with REFERENCE
#   EXPECT_ERROR  (optional) pumgen has to fail with a message matching this regular expression
#   EXPECT_OUTPUT (optional) the output of a successful run has to match this regular expression
#   EXPECT_XDMF   (optional) the XDMF file next to the output has to match this regular expression

string(REPLACE "|" ";" command "${COMMAND}")
file(REMOVE "${OUTPUT}")

execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_VARIABLE log
  ERROR_VARIABLE log
  TIMEOUT 100)
message("${log}")

if(DEFINED EXPECT_ERROR)
  if(result EQUAL 0)
    message(FATAL_ERROR "pumgen succeeded, but was expected to fail")
  endif()
  if(NOT log MATCHES "${EXPECT_ERROR}")
    message(FATAL_ERROR "pumgen failed without a message matching \"${EXPECT_ERROR}\"")
  endif()
  return()
endif()

if(NOT result EQUAL 0)
  message(FATAL_ERROR "pumgen failed: ${result}")
endif()

if(DEFINED EXPECT_OUTPUT AND NOT log MATCHES "${EXPECT_OUTPUT}")
  message(FATAL_ERROR "The output does not match \"${EXPECT_OUTPUT}\"")
endif()

if(DEFINED EXPECT_XDMF)
  if(OUTPUT MATCHES "[.]puml[.]h5$")
    string(REGEX REPLACE "[.]puml[.]h5$" ".xdmf" xdmf "${OUTPUT}")
  elseif(OUTPUT MATCHES "[.]h5$")
    string(REGEX REPLACE "[.]h5$" ".xdmf" xdmf "${OUTPUT}")
  else()
    set(xdmf "${OUTPUT}.xdmf")
  endif()
  if(NOT EXISTS "${xdmf}")
    message(FATAL_ERROR "${xdmf} was not written")
  endif()
  file(READ "${xdmf}" content)
  if(NOT content MATCHES "${EXPECT_XDMF}")
    message(FATAL_ERROR "${xdmf} does not match \"${EXPECT_XDMF}\":\n${content}")
  endif()
endif()

if(DEFINED REFERENCE)
  set(failures "")
  set(datasets connect geometry group boundary)
  if(IDENTIFY)
    list(APPEND datasets identify)
  endif()
  if(HIGHORDER)
    list(APPEND datasets geometry_ho geometry_ho_offsets order)
  endif()
  if(VTKHDF)
    # the view of a rectangular connect: a virtual connectivity and the offsets and kinds
    list(APPEND datasets VTKHDF/Points=geometry VTKHDF/Connectivity VTKHDF/Offsets VTKHDF/Types)
  endif()
  if(MIXED)
    # the VTKHDF view links the datasets under the names of VTK
    list(APPEND datasets connect_offsets cell_type VTKHDF/Points=geometry
         VTKHDF/Connectivity=connect VTKHDF/Offsets=connect_offsets VTKHDF/Types=cell_type)
  endif()
  foreach(entry ${datasets})
    # an entry is a dataset, or a dataset of the output and its dataset in the reference
    string(REPLACE "=" ";" paths ${entry})
    list(GET paths 0 dataset)
    list(GET paths -1 referenceDataset)
    set(options "")
    if(referenceDataset MATCHES "^geometry(_ho)?$" AND DEFINED TOLERANCE)
      set(options -d ${TOLERANCE})
    endif()
    execute_process(
      COMMAND ${H5DIFF} -n 10 ${options} ${OUTPUT} ${REFERENCE} /${dataset} /${referenceDataset}
      RESULT_VARIABLE differs
      OUTPUT_VARIABLE difflog
      ERROR_VARIABLE difflog)
    if(NOT differs EQUAL 0)
      string(APPEND failures "\n/${dataset} differs from the reference:\n${difflog}")
    endif()
  endforeach()
  if(NOT failures STREQUAL "")
    message(FATAL_ERROR "${failures}")
  endif()
endif()
