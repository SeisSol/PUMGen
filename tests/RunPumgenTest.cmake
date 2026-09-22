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
#   H5DIFF        h5diff executable, required together with REFERENCE
#   EXPECT_ERROR  (optional) pumgen has to fail with a message matching this regular expression

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

if(DEFINED REFERENCE)
  set(failures "")
  foreach(dataset connect geometry group boundary)
    set(options "")
    if(dataset STREQUAL "geometry" AND DEFINED TOLERANCE)
      set(options -d ${TOLERANCE})
    endif()
    execute_process(
      COMMAND ${H5DIFF} -n 10 ${options} ${OUTPUT} ${REFERENCE} /${dataset} /${dataset}
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
