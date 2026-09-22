# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

# Runs pumgen-sizefield and compares the written grids with the expected ones.
#
# Variables:
#   COMMAND   pumgen-sizefield command line, with '|' as the argument separator
#   PREFIX    output prefix of the command
#   EXPECTED  prefix of the expected grids
#   GRIDS     number of grids

string(REPLACE "|" ";" command "${COMMAND}")
file(REMOVE "${PREFIX}.geo")
execute_process(
  COMMAND ${command}
  RESULT_VARIABLE result
  OUTPUT_VARIABLE log
  ERROR_VARIABLE log
  TIMEOUT 100)
message("${log}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "pumgen-sizefield failed: ${result}")
endif()

math(EXPR last "${GRIDS} - 1")
foreach(grid RANGE ${last})
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${PREFIX}-${grid}.bin"
                          "${EXPECTED}-${grid}.bin" RESULT_VARIABLE differs)
  if(NOT differs EQUAL 0)
    message(FATAL_ERROR "${PREFIX}-${grid}.bin differs from ${EXPECTED}-${grid}.bin")
  endif()
endforeach()
if(NOT EXISTS "${PREFIX}.geo")
  message(FATAL_ERROR "${PREFIX}.geo was not written")
endif()
