cmake_minimum_required(VERSION 3.9)

project(pumgen LANGUAGES C CXX)

option(SIMMETRIX "Use simmetrix libraries" OFF)
set(SIMMETRIX_ROOT "" CACHE STRING "Root directory of simmetrix dev files")
set(BEFORE_SIM_11 OFF CACHE BOOL "Use SimModSuite older than version 11.")
set(BEFORE_SIM_14 OFF CACHE BOOL "Use SimModSuite older than version 14.")
option(NETCDF "Use netcdf" OFF)
set(LOG_LEVEL "info" CACHE STRING "Log level for the code")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS "debug" "info" "warning" "error")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()

add_executable(pumgen
src/input/AnalysisAttributes.cpp
src/input/MeshAttributes.cpp
src/input/NetCDFPartition.cpp
src/input/ParallelVertexFilter.cpp
src/input/SimModSuite.cpp
src/input/split.cpp
src/meshreader/FidapReader.cpp
src/meshreader/GambitReader.cpp
src/pumgen.cpp
)

target_include_directories(pumgen PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

if (NETCDF)
  find_package(NetCDF REQUIRED)
  target_sources(pumgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/input/NetCDFPartition.cpp)
  target_include_directories(pumgen PUBLIC ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(pumgen PUBLIC ${NetCDF_LIBRARY})
  target_compile_definitions(pumgen PUBLIC USE_NETCDF)
endif()

set(HDF5_PREFER_PARALLEL True)
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(pumgen PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(pumgen PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
target_compile_definitions(pumgen PUBLIC USE_HDF)


find_package(OpenMP REQUIRED)
target_link_libraries(pumgen PUBLIC OpenMP::OpenMP_CXX)

find_package(MPI REQUIRED)
target_include_directories(pumgen SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(pumgen PUBLIC MPI::MPI_C)

#target_compile_definitions(pumgen PUBLIC USE_MPI PARALLEL)

find_package(Apf REQUIRED)
target_include_directories(pumgen PUBLIC ${APF_INCLUDE_DIR})
target_link_libraries(pumgen PUBLIC ${APF_LIBRARIES})

add_library(tinyxml2 submodules/tinyxml2/tinyxml2.cpp)
target_link_libraries(pumgen PUBLIC tinyxml2)

if (SIMMETRIX)
    find_package(Simmetrix REQUIRED)
    target_sources(pumgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/input/SimModSuite.cpp)
    target_compile_definitions(pumgen PUBLIC USE_SIMMOD)
    target_include_directories(pumgen PUBLIC ${SIMMETRIX_INCLUDE_DIR})
    target_link_libraries(pumgen PUBLIC ${SIMMETRIX_LIBRARIES})
endif()

if (BEFORE_SIM_11)
  target_compile_definitions(pumgen PUBLIC BEFORE_SIM_11)
  set(BEFORE_SIM_14 ON)
endif()
if (BEFORE_SIM_14)
  target_compile_definitions(pumgen PUBLIC BEFORE_SIM_14)
endif()

function(cast_log_level_to_int log_level_str log_level_int)
  if (${log_level_str} STREQUAL "debug")
    set(${log_level_int} 3 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "info")
    set(${log_level_int} 2 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "warning")
    set(${log_level_int} 1 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "error")
    set(${log_level_int} 0 PARENT_SCOPE)
  endif()
endfunction()

cast_log_level_to_int(LOG_LEVEL LOG_LEVEL)
target_compile_definitions(pumgen PUBLIC LOG_LEVEL=${LOG_LEVEL})
set_target_properties(pumgen PROPERTIES INSTALL_RPATH ${SIM_PS_KRNL_LIB_DIR})
install(TARGETS pumgen RUNTIME DESTINATION bin DESTINATION ${CMAKE_INSTALL_BINDIR})