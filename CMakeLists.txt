cmake_minimum_required(VERSION 3.9)

project(pumgen LANGUAGES C CXX)

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
  src/input/NetCDFPartition.cpp
  src/input/ParallelVertexFilter.cpp
  #src/input/SimModSuite.cpp
  src/pumgen.cpp
  src/meshreader/GambitReader.cpp
  src/meshreader/FidapReader.cpp
)

target_include_directories(pumgen PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src" # SeisSol :(
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

find_package(NetCDF REQUIRED)
target_include_directories(pumgen PUBLIC ${NetCDF_INCLUDE_DIRS})
target_link_libraries(pumgen PUBLIC ${NetCDF_LIBRARY})
target_compile_definitions(pumgen PUBLIC USE_NETCDF)

#if (MPI)
    set(HDF5_PREFER_PARALLEL True)
#endif()
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(pumgen PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(pumgen PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
target_compile_definitions(pumgen PUBLIC USE_HDF)


find_package(OpenMP REQUIRED)
target_link_libraries(pumgen PUBLIC OpenMP::OpenMP_CXX)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#target_compile_definitions(pumgen PUBLIC OMP OMPI_SKIP_MPICXX)

find_package(MPI REQUIRED)

target_include_directories(pumgen SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(pumgen PUBLIC MPI::MPI_C)

#target_compile_definitions(pumgen PUBLIC USE_MPI PARALLEL)

find_package(Apf REQUIRED)
message(${APF_INCLUDE_DIR})
message(${APF_LIBRARIES})
target_include_directories(pumgen PUBLIC ${APF_INCLUDE_DIR})
target_link_libraries(pumgen PUBLIC ${APF_LIBRARIES})
