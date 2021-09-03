cmake_minimum_required(VERSION 3.9)

project(pumgen LANGUAGES C CXX)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

option(SIMMETRIX "Use simmetrix libraries" OFF)
set(SIMMETRIX_ROOT "" CACHE STRING "Root directory of simmetrix dev files")
set(SIM_MPI "mpich3" CACHE STRING "MPI version used by simmetrix")
option(NETCDF "Use netcdf" OFF)
option(PARASOLID "enable support for parasolid files" OFF)
set(LOG_LEVEL "info" CACHE STRING "Log level for the code")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS "debug" "info" "warning" "error")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()

add_executable(pumgen
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pumgen.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/meshreader/FidapReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/meshreader/GambitReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/meshreader/ParallelGMSHReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/GMSHLexer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/GMSHParser.cpp
  )

target_include_directories(pumgen PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

if (PARASOLID)
  target_compile_definitions(pumgen PUBLIC PARASOLID)
endif()

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
target_link_libraries(pumgen PUBLIC MPI::MPI_CXX)

target_compile_definitions(pumgen PUBLIC PARALLEL)

find_package(Apf REQUIRED)
target_include_directories(pumgen PUBLIC ${APF_INCLUDE_DIR})
target_link_libraries(pumgen PUBLIC ${APF_LIBRARIES})

add_library(tinyxml2 submodules/tinyxml2/tinyxml2.cpp)
target_link_libraries(pumgen PUBLIC tinyxml2)

if (SIMMETRIX)
    # Libs
    include(ExternalProject)

    find_package(YAML-CPP 0.6.3 QUIET)
    if (YAML-CPP_FOUND)
        target_link_libraries(SeisSol-lib PUBLIC ${YAML_CPP_LIBRARIES})
        target_include_directories(SeisSol-lib PUBLIC ${YAML_CPP_INCLUDE_DIR})
    else()
        ExternalProject_Add(yaml-cpp-build
                SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/submodules/yaml-cpp/
                INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install
                PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
                CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                -DYAML_CPP_BUILD_TOOLS=OFF
                -DYAML_CPP_BUILD_TESTS=OFF
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_CXX_FLAGS=${EXTRA_CXX_FLAGS}
                )
        add_dependencies(pumgen yaml-cpp-build)
        target_link_libraries(pumgen PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install/lib/libyaml-cpp.a")
        target_include_directories(pumgen PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp-install/include/)
    endif()

    ExternalProject_Add(ImpalaJIT-build
            SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/submodules/ImpalaJIT/
            INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_CXX_FLAGS=${EXTRA_CXX_FLAGS}
            )
    add_dependencies(pumgen ImpalaJIT-build)
    target_link_libraries(pumgen PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install/lib/libimpalajit.a")
    target_include_directories(pumgen PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ImpalaJIT-install/include/)

    # asagi
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(ASAGI REQUIRED asagi)
    target_compile_definitions(pumgen PUBLIC USE_ASAGI)
    target_link_libraries(pumgen PUBLIC ${ASAGI_STATIC_LDFLAGS})
    target_include_directories(pumgen PUBLIC ${ASAGI_INCLUDE_DIRS})
    target_compile_options(pumgen PUBLIC ${ASAGI_CFLAGS} ${ASAGI_CFLAGS_OTHER})

    # easi
    target_include_directories(pumgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/easi/include/)

find_package(Simmetrix REQUIRED)
    target_sources(pumgen PUBLIC 
      ${CMAKE_CURRENT_SOURCE_DIR}/src/input/AnalysisAttributes.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/input/EasiMeshSize.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/input/MeshAttributes.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/input/ParallelVertexFilter.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/src/input/split.cpp
    )
    target_compile_definitions(pumgen PUBLIC USE_SIMMOD)
    if(${SIM_MAJOR_VER} LESS 11)
      target_compile_definitions(pumgen PUBLIC BEFORE_SIM_11)
    elseif(${SIM_MAJOR_VER} LESS 15)
      target_compile_definitions(pumgen PUBLIC BEFORE_SIM_15)
    endif()
    target_include_directories(pumgen PUBLIC ${SIMMETRIX_INCLUDE_DIR})
    target_link_libraries(pumgen PUBLIC ${SIMMETRIX_LIBRARIES})

    target_link_libraries(pumgen PUBLIC tirpc) # Only for Fedora
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
if(SIMMETRIX)
  set_target_properties(pumgen PROPERTIES INSTALL_RPATH ${SIM_PS_KRNL_LIB_DIR})
endif()
install(TARGETS pumgen RUNTIME DESTINATION bin DESTINATION ${CMAKE_INSTALL_BINDIR})
