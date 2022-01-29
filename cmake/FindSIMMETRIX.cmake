include(FindPackageHandleStandardArgs)

find_path(GMI_SIM_INCLUDE_DIR gmi_sim.h)
find_path(APF_SIM_INCLUDE_DIR apfSIM.h)
find_path(MESH_SIM_INCLUDE_DIR MeshSim.h)

list(APPEND SIMMETRIX_INCLUDE_DIR
  ${GMI_SIM_INCLUDE_DIR}
  ${APF_SIM_INCLUDE_DIR}
  ${MESH_SIM_INCLUDE_DIR}
)

find_library(GMI_SIM_LIB gmi_sim)
find_library(APF_SIM_LIB apf_sim)

set(SIM_LIB_HINT ${SIMMETRIX_ROOT}/lib/x64_rhel8_gcc83)

find_library(SIM_DISCRETE_LIB SimDiscrete ${SIM_LIB_HINT})
find_library(SIM_MESHING_LIB SimMeshing ${SIM_LIB_HINT})
find_library(SIM_MESH_TOOLS_LIB SimMeshTools ${SIM_LIB_HINT})
find_library(SIM_MODEL_LIB SimModel ${SIM_LIB_HINT})
if (PARASOLID)
find_library(SIM_PARASOLID_LIB NAMES SimParasolid310 SimParasolid320 SimParasolid330 PATHS ${SIM_LIB_HINT})
endif()
find_library(SIM_PARTITIONED_MESH_LIB SimPartitionedMesh ${SIM_LIB_HINT})
find_library(SIM_PARTITIONED_MESH_MPI_LIB SimPartitionedMesh-mpi ${SIM_LIB_HINT})
find_library(SIM_PARTITINED_WRAPPER_LIB SimPartitionWrapper-${SIM_MPI} ${SIM_LIB_HINT})
find_library(SIM_PS_KRNL_LIB pskernel ${SIM_LIB_HINT}/psKrnl)
  
get_filename_component(SIM_PS_KRNL_LIB_DIR ${SIM_PS_KRNL_LIB} DIRECTORY)

list(APPEND SIMMETRIX_LIBRARIES
  "${GMI_SIM_LIB}"
  "${APF_SIM_LIB}"
  "${SIM_DISCRETE_LIB}"
  "${SIM_EXPORT_LIB}"
  "${SIM_MESHING_LIB}"
  "${SIM_MESH_TOOLS_LIB}"
  "${SIM_PARTITIONED_MESH_MPI_LIB}"
  "${SIM_PARTITINED_WRAPPER_LIB}"
  "${SIM_MODEL_LIB}"
  "${SIM_PS_KRNL_LIB}"
)

if (PARASOLID)
list(APPEND SIMMETRIX_LIBRARIES
  "${SIM_PARASOLID_LIB}"
)
endif()

message(SIM_MODEL_LIB : ${SIM_MODEL_LIB})

string(REGEX REPLACE ".*/.*([0-9]+).[0-9]-.*" "\\1" SIM_MAJOR_VER ${SIM_MODEL_LIB})

message(SIM_MAJOR_VER : ${SIM_MAJOR_VER})

find_package_handle_standard_args(SIMMETRIX DEFAULT_MSG
  SIMMETRIX_INCLUDE_DIR SIMMETRIX_LIBRARIES SIM_MAJOR_VER)

mark_as_advanced(SIMMETRIX_INCLUDE_DIR SIMMETRIX_LIBRARIES SIM_MAJOR_VER)
