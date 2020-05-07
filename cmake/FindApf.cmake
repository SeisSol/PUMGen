include(FindPackageHandleStandardArgs)

# simmetrixLibs = [('gmi_sim', 'gmi_sim.h'), ('apf_sim', 'apfSIM.h')]
# # Libs required to couple Zoltan
# zoltanLibs = [('apf_zoltan', 'apfZoltan.h')]
# # Other APF libs
# libs = [('gmi', 'gmi.h'), ('mds', 'apfMDS.h'), ('ma', 'ma.h'),
#     ('apf', 'apf.h'), ('pcu', 'PCU.h'),
#     ('lion','lionCompress.h'), ('mth','mth.h')]


find_path(GMI_INCLUDE_DIR gmi.h)
find_path(MDS_INCLUDE_DIR apfMDS.h)
find_path(MA_INCLUDE_DIR ma.h)
find_path(APF_INCLUDE_DIR apf.h)
find_path(PCU_INCLUDE_DIR PCU.h)
find_path(LION_INCLUDE_DIR lionCompress.h)
find_path(MTH_INCLUDE_DIR mth.h)
find_path(APF_ZOLTAN_INCLUDE_DIR apfZoltan.h)

find_path(GMI_SIM_INCLUDE_DIR gmi_sim.h)
find_path(APF_SIM_INCLUDE_DIR apf_sim.h)

set(APF_INCLUDE_DIR
  ${GMI_INCLUDE_DIR}
  ${MDS_INCLUDE_DIR}
  ${MA_INCLUDE_DIR}
  ${APF_INCLUDE_DIR}
  ${PCU_INCLUDE_DIR}
  ${LION_INCLUDE_DIR}
  ${MTH_INCLUDE_DIR}
  ${APF_ZOLTAN_INCLUDE_DIR}
  ${GMI_SIM_INCLUDE_DIR}
  ${APF_SIM_INCLUDE_DIR}
)

find_library(GMI_LIB gmi)
find_library(MDS_LIB mds)
find_library(MA_LIB ma)
find_library(APF_LIB apf)
find_library(PCU_LIB pcu)
find_library(LION_LIB lion)
find_library(MTH_LIB mth)
find_library(APF_ZOLTAN_LIB apf_zoltan)

find_library(GMI_SIM_LIB gmi_sim)
find_library(APF_SIM_LIB apf_sim)

set(APF_LIBRARIES
  ${GMI_LIB}
  ${MDS_LIB}
  ${MA_LIB}
  ${APF_LIB}
  ${PCU_LIB}
  ${LION_LIB}
  ${MTH_LIB}
  ${APF_ZOLTAN_LIB}
  ${GMI_SIM_LIB}
  ${APF_SIM_LIB}
)

find_package_handle_standard_args(APF DEFAULT_MSG
        APF_INCLUDE_DIR APF_LIBRARIES)

mark_as_advanced(APF_INCLUDE_DIR APF_LIBRARIES)
