# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Option not to run with large tests
# ----------------------------------
option(SKIP_LARGE_TESTS "Skip tests with more than 6 processors" ON)

# Create the Data/ directory
# --------------------------
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# FV3 config files
# ----------------
list( APPEND fv3jedi_test_fv3files
Data/fv3files/akbk64.nc4
Data/fv3files/akbk72.nc4
Data/fv3files/akbk127.nc4
Data/fv3files/field_table
Data/fv3files/field_table_cold_starts
Data/fv3files/field_table_lam_cmaq
Data/fv3files/fmsmpp.nml
Data/fv3files/input_gfs_c12.nml
Data/fv3files/input_gfs_c12_p12.nml
Data/fv3files/input_gfs_c24.nml
Data/fv3files/input_gfs_c24_p12.nml
Data/fv3files/input_geos_c12.nml
Data/fv3files/input_geos_c12_p12.nml
Data/fv3files/input_geos_c24.nml
Data/fv3files/input_lam_cmaq.nml
Data/fv3files/inputpert_4dvar.nml
Data/fv3files/inputpert_fsoi.nml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fv3files)
LINK_FILES( "${fv3jedi_test_fv3files}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Test input scripts
# ------------------
list( APPEND fv3jedi_test_input
  testinput/3denvar_geos_aero.yaml
  testinput/3dvar_gfs_bal.yaml
  testinput/3dvar_lam_cmaq.yaml
  testinput/4denvar.yaml
  testinput/addincrement_geos.yaml
  testinput/addincrement_gfs.yaml
  testinput/bumpparameters_nicas_geos_p12.yaml
  testinput/bumpparameters_nicas_geos.yaml
  testinput/bumpparameters_nicas_gfs_p12.yaml
  testinput/bumpparameters_nicas_gfs.yaml
  testinput/bumpparameters_nicas_lam_cmaq.yaml
  testinput/convertstate_geos_crtm_surface.yaml
  testinput/convertstate_geos.yaml
  testinput/convertstate_gfs_cold2fv3.yaml
  testinput/convertstate_gfs.yaml
  testinput/convertstate_gfs_c2ll.yaml
  testinput/diffstates_geos.yaml
  testinput/diffstates_gfs.yaml
  testinput/diffstates_lam_cmaq.yaml
  testinput/dirac_geos.yaml
  testinput/dirac_gfs_bal.yaml
  testinput/dirac_gfs_c2a.yaml
  testinput/dirac_gfs.yaml
  testinput/eda_3dvar_1.yaml
  testinput/eda_3dvar_2.yaml
  testinput/eda_3dvar.yaml
  testinput/eda_4dvar_fv3lm_1.yaml
  testinput/eda_4dvar_fv3lm_2.yaml
  testinput/eda_4dvar_fv3lm.yaml
  testinput/enshofx_fv3lm_1.yaml
  testinput/enshofx_fv3lm_2.yaml
  testinput/enshofx_fv3lm.yaml
  testinput/ensvariance_geos.yaml
  testinput/ensvariance_gfs.yaml
  testinput/errorcovariance_bump.yaml
  testinput/errorcovariance.yaml
  testinput/forecast_fv3lm.yaml
  testinput/forecast_geos.yaml
  testinput/forecast_ufs.yaml
  testinput/geometry_2d.yaml
  testinput/geometry_geos.yaml
  testinput/geometry_gfs127.yaml
  testinput/geometry_gfs.yaml
  testinput/geometry_iterator_geos.yaml
  testinput/geometry_iterator_gfs.yaml
  testinput/geometry_iterator_lam_cmaq.yaml
  testinput/geometry_lam_cmaq.yaml
  testinput/geometry_nested.yaml
  testinput/getvalues_geos.yaml
  testinput/getvalues_gfs.yaml
  testinput/hofx_fv3-geos.yaml
  testinput/hofx_fv3-gfs.yaml
  testinput/hofx_fv3lm.yaml
  testinput/hofx_gfs_aero.yaml
  testinput/hofx_gnssro_ropp.yaml
  testinput/hofx_nomodel.yaml
  testinput/hofx_nomodel_lam_cmaq.yaml
  testinput/hofx_save_geovals.yaml
  testinput/hyb-3dvar_geos.yaml
  testinput/hyb-3dvar_gfs_0obs.yaml
  testinput/hyb-3dvar_gfs_1obs.yaml
  testinput/hyb-3dvar_gfs_aero.yaml
  testinput/hyb-3dvar_gfs_gnssro_ropp.yaml
  testinput/hyb-3dvar_gfs.yaml
  testinput/hyb-3dvar_p12.yaml
  testinput/hyb-3dvar.yaml
  testinput/hyb-4dvar_pseudo-geos.yaml
  testinput/hyb-fgat_fv3lm.yaml
  testinput/hyb-fgat_gfs_aero.yaml
  testinput/increment_geos.yaml
  testinput/increment_gfs.yaml
  testinput/letkf.yaml
  testinput/lgetkf.yaml
  testinput/lineargetvalues_geos.yaml
  testinput/lineargetvalues_gfs.yaml
  testinput/linearmodel_physics.yaml
  testinput/linearmodel.yaml
  testinput/linearvariablechange_geos.yaml
  testinput/linearvariablechange_gfs.yaml
  testinput/localization_bump.yaml
  testinput/model_fv3lm.yaml
  testinput/model_geos.yaml
  testinput/model_pseudo-geos.yaml
  testinput/model_pseudo-gfs.yaml
  testinput/state_geos_aerosol_ext.yaml
  testinput/state_geos_aerosol.yaml
  testinput/state_geos.yaml
  testinput/state_gfs.yaml
  testinput/state_lam_cmaq.yaml
  testinput/state_lam_atmos.yaml
  testinput/test_time.yaml
  testinput/variablechange_geos.yaml
  testinput/variablechange_gfs.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
LINK_FILES( "${fv3jedi_test_input}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Test output scripts
# -------------------
list( APPEND fv3jedi_testoutput
  testoutput/3dvar_lam_cmaq.ref
  testoutput/3denvar_geos_aero.ref
  testoutput/4denvar.ref
  testoutput/addincrement_geos.ref
  testoutput/addincrement_gfs.ref
  testoutput/bumpparameters_nicas_geos.ref
  testoutput/bumpparameters_nicas_gfs.ref
  testoutput/bumpparameters_nicas_lam_cmaq.ref
  testoutput/convertstate_geos_crtm_surface.ref
  testoutput/convertstate_geos.ref
  testoutput/convertstate_gfs_cold2fv3.ref
  testoutput/convertstate_gfs.ref
  testoutput/convertstate_gfs_c2ll.ref
  testoutput/CTestCostData.txt.awsgnu.test
  testoutput/CTestCostData.txt.awsintel.test
  testoutput/diffstates_geos.ref
  testoutput/diffstates_gfs.ref
  testoutput/diffstates_lam_cmaq.ref
  testoutput/dirac_geos.ref
  testoutput/dirac_gfs_bal.ref
  testoutput/dirac_gfs_c2a.ref
  testoutput/dirac_gfs.ref
  testoutput/ensvariance_geos.ref
  testoutput/ensvariance_gfs.ref
  testoutput/forecast_fv3-geos.ref
  testoutput/forecast_fv3-gfs.ref
  testoutput/forecast_fv3lm.ref
  testoutput/hofx_fv3-geos.ref
  testoutput/hofx_fv3-gfs_gnssro_ropp.ref
  testoutput/hofx_fv3-gfs.ref
  testoutput/hofx_fv3lm.ref
  testoutput/hofx_gfs_aero.ref
  testoutput/hofx_gnssro_ropp.ref
  testoutput/hofx_nomodel_lam_cmaq.ref
  testoutput/hofx_nomodel_gfs.ref
  testoutput/hofx_nomodel.ref
  testoutput/hofx_save_geovals.ref
  testoutput/hyb-3dvar_geos.ref
  testoutput/hyb-3dvar_gfs_0obs.ref
  testoutput/hyb-3dvar_gfs_1obs.ref
  testoutput/hyb-3dvar_gfs_aero.ref
  testoutput/hyb-3dvar_gfs_gnssro_ropp.ref
  testoutput/hyb-3dvar_gfs.ref
  testoutput/hyb-3dvar.ref
  testoutput/hyb-4dvar_fv3-gfs.ref
  testoutput/hyb-4dvar_pseudo-geos.ref
  testoutput/hyb-fgat_fv3-geos.ref
  testoutput/hyb-fgat_fv3-gfs.ref
  testoutput/hyb-fgat_fv3lm.ref
  testoutput/hyb-fgat_gfs_aero.ref
  testoutput/letkf_geos.ref
  testoutput/letkf_gfs.ref
  testoutput/letkf.ref
  testoutput/lgetkf_geos.ref
  testoutput/lgetkf_gfs.ref
  testoutput/lgetkf.ref
  testoutput/update_output.py
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
LINK_FILES( "${fv3jedi_testoutput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Grid data for femps
# -------------------
list( APPEND fv3jedi_femps_grid
Data/femps/fv3grid_c0006.nc4
Data/femps/fv3grid_c0012.nc4
Data/femps/fv3grid_c0024.nc4
Data/femps/fv3grid_c0048.nc4
Data/femps/fv3grid_c0096.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/femps)
LINK_FILES( "${fv3jedi_femps_grid}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Field sets
# ----------
list( APPEND fv3jedi_field_sets
Data/fieldsets/aerosols_geos.yaml
Data/fieldsets/aerosols_gfs.yaml
Data/fieldsets/cold_start_64.yaml
Data/fieldsets/cold_start_127.yaml
Data/fieldsets/dynamics_lam_cmaq.yaml
Data/fieldsets/dynamics.yaml
Data/fieldsets/lam_cmaq.yaml
Data/fieldsets/oro_data.yaml
Data/fieldsets/ufo.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fieldsets)
LINK_FILES( "${fv3jedi_field_sets}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Resources from other respositores
# ---------------------------------

# CRTM coefficient files
list( APPEND crtm_coeff_files
AerosolCoeff/Little_Endian/AerosolCoeff.bin
CloudCoeff/Little_Endian/CloudCoeff.bin
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/NPOESS.IRland.EmisCoeff.bin
EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin
EmisCoeff/VIS_Ice/SEcategory/Little_Endian/NPOESS.VISice.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/NPOESS.VISland.EmisCoeff.bin
EmisCoeff/VIS_Snow/SEcategory/Little_Endian/NPOESS.VISsnow.EmisCoeff.bin
EmisCoeff/VIS_Water/SEcategory/Little_Endian/NPOESS.VISwater.EmisCoeff.bin
EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/Little_Endian/USGS.IRland.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/USGS.VISland.EmisCoeff.bin
SpcCoeff/Little_Endian/hirs4_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/hirs4_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n15.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n18.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
SpcCoeff/Little_Endian/atms_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_npp.TauCoeff.bin
SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n18.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n19.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD1_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD1_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD2_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD2_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD3_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD3_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD4_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD4_g15.TauCoeff.bin
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g16.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g16.TauCoeff.bin
SpcCoeff/Little_Endian/ahi_himawari8.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ahi_himawari8.TauCoeff.bin
SpcCoeff/Little_Endian/radiometer_smap.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/radiometer_smap.TauCoeff.bin
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)

# Files needed for GEOS aerosol testing
if( ${GEOS-AERO_FOUND} )
  list( APPEND geos-aero_test_data
  Data/optics_BC.v1_3_.nc
  Data/optics_BRC.v1_5_.nc
  Data/optics_DU.v15_3_.nc
  Data/optics_NI.v2_5_.nc
  Data/optics_OC.v1_3_.nc
  Data/optics_SS.v3_3_.nc
  Data/optics_SU.v1_3_.nc
  )

  CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                           ${CMAKE_CURRENT_BINARY_DIR}/Data
                           ${geos-aero_test_data} )

  list( APPEND geos-aero_test_input
  testinput/geosaod.rc
  testinput/Chem_MieRegistry.rc)

  CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${geos-aero_test_input} )
endif()

# --------------------------------------------------------------------------------------------------

# Set the test tiering
# --------------------
set( FV3JEDI_TEST_TIER 1 )

# Override FV3JEDI_TEST_TIER using environment variable
if( DEFINED ENV{FV3JEDI_TEST_TIER} )
    set( FV3JEDI_TEST_TIER "$ENV{FV3JEDI_TEST_TIER}" )
endif()

message(STATUS "Running fv3jedi tests for tiers up to and including tier " ${FV3JEDI_TEST_TIER})

# --------------------------------------------------------------------------------------------------

# Repo wide configuration for getting test data
# ---------------------------------------------

# Get the current git branch/tag used by fv3-jedi source code
execute_process( COMMAND bash "-c" "git symbolic-ref -q --short HEAD || git describe --tags --exact-match"
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                 OUTPUT_VARIABLE GIT_BRANCH_FV3JEDI
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

message(STATUS "fv3-jedi is on branch/tag ${GIT_BRANCH_FV3JEDI}")

# For releases data is obtained for UCAR dash
set(ECBUILD_DOWNLOAD_BASE_URL_DEFAULT https://jedi-test-files.s3.amazonaws.com)

if(GIT_BRANCH_FV3JEDI STREQUAL "1.0.0")
  set(ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/dataset/147_miesch/file/)
  set(BRANCHPATH "/")
  set(MD5CHECK 0)
else()
  set(ECBUILD_DOWNLOAD_BASE_URL ${ECBUILD_DOWNLOAD_BASE_URL_DEFAULT})
  set(BRANCHPATH ${GIT_BRANCH_FV3JEDI})
  set(MD5CHECK 1)
endif()

message(STATUS "Git branch/tag path to use to obtain data " ${BRANCHPATH})

# Copy test data downloader script to bin
set ( DOWNLOADER fv3jedi_testdata_downloader.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${DOWNLOADER} ${CMAKE_BINARY_DIR}/bin/${DOWNLOADER} COPYONLY)

# Check for local path to data
if (DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES})
  set(LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}")
endif()

# Function for getting test data for different repos
# --------------------------------------------------

function(fv3jedi_data_download_test REPO DEFAULTBRANCH test_data_file LINKNG_PATH)

  # Special case of data structure for UCAR DASH release data
  if(GIT_BRANCH_FV3JEDI STREQUAL "1.0.0")
    set(REPOPATH "/")
  else()
    set(REPOPATH ${REPO})
  endif()

  message(STATUS "Searching for test data for ${REPO}")

  if( NOT DEFINED LOCAL_PATH_JEDI_TESTFILES )

    message(STATUS " LOCAL_PATH_JEDI_TESTFILES is not defined, download test files")

    # Check whether the URL exists
    ecbuild_check_urls(NAMES ${REPOPATH}/${BRANCHPATH}/${test_data_file}
                       RESULT SPECIFIC_TEST_FILES)

    # Check whether data exists at the remote path, if not use defaults
    message(STATUS " Looking for data at ${ECBUILD_DOWNLOAD_BASE_URL}/${REPOPATH}/${BRANCHPATH}/${test_data_file} ... ")
    if(${SPECIFIC_TEST_FILES} MATCHES 0)
      set(REPOPATHREMOTE   ${REPOPATH})
      set(BRANCHPATHREMOTE ${BRANCHPATH})
      message(STATUS "  ... data found, will download this" )
    else()
      set(ECBUILD_DOWNLOAD_BASE_URL  ${ECBUILD_DOWNLOAD_BASE_URL_DEFAULT})
      set(REPOPATHREMOTE             ${REPO})
      set(BRANCHPATHREMOTE           ${DEFAULTBRANCH})
      message(STATUS "  ... data not found, will download from defaults: ${ECBUILD_DOWNLOAD_BASE_URL}/${REPOPATHREMOTE}/${BRANCHPATHREMOTE}/${test_data_file}" )
    endif()

    # Note: do we at some point need to force the path using an ENV variable? Could be useful when
    # building where internet access is not available but tests run where it is. E.g. in SLURM.
    #message(STATUS " ${REPO} test data will be downloaded from: " ${ECBUILD_DOWNLOAD_BASE_URL}/${REPOPATHREMOTE}/${BRANCHPATHREMOTE}/${test_data_file})

    # Set target directory for the downloaded data
    set(OUTPUT_PATH ${CMAKE_BINARY_DIR}/test_data/${REPO}/${BRANCHPATHREMOTE})
    message(STATUS " will save ${REPO} test data to " ${OUTPUT_PATH})
    file(MAKE_DIRECTORY ${OUTPUT_PATH})

    # Add test for downloading IODA data
    ecbuild_add_test( TARGET    fv3jedi_test_tier1_get_test_data_${REPO}
                      TYPE      SCRIPT
                      COMMAND   ${CMAKE_BINARY_DIR}/bin/${DOWNLOADER}
                      ARGS      ${REPOPATHREMOTE}/${BRANCHPATHREMOTE}/${test_data_file}  # Part after URL
                                ${test_data_file}                                        # Filename
                                ${OUTPUT_PATH}                                           # Output path
                                ${ECBUILD_DOWNLOAD_BASE_URL}                             # Url
                                ${MD5CHECK}                                              # Check for md5 match
                                testoutput/fv3jedi_test_tier1_get_test_data_${REPO}.log )
  else()
    if(EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/${REPO}/${GIT_BRANCH_FV3JEDI}/${test_data_file})
      set(INPUT_PATH ${LOCAL_PATH_JEDI_TESTFILES}/${REPO}/${GIT_BRANCH_FV3JEDI})
      set(OUTPUT_PATH ${CMAKE_BINARY_DIR}/test_data/${REPO}/${GIT_BRANCH_FV3JEDI})
    else()
      set(INPUT_PATH ${LOCAL_PATH_JEDI_TESTFILES}/${REPO}/develop)
      set(OUTPUT_PATH ${CMAKE_BINARY_DIR}/test_data/${REPO}/develop)
    endif()
    message(STATUS " Use local test files from ${INPUT_PATH}")
    message(STATUS "  ... copy test files to ${OUTPUT_PATH}")
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy_directory ${INPUT_PATH} ${OUTPUT_PATH} )
  endif()

  set(${LINKNG_PATH} ${OUTPUT_PATH} PARENT_SCOPE)

endfunction()

# Perform fv3-jedi test data download
# -----------------------------------

fv3jedi_data_download_test( fv3-jedi develop fv3jedi_testdatainput_tier_1.tar.gz LINKNG_PATH)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${LINKNG_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs )

# Perform ioda test data download
# -------------------------------

fv3jedi_data_download_test( ioda develop ioda_testinput_tier_1.tar.gz LINKNG_PATH)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${LINKNG_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/obs )

# Perform crtm test data download
# -------------------------------

fv3jedi_data_download_test( crtm 2.3.0 crtm_coefficients.tar.gz LINKNG_PATH)

create_symlink_filename( ${LINKNG_PATH}
                         ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm
                         ${crtm_coeff_files} )

# --------------------------------------------------------------------------------------------------

# The following is a wrapper to simplify the generation of tests.
# There are two types of tests:
#  1) fv3-jedi executables (EXE must be given)
#  2) interface tests  (SRC must be given)
# fv3-jedi executables by default will also run a second step where the output test
# log stream results are compared with reference values. The interface tests
# do not do this.
#
# Arguments:
#  NOCOMPARE - if flag is present, no compare test is performed
#  NAME      - the name of the test (after test_fv3jedi is prepended), as well
#              as the base name of the input confg file
#  EXE       - name of executable to be used.
#  SRC       - The name of the unit test source file to compile
#  CFG       - The name of the yaml file, if the default testinput/${NAME}.yaml
#              is not to be used
#  TOL       - The tolerances to use for the test "float_tolerance int_tolerance"
#              If none are given, ${FV3JEDI_TESTS_DEFAULT_TOL} is used. Only used
#              if an EXE with no NOCOMPARE flag set.
#  MPI       - The number of MPI PEs to use. If not specified, ${FV3JEDI_TESTS_DEFAULT_MPI}
#              will be used
# TEST_DEPENDS - list of tests this test depends on

set( FV3JEDI_TESTS_DEFAULT_TOL "0.0;0" CACHE STRING
     "Default test tolerances for compare.py based tests. (float_tol;int_tol) ")

set( FV3JEDI_TESTS_DEFAULT_MPI 6 CACHE STRING
     "Default number of PEs to use for fv3-jedi MPI based tests")

set( FV3JEDI_TESTS_FORC_DEFAULT_TOL OFF CACHE BOOL
     "If true, the specified test tolerances are ignored and FV3JEDI_TESTS_DEFAULT_TOL is used" )

function(fv3jedi_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOCOMPARE TRAPFPE INTERFACETEST)
  set(singlevals NAME EXE CFG MPI REF TIER)
  set(multivals  TEST_DEPENDS TOL)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # set the tolerances to defaults if none given, or if we want to override
  IF(NOT ARG_TOL OR FV3JEDI_TESTS_FORC_DEFAULT_TOL)
    SET(ARG_TOL ${FV3JEDI_TESTS_DEFAULT_TOL})
  ENDIF()

  # determine if floating point error trapping should be set
  if ( ARG_TRAPFPE )
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=1")
  else()
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=0")
  endif()

  # determine the default config file name
  if ( ARG_CFG )
    set ( CONFIG_FILE testinput/${ARG_CFG} )
  else()
    set ( CONFIG_FILE testinput/${ARG_NAME}.yaml )
  endif()

  # MPI PEs
  set( MPI ${FV3JEDI_TESTS_DEFAULT_MPI} )
  if ( ARG_MPI )
    set( MPI ${ARG_MPI})
  endif()

  # Test executable must be provided
  if (NOT ARG_EXE)
    message( FATAL_ERROR "fv3jedi_add_test required EXE argument." )
  endif()

  # Set the tier for the test
  if (NOT ARG_TIER)
    set( TESTTIER 1 )
  else()
    set( TESTTIER ${ARG_TIER} )
  endif()

  # Default reference output is same as input name
  if ( ARG_REF)
    set( TESTREF ${ARG_REF})
  else()
    set( TESTREF ${ARG_NAME})
  endif()

  # Are we building a unit test / or running an executable?
  if ( TESTTIER LESS_EQUAL FV3JEDI_TEST_TIER )
    if ( ARG_INTERFACETEST )

      # Building a unit test, therfore also assume no compare step
      ecbuild_add_test( TARGET  fv3jedi_test_tier${TESTTIER}_${ARG_NAME}
                        TYPE    EXE
                        COMMAND ${CMAKE_BINARY_DIR}/bin/${ARG_EXE}
                        ARGS    ${CONFIG_FILE}
                        MPI     ${MPI}
                        ENVIRONMENT ${TRAPFPE_ENV}
                        TEST_DEPENDS test_fv3jedi_${ARG_TEST_DEPENDS} )

    else()

      # Get tolerances to run the test with
      LIST(GET ARG_TOL 0 TOL_F)
      LIST(GET ARG_TOL 1 TOL_I)

      # Find the MPI command
      set(MPI_CMD "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI}")

      # Run test executable and the compare script
      ecbuild_add_test( TARGET  fv3jedi_test_tier${TESTTIER}_${ARG_NAME}
                        TYPE    SCRIPT
                        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/test_wrapper.sh"
                        ARGS    "${CMAKE_BINARY_DIR}/bin/${ARG_EXE}"
                                ${CONFIG_FILE}
                        ENVIRONMENT
                                ${TRAPFPE_ENV}
                                COMPARE_SCRIPT=${CMAKE_BINARY_DIR}/bin/compare.py
                                COMPARE_TESTNAME=${ARG_NAME}
                                COMPARE_REFERENCE=${TESTREF}
                                COMPARE_TOL_F=${TOL_F}
                                COMPARE_TOL_I=${TOL_I}
                                MPI_CMD=${MPI_CMD}
                                SKIP_COMPARE=${ARG_NOCOMPARE}
                        DEPENDS ${ARG_EXE}
                        TEST_DEPENDS test_fv3jedi_${ARG_TEST_DEPENDS})
    endif()
  endif()

endfunction()

# --------------------------------------------------------------------------------------------------


# Create tests
# ------------

# Executables
ecbuild_add_executable( TARGET  test_fv3jedi_errorcovariance.x
                        SOURCES mains/TestErrorCovariance.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_geometry.x
                        SOURCES mains/TestGeometry.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_getvalues.x
                        SOURCES mains/TestGetValues.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_geometryiterator.x
                        SOURCES mains/TestGeometryIterator.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_increment.x
                        SOURCES mains/TestIncrement.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_lineargetvalues.x
                        SOURCES mains/TestLinearGetValues.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_linearmodel.x
                        SOURCES mains/TestLinearModel.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_linearvariablechange.x
                        SOURCES mains/TestLinVarCha.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_localization.x
                        SOURCES mains/TestLocalization.cc
                        LIBS    fv3jedi )

#ecbuild_add_executable( TARGET  test_fv3jedi_modelauxcontrol.x
#                        SOURCES mains/TestModelAuxControl.cc
#                        LIBS    fv3jedi )
#
#ecbuild_add_executable( TARGET  test_fv3jedi_modelauxcovariance.x
#                        SOURCES mains/TestModelAuxCovariance.cc
#                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_model.x
                        SOURCES mains/TestModel.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_state.x
                        SOURCES mains/TestState.cc
                        LIBS    fv3jedi )

ecbuild_add_executable( TARGET  test_fv3jedi_variablechange.x
                        SOURCES mains/TestVarCha.cc
                        LIBS    fv3jedi )

#Make some output directories for test data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/gsirf)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/forecast)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/increment)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/ensvariance)

set(members mem001;mem002;mem003;mem004)
foreach(MEM IN LISTS members)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/geos/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/gfs/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/geos/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/gfs/${MEM})
endforeach()

set(members mem000;mem001;mem002;mem003;mem004;mem005;mem006;mem007;mem008;mem009;mem010)
foreach(MEM IN LISTS members)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/geos/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/gfs/${MEM})
endforeach()

# Stage directory with directories for running the models
# -------------------------------------------------------
if(FV3_FORECAST_MODEL MATCHES "GEOS" OR FV3_FORECAST_MODEL MATCHES "UFS")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/ModelRunDirs/${FV3_FORECAST_MODEL})
  file(GLOB model_files "${FV3_FORECAST_MODEL_RUNDIR}/*")
  LINK_FILES_DIR("${model_files}" ${CMAKE_CURRENT_BINARY_DIR}/Data/ModelRunDirs/${FV3_FORECAST_MODEL}/)
  if (FV3_FORECAST_MODEL MATCHES "GEOS")
    execute_process( COMMAND rm ${CMAKE_CURRENT_BINARY_DIR}/Data/ModelRunDirs/${FV3_FORECAST_MODEL}/cap_restart )
    execute_process( COMMAND cp ${FV3_FORECAST_MODEL_RUNDIR}/cap_restart ${CMAKE_CURRENT_BINARY_DIR}/Data/ModelRunDirs/${FV3_FORECAST_MODEL}/cap_restart )
    execute_process( COMMAND chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/Data/ModelRunDirs/${FV3_FORECAST_MODEL}/cap_restart )
  endif()
endif()

#Tolerances for compare script
set (ctol 0.0) #Max relative difference
set (idif 0)   #Max difference

ecbuild_add_resources( TARGET   fv3jedi_test_scripts
                       SOURCES_PACK
                       ${fv3jedi_test_input}
                     )

# Tier 1 tests
# ------------

# Interface tests
# ---------------
fv3jedi_add_test( NAME    geometry_gfs
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_gfs127
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_geos
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_lam_cmaq
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_2d
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_nested
                  EXE     test_fv3jedi_geometry.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_iterator_gfs
                  EXE     test_fv3jedi_geometryiterator.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_iterator_geos
                  EXE     test_fv3jedi_geometryiterator.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    geometry_iterator_lam_cmaq
                  EXE     test_fv3jedi_geometryiterator.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_gfs
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_geos
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_lam_cmaq
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_lam_atmos
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_geos_aerosol
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    state_geos_aerosol_ext
                  EXE     test_fv3jedi_state.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    increment_gfs
                  EXE     test_fv3jedi_increment.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    increment_geos
                  EXE     test_fv3jedi_increment.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    getvalues_gfs
                  EXE     test_fv3jedi_getvalues.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    getvalues_geos
                  EXE     test_fv3jedi_getvalues.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    lineargetvalues_gfs
                  EXE     test_fv3jedi_lineargetvalues.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    lineargetvalues_geos
                  EXE     test_fv3jedi_lineargetvalues.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    errorcovariance
                  EXE     test_fv3jedi_errorcovariance.x
                  INTERFACETEST )

# Variable changes tests (some required for latter tests)
# -------------------------------------------------------
fv3jedi_add_test( NAME    convertstate_gfs
                  EXE     fv3jedi_convertstate.x
                  TOL     1.0e-3 0)

fv3jedi_add_test( NAME    convertstate_geos
                  EXE     fv3jedi_convertstate.x )

fv3jedi_add_test( NAME    convertstate_geos_crtm_surface
                  EXE     fv3jedi_convertstate.x )

if( FV3_FORECAST_MODEL MATCHES "FV3CORE" )
  fv3jedi_add_test( NAME    convertstate_gfs_cold2fv3
                    EXE     fv3jedi_convertstate.x )
endif()

# Remaining interface tests
# -------------------------
fv3jedi_add_test( NAME    variablechange_gfs
                  EXE     test_fv3jedi_variablechange.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    variablechange_geos
                  EXE     test_fv3jedi_variablechange.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    linearvariablechange_gfs
                  EXE     test_fv3jedi_linearvariablechange.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    linearvariablechange_geos
                  EXE     test_fv3jedi_linearvariablechange.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    model_pseudo-gfs
                  EXE     test_fv3jedi_model.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    model_pseudo-geos
                  EXE     test_fv3jedi_model.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    model_fv3lm
                  EXE     test_fv3jedi_model.x
                  INTERFACETEST )

if( FV3_FORECAST_MODEL MATCHES "UFS" AND NOT SKIP_LARGE_TESTS)
fv3jedi_add_test( NAME    model_ufs
                  EXE     test_fv3jedi_model.x
                  INTERFACETEST )
endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT SKIP_LARGE_TESTS)
fv3jedi_add_test( NAME    model_geos
                  EXE     test_fv3jedi_model.x
                  INTERFACETEST )
endif()

fv3jedi_add_test( NAME    linearmodel
                  EXE     test_fv3jedi_linearmodel.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    linearmodel_physics
                  EXE     test_fv3jedi_linearmodel.x
                  INTERFACETEST )

# Forecast tests
# --------------

fv3jedi_add_test( NAME    forecast_fv3lm
                  EXE     fv3jedi_forecast.x
                  TOL     1e-3 0)

if( FV3_FORECAST_MODEL MATCHES "UFS" AND NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME    forecast_ufs
                    EXE     fv3jedi_forecast.x
                    TOL     0.0 0)
endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME    forecast_geos
                    EXE     fv3jedi_forecast.x)
endif()

# H(x)
# ----
fv3jedi_add_test( NAME    hofx_fv3lm
                  EXE     fv3jedi_hofx.x
                  TOL     1e-3 0)

# H(x) when ROPP is enabled
# -------------------------
if( ${ROPP-UFO_FOUND} )
fv3jedi_add_test( NAME    hofx_gnssro_ropp
                  EXE     fv3jedi_hofx.x
                  TOL     1e-1 0)
endif()

# H(x) no model
# -------------
fv3jedi_add_test( NAME    hofx_nomodel
                  EXE     fv3jedi_hofx_nomodel.x
                  TOL     1e-3 0)

fv3jedi_add_test( NAME    hofx_gfs_aero
                  EXE     fv3jedi_hofx_nomodel.x
                  TOL     1e-3 0)

fv3jedi_add_test( NAME    hofx_nomodel_lam_cmaq
                  EXE     fv3jedi_hofx_nomodel.x
                  TOL     1e-3 0)

# Test to save GeoVaLs to file
# ----------------------------
fv3jedi_add_test( NAME    hofx_save_geovals
                  EXE     fv3jedi_hofx_nomodel.x
                  TOL     1e-3 0)

# B Matrix tests
# --------------
fv3jedi_add_test( NAME    bumpparameters_nicas_geos
                  EXE     fv3jedi_parameters.x)

fv3jedi_add_test( NAME    bumpparameters_nicas_gfs
                  EXE     fv3jedi_parameters.x)

fv3jedi_add_test( NAME    bumpparameters_nicas_lam_cmaq
                  EXE     fv3jedi_parameters.x)

if (NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME  bumpparameters_nicas_geos_p12
                    EXE   fv3jedi_parameters.x
                    MPI   12
                    REF   bumpparameters_nicas_geos.run )

  fv3jedi_add_test( NAME  bumpparameters_nicas_gfs_p12
                    EXE   fv3jedi_parameters.x
                    MPI   12
                    REF   bumpparameters_nicas_gfs.run )
endif()

fv3jedi_add_test( NAME    errorcovariance_bump
                  EXE     test_fv3jedi_errorcovariance.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    localization_bump
                  EXE     test_fv3jedi_localization.x
                  INTERFACETEST )

fv3jedi_add_test( NAME    dirac_geos
                  EXE     fv3jedi_dirac.x)

fv3jedi_add_test( NAME    dirac_gfs
                  EXE     fv3jedi_dirac.x)

fv3jedi_add_test( NAME    dirac_gfs_bal
                  EXE     fv3jedi_dirac.x)

fv3jedi_add_test( NAME    dirac_gfs_c2a
                  EXE     fv3jedi_dirac.x)

# Ensemble variance
# -----------------
fv3jedi_add_test( NAME    ensvariance_gfs
                  EXE     fv3jedi_ensvariance.x)

fv3jedi_add_test( NAME    ensvariance_geos
                  EXE     fv3jedi_ensvariance.x)

# 3DVar, 3DEnVar and hybrid 3DVar
# -------------------------------
fv3jedi_add_test( NAME         hyb-3dvar
                  EXE          fv3jedi_var.x
                  TOL          1.0e-1 0
                  TEST_DEPENDS bumpparameters_nicas_gfs )

if (NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME  hyb-3dvar_p12
                    EXE   fv3jedi_var.x
                    MPI   12
                    REF   hyb-3dvar.run
                    TEST_DEPENDS bumpparameters_nicas_gfs_p12 hyb-3dvar)
endif()

# 3DVar with GSI balance operator
#fv3jedi_add_test( NAME    3dvar_gfs_bal
#                  EXE     fv3jedi_var.x )

# 3DVar with zero obs is not working yet, so these 2 tests are are commented out. The YAML,
# reference and obs files are are already in place.
#fv3jedi_add_test( NAME         hyb-3dvar_gfs_0obs
#                  EXE          fv3jedi_var.x
#                  TEST_DEPENDS bumpparameters_nicas_gfs)

if( ${ROPP-UFO_FOUND} )
  fv3jedi_add_test( NAME         hyb-3dvar_gfs_gnssro_ropp
                    EXE          fv3jedi_var.x
                    TOL          1.0e-1 0
                    TEST_DEPENDS bumpparameters_nicas_gfs )
endif()

# Chemistry / AQ Variational Tests
# -----------
fv3jedi_add_test( NAME         3dvar_lam_cmaq
                  EXE          fv3jedi_var.x
                  TOL          1.0e-1 0
                  TEST_DEPENDS bumpparameters_nicas_lam_cmaq )

fv3jedi_add_test( NAME         hyb-3dvar_gfs_aero
                  EXE          fv3jedi_var.x
                  TOL          1.0e-1 0
                  TEST_DEPENDS bumpparameters_nicas_gfs )

if( ${GEOS-AERO_FOUND} )
  fv3jedi_add_test( NAME         3denvar_geos_aero
                    EXE          fv3jedi_var.x
                    NOCOMPARE
                    TEST_DEPENDS test_fv3jedi_bumpparameters_nicas_geos )
endif()

# Hybrid FGAT
# -----------
fv3jedi_add_test( NAME         hyb-fgat_fv3lm
                  EXE          fv3jedi_var.x
                  TEST_DEPENDS bumpparameters_nicas_geos )

fv3jedi_add_test( NAME         hyb-fgat_gfs_aero
                  EXE          fv3jedi_var.x
                  TEST_DEPENDS bumpparameters_nicas_gfs )

# 4DEnVar
# -------
if (NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME         4denvar
                    EXE          fv3jedi_var.x
                    MPI          18
                    TEST_DEPENDS bumpparameters_nicas_gfs )
endif()

# Hybrid 4DVar
# ------------
fv3jedi_add_test( NAME         hyb-4dvar_pseudo-geos
                  EXE          fv3jedi_var.x
                  TEST_DEPENDS bumpparameters_nicas_geos )

# Diff states
# -----------
fv3jedi_add_test( NAME         diffstates_gfs
                  EXE          fv3jedi_diffstates.x
                  TOL          1.0e-3 0
                  TEST_DEPENDS hyb-3dvar )

fv3jedi_add_test( NAME         diffstates_geos
                  EXE          fv3jedi_diffstates.x
                  TOL          1.0e-4 0
                  TEST_DEPENDS hyb-4dvar_pseudo-geos )

fv3jedi_add_test( NAME         diffstates_lam_cmaq
                  EXE          fv3jedi_diffstates.x
                  TOL          1.0e-3 0
                  TEST_DEPENDS 3dvar_lam_cmaq )

# Add increment
# -------------
fv3jedi_add_test( NAME         addincrement_gfs
                  EXE          fv3jedi_addincrement.x
                  TOL          1.0e-3 0
                  TEST_DEPENDS diffstates_gfs )

if (NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME         addincrement_geos
                    EXE          fv3jedi_addincrement.x
                    MPI          12
                    TOL          1.0e-4 0
                    TEST_DEPENDS diffstates_geos )
endif()

# LETKF / LGETKF
# --------------
fv3jedi_add_test( NAME         letkf
                  EXE          fv3jedi_letkf.x
                  TOL          1.0e-2 0 )

fv3jedi_add_test( NAME         lgetkf
                  EXE          fv3jedi_letkf.x
                  TOL          1.0e-2 0 )

# EnsHofx tests
# -------------
if (NOT SKIP_LARGE_TESTS)
  fv3jedi_add_test( NAME         enshofx_fv3lm
                    EXE          fv3jedi_enshofx.x
                    REF          hofx_fv3lm.run
                    MPI          12
                    TOL          1.0e-2 0 )
endif()

# Ensemble data assimilation tests
# --------------------------------
if (NOT SKIP_LARGE_TESTS)

  # 3DVar
  fv3jedi_add_test( NAME         eda_3dvar
                    EXE          fv3jedi_eda.x
                    MPI          12
                    NOCOMPARE
                    TEST_DEPENDS bumpparameters_nicas_gfs )

endif()

# Test of time taken
# ------------------
if ( DEFINED TEST_TIME_REF_FV3 )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/testoutput/${TEST_TIME_REF_FV3}
          ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${TEST_TIME_REF_FV3} )

  ecbuild_add_test( TARGET fv3jedi_test_tier1_run_times
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/run_time_test.py
                    ARGS testoutput/${TEST_TIME_REF_FV3}
                    ${CMAKE_BINARY_DIR}/fv3-jedi/Testing/Temporary/CTestCostData.txt
                    testinput/test_time.yaml)
else()
  message ( STATUS "TEST_TIME not built")
endif()



# Tier 2 tests
# ------------

# Full observation h(x) for each model
# ------------------------------------
fv3jedi_add_test( NAME    hofx_fv3-gfs
                  EXE     fv3jedi_hofx.x
                  TOL     1e-3 0
                  TIER    2)

fv3jedi_add_test( NAME    hofx_fv3-geos
                  EXE     fv3jedi_hofx.x
                  TOL     1e-2 0
                  TIER    2)

# Full observation hyb-3DVar for each model
# -----------------------------------------
fv3jedi_add_test( NAME         hyb-3dvar_gfs
                  EXE          fv3jedi_var.x
                  TOL          1.0e-1 0
                  TEST_DEPENDS bumpparameters_nicas_gfs
                  TIER    2 )

fv3jedi_add_test( NAME         hyb-3dvar_geos
                  EXE          fv3jedi_var.x
                  TOL          1.0e-1 0
                  TEST_DEPENDS bumpparameters_nicas_geos
                  TIER    2 )

# EDA 4DVar
# ---------
fv3jedi_add_test( NAME         eda_4dvar_fv3lm
                  EXE          fv3jedi_eda.x
                  MPI          12
                  NOCOMPARE
                  TEST_DEPENDS bumpparameters_nicas_gfs
                  TIER    2 )
