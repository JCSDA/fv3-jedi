# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include( ${CMAKE_SOURCE_DIR}/jedicmake/cmake/Functions/git_functions.cmake )

# Option not to run with large tests
# ----------------------------------
option(BUILD_LARGE_TESTS "build tests with more than 6 processors" ON)

# Create the Data/ directory
# --------------------------
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# FV3 config files
# ----------------
list( APPEND fv3jedi_test_fv3files
Data/fv3files/akbk64.nc4
Data/fv3files/akbk72.nc4
Data/fv3files/akbk127.nc4
Data/fv3files/field_table_gmao
Data/fv3files/field_table_gfdl
Data/fv3files/field_table_lam_cmaq
Data/fv3files/field_table_lam_rrfs
Data/fv3files/fmsmpp.nml
Data/fv3files/input_gfs_c12.nml
Data/fv3files/input_gfs_c12_p12.nml
Data/fv3files/input_gfs_c24.nml
Data/fv3files/input_gfs_c24_p12.nml
Data/fv3files/input_geos_c12.nml
Data/fv3files/input_geos_c12_p12.nml
Data/fv3files/input_geos_c24.nml
Data/fv3files/input_lam_cmaq.nml
Data/fv3files/input_lam_rrfs.nml
Data/fv3files/inputpert_4dvar.nml
Data/fv3files/inputpert_fsoi.nml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fv3files)
LINK_FILES( "${fv3jedi_test_fv3files}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Static B-specific dates
# -----------------------
set(staticb_dates 2018041421;2018041500;2018041503)

# Test input scripts
# ------------------
list( APPEND fv3jedi_testinput
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
  testinput/bumpparameters_nicas_gfs_aero.yaml
  testinput/bumpparameters_nicas_gfs_aero_logp.yaml
  testinput/bumpparameters_nicas_gfs_aero_p12.yaml
  testinput/bumpparameters_nicas_lam_cmaq.yaml
  testinput/convertincrement_geos_model2geovals.yaml
  testinput/convertstate_geos_crtm_surface.yaml
  testinput/convertstate_geos_model2geovals.yaml
  testinput/convertstate_geos.yaml
  testinput/convertstate_gfs_coldstartwinds.yaml
  testinput/convertstate_gfs_vertremap.yaml
  testinput/convertstate_gfs_model2geovals.yaml
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
  testinput/eda_block_3dvar_1.yaml
  testinput/eda_block_3dvar_2.yaml
  testinput/eda_block_3dvar.yaml
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
  testinput/hofx_nomodel_lam_rrfs.yaml
  testinput/hofx_save_geovals.yaml
  testinput/hyb-3dvar_bc.yaml
  testinput/hyb-3dvar_fsoi_backward.yaml
  testinput/hyb-3dvar_fsoi_forward.yaml
  testinput/hyb-3dvar_geos.yaml
  testinput/hyb-3dvar_gfs_0obs.yaml
  testinput/hyb-3dvar_gfs_1obs.yaml
  testinput/hyb-3dvar_gfs_aero.yaml
  testinput/hyb-3dvar_gfs_aero_bump.yaml
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
  testinput/letkf_snow.yaml
  testinput/lineargetvalues_geos.yaml
  testinput/lineargetvalues_gfs.yaml
  testinput/linearmodel_physics.yaml
  testinput/linearmodel.yaml
  testinput/linearvariablechange_geos.yaml
  testinput/linearvariablechange_gfs.yaml
  testinput/localization_bump.yaml
  testinput/model_fv3lm.yaml
  testinput/model_geos.yaml
  testinput/model_ufs.yaml
  testinput/model_pseudo-geos.yaml
  testinput/model_pseudo-gfs.yaml
  testinput/obslocalization.yaml
  testinput/state_geos_aerosol_ext.yaml
  testinput/state_geos_aerosol.yaml
  testinput/state_geos.yaml
  testinput/state_gfs.yaml
  testinput/state_lam_cmaq.yaml
  testinput/state_lam_atmos.yaml
  testinput/staticb_dirac_global_gfs_6pe.yaml
  testinput/staticb_dirac_global_gfs_12pe.yaml
  testinput/staticb_dirac_local_gfs_6pe.yaml
  testinput/staticb_dirac_local_aero_12pe.yaml
  testinput/staticb_dirac_local_gfs_12pe.yaml
  testinput/staticb_cor_aero.yaml
  testinput/staticb_cor_gfs.yaml
  testinput/staticb_nicas_aero.yaml
  testinput/staticb_nicas_gfs.yaml
  testinput/staticb_psichitouv_gfs.yaml
  testinput/staticb_split_nicas_gfs.yaml
  testinput/staticb_split_psichitouv_gfs.yaml
  testinput/staticb_split_vbal_gfs.yaml
  testinput/staticb_var_gfs.yaml
  testinput/staticb_var_aero.yaml
  testinput/staticb_vbal_gfs.yaml
  testinput/test_time.yaml
  testinput/variablechange_geos.yaml
  testinput/variablechange_gfs.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
LINK_FILES( "${fv3jedi_testinput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

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
  testoutput/bumpparameters_nicas_gfs_aero.ref
  testoutput/bumpparameters_nicas_gfs_aero_logp.ref
  testoutput/bumpparameters_nicas_lam_cmaq.ref
  testoutput/convertincrement_geos_model2geovals.ref
  testoutput/convertstate_geos_crtm_surface.ref
  testoutput/convertstate_geos_model2geovals.ref
  testoutput/convertstate_geos.ref
  testoutput/convertstate_gfs_coldstartwinds.ref
  testoutput/convertstate_gfs_vertremap.ref
  testoutput/convertstate_gfs_model2geovals.ref
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
  testoutput/eda_3dvar.ref
  testoutput/eda_4dvar_fv3lm.ref
  testoutput/eda_block_3dvar.ref
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
  testoutput/hofx_nomodel_lam_rrfs.ref
  testoutput/hofx_nomodel_gfs.ref
  testoutput/hofx_nomodel.ref
  testoutput/hofx_save_geovals.ref
  testoutput/hyb-3dvar_bc.ref
  testoutput/hyb-3dvar_fsoi_backward.ref
  testoutput/hyb-3dvar_fsoi_forward.ref
  testoutput/hyb-3dvar_geos.ref
  testoutput/hyb-3dvar_gfs_0obs.ref
  testoutput/hyb-3dvar_gfs_1obs.ref
  testoutput/hyb-3dvar_gfs_aero.ref
  testoutput/hyb-3dvar_gfs_aero_bump.ref
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
  testoutput/letkf_snow.ref
  testoutput/lgetkf_geos.ref
  testoutput/lgetkf_gfs.ref
  testoutput/lgetkf.ref
  testoutput/staticb_dirac_global_gfs_6pe.ref
  testoutput/staticb_dirac_global_gfs_12pe.ref
  testoutput/staticb_dirac_local_gfs_6pe.ref
  testoutput/staticb_dirac_local_aero_12pe.ref
  testoutput/staticb_dirac_local_gfs_12pe.ref
  testoutput/staticb_cor_aero.ref
  testoutput/staticb_cor_gfs.ref
  testoutput/staticb_nicas_aero.ref
  testoutput/staticb_nicas_gfs.ref
  testoutput/staticb_psichitouv_gfs.ref
  testoutput/staticb_split_vbal_gfs.ref
  testoutput/staticb_split_psichitouv_gfs.ref
  testoutput/staticb_split_nicas_gfs.ref
  testoutput/staticb_var_aero.ref
  testoutput/staticb_var_gfs.ref
  testoutput/staticb_vbal_gfs.ref
  testoutput/update_output.py
)

foreach(DATE ${staticb_dates})
  list( APPEND fv3jedi_testoutput
    testoutput/staticb_ens-to-psichi_${DATE}_gfs.ref
    testoutput/staticb_unbal_${DATE}_gfs.ref
    testoutput/staticb_vbal_${DATE}_gfs.ref
    testoutput/staticb_var-mom_${DATE}_aero.ref
    testoutput/staticb_var-mom_${DATE}_gfs.ref
  )
endforeach()

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
Data/fieldsets/land.yaml
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
SpcCoeff/Little_Endian/amsr2_gcom-w1.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsr2_gcom-w1.TauCoeff.bin
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

# Files needed for GEOS aerosol testing
if( ${geos-aero_FOUND} )
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

  list( APPEND geos-aero_testinput
  testinput/geosaod.rc
  testinput/geosaod_aeronet.rc
  testinput/Chem_MieRegistry.rc)

  CREATE_SYMLINK_FILENAME( ${geos-aero_SOURCE_DIR}/test/
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${geos-aero_testinput} )
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
# Currently only use 2.3.0 CRTM coefficients
set( CRTM_COEFFS_BRANCH "2.3.0" )

function(ADD_DOWNLOAD_TEST repo test_files_dirname output_data_path)
  # If REPO branch is being built set GIT_BRANCH_FUNC to REPO's current branch.
  # If a tagged version of REPO is being built set GIT_TAG_FUNC to REPO's current tag.
  find_branch_name(REPO_DIR_NAME ${repo})

  if ( ${repo} STREQUAL "crtm" )
    set( GIT_TAG_FUNC ${CRTM_COEFFS_BRANCH} )
  endif()

  if( DEFINED GIT_BRANCH_FUNC )
    set( REPO_GIT_BRANCH ${GIT_BRANCH_FUNC} )
  elseif( DEFINED GIT_TAG_FUNC )
    set( REPO_GIT_BRANCH ${GIT_TAG_FUNC} )
  endif()



  # When env veriable LOCAL_PATH_JEDI_TESTFILES is set, simply link test files
  # to build directory. get_${repo}_test_data checks the existence of test file directory.
  # ${repo} test data must be stored in
  # ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${BRANCH}/testinput_tier_1
  if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} )
    set( LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}" )
    message( STATUS "use LOCAL_PATH_JEDI_TESTFILES: ${LOCAL_PATH_JEDI_TESTFILES}" )
    # If test data specific to testing branch exists locally use it.
    # If not use test data specific to develop branch.
    if( EXISTS ${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH} )
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/${REPO_GIT_BRANCH}" )
    else()
      set( TESTFILE_DIR_REPO "${LOCAL_PATH_JEDI_TESTFILES}/${repo}/develop" )
    endif()

    message( STATUS "Test data in ${TESTFILE_DIR_REPO} is linked to build directory" )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
    ${TESTFILE_DIR_REPO} )
    set ( REPO_DATA_DOWNLOADER fv3jedi_data_checker.py )

  # Any tagged version of REPO is being built.
  # Build fv3jedi_get_{repo}_test_data test to download test data from DASH.
  elseif( DEFINED GIT_TAG_FUNC)
    message( STATUS "Tagged version of ${repo} is used" )
    # set ARGS for get_${repo}_test_data

    # ECBUILD_DOWNLOAD_BASE_URL env var can be used to force test files and
    # crtm coef to be downloaded from other databases such as S3 instead of DASH
    # example ECBUILD_DOWNLOAD_BASE_URL=https://jedi-test-files.s3.amazonaws.com

    if( DEFINED ENV{ECBUILD_DOWNLOAD_BASE_URL} )
      set( ECBUILD_DOWNLOAD_BASE_URL "$ENV{ECBUILD_DOWNLOAD_BASE_URL}/${repo}/${GIT_TAG_FUNC}" )
    else()
      set( ECBUILD_DOWNLOAD_BASE_URL https://dashrepo.ucar.edu/api/v1/dataset/147_miesch/version/1.1.0/file )
    endif()

    set( DIRNAME ${GIT_TAG_FUNC} )
    set( checksum "0" )
    set( TESTFILE_DIR_REPO "${CMAKE_SOURCE_DIR}/test-data-release/${repo}/${DIRNAME}" )

    # Create test-data-release in source directory
    file( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/test-data-release )
    # temporary solution for release 1.1.0
    if ( ${repo} STREQUAL "crtm" )
      list( APPEND REPO_DATA_DOWNLOADER_ARGS
	          ${ECBUILD_DOWNLOAD_BASE_URL}
            ${TESTFILE_DIR_REPO}
            ${test_files_dirname}
            ${checksum} )
      file( MAKE_DIRECTORY ${TESTFILE_DIR_REPO} )
    else()
      list( APPEND REPO_DATA_DOWNLOADER_ARGS
            ${ECBUILD_DOWNLOAD_BASE_URL}
            ${CMAKE_SOURCE_DIR}/test-data-release
            ${test_files_dirname}
            ${checksum} )
    endif()
    message( STATUS "ECBUILD_DOWNLOAD_BASE_URL is set, ${repo} test files will
    be downloaded from ${ECBUILD_DOWNLOAD_BASE_URL}/test-data-release")
    set ( REPO_DATA_DOWNLOADER fv3jedi_testdata_downloader.py )

  # Any branch of REPO is being built.
  # repo-data repository is already cloned by bundle/CMakeLists.txt.
  # Link test files in repo-data repo to build directory.
  # get_repo_test_data checks the existence of test file directory
  else()
    set( TESTFILE_DIR_REPO ${CMAKE_SOURCE_DIR}/${repo}-data )
    list( APPEND REPO_DATA_DOWNLOADER_ARGS
          ${TESTFILE_DIR_REPO} )
    set (REPO_DATA_DOWNLOADER fv3jedi_data_checker.py )
  endif()

  # add test for downloading data
  ecbuild_add_test( TARGET    fv3jedi_get_${repo}_test_data
                    TYPE      SCRIPT
                    COMMAND   ${CMAKE_BINARY_DIR}/bin/${REPO_DATA_DOWNLOADER}
                    ARGS  ${REPO_DATA_DOWNLOADER_ARGS} )

  # pass output data path back to calling code
  set(${output_data_path} ${TESTFILE_DIR_REPO} PARENT_SCOPE)

  # Create download script for fv3jedi_get_*_test_data test
  set ( FILENAME ${REPO_DATA_DOWNLOADER})
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()


  add_custom_target( bin_fv3jedi_${repo}_test_download_scripts ALL
    COMMAND chmod +x ${DEST_FILE}
    DEPENDS ${DEST_FILE} )

endfunction()
#############################################
# add tests to download fv3-jedi & crtm data

ADD_DOWNLOAD_TEST( "fv3-jedi" fv3-jedi_testinput_tier_1.tar.gz FV3JEDI_FV3JEDI_TESTFILES_PATH )
ADD_DOWNLOAD_TEST( "crtm" crtm_coefficients.tar.gz CRTM_FV3JEDI_TESTFILES_PATH )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/obs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${FV3JEDI_FV3JEDI_TESTFILES_PATH}/testinput_tier_1/inputs
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${FV3JEDI_FV3JEDI_TESTFILES_PATH}/testinput_tier_1/obs
               ${CMAKE_CURRENT_BINARY_DIR}/Data/obs/testinput_tier_1)

create_symlink_filename( ${CRTM_FV3JEDI_TESTFILES_PATH}
                         ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm
                         ${crtm_coeff_files} )

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

ecbuild_add_executable( TARGET  test_fv3jedi_obslocalization.x
                        SOURCES mains/TestObsLocalization.cc
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
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/gfs_land)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/increment)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/ensvariance)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/convertstate)

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
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/gfs_land/${MEM})
endforeach()

# Static B directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb_aero)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb/balanced)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb/unbalanced)
set(members mem001;mem002;mem003;mem004;mem005;mem006;mem007;mem008;mem009;mem010)
foreach(MEM IN LISTS members)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb/balanced/${MEM})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb/unbalanced/${MEM})
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

ecbuild_add_resources( TARGET   fv3jedi_test_scripts
                       SOURCES_PACK
                       ${fv3jedi_testinput}
                     )

# Tier 1 tests
# ------------

# Interface tests
# ---------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_gfs
                  MPI      6
                  ARGS     testinput/geometry_gfs.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_gfs127
                  MPI      6
                  ARGS     testinput/geometry_gfs127.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_geos
                  MPI      6
                  ARGS     testinput/geometry_geos.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_lam_cmaq
                  MPI      6
                  ARGS     testinput/geometry_lam_cmaq.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_2d
                  MPI      6
                  ARGS     testinput/geometry_2d.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_nested
                  MPI      6
                  ARGS     testinput/geometry_nested.yaml
                  COMMAND  test_fv3jedi_geometry.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_iterator_gfs
                  MPI      6
                  ARGS     testinput/geometry_iterator_gfs.yaml
                  COMMAND  test_fv3jedi_geometryiterator.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_iterator_geos
                  MPI      6
                  ARGS     testinput/geometry_iterator_geos.yaml
                  COMMAND  test_fv3jedi_geometryiterator.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_geometry_iterator_lam_cmaq
                  MPI      6
                  ARGS     testinput/geometry_iterator_lam_cmaq.yaml
                  COMMAND  test_fv3jedi_geometryiterator.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_gfs
                  MPI      6
                  ARGS     testinput/state_gfs.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_geos
                  MPI      6
                  ARGS     testinput/state_geos.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_lam_cmaq
                  MPI      6
                  ARGS     testinput/state_lam_cmaq.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_lam_atmos
                  MPI      6
                  ARGS     testinput/state_lam_atmos.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_geos_aerosol
                  MPI      6
                  ARGS     testinput/state_geos_aerosol.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_state_geos_aerosol_ext
                  MPI      6
                  ARGS     testinput/state_geos_aerosol_ext.yaml
                  COMMAND  test_fv3jedi_state.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_increment_gfs
                  MPI      6
                  ARGS     testinput/increment_gfs.yaml
                  COMMAND  test_fv3jedi_increment.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_increment_geos
                  MPI      6
                  ARGS     testinput/increment_geos.yaml
                  COMMAND  test_fv3jedi_increment.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_getvalues_gfs
                  MPI      6
                  ARGS     testinput/getvalues_gfs.yaml
                  COMMAND  test_fv3jedi_getvalues.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_getvalues_geos
                  MPI      6
                  ARGS     testinput/getvalues_geos.yaml
                  COMMAND  test_fv3jedi_getvalues.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_lineargetvalues_gfs
                  MPI      6
                  ARGS     testinput/lineargetvalues_gfs.yaml
                  COMMAND  test_fv3jedi_lineargetvalues.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_lineargetvalues_geos
                  MPI      6
                  ARGS     testinput/lineargetvalues_geos.yaml
                  COMMAND  test_fv3jedi_lineargetvalues.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_errorcovariance
                  MPI      6
                  ARGS     testinput/errorcovariance.yaml
                  COMMAND  test_fv3jedi_errorcovariance.x )


# Variable changes tests (some required for latter tests)
# -------------------------------------------------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_gfs
                  MPI      6
                  ARGS     testinput/convertstate_gfs.yaml
                  COMMAND  fv3jedi_convertstate.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_geos
                  MPI      6
                  ARGS     testinput/convertstate_geos.yaml
                  COMMAND  fv3jedi_convertstate.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_geos_crtm_surface
                  MPI      6
                  ARGS     testinput/convertstate_geos_crtm_surface.yaml
                  COMMAND  fv3jedi_convertstate.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_gfs_c2ll
                  MPI      6
                  ARGS     testinput/convertstate_gfs_c2ll.yaml
                  COMMAND  fv3jedi_convertstate.x )

if( FV3_FORECAST_MODEL MATCHES "FV3CORE" )
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_gfs_coldstartwinds
                    MPI      6
                    ARGS     testinput/convertstate_gfs_coldstartwinds.yaml
                    COMMAND  fv3jedi_convertstate.x )

  ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_gfs_vertremap
                    MPI      6
                    ARGS     testinput/convertstate_gfs_vertremap.yaml
                    COMMAND  fv3jedi_convertstate.x
                    TEST_DEPENDS fv3jedi_test_tier1_convertstate_gfs_coldstartwinds)
endif()

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_gfs_model2geovals
                  MPI      6
                  ARGS     testinput/convertstate_gfs_model2geovals.yaml
                  COMMAND  fv3jedi_convertstate.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertstate_geos_model2geovals
                  MPI      6
                  ARGS     testinput/convertstate_geos_model2geovals.yaml
                  COMMAND  fv3jedi_convertstate.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_convertincrement_geos_model2geovals
                  MPI      6
                  ARGS     testinput/convertincrement_geos_model2geovals.yaml
                COMMAND  fv3jedi_convertincrement.x )

# Remaining interface tests
# -------------------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_variablechange_gfs
                  MPI      6
                  ARGS     testinput/variablechange_gfs.yaml
                  COMMAND  test_fv3jedi_variablechange.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_variablechange_geos
                  MPI      6
                  ARGS     testinput/variablechange_geos.yaml
                  COMMAND  test_fv3jedi_variablechange.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_linearvariablechange_gfs
                  MPI      6
                  ARGS     testinput/linearvariablechange_gfs.yaml
                  COMMAND  test_fv3jedi_linearvariablechange.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_linearvariablechange_geos
                  MPI      6
                  ARGS     testinput/linearvariablechange_geos.yaml
                  COMMAND  test_fv3jedi_linearvariablechange.x )


ecbuild_add_test( TARGET   fv3jedi_test_tier1_model_pseudo-gfs
                  MPI      6
                  ARGS     testinput/model_pseudo-gfs.yaml
                  COMMAND  test_fv3jedi_model.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_model_pseudo-geos
                  MPI      6
                  ARGS     testinput/model_pseudo-geos.yaml
                  COMMAND  test_fv3jedi_model.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_model_fv3lm
                  MPI      6
                  ARGS     testinput/model_fv3lm.yaml
                  COMMAND  test_fv3jedi_model.x )


#if( FV3_FORECAST_MODEL MATCHES "UFS" AND BUILD_LARGE_TESTS)
#ecbuild_add_test( TARGET   fv3jedi_test_tier1_model_ufs
#                  MPI      6
#                  ARGS     testinput/model_ufs.yaml
#                  COMMAND  test_fv3jedi_model.x )
#endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND BUILD_LARGE_TESTS)
ecbuild_add_test( TARGET   fv3jedi_test_tier1_model_geos
                  MPI      6
                  ARGS     testinput/model_geos.yaml
                  COMMAND  test_fv3jedi_model.x )
endif()

ecbuild_add_test( TARGET   fv3jedi_test_tier1_linearmodel
                  MPI      6
                  ARGS     testinput/linearmodel.yaml
                  COMMAND  test_fv3jedi_linearmodel.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_linearmodel_physics
                  MPI      6
                  ARGS     testinput/linearmodel_physics.yaml
                  COMMAND  test_fv3jedi_linearmodel.x )


# Forecast tests
# --------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_forecast_fv3lm
                  MPI      6
                  ARGS     testinput/forecast_fv3lm.yaml
                  COMMAND  fv3jedi_forecast.x )

if( FV3_FORECAST_MODEL MATCHES "UFS" AND BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_forecast_ufs
                    MPI      6
                    ARGS     testinput/forecast_ufs.yaml
                    COMMAND  fv3jedi_forecast.x )
endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_forecast_geos
                    MPI      6
                    ARGS     testinput/forecast_geos.yaml
                    COMMAND  fv3jedi_forecast.x )
endif()

# H(x)
# ----
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_fv3lm
                  MPI      6
                  ARGS     testinput/hofx_fv3lm.yaml
                  COMMAND  fv3jedi_hofx.x )

# H(x) when ROPP is enabled
# -------------------------
if( ${ropp-ufo_FOUND} )
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_gnssro_ropp
                    MPI      6
                    ARGS     testinput/hofx_gnssro_ropp.yaml
                    COMMAND  fv3jedi_hofx.x )
endif()

# H(x) no model
# -------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_nomodel
                  MPI      6
                  ARGS     testinput/hofx_nomodel.yaml
                  COMMAND  fv3jedi_hofx_nomodel.x )

if( ${geos-aero_FOUND} )
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_gfs_aero
                    MPI      6
                    ARGS     testinput/hofx_gfs_aero.yaml
                    COMMAND  fv3jedi_hofx.x )
endif()

ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_nomodel_lam_cmaq
                  MPI      6
                  ARGS     testinput/hofx_nomodel_lam_cmaq.yaml
                  COMMAND  fv3jedi_hofx_nomodel.x )

# Test to save GeoVaLs to file
# ----------------------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hofx_save_geovals
                  MPI      6
                  ARGS     testinput/hofx_save_geovals.yaml
                  COMMAND  fv3jedi_hofx_nomodel.x )

# B Matrix tests
# --------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_geos
                  MPI      6
                  ARGS     testinput/bumpparameters_nicas_geos.yaml
                  COMMAND  fv3jedi_error_covariance_training.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_gfs
                  MPI      6
                  ARGS     testinput/bumpparameters_nicas_gfs.yaml
                  COMMAND  fv3jedi_error_covariance_training.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero
                  MPI      6
                  ARGS     testinput/bumpparameters_nicas_gfs_aero.yaml
                  COMMAND  fv3jedi_error_covariance_training.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero_logp
                  MPI      6
                  ARGS     testinput/bumpparameters_nicas_gfs_aero_logp.yaml
                  COMMAND  fv3jedi_error_covariance_training.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_lam_cmaq
                  MPI      6
                  ARGS     testinput/bumpparameters_nicas_lam_cmaq.yaml
                  COMMAND  fv3jedi_error_covariance_training.x )

if (BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_geos_p12
                    MPI      12
                    ARGS     testinput/bumpparameters_nicas_geos_p12.yaml
                    COMMAND  fv3jedi_error_covariance_training.x )

  ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_gfs_p12
                    MPI      12
                    ARGS     testinput/bumpparameters_nicas_gfs_p12.yaml
                    COMMAND  fv3jedi_error_covariance_training.x )

  ecbuild_add_test( TARGET   fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero_p12
                    MPI      12
                    ARGS     testinput/bumpparameters_nicas_gfs_aero_p12.yaml
                    COMMAND  fv3jedi_error_covariance_training.x )
endif()

ecbuild_add_test( TARGET   fv3jedi_test_tier1_errorcovariance_bump
                  MPI      6
                  ARGS     testinput/errorcovariance_bump.yaml
                  COMMAND  test_fv3jedi_errorcovariance.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_localization_bump
                  MPI      6
                  ARGS     testinput/localization_bump.yaml
                  COMMAND  test_fv3jedi_localization.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_dirac_geos
                  MPI      6
                  ARGS     testinput/dirac_geos.yaml
                  COMMAND  fv3jedi_dirac.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_dirac_gfs
                  MPI      6
                  ARGS     testinput/dirac_gfs.yaml
                  COMMAND  fv3jedi_dirac.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_dirac_gfs_bal
                  MPI      6
                  ARGS     testinput/dirac_gfs_bal.yaml
                  COMMAND  fv3jedi_dirac.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_dirac_gfs_c2a
                  MPI      6
                  ARGS     testinput/dirac_gfs_c2a.yaml
                  COMMAND  fv3jedi_dirac.x )

# Ensemble variance
# -----------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_ensvariance_gfs
                  MPI      6
                  ARGS     testinput/ensvariance_gfs.yaml
                  COMMAND  fv3jedi_ensvariance.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_ensvariance_geos
                  MPI      6
                  ARGS     testinput/ensvariance_geos.yaml
                  COMMAND  fv3jedi_ensvariance.x )

# 3DVar, 3DEnVar and hybrid 3DVar
# -------------------------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar
                  MPI      6
                  ARGS     testinput/hyb-3dvar.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

if (BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_p12
                    MPI      12
                    ARGS     testinput/hyb-3dvar_p12.yaml
                    COMMAND  fv3jedi_var.x
                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs_p12 fv3jedi_test_tier1_hyb-3dvar )
endif()

# 3DVar with GSI balance operator
#ecbuild_add_test( TARGET   fv3jedi_test_tier1_3dvar_gfs_bal
#                  MPI      6
#                  ARGS     testinput/3dvar_gfs_bal.yaml
#                  COMMAND  fv3jedi_var.x )

# 3DVar with zero obs is not working yet, so these 2 tests are are commented out. The YAML,
# reference and obs files are are already in place.
#ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_gfs_0obs
#                  MPI      6
#                  ARGS     testinput/hyb-3dvar_gfs_0obs.yaml
#                  COMMAND  fv3jedi_var.x
#                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

# 3DVar with bias correction on for one of the satellites
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_bc
                  MPI      6
                  ARGS     testinput/hyb-3dvar_bc.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

if( ${ropp-ufo_FOUND} )
#  ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_gfs_gnssro_ropp
#                    MPI      6
#                    ARGS     testinput/hyb-3dvar_gfs_gnssro_ropp.yaml
#                    COMMAND  fv3jedi_var.x
#                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )
endif()

# Chemistry / AQ Variational Tests
# -----------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_3dvar_lam_cmaq
                  MPI      6
                  ARGS     testinput/3dvar_lam_cmaq.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_lam_cmaq )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_gfs_aero
                  MPI      6
                  ARGS     testinput/hyb-3dvar_gfs_aero.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero )

if( ${geos-aero_FOUND} )
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_3denvar_geos_aero
                    MPI      6
                    ARGS     testinput/3denvar_geos_aero.yaml
                    COMMAND  fv3jedi_var.x
                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_geos )
endif()

# Hybrid FGAT
# -----------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-fgat_fv3lm
                  MPI      6
                  ARGS     testinput/hyb-fgat_fv3lm.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_geos )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-fgat_gfs_aero
                  MPI      6
                  ARGS     testinput/hyb-fgat_gfs_aero.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero )

# 4DEnVar
# -------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_4denvar
                  MPI      18
                  ARGS     testinput/4denvar.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

# Hybrid 4DVar
# ------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-4dvar_pseudo-geos
                  MPI      6
                  ARGS     testinput/hyb-4dvar_pseudo-geos.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_geos )

# FSOI
# ----
ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_fsoi_forward
                  MPI      6
                  ARGS     testinput/hyb-3dvar_fsoi_forward.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_hyb-3dvar_fsoi_backward
                  MPI      6
                  ARGS     testinput/hyb-3dvar_fsoi_backward.yaml
                  COMMAND  fv3jedi_var.x
                  TEST_DEPENDS fv3jedi_test_tier1_hyb-3dvar_fsoi_forward )

# Diff states
# -----------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_diffstates_gfs
                  MPI      6
                  ARGS     testinput/diffstates_gfs.yaml
                  COMMAND  fv3jedi_diffstates.x
                  TEST_DEPENDS fv3jedi_test_tier1_hyb-3dvar )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_diffstates_geos
                  MPI      6
                  ARGS     testinput/diffstates_geos.yaml
                  COMMAND  fv3jedi_diffstates.x
                  TEST_DEPENDS fv3jedi_test_tier1_hyb-4dvar_pseudo-geos )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_diffstates_lam_cmaq
                  MPI      6
                  ARGS     testinput/diffstates_lam_cmaq.yaml
                  COMMAND  fv3jedi_diffstates.x
                  TEST_DEPENDS fv3jedi_test_tier1_3dvar_lam_cmaq )

# Add increment
# -------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_addincrement_gfs
                  MPI      6
                  ARGS     testinput/addincrement_gfs.yaml
                  COMMAND  fv3jedi_addincrement.x
                  TEST_DEPENDS fv3jedi_test_tier1_diffstates_gfs )

if (BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_addincrement_geos
                    MPI      12
                    ARGS     testinput/addincrement_geos.yaml
                    COMMAND  fv3jedi_addincrement.x
                    TEST_DEPENDS fv3jedi_test_tier1_diffstates_geos )
endif()

# LETKF / LGETKF
# --------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_letkf
                  MPI      6
                  ARGS     testinput/letkf.yaml
                  COMMAND  fv3jedi_letkf.x )

ecbuild_add_test( TARGET   fv3jedi_test_tier1_lgetkf
                  MPI      6
                  ARGS     testinput/lgetkf.yaml
                  COMMAND  fv3jedi_letkf.x )

# Land DA
# --------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_letkf_snow
                  MPI      6
                  ARGS     testinput/letkf_snow.yaml
                  COMMAND  fv3jedi_letkf.x )

# Obs Localizaton Test
# --------------------
ecbuild_add_test( TARGET   fv3jedi_test_tier1_obslocalization
                  MPI      6
                  ARGS     testinput/obslocalization.yaml
                  COMMAND  test_fv3jedi_obslocalization.x )

# EnsHofx tests
# -------------
if (BUILD_LARGE_TESTS)
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_enshofx_fv3lm
                    MPI      12
                    ARGS     testinput/enshofx_fv3lm.yaml
                    COMMAND  fv3jedi_enshofx.x )
endif()

# Ensemble data assimilation tests
# --------------------------------
if (BUILD_LARGE_TESTS)
  # 3DVar
  ecbuild_add_test( TARGET   fv3jedi_test_tier1_eda_3dvar
                    MPI      12
                    ARGS     testinput/eda_3dvar.yaml
                    COMMAND  fv3jedi_eda.x
                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )
endif()

ecbuild_add_test( TARGET   fv3jedi_test_tier1_eda_block_3dvar
                  MPI      12
                  ARGS     testinput/eda_block_3dvar.yaml
                  COMMAND  fv3jedi_eda.x
                  TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

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

if (FV3JEDI_TEST_TIER GREATER_EQUAL 2 AND BUILD_LARGE_TESTS)

  # Full BUMP static covariance setup
  # ----------------------------------

  # Loop over dates
  foreach(DATE ${staticb_dates})
    string(SUBSTRING ${DATE} 0 4 YYYY)
    string(SUBSTRING ${DATE} 4 2 MM)
    string(SUBSTRING ${DATE} 6 2 DD)
    string(SUBSTRING ${DATE} 8 2 HH)
    message(STATUS "Creating static B yaml files for cycle: ${YYYY}-${MM}-${DD}:${HH}")

    # Apply sed on template yaml file
    execute_process( COMMAND     sed "-e s/_YYYY_/${YYYY}/g;s/_MM_/${MM}/g;s/_DD_/${DD}/g;s/_HH_/${HH}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/staticb_ens-to-psichi_DATE_gfs.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/staticb_ens-to-psichi_${DATE}_gfs.yaml )
    execute_process( COMMAND     sed "-e s/_YYYY_/${YYYY}/g;s/_MM_/${MM}/g;s/_DD_/${DD}/g;s/_HH_/${HH}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/staticb_vbal_DATE_gfs.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/staticb_vbal_${DATE}_gfs.yaml )
    execute_process( COMMAND     sed "-e s/_YYYY_/${YYYY}/g;s/_MM_/${MM}/g;s/_DD_/${DD}/g;s/_HH_/${HH}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/staticb_unbal_DATE_gfs.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/staticb_unbal_${DATE}_gfs.yaml )
    execute_process( COMMAND     sed "-e s/_YYYY_/${YYYY}/g;s/_MM_/${MM}/g;s/_DD_/${DD}/g;s/_HH_/${HH}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/staticb_var-mom_DATE_aero.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/staticb_var-mom_${DATE}_aero.yaml )
    execute_process( COMMAND     sed "-e s/_YYYY_/${YYYY}/g;s/_MM_/${MM}/g;s/_DD_/${DD}/g;s/_HH_/${HH}/g"
                     INPUT_FILE  ${CMAKE_CURRENT_SOURCE_DIR}/testinput/staticb_var-mom_DATE_gfs.yaml
                     OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/testinput/staticb_var-mom_${DATE}_gfs.yaml )

    # Change of variable from u/v to psi/chi for the ensemble members
    ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_ens-to-psichi_${DATE}_gfs
                      MPI          12
                      ARGS         testinput/staticb_ens-to-psichi_${DATE}_gfs.yaml
                      COMMAND      fv3jedi_convertstate.x )

    # Compute vertical covariances and vertical balance using ensemble members
    ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_vbal_${DATE}_gfs
                      MPI          12
                      ARGS         testinput/staticb_vbal_${DATE}_gfs.yaml
                      COMMAND      fv3jedi_error_covariance_training.x
                      TEST_DEPENDS fv3jedi_test_tier2_staticb_ens-to-psichi_${DATE}_gfs )

    # Unbalance ensemble members
    ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_unbal_${DATE}_gfs
                      MPI          12
                      ARGS         testinput/staticb_unbal_${DATE}_gfs.yaml
                      COMMAND      fv3jedi_error_covariance_training.x
                      TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_${DATE}_gfs )

    # Compute variance and moments using unbalanced ensemble members
    ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_var-mom_${DATE}_aero
                      MPI          12
                      ARGS         testinput/staticb_var-mom_${DATE}_aero.yaml
                      COMMAND      fv3jedi_error_covariance_training.x )

    ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_var-mom_${DATE}_gfs
                      MPI          12
                      ARGS         testinput/staticb_var-mom_${DATE}_gfs.yaml
                      COMMAND      fv3jedi_error_covariance_training.x
                      TEST_DEPENDS fv3jedi_test_tier2_staticb_unbal_${DATE}_gfs )
  endforeach()

  # Compute vertical balance using previously pre-computed vertical covariances, and unbalance ensemble members
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_vbal_gfs
                    MPI          12
                    ARGS         testinput/staticb_vbal_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_2018041421_gfs
                                 fv3jedi_test_tier2_staticb_vbal_2018041500_gfs
                                 fv3jedi_test_tier2_staticb_vbal_2018041503_gfs )

  # Compute variance using pre-computed variance
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_var_aero
                    MPI          12
                    ARGS         testinput/staticb_var_aero.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var-mom_2018041421_aero
                                 fv3jedi_test_tier2_staticb_var-mom_2018041500_aero
                                 fv3jedi_test_tier2_staticb_var-mom_2018041503_aero )

  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_var_gfs
                    MPI          12
                    ARGS         testinput/staticb_var_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var-mom_2018041421_gfs
                                 fv3jedi_test_tier2_staticb_var-mom_2018041500_gfs
                                 fv3jedi_test_tier2_staticb_var-mom_2018041503_gfs )

  # Compute correlation length-scales using pre-computed moments
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_cor_aero
                    MPI          12
                    ARGS         testinput/staticb_cor_aero.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var-mom_2018041421_aero
                                 fv3jedi_test_tier2_staticb_var-mom_2018041500_aero
                                 fv3jedi_test_tier2_staticb_var-mom_2018041503_aero )

  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_cor_gfs
                    MPI          12
                    ARGS         testinput/staticb_cor_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var-mom_2018041421_gfs
                                 fv3jedi_test_tier2_staticb_var-mom_2018041500_gfs
                                 fv3jedi_test_tier2_staticb_var-mom_2018041503_gfs )

  # Compute NICAS operator using the length-scales
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_nicas_aero
                    MPI          12
                    ARGS         testinput/staticb_nicas_aero.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_cor_aero )

  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_nicas_gfs
                    MPI          12
                    ARGS         testinput/staticb_nicas_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_cor_gfs )

  # Compute psi/chi to u/v operator
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_psichitouv_gfs
                    MPI          12
                    ARGS         testinput/staticb_psichitouv_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x )

  # Split VBAL operator to get local files on 6 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_split_vbal_gfs
                    MPI          6
                    ARGS         testinput/staticb_split_vbal_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_gfs )

  # Split PSICHITOUV operator to get local files on 6 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_split_psichitouv_gfs
                    MPI          6
                    ARGS         testinput/staticb_split_psichitouv_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_psichitouv_gfs )

  # Split NICAS operator to get local files on 6 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_split_nicas_gfs
                    MPI          6
                    ARGS         testinput/staticb_split_nicas_gfs.yaml
                    COMMAND      fv3jedi_error_covariance_training.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_nicas_gfs )

  # Dirac test using the full BUMP static covariance, local files, 6 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_dirac_local_gfs_6pe
                    MPI          6
                    ARGS         testinput/staticb_dirac_local_gfs_6pe.yaml
                    COMMAND      fv3jedi_dirac.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_split_vbal_gfs
                                 fv3jedi_test_tier2_staticb_var_gfs
                                 fv3jedi_test_tier2_staticb_split_nicas_gfs
                                 fv3jedi_test_tier2_staticb_split_psichitouv_gfs )

  # Dirac test using the full BUMP static covariance, local files 12 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_dirac_local_aero_12pe
                    MPI          12
                    ARGS         testinput/staticb_dirac_local_aero_12pe.yaml
                    COMMAND      fv3jedi_dirac.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var_aero
                                 fv3jedi_test_tier2_staticb_nicas_aero )

  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_dirac_local_gfs_12pe
                    MPI          12
                    ARGS         testinput/staticb_dirac_local_gfs_12pe.yaml
                    COMMAND      fv3jedi_dirac.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_gfs
                                 fv3jedi_test_tier2_staticb_var_gfs
                                 fv3jedi_test_tier2_staticb_nicas_gfs
                                 fv3jedi_test_tier2_staticb_psichitouv_gfs )

  # Dirac test using the full BUMP static covariance, global files, 6 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_dirac_global_gfs_6pe
                    MPI          6
                    ARGS         testinput/staticb_dirac_global_gfs_6pe.yaml
                    COMMAND      fv3jedi_dirac.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_gfs
                                 fv3jedi_test_tier2_staticb_var_gfs
                                 fv3jedi_test_tier2_staticb_nicas_gfs )

  # Dirac test using the full BUMP static covariance, global files 12 PEs
  ecbuild_add_test( TARGET       fv3jedi_test_tier2_staticb_dirac_global_gfs_12pe
                    MPI          12
                    ARGS         testinput/staticb_dirac_global_gfs_12pe.yaml
                    COMMAND      fv3jedi_dirac.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_gfs
                                 fv3jedi_test_tier2_staticb_var_gfs
                                 fv3jedi_test_tier2_staticb_nicas_gfs )

  # Full observation h(x) for each model
  # ------------------------------------

  ecbuild_add_test( TARGET   fv3jedi_test_tier2_hofx_fv3-gfs
                    MPI      6
                    ARGS     testinput/hofx_fv3-gfs.yaml
                    COMMAND  fv3jedi_hofx.x )

  ecbuild_add_test( TARGET   fv3jedi_test_tier2_hofx_fv3-geos
                    MPI      6
                    ARGS     testinput/hofx_fv3-geos.yaml
                    COMMAND  fv3jedi_hofx.x )

  ecbuild_add_test( TARGET   fv3jedi_test_tier2_hofx_nomodel_lam_rrfs
                    MPI      12
                    ARGS     testinput/hofx_nomodel_lam_rrfs.yaml
                    COMMAND  fv3jedi_hofx_nomodel.x )


  # Full observation hyb-3DVar for each model
  # -----------------------------------------

  ecbuild_add_test( TARGET  fv3jedi_test_tier2_hyb-3dvar_gfs_aero_bump
                    MPI     12
                    ARGS    testinput/hyb-3dvar_gfs_aero_bump.yaml
                    COMMAND fv3jedi_var.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_var_aero
                                 fv3jedi_test_tier2_staticb_split_nicas_aero
                                 fv3jedi_test_tier1_bumpparameters_nicas_gfs_aero )

  ecbuild_add_test( TARGET   fv3jedi_test_tier2_hyb-3dvar_gfs
                    MPI      12
                    ARGS     testinput/hyb-3dvar_gfs.yaml
                    COMMAND  fv3jedi_var.x
                    TEST_DEPENDS fv3jedi_test_tier2_staticb_vbal_gfs
                                 fv3jedi_test_tier2_staticb_var_gfs
                                 fv3jedi_test_tier2_staticb_split_nicas_gfs
                                 fv3jedi_test_tier1_bumpparameters_nicas_gfs )

  ecbuild_add_test( TARGET   fv3jedi_test_tier2_hyb-3dvar_geos
                    MPI      6
                    ARGS     testinput/hyb-3dvar_geos.yaml
                    COMMAND  fv3jedi_var.x
                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_geos )

  # EDA 4DVar
  # ---------
  ecbuild_add_test( TARGET   fv3jedi_test_tier2_eda_4dvar_fv3lm
                    MPI      12
                    ARGS     testinput/eda_4dvar_fv3lm.yaml
                    COMMAND  fv3jedi_eda.x
                    TEST_DEPENDS fv3jedi_test_tier1_bumpparameters_nicas_gfs )

endif()
