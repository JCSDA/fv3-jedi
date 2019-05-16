# Create the Data/ and Data/inputs/ directories
# ---------------------------------------------
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs)


# GEOS test files
# ---------------
list( APPEND fv3jedi_geos_test_data
Data/inputs/geos_c12/akbk72.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.nc4
Data/inputs/geos_c12/geos.bkg.20180414_210000z.nc4
Data/inputs/geos_c12/geos.bkg.20180414_220000z.nc4
Data/inputs/geos_c12/geos.bkg.20180414_230000z.nc4
Data/inputs/geos_c12/geos.bkg.20180415_000000z.nc4
Data/inputs/geos_c12/geos.bkg.20180415_010000z.nc4
Data/inputs/geos_c12/geos.bkg.20180415_020000z.nc4
Data/inputs/geos_c12/geos.bkg.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem001.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem001.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem001.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem002.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem002.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem002.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem003.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem003.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem003.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem004.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem004.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem004.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem005.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem005.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem005.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem006.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem006.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem006.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem007.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem007.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem007.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem008.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem008.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem008.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem009.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem009.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem009.20180415_030000z.nc4
Data/inputs/geos_c12/geos.mem010.20180414_210000z.nc4
Data/inputs/geos_c12/geos.mem010.20180415_000000z.nc4
Data/inputs/geos_c12/geos.mem010.20180415_030000z.nc4
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/geos_c12)
LINK_FILES( "${fv3jedi_geos_test_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# GFS test files
# --------------
list( APPEND fv3jedi_gfs_test_data Data/inputs/gfs_c12/akbk.nc)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180414.210000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180414.220000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180414.230000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180415.000000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180415.010000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180415.020000)
GFS_FILES_BKG(Data/inputs/gfs_c12/bkg 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem001 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem001 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem001 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem002 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem002 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem002 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem003 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem003 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem003 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem004 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem004 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem004 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem005 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem005 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem005 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem006 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem006 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem006 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem007 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem007 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem007 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem008 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem008 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem008 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem009 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem009 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem009 20180415.030000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem010 20180414.210000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem010 20180415.000000)
GFS_FILES_ENS(Data/inputs/gfs_c12/mem010 20180415.030000)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/bkg)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem001)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem002)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem003)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem004)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem005)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem006)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem007)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem008)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem009)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_c12/mem010)
LINK_FILES( "${fv3jedi_gfs_test_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# FV3 config files
# ----------------

list( APPEND fv3jedi_test_fv3files
Data/fv3files/field_table
Data/fv3files/fmsmpp.nml
Data/fv3files/input_gfs_c12.nml
Data/fv3files/input_gfs_c24.nml
Data/fv3files/input_gfs_c48.nml
Data/fv3files/input_gfs_c96.nml
Data/fv3files/input_gfs_c192.nml
Data/fv3files/input_gfs_c384.nml
Data/fv3files/input_gfs_c768.nml
Data/fv3files/input_geos_c12.nml
Data/fv3files/input_geos_c24.nml
Data/fv3files/input_geos_c48.nml
Data/fv3files/input_geos_c90.nml
Data/fv3files/input_geos_c180.nml
Data/fv3files/input_geos_c360.nml
Data/fv3files/input_geos_c720.nml
Data/fv3files/inputpert_4dvar.nml
Data/fv3files/inputpert_fsoi.nml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fv3files)
LINK_FILES( "${fv3jedi_test_fv3files}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Test input scripts
# ------------------
list( APPEND fv3jedi_test_input
testinput/3dvar_gfs_gnssrobndropp1d.yaml
testinput/4denvar_geos.yaml
testinput/4denvar_gfs.yaml
testinput/bumpparameters_cor_geos.yaml
testinput/bumpparameters_cor_gfs.yaml
testinput/bumpparameters_loc_geos.yaml
testinput/bumpparameters_loc_gfs.yaml
testinput/dirac_geos.yaml
testinput/dirac_gfs.yaml
testinput/errorcov.yaml
testinput/forecast_geos-fv3.yaml
testinput/forecast_geos.yaml
testinput/forecast_gfs-fv3.yaml
testinput/forecast_gfs.yaml
testinput/geometry_geos.yaml
testinput/geometry_gfs.yaml
testinput/hofx_geos.yaml
testinput/hofx_gfs.yaml
testinput/hyb-3dvar_geos.yaml
testinput/hyb-3dvar_gfs.yaml
testinput/hyb-4dvar_geos.yaml
testinput/hyb-4dvar_gfs.yaml
testinput/hyb-fgat_geos.yaml
testinput/hyb-fgat_gfs.yaml
testinput/increment_geos.yaml
testinput/increment_gfs.yaml
testinput/linearmodel_geos.yaml
testinput/linearmodel_gfs.yaml
testinput/linearmodel_pseudogeos.yaml
testinput/linvarcha_geos.yaml
testinput/linvarcha_gfs.yaml
testinput/localization.yaml
testinput/model_geos-fv3.yaml
testinput/model_geos-pseudo.yaml
testinput/model_gfs-fv3.yaml
testinput/state_geosaero.yaml
testinput/state_geos.yaml
testinput/state_gfs.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
LINK_FILES( "${fv3jedi_test_input}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Test output scripts
# -------------------
list( APPEND fv3jedi_testoutput
testoutput/3dvar_geos.test
testoutput/3dvar_gfs.test
testoutput/4denvar_geos.test
testoutput/4denvar_gfs.test
testoutput/dirac_cor.test
testoutput/dirac_loc.test
testoutput/hofx_geos.test
testoutput/hofx_gfs.test
testoutput/hyb-4dvar_geos.test
testoutput/hyb-4dvar_gfs.test
testoutput/hyb-fgat_geos.test
testoutput/hyb-fgat_gfs.test
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
LINK_FILES( "${fv3jedi_testoutput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Resources from other respositores
# ---------------------------------

# CRTM coefficient files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
file(GLOB crtm_spc_coeffs "${crtm_SOURCE_DIR}/fix/SpcCoeff/Little_Endian/*bin")
file(GLOB crtm_taa_coeffs "${crtm_SOURCE_DIR}/fix/TauCoeff/ODAS/Little_Endian/*bin")
file(GLOB crtm_tap_coeffs "${crtm_SOURCE_DIR}/fix/TauCoeff/ODPS/Little_Endian/*bin")
file(GLOB crtm_cld_coeffs "${crtm_SOURCE_DIR}/fix/CloudCoeff/Little_Endian/*bin")
file(GLOB crtm_aer_coeffs "${crtm_SOURCE_DIR}/fix/AerosolCoeff/Little_Endian/*bin")
file(GLOB crtm_iri_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/IR_Ice/SEcategory/Little_Endian/*bin")
file(GLOB crtm_irl_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/IR_Land/SEcategory/Little_Endian/*bin")
file(GLOB crtm_irs_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/IR_Snow/SEcategory/Little_Endian/*bin")
file(GLOB crtm_irw_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/IR_Water/Little_Endian/*bin")
file(GLOB crtm_mww_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/MW_Water/Little_Endian/*bin")
file(GLOB crtm_vii_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/VIS_Ice/SEcategory/Little_Endian/*bin")
file(GLOB crtm_vil_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/VIS_Land/SEcategory/Little_Endian/*bin")
file(GLOB crtm_vis_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/VIS_Snow/SEcategory/Little_Endian/*bin")
file(GLOB crtm_viw_coeffs "${crtm_SOURCE_DIR}/fix/EmisCoeff/VIS_Water/SEcategory/Little_Endian/*bin")
LINK_FILES_DIR("${crtm_spc_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_tau_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_tap_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_cld_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_aer_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_iri_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_irl_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_irs_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_irw_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_mww_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_vii_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_vil_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_vis_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
LINK_FILES_DIR("${crtm_viw_coeffs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)

# IODA observation files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/obs/)
file(GLOB ioda_obs "${ioda_SOURCE_DIR}/test/testinput/atmosphere/*nc4")
LINK_FILES_DIR("${ioda_obs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/obs/)


# Create tests
# ------------

#Make some output directories for test data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis)


ecbuild_add_resources( TARGET   fv3jedi_test_scripts
                       SOURCES_PACK
                       ${fv3jedi_test_input}
                     )

# Test executables

ecbuild_add_test( TARGET  test_fv3jedi_geometry_gfs
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_geos
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_gfs
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_geos
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_geos_aerosol
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_geosaero.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_gfs-fv3
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_gfs-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_geos-fv3
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_geos-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_geos-pseudo
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_geos-pseudo.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_gfs
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_geos
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelaux_gfs-fv3
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model_gfs-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelaux_geos-fv3
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model_geos-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_increment_gfs
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_increment_geos
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_gfs
                  SOURCES executables/TestLinVarCha.cc
                  ARGS    "testinput/linvarcha_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_geos
                  SOURCES executables/TestLinVarCha.cc
                  ARGS    "testinput/linvarcha_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_gfs
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/model_gfs-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_geos
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/model_geos-fv3.yaml"
                  MPI     6
                  LIBS    fv3jedi )


# H(x)
# ----

ecbuild_add_test( TARGET test_fv3jedi_hofx_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediHofX.cc
                  ARGS "testinput/hofx_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_hofx_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediHofX.cc
                  ARGS "testinput/hofx_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )


# B Matrix tests
# --------------

#Parameter generation
ecbuild_add_test( TARGET test_fv3jedi_bumpparameters_cor_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediEstimateParams.cc
                  ARGS "testinput/bumpparameters_cor_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_bumpparameters_cor_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediEstimateParams.cc
                  ARGS "testinput/bumpparameters_cor_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_bumpparameters_loc_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediEstimateParams.cc
                  ARGS "testinput/bumpparameters_loc_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_bumpparameters_loc_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediEstimateParams.cc
                  ARGS "testinput/bumpparameters_loc_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

# Error covariance
ecbuild_add_test( TARGET  test_fv3jedi_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcov.yaml"
                  MPI     6
                  LIBS    fv3jedi )


# Localization
ecbuild_add_test( TARGET  test_fv3jedi_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/localization.yaml"
                  MPI     6
                  LIBS    fv3jedi )


# Diracs
ecbuild_add_test( TARGET test_fv3jedi_dirac_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediDirac.cc
                  ARGS "testinput/dirac_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_dirac_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediDirac.cc
                  ARGS "testinput/dirac_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )


# Data assimilation runs
# ----------------------

#3DVar
ecbuild_add_test( TARGET test_fv3jedi_hyb3dvar_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-3dvar_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_hyb3dvar_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-3dvar_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

#Hybrid FGAT
ecbuild_add_test( TARGET test_fv3jedi_hybfgat_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-fgat_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_hybfgat_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-fgat_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

#4DEnVar
ecbuild_add_test( TARGET test_fv3jedi_4denvar_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/4denvar_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_4denvar_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/4denvar_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )

#Hybrid 4DVar
ecbuild_add_test( TARGET test_fv3jedi_hyb4dvar_gfs
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-4dvar_gfs.yaml"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET test_fv3jedi_hyb4dvar_geos
                  TYPE EXE
                  SOURCES ../mains/fv3jediVar.cc
                  ARGS "testinput/hyb-4dvar_geos.yaml"
                  MPI     6
                  LIBS    fv3jedi )
