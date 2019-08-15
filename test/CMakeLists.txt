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

# GFS AERO test files
# -------------------
list( APPEND fv3jedi_gfs_test_data Data/inputs/gfs_aero_c12/akbk.nc)
GFS_AERO_FILES_BKG(Data/inputs/gfs_aero_c12/ensmean 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem001 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem002 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem003 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem004 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem005 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem006 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem007 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem008 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem009 20180415.000000)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/ensmean)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem001)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem002)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem003)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem004)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem005)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem006)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem007)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem008)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem009)
LINK_FILES( "${fv3jedi_gfs_test_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# FV3 config files
# ----------------

list( APPEND fv3jedi_test_fv3files
Data/fv3files/field_table
Data/fv3files/fmsmpp.nml
Data/fv3files/input_gfs_c12.nml
Data/fv3files/input_gfs_c12_p12.nml
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
testinput/bumpparameters_cor_gfs_p12.yaml
testinput/bumpparameters_loc_geos.yaml
testinput/bumpparameters_loc_gfs.yaml
testinput/bumpparameters_cor_gfs_aero.yaml
testinput/bumpparameters_loc_gfs_aero.yaml
testinput/bumpparameters_loc_gfs_p12.yaml
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
testinput/hyb-3dvar_gfs_p12.yaml
testinput/hyb-3dvar_gfs.yaml
testinput/hyb-3dvar_gfs_aero.yaml
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
testoutput/4denvar_geos.ref
testoutput/4denvar_gfs.ref
testoutput/bumpparameters_cor_geos.ref
testoutput/bumpparameters_cor_gfs.ref
testoutput/bumpparameters_loc_geos.ref
testoutput/bumpparameters_loc_gfs.ref
testoutput/dirac_geos.ref
testoutput/dirac_gfs.ref
#testoutput/forecast_geos-fv3.ref
#testoutput/forecast_geos.ref
#testoutput/forecast_gfs-fv3.ref
#testoutput/forecast_gfs.ref
testoutput/hofx_geos.ref
testoutput/hofx_gfs.ref
testoutput/hyb-3dvar_geos.ref
testoutput/hyb-3dvar_gfs.ref
testoutput/hyb-4dvar_geos.ref
testoutput/hyb-4dvar_gfs.ref
testoutput/hyb-fgat_geos.ref
testoutput/hyb-fgat_gfs.ref
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
LINK_FILES( "${fv3jedi_testoutput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Resources from other respositores
# ---------------------------------

# CRTM coefficient files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/crtm/)
file(GLOB crtm_spc_coeffs "${crtm_SOURCE_DIR}/fix/SpcCoeff/Little_Endian/*bin")
file(GLOB crtm_tau_coeffs "${crtm_SOURCE_DIR}/fix/TauCoeff/ODAS/Little_Endian/*bin")
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
file(GLOB ioda_atm_obs "${ioda_SOURCE_DIR}/test/testinput/atmosphere/*nc4")
file(GLOB ioda_mar_obs "${ioda_SOURCE_DIR}/test/testinput/marine/*nc4")
LINK_FILES_DIR("${ioda_atm_obs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/obs/)
LINK_FILES_DIR("${ioda_mar_obs}" ${CMAKE_CURRENT_BINARY_DIR}/Data/obs/)


# Create tests
# ------------

#Make some output directories for test data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis)

#Tolernaces for compare script
set (ctol 0.0) #Max relative difference
set (idif 0)      #Max difference

ecbuild_add_resources( TARGET   fv3jedi_test_scripts
                       SOURCES_PACK
                       ${fv3jedi_test_input}
                     )

# Test executables

ecbuild_add_test( TARGET  test_fv3jedi_geometry_gfs
                  SOURCES executables/TestGeometry.cc
                  ARGS    testinput/geometry_gfs.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_geos
                  SOURCES executables/TestGeometry.cc
                  ARGS    testinput/geometry_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_gfs
                  SOURCES executables/TestState.cc
                  ARGS    testinput/state_gfs.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_geos
                  SOURCES executables/TestState.cc
                  ARGS    testinput/state_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_geos_aerosol
                  SOURCES executables/TestState.cc
                  ARGS    testinput/state_geosaero.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_gfs-fv3
                  SOURCES executables/TestModel.cc
                  ARGS    testinput/model_gfs-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_geos-fv3
                  SOURCES executables/TestModel.cc
                  ARGS    testinput/model_geos-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model_geos-pseudo
                  SOURCES executables/TestModel.cc
                  ARGS    testinput/model_geos-pseudo.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_gfs
                  SOURCES executables/TestLinearModel.cc
                  ARGS    testinput/linearmodel_gfs.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_geos
                  SOURCES executables/TestLinearModel.cc
                  ARGS    testinput/linearmodel_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelaux_gfs-fv3
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    testinput/model_gfs-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelaux_geos-fv3
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    testinput/model_geos-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_increment_gfs
                  SOURCES executables/TestIncrement.cc
                  ARGS    testinput/increment_gfs.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_increment_geos
                  SOURCES executables/TestIncrement.cc
                  ARGS    testinput/increment_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_gfs
                  SOURCES executables/TestLinVarCha.cc
                  ARGS    testinput/linvarcha_gfs.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_geos
                  SOURCES executables/TestLinVarCha.cc
                  ARGS    testinput/linvarcha_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_gfs
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    testinput/model_gfs-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_geos
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    testinput/model_geos-fv3.yaml
                  MPI     6
                  LIBS    fv3jedi )


# H(x)
# ----

ecbuild_add_test( TARGET  test_fv3jedi_hofx_gfs
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x
                  ARGS    testinput/hofx_gfs.yaml testoutput/hofx_gfs.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_gfs.run testoutput/hofx_gfs.ref 1.0e-3 ${idif}
                  TEST_DEPENDS test_fv3jedi_hofx_gfs )

ecbuild_add_test( TARGET  test_fv3jedi_hofx_geos
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x
                  ARGS    testinput/hofx_geos.yaml testoutput/hofx_geos.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_geos.run testoutput/hofx_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hofx_geos )

# B Matrix tests
# --------------

#Parameter generation
ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_cor_geos.yaml testoutput/bumpparameters_cor_geos.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_cor_geos.run testoutput/bumpparameters_cor_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_geos )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_cor_gfs.yaml testoutput/bumpparameters_cor_gfs.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_cor_gfs.run testoutput/bumpparameters_cor_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_geos.yaml testoutput/bumpparameters_loc_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_geos )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_loc_geos.run testoutput/bumpparameters_loc_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_gfs.yaml testoutput/bumpparameters_loc_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_gfs_aero.yaml testoutput/bumpparameters_loc_gfs_aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs_aero )


ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_loc_gfs.run testoutput/bumpparameters_loc_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

# Error covariance
ecbuild_add_test( TARGET       test_fv3jedi_errorcovariance
                  SOURCES      executables/TestErrorCovariance.cc
                  ARGS         testinput/errorcov.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs )

# Localization
ecbuild_add_test( TARGET       test_fv3jedi_localization
                  SOURCES      executables/TestLocalization.cc
                  ARGS         testinput/localization.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

# Diracs
ecbuild_add_test( TARGET       test_fv3jedi_dirac_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_dirac.x
                  ARGS         testinput/dirac_geos.yaml testoutput/dirac_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/dirac_geos.run testoutput/dirac_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_dirac_geos )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_dirac.x
                  ARGS         testinput/dirac_gfs.yaml testoutput/dirac_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/dirac_gfs.run testoutput/dirac_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_dirac_gfs )

# Data assimilation runs
# ----------------------

#3DVar
ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs.yaml testoutput/hyb-3dvar_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs_aero.yaml testoutput/hyb-3dvar_gfs-aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_aero )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs.run testoutput/hyb-3dvar_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_geos.yaml testoutput/hyb-3dvar_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_geos.run testoutput/hyb-3dvar_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_geos )

#Hybrid FGAT
ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-fgat_gfs.yaml testoutput/hyb-fgat_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-fgat_gfs.run testoutput/hyb-fgat_gfs.ref 1.0e-4 ${idif}
                  TEST_DEPENDS test_fv3jedi_hybfgat_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-fgat_geos.yaml testoutput/hyb-fgat_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-fgat_geos.run testoutput/hyb-fgat_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hybfgat_geos )

#4DEnVar
ecbuild_add_test( TARGET       test_fv3jedi_4denvar_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/4denvar_gfs.yaml testoutput/4denvar_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/4denvar_gfs.run testoutput/4denvar_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_4denvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/4denvar_geos.yaml testoutput/4denvar_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/4denvar_geos.run testoutput/4denvar_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_4denvar_geos )

#Hybrid 4DVar
ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-4dvar_gfs.yaml testoutput/hyb-4dvar_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-4dvar_gfs.run testoutput/hyb-4dvar_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb4dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-4dvar_geos.yaml testoutput/hyb-4dvar_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-4dvar_geos.run testoutput/hyb-4dvar_geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb4dvar_geos )

# Test wit 12 (instead of 6) processors

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_cor_gfs_p12.yaml testoutput/bumpparameters_cor_gfs_p12.run
                  MPI          12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_cor_gfs_p12.run testoutput/bumpparameters_cor_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_gfs_p12.yaml testoutput/bumpparameters_loc_gfs_p12.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_compare_p12
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_loc_gfs_p12.run testoutput/bumpparameters_loc_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs_p12.yaml testoutput/hyb-3dvar_gfs_p12.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs_p12.run testoutput/hyb-3dvar_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs_p12 )
