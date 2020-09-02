# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Option not to run with large tests
# ----------------------------------
option(SKIP_LARGE_TESTS "Skip tests with more than 6 processors" OFF)

# Create the Data/ and Data/inputs/ directories
# ---------------------------------------------
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs)


# GEOS test files
# ---------------
list( APPEND fv3jedi_geos_test_data
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.mem001.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.mem002.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.mem003.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.mem004.nc4
Data/inputs/geos_c12/geos.aero.bkg.20180415_000000z.mem005.nc4
Data/inputs/geos_c12/geos.aero.bkg.ext.20180415_000000z.nc4
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
GFS_AERO_FILES_BKG(Data/inputs/gfs_aero_c12/bkg 20180414.210000)
GFS_AERO_FILES_BKG(Data/inputs/gfs_aero_c12/bkg 20180415.000000)
GFS_AERO_FILES_BKG(Data/inputs/gfs_aero_c12/bkg 20180415.030000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem001 20180414.210000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem001 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem001 20180415.030000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem002 20180414.210000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem002 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem002 20180415.030000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem003 20180414.210000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem003 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem003 20180415.030000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem004 20180414.210000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem004 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem004 20180415.030000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem005 20180414.210000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem005 20180415.000000)
GFS_AERO_FILES_ENS(Data/inputs/gfs_aero_c12/mem005 20180415.030000)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/bkg)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem001)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem002)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem003)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem004)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/gfs_aero_c12/mem005)
LINK_FILES( "${fv3jedi_gfs_test_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# FV3-LAM-CMAQ test files
LAM_CMAQ_FILES_BKG(Data/inputs/lam_cmaq/bkg 20190917.000000)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/lam_cmaq)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/lam_cmaq/bkg)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/lam_cmaq/INPUT)
LINK_FILES( "${fv3jedi_lam_test_data}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

list( APPEND fv3jedi_test_lamfiles
Data/inputs/lam_cmaq/INPUT/grid_spec.nc
Data/inputs/lam_cmaq/INPUT/C48_grid.tile7.nc
)
LINK_FILES( "${fv3jedi_test_lamfiles}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# FV3 config files
# ----------------

list( APPEND fv3jedi_test_fv3files
Data/fv3files/akbk64.nc4
Data/fv3files/akbk72.nc4
Data/fv3files/field_table
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

# NMC balance coefficient files
# -----------------------------

list( APPEND fv3jedi_test_nmcbalance
Data/inputs/nmcbalance/global_berror.l64y192.nc
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/inputs/nmcbalance)
LINK_FILES( "${fv3jedi_test_nmcbalance}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Test input scripts
# ------------------
list( APPEND fv3jedi_test_input
  testinput/3denvar_geos_aero.yaml
  testinput/3dvar_gfs_bal.yaml
  testinput/4denvar_pseudo-geos.yaml
  testinput/4denvar_pseudo-gfs.yaml
  testinput/addincrement_geos.yaml
  testinput/addincrement_gfs.yaml
  testinput/bumpparameters_cor_geos_p12.yaml
  testinput/bumpparameters_cor_geos.yaml
  testinput/bumpparameters_cor_gfs_aero.yaml
  testinput/bumpparameters_cor_gfs_p12.yaml
  testinput/bumpparameters_cor_gfs.yaml
  testinput/bumpparameters_loc_geos_aero.yaml
  testinput/bumpparameters_loc_geos_p12.yaml
  testinput/bumpparameters_loc_geos.yaml
  testinput/bumpparameters_loc_gfs_aero.yaml
  testinput/bumpparameters_loc_gfs_p12.yaml
  testinput/bumpparameters_loc_gfs.yaml
  testinput/convertstate_geos_bkg2rst.yaml
  testinput/convertstate_geos.yaml
  testinput/convertstate_gfs.yaml
  testinput/crtm_surface_for_geos.yaml
  testinput/diffstates_geos.yaml
  testinput/diffstates_gfs.yaml
  testinput/dirac_geos.yaml
  testinput/dirac_gfs_bal.yaml
  testinput/dirac_gfs_c2a.yaml
  testinput/dirac_gfs.yaml
  testinput/eda_3dvar_geos.yaml
  testinput/eda_3dvar_gfs_member2.yaml
  testinput/eda_3dvar_gfs.yaml
  testinput/eda_4dvar_fv3-geos.yaml
  testinput/eda_4dvar_fv3-gfs.yaml
  testinput/enshofx_fv3-geos.yaml
  testinput/ensvariance_geos.yaml
  testinput/ensvariance_gfs.yaml
  testinput/errorcov.yaml
  testinput/errorcov_id.yaml
  testinput/forecast_fv3-geos.yaml
  testinput/forecast_fv3-gfs.yaml
  testinput/forecast_geos.yaml
  testinput/forecast_gfs.yaml
  testinput/geometry_geos.yaml
  testinput/geometry_gfs.yaml
  testinput/geometry_lam_cmaq.yaml
  testinput/geometry_iterator_geos.yaml
  testinput/geometry_iterator_gfs.yaml
  testinput/geometry_iterator_lam_cmaq.yaml
  testinput/getvalues_geos.yaml
  testinput/getvalues_gfs.yaml
  testinput/hofx_fv3-geos.yaml
  testinput/hofx_fv3-gfs.yaml
  testinput/hofx_gfs_aero.yaml
  testinput/hofx_fv3-gfs_gnssro_ropp.yaml
  testinput/hofx_nomodel_geos.yaml
  testinput/hofx_nomodel_gfs.yaml
  testinput/hyb-3dvar_geos_p12.yaml
  testinput/hyb-3dvar_geos.yaml
  testinput/hyb-3dvar_gfs_0obs.yaml
  testinput/hyb-3dvar_gfs_1obs.yaml
  testinput/hyb-3dvar_gfs_aero.yaml
  testinput/hyb-3dvar_gfs_p12.yaml
  testinput/hyb-3dvar_gfs.yaml
  testinput/hyb-3dvar_gfs_gnssro_ropp.yaml
  testinput/hyb-4dvar_pseudo-geos.yaml
  testinput/hyb-4dvar_fv3-gfs.yaml
  testinput/hyb-fgat_fv3-geos.yaml
  testinput/hyb-fgat_fv3-gfs.yaml
  testinput/hyb-fgat_gfs_aero.yaml
  testinput/increment_geos.yaml
  testinput/increment_gfs.yaml
  testinput/letkf_geos.yaml
  testinput/letkf_gfs.yaml
  testinput/lgetkf_geos.yaml
  testinput/lgetkf_gfs.yaml
  testinput/lineargetvalues_geos.yaml
  testinput/lineargetvalues_gfs.yaml
  testinput/linearmodel_geos.yaml
  testinput/linearmodel_gfs.yaml
  testinput/linearmodel-physics_geos.yaml
  testinput/linvarcha_geos.yaml
  testinput/linvarcha_gfs.yaml
  testinput/localization.yaml
  testinput/model_fv3-gfs.yaml
  testinput/model_fv3-geos.yaml
  testinput/model_geos.yaml
  testinput/model_pseudo-geos.yaml
  testinput/model_pseudo-gfs.yaml
  testinput/state_geosaero_ext.yaml
  testinput/state_geosaero.yaml
  testinput/state_geos.yaml
  testinput/state_gfs.yaml
  testinput/state_lam_cmaq.yaml
  testinput/test_time.yaml
  testinput/varcha_geos.yaml
  testinput/varcha_gfs.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
LINK_FILES( "${fv3jedi_test_input}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


# Test output scripts
# -------------------
list( APPEND fv3jedi_testoutput
  testoutput/4denvar_pseudo-geos.ref
  testoutput/4denvar_pseudo-gfs.ref
  testoutput/bumpparameters_cor_geos.ref
  testoutput/bumpparameters_cor_gfs.ref
  testoutput/bumpparameters_loc_geos.ref
  testoutput/bumpparameters_loc_gfs_aero.ref
  testoutput/bumpparameters_loc_gfs.ref
  testoutput/CTestCostData.txt.awsgnu.test
  testoutput/CTestCostData.txt.awsintel.test
  testoutput/dirac_geos.ref
  testoutput/dirac_gfs_bal.ref
  testoutput/dirac_gfs_c2a.ref
  testoutput/dirac_gfs.ref
  testoutput/forecast_fv3-geos.ref
  testoutput/forecast_fv3-gfs.ref
  testoutput/hofx_fv3-geos.ref
  testoutput/hofx_fv3-gfs.ref
  testoutput/hofx_nomodel_geos.ref
  testoutput/hofx_nomodel_gfs.ref
  testoutput/hyb-3dvar_geos.ref
  testoutput/hyb-3dvar_gfs_0obs.ref
  testoutput/hyb-3dvar_gfs_1obs.ref
  testoutput/hyb-3dvar_gfs-aero.ref
  testoutput/hyb-3dvar_gfs.ref
  testoutput/hyb-4dvar_fv3-gfs.ref
  testoutput/hyb-4dvar_pseudo-geos.ref
  testoutput/hyb-fgat_fv3-geos.ref
  testoutput/hyb-fgat_fv3-gfs.ref
  testoutput/letkf_geos.ref
  testoutput/letkf_gfs.ref
  testoutput/lgetkf_geos.ref
  testoutput/lgetkf_gfs.ref
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
Data/fieldsets/aerosols_gfs.yaml
Data/fieldsets/aerosols_geos.yaml
Data/fieldsets/lam_cmaq.yaml
Data/fieldsets/dynamics.yaml
Data/fieldsets/dynamics_lam_cmaq.yaml
Data/fieldsets/ufo.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fieldsets)
LINK_FILES( "${fv3jedi_field_sets}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

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
list( APPEND fv3_ioda_test_data
      testinput_tier_1.tar.gz
    )

# Set URL for IODA test files
set(FV3_IODA_DOWNLOAD_BASE_URL https://jedi-test-files.s3.amazonaws.com)


# If local path to testfiles is defined don't download
if (DEFINED ENV{LOCAL_PATH_TESTFILES_IODA})
        set(LOCAL_PATH_TESTFILES_IODA "$ENV{LOCAL_PATH_TESTFILES_IODA}")
endif()

if( NOT DEFINED LOCAL_PATH_TESTFILES_IODA )

  message(STATUS "LOCAL_PATH_TESTFILES_IODA is not defined, download test files")
  if ( NOT DEFINED FV3_IODA_TESTFILES_BRANCH)

    # Get the current git branch in fv3-jedi
    execute_process(
      COMMAND git rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH_FV3
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "fv3-jedi is in branch " ${GIT_BRANCH_FV3})

  else()
    message(STATUS "Branch name provided by user")
    set(GIT_BRANCH_FV3 ${FV3_IODA_TESTFILES_BRANCH})
  endif()

  # Check whether the URL exists or not
  set( ECBUILD_DOWNLOAD_BASE_URL ${FV3_IODA_DOWNLOAD_BASE_URL}/ioda)

  ecbuild_check_urls(NAMES ${GIT_BRANCH_FV3}/${fv3_ioda_test_data}
                     RESULT  FV3_SPECIFIC_TEST_FILES)

  # Set distant directory
  if(${FV3_SPECIFIC_TEST_FILES} MATCHES 0)
    # Download and extract new test files (distant directory = git branch)
    set(DIRNAME ${GIT_BRANCH_FV3})
    message(STATUS "GIT_BRANCH_FV3 found, will download " ${GIT_BRANCH_FV3})
  else()
    # Download and extract develop test files (distant directory = develop)
    message(STATUS "GIT_BRANCH_FV3 not found on S3, will download develop")
    set(DIRNAME "develop")
  endif()

  message(STATUS "Test data will be downloaded from: " ${ECBUILD_DOWNLOAD_BASE_URL}/${DIRNAME})

  set(FV3_IODA_REP_NAME ioda)
  set(FV3_IODA_TESTFILES_NAME ${fv3_ioda_test_data})
  set(FV3_IODA_BRANCH_NAME ${DIRNAME})
  set(FV3_IODA_TESTFILES_PATH ${CMAKE_BINARY_DIR}/test_data/ioda/${DIRNAME})
  message(STATUS "Save data to " ${FV3_IODA_TESTFILES_PATH})
  file(MAKE_DIRECTORY ${FV3_IODA_TESTFILES_PATH})

  # Create download script for get_ioda_test_data test
  set ( FILENAME fv3_ioda_data_downloader.py)
  set ( SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME} )
  set ( DEST_FILE ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
  list( APPEND bin_fv3_ioda_test_download_scripts_deps ${DEST_FILE} )

  if( EXISTS "${SOURCE_FILE}.in" )
    configure_file( ${SOURCE_FILE}.in ${DEST_FILE} @ONLY )
  else()
    configure_file( ${SOURCE_FILE}    ${DEST_FILE} @ONLY )
  endif()

  add_custom_target( bin_fv3_ioda_test_download_scripts ALL
      COMMAND chmod +x ${bin_fv3_ioda_test_download_scripts_deps}
      DEPENDS ${bin_fv3_ioda_test_download_scripts_deps})

  ecbuild_add_test( TARGET    fv3_get_ioda_test_data
                    TYPE      SCRIPT
                    COMMAND    ${CMAKE_BINARY_DIR}/bin/${FILENAME}
                    ARGS      testoutput/download_test_fv3.log )
else()
  set(FV3_IODA_TESTFILES_PATH ${LOCAL_PATH_TESTFILES_IODA})
  message(STATUS "use LOCAL_PATH_TESTFILES_IODA : " ${FV3_IODA_TESTFILES_PATH})
endif()

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${FV3_IODA_TESTFILES_PATH}
                 ${CMAKE_CURRENT_BINARY_DIR}/Data/obs )

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

ecbuild_add_executable( TARGET  test_fv3jedi_linvarcha.x
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

ecbuild_add_executable( TARGET  test_fv3jedi_varcha.x
                        SOURCES mains/TestVarCha.cc
                        LIBS    fv3jedi )

#Make some output directories for test data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/gsirf)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx/eda/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/eda/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/geos)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/analysis/letkf/gfs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/increment)

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

# Basic tests
# -----------

ecbuild_add_test( TARGET  test_fv3jedi_geometry_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometry.x
                  ARGS    testinput/geometry_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometry.x
                  ARGS    testinput/geometry_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_lam_cmaq
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometry.x
                  ARGS    testinput/geometry_lam_cmaq.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_iterator_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometryiterator.x
                  ARGS    testinput/geometry_iterator_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_iterator_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometryiterator.x
                  ARGS    testinput/geometry_iterator_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_geometry_iterator_lam_cmaq
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_geometryiterator.x
                  ARGS    testinput/geometry_iterator_lam_cmaq.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_crtm_surface_for_geos
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_convertstate.x
                  ARGS    testinput/crtm_surface_for_geos.yaml
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_state.x
                  ARGS    testinput/state_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_state_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_state.x
                  ARGS    testinput/state_geos.yaml
                  MPI     6
                  LIBS    fv3jedi
                  TEST_DEPENDS test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET  test_fv3jedi_state_lam_cmaq
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_state.x
                  ARGS    testinput/state_lam_cmaq.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_increment_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_increment.x
                  ARGS    testinput/increment_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_increment_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_increment.x
                  ARGS    testinput/increment_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_getvalues_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_getvalues.x
                  ARGS    testinput/getvalues_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_getvalues_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_getvalues.x
                  ARGS    testinput/getvalues_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_lineargetvalues_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_lineargetvalues.x
                  ARGS    testinput/lineargetvalues_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_lineargetvalues_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_lineargetvalues.x
                  ARGS    testinput/lineargetvalues_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_convertstate_gfs
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_convertstate.x
                  ARGS    testinput/convertstate_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_convertstate_geos_bkg2rst
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_convertstate.x
                  ARGS    testinput/convertstate_geos_bkg2rst.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_convertstate_geos
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_convertstate.x
                  ARGS    testinput/convertstate_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_varcha_gfs
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_varcha.x
                  ARGS         testinput/varcha_gfs.yaml
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_convertstate_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_varcha_geos
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_varcha.x
                  ARGS         testinput/varcha_geos.yaml
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_convertstate_geos )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_linvarcha.x
                  ARGS    testinput/linvarcha_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_linvarcha_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_linvarcha.x
                  ARGS    testinput/linvarcha_geos.yaml
                  MPI     6 )

# Model tests
# -----------

ecbuild_add_test( TARGET  test_fv3jedi_model_pseudo-gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_model.x
                  ARGS    testinput/model_pseudo-gfs.yaml
                  MPI     6
                  TEST_DEPENDS test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET  test_fv3jedi_model_pseudo-geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_model.x
                  ARGS    testinput/model_pseudo-geos.yaml
                  MPI     6
                  TEST_DEPENDS test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET  test_fv3jedi_model_fv3-geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_model.x
                  ARGS    testinput/model_fv3-geos.yaml
                  MPI     6 )

if( FV3_FORECAST_MODEL MATCHES "UFS" AND NOT SKIP_LARGE_TESTS)
  ecbuild_add_test( TARGET  test_fv3jedi_model_gfs
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_model.x
                    ARGS    testinput/model_gfs.yaml
                    MPI     6 )
endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT SKIP_LARGE_TESTS)
  ecbuild_add_test( TARGET  test_fv3jedi_model_geos
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_model.x
                    ARGS    testinput/model_geos.yaml
                    MPI     24 )
endif()

#ecbuild_add_test( TARGET  test_fv3jedi_modelaux_gfs-fv3
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_modelauxcontrol.x
#                  ARGS    testinput/model_gfs-fv3.yaml
#                  MPI     6 )
#
#ecbuild_add_test( TARGET  test_fv3jedi_modelaux_geos-fv3
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_modelauxcontrol.x
#                  ARGS    testinput/model_geos-fv3.yaml
#                  MPI     6 )
#
#ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_gfs
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_modelauxcovariance.x
#                  ARGS    testinput/model_gfs-fv3.yaml
#                  MPI     6)
#
#ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance_geos
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_modelauxcovariance.x
#                  ARGS    testinput/model_geos-fv3.yaml
#                  MPI     6 )

# Linear model tests
# ------------------
ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_linearmodel.x
                  ARGS    testinput/linearmodel_gfs.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_linearmodel.x
                  ARGS    testinput/linearmodel_geos.yaml
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel-physics_geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_linearmodel.x
                  ARGS    testinput/linearmodel-physics_geos.yaml
                  MPI     6 )

# Forecast tests
# --------------
ecbuild_add_test( TARGET  test_fv3jedi_forecast_fv3-gfs
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_forecast.x
                  ARGS    testinput/forecast_fv3-gfs.yaml testoutput/forecast_fv3-gfs.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_forecast_fv3-gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/forecast_fv3-gfs.run testoutput/forecast_fv3-gfs.ref 1.0e-3 ${idif}
                  TEST_DEPENDS test_fv3jedi_forecast_fv3-gfs )

ecbuild_add_test( TARGET  test_fv3jedi_forecast_fv3-geos
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_forecast.x
                  ARGS    testinput/forecast_fv3-geos.yaml testoutput/forecast_fv3-geos.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_forecast_fv3-geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/forecast_fv3-geos.run testoutput/forecast_fv3-geos.ref 1.0e-3 ${idif}
                TEST_DEPENDS test_fv3jedi_forecast_fv3-geos )

if( FV3_FORECAST_MODEL MATCHES "UFS" AND NOT SKIP_LARGE_TESTS)
  ecbuild_add_test( TARGET  test_fv3jedi_forecast_gfs
                    COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_forecast.x
                    ARGS    testinput/forecast_gfs.yaml testoutput/forecast_gfs.run
                    MPI     6 )

  ecbuild_add_test( TARGET       test_fv3jedi_forecast_gfs_compare
                    TYPE         SCRIPT
                    COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                    ARGS         testoutput/forecast_gfs.run testoutput/forecast_gfs.ref 1.0e-3 ${idif}
                    TEST_DEPENDS test_fv3jedi_forecast_gfs )
endif()

if( FV3_FORECAST_MODEL MATCHES "GEOS" AND NOT SKIP_LARGE_TESTS)
  ecbuild_add_test( TARGET  test_fv3jedi_forecast_geos
                    COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_forecast.x
                    ARGS    testinput/forecast_geos.yaml testoutput/forecast_geos.run
                    MPI     24 )

  ecbuild_add_test( TARGET       test_fv3jedi_forecast_geos_compare
                    TYPE         SCRIPT
                    COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                    ARGS         testoutput/forecast_geos.run testoutput/forecast_geos.ref 1.0e-3 ${idif}
                    TEST_DEPENDS test_fv3jedi_forecast_geos )
endif()

# H(x)
# ----

ecbuild_add_test( TARGET  test_fv3jedi_hofx_fv3-gfs
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x
                  ARGS    testinput/hofx_fv3-gfs.yaml testoutput/hofx_fv3-gfs.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_fv3-gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_fv3-gfs.run testoutput/hofx_fv3-gfs.ref 1.0e-3 ${idif}
                  TEST_DEPENDS test_fv3jedi_hofx_fv3-gfs )

ecbuild_add_test( TARGET  test_fv3jedi_hofx_fv3-geos
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x
                  ARGS    testinput/hofx_fv3-geos.yaml testoutput/hofx_fv3-geos.run
                  MPI     6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_fv3-geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_fv3-geos.run testoutput/hofx_fv3-geos.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_hofx_fv3-geos )

# H(x) when ROPP is enabled
# -------------------------

if ( NOT BUNDLE_SKIP_ROPP-UFO )
list( APPEND ropp_hofx_testoutput
  testoutput/hofx_fv3-gfs_gnssro_ropp.ref
)
LINK_FILES( "${ropp_hofx_testoutput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

ecbuild_add_test( TARGET  test_fv3jedi_hofx_fv3-gfs_gnssro_ropp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x
                  ARGS    testinput/hofx_fv3-gfs_gnssro_ropp.yaml testoutput/hofx_fv3-gfs_gnssro_ropp.run
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_hofx_fv3-gfs_gnssro_ropp_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_fv3-gfs_gnssro_ropp.run  testoutput/hofx_fv3-gfs_gnssro_ropp.ref 1.0e-1 ${idif}
                  TEST_DEPENDS test_fv3jedi_hofx_fv3-gfs_gnssro_ropp )
endif()

# H(x) no model
# -------------

ecbuild_add_test( TARGET       test_fv3jedi_hofx_nomodel_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx_nomodel.x
                  ARGS         testinput/hofx_nomodel_gfs.yaml testoutput/hofx_nomodel_gfs.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_nomodel_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_nomodel_gfs.run testoutput/hofx_nomodel_gfs.ref 1.0e-3 ${idif}
                TEST_DEPENDS test_fv3jedi_hofx_nomodel_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_nomodel_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx_nomodel.x
                  ARGS         testinput/hofx_nomodel_geos.yaml testoutput/hofx_nomodel_geos.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_nomodel_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hofx_nomodel_geos.run testoutput/hofx_nomodel_geos.ref 1.0e-2 ${idif}
                TEST_DEPENDS test_fv3jedi_hofx_nomodel_geos )


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

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_loc_gfs.run testoutput/bumpparameters_loc_gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

# Error covariance
ecbuild_add_test( TARGET       test_fv3jedi_errorcovariance
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_errorcovariance.x
                  ARGS         testinput/errorcov.yaml
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_errorcovariance_id
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_errorcovariance.x
                  ARGS         testinput/errorcov_id.yaml
                  MPI          6 )

# Localization
ecbuild_add_test( TARGET       test_fv3jedi_localization
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_localization.x
                  ARGS         testinput/localization.yaml
                  MPI          6
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

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs_bal
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_dirac.x
                  ARGS         testinput/dirac_gfs_bal.yaml testoutput/dirac_gfs_bal.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs_bal_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/dirac_gfs_bal.run testoutput/dirac_gfs_bal.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_dirac_gfs_bal )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs_c2a
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_dirac.x
                  ARGS         testinput/dirac_gfs_c2a.yaml testoutput/dirac_gfs_c2a.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_dirac_gfs_c2a_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/dirac_gfs_c2a.run testoutput/dirac_gfs_c2a.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_dirac_gfs_c2a )

# Ensemble variance
ecbuild_add_test( TARGET       test_fv3jedi_ensvariance_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_ensvariance.x
                  ARGS         testinput/ensvariance_gfs.yaml testoutput/ensvariance_gfs.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_ensvariance_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_ensvariance.x
                  ARGS         testinput/ensvariance_geos.yaml testoutput/ensvariance_geos.run
                  MPI          6 )

# Data assimilation runs
# ----------------------

#3DVar
ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs.yaml testoutput/hyb-3dvar_gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs.run testoutput/hyb-3dvar_gfs.ref 1.0e-1 ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs )

# 3DVar with GSI balance operator
#ecbuild_add_test( TARGET       test_fv3jedi_3dvar_gfs_bal
#                  TYPE         EXE
#                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
#                  ARGS         testinput/3dvar_gfs_bal.yaml testoutput/3dvar_gfs_bal.run
#                  MPI          6
#                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

# 3DVar with zero obs is not working yet, so these 2 tests are are commented out.
# When you want to debug this case, all you need to do is uncomment these tests.
# The YAML and obs files are are already in place. The appropriate ref file for checking
# the output is also in place, but needs to be updated.
#ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_0obs
#                  TYPE         EXE
#                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
#                  ARGS         testinput/hyb-3dvar_gfs_0obs.yaml testoutput/hyb-3dvar_gfs_0obs.run
#                  MPI          6
#                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )
#
#ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_0obs_compare
#                  TYPE         SCRIPT
#                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
#                  ARGS         testoutput/hyb-3dvar_gfs_0obs.run testoutput/hyb-3dvar_gfs_0obs.ref 1.0e-1 ${idif}
#                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_1obs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs_1obs.yaml testoutput/hyb-3dvar_gfs_1obs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_1obs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs_1obs.run testoutput/hyb-3dvar_gfs_1obs.ref 1.0e-1 ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_geos.yaml testoutput/hyb-3dvar_geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_geos.run testoutput/hyb-3dvar_geos.ref 1.0e-1 1
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_geos )

# 3DVar when ROPP is enabled
# -------------------------
if ( NOT BUNDLE_SKIP_ROPP-UFO )
list( APPEND ropp_var_testoutput
  testoutput/hyb-3dvar_gfs_gnssro_ropp.ref
)
LINK_FILES( "${ropp_var_testoutput}" ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

ecbuild_add_test( TARGET  test_fv3jedi_hyb-3dvar_gfs_gnssro_ropp
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS    testinput/hyb-3dvar_gfs_gnssro_ropp.yaml testoutput/hyb-3dvar_gfs_gnssro_ropp.run
                  MPI     6 )

ecbuild_add_test( TARGET  test_fv3jedi_hyb-3dvar_gfs_gnssro_ropp_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs_gnssro_ropp.run testoutput/hyb-3dvar_gfs_gnssro_ropp.ref 1.0e-1 ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb-3dvar_gfs_gnssro_ropp )
endif()

# Diff states
ecbuild_add_test( TARGET       test_fv3jedi_diffstates_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_diffstates.x
                  ARGS         testinput/diffstates_gfs.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_diffstates_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_diffstates.x
                  ARGS         testinput/diffstates_geos.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_geos )

# Add increment
ecbuild_add_test( TARGET       test_fv3jedi_addincrement_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_addincrement.x
                  ARGS         testinput/addincrement_gfs.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs test_fv3jedi_diffstates_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_addincrement_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_addincrement.x
                  ARGS         testinput/addincrement_geos.yaml
                  MPI          6
                  LIBS         fv3jedi
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_geos test_fv3jedi_diffstates_geos )



#Hybrid FGAT
ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_fv3-gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-fgat_fv3-gfs.yaml testoutput/hyb-fgat_fv3-gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_fv3-gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-fgat_fv3-gfs.run testoutput/hyb-fgat_fv3-gfs.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_hybfgat_fv3-gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_fv3-geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-fgat_fv3-geos.yaml testoutput/hyb-fgat_fv3-geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_fv3-geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-fgat_fv3-geos.run testoutput/hyb-fgat_fv3-geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hybfgat_fv3-geos )

#4DEnVar
ecbuild_add_test( TARGET       test_fv3jedi_4denvar_pseudo-gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/4denvar_pseudo-gfs.yaml testoutput/4denvar_pseudo-gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_pseudo-gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/4denvar_pseudo-gfs.run testoutput/4denvar_pseudo-gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_4denvar_pseudo-gfs )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_pseudo-geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/4denvar_pseudo-geos.yaml testoutput/4denvar_pseudo-geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET       test_fv3jedi_4denvar_pseudo-geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/4denvar_pseudo-geos.run testoutput/4denvar_pseudo-geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_4denvar_pseudo-geos )

#Hybrid 4DVar
ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_fv3-gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-4dvar_fv3-gfs.yaml testoutput/hyb-4dvar_fv3-gfs.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_fv3-gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-4dvar_fv3-gfs.run testoutput/hyb-4dvar_fv3-gfs.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb4dvar_fv3-gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_pseudo-geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-4dvar_pseudo-geos.yaml testoutput/hyb-4dvar_pseudo-geos.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos test_fv3jedi_crtm_surface_for_geos )

ecbuild_add_test( TARGET       test_fv3jedi_hyb4dvar_pseudo-geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-4dvar_pseudo-geos.run testoutput/hyb-4dvar_pseudo-geos.ref ${ctol} ${idif}
                  TEST_DEPENDS test_fv3jedi_hyb4dvar_geos )

#LETKF
ecbuild_add_test( TARGET       test_fv3jedi_letkf_geos_run
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_letkf.x
                  ARGS         testinput/letkf_geos.yaml testoutput/letkf_geos.run
                  MPI          6)

ecbuild_add_test( TARGET       test_fv3jedi_letkf_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/letkf_geos.run testoutput/letkf_geos.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_letkf_geos_run )

ecbuild_add_test( TARGET       test_fv3jedi_letkf_gfs_run
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_letkf.x
                  ARGS         testinput/letkf_gfs.yaml testoutput/letkf_gfs.run
                  MPI          6)

ecbuild_add_test( TARGET       test_fv3jedi_letkf_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/letkf_gfs.run testoutput/letkf_gfs.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_letkf_gfs_run )

#LGETKF
ecbuild_add_test( TARGET       test_fv3jedi_lgetkf_gfs_run
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_letkf.x
                  ARGS         testinput/lgetkf_gfs.yaml testoutput/lgetkf_gfs.run
                  MPI          6)

ecbuild_add_test( TARGET       test_fv3jedi_lgetkf_gfs_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/lgetkf_gfs.run testoutput/lgetkf_gfs.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_lgetkf_gfs_run )

ecbuild_add_test( TARGET       test_fv3jedi_lgetkf_geos_run
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_letkf.x
                  ARGS         testinput/lgetkf_geos.yaml testoutput/lgetkf_geos.run
                  MPI          6)

ecbuild_add_test( TARGET       test_fv3jedi_lgetkf_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/lgetkf_geos.run testoutput/lgetkf_geos.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_lgetkf_geos_run )


# Test with 12 (instead of 6) processors

if (NOT SKIP_LARGE_TESTS)

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_geos_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_cor_geos_p12.yaml testoutput/bumpparameters_cor_geos_p12.run
                  MPI          12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_geos_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_cor_geos_p12.run testoutput/bumpparameters_cor_geos.run.ref 0.0 0
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_geos_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_geos_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_geos_p12.yaml testoutput/bumpparameters_loc_geos_p12.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_geos_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_geos_compare_p12
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_loc_geos_p12.run testoutput/bumpparameters_loc_geos.run.ref 0.0 0
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_geos_p12.yaml testoutput/hyb-3dvar_geos_p12.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos_p12 test_fv3jedi_crtm_surface_for_geos test_fv3jedi_hyb3dvar_geos)

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_geos_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_geos_p12.run testoutput/hyb-3dvar_geos.run.ref 0.0 0
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_geos_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_cor_gfs_p12.yaml testoutput/bumpparameters_cor_gfs_p12.run
                  MPI          12 )

ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_cor_gfs_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/bumpparameters_cor_gfs_p12.run testoutput/bumpparameters_cor_gfs.run.ref 0.0 0
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
                  ARGS         testoutput/bumpparameters_loc_gfs_p12.run testoutput/bumpparameters_loc_gfs.run.ref 0.0 0
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_p12 )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_p12
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs_p12.yaml testoutput/hyb-3dvar_gfs_p12.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_p12 test_fv3jedi_hyb3dvar_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_p12_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/hyb-3dvar_gfs_p12.run testoutput/hyb-3dvar_gfs.run.ref 0.0 0
                  TEST_DEPENDS test_fv3jedi_hyb3dvar_gfs_p12 )

endif()

# Aerosol tests

ecbuild_add_test( TARGET  test_fv3jedi_state_geos_aerosol
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_state.x
                  ARGS    testinput/state_geosaero.yaml
                  MPI     6 )


ecbuild_add_test( TARGET  test_fv3jedi_state_geos_aerosol_ext
                  COMMAND ${CMAKE_BINARY_DIR}/bin/test_fv3jedi_state.x
                  ARGS    testinput/state_geosaero_ext.yaml
                  MPI     6 )


ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_gfs_aero.yaml testoutput/bumpparameters_loc_gfs_aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_cor_gfs_aero )

#ecbuild_add_test( TARGET       test_fv3jedi_bumpparameters_loc_gfs_aero_compare
#                  TYPE         SCRIPT
#                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
#                  ARGS         testoutput/bumpparameters_loc_gfs_aero.run# testoutput/bumpparameters_loc_gfs_aero.run.ref 0.0 0
#                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_aero )

ecbuild_add_test( TARGET       test_fv3jedi_hyb3dvar_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-3dvar_gfs_aero.yaml testoutput/hyb-3dvar_gfs-aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_aero )

ecbuild_add_test( TARGET       test_fv3jedi_hofx_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx_nomodel.x
                  ARGS         testinput/hofx_gfs_aero.yaml testoutput/hofx_gfs_aero.run
                  MPI          6 )

ecbuild_add_test( TARGET       test_fv3jedi_hybfgat_gfs_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/hyb-fgat_gfs_aero.yaml testoutput/hyb-fgat_gfs_aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs_aero )


if( NOT ${BUNDLE_SKIP_GEOS-AERO} )

###
#files needed
###

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

ecbuild_add_test( TARGET  test_fv3jedi_bumpparameters_loc_geos_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_parameters.x
                  ARGS         testinput/bumpparameters_loc_geos_aero.yaml testoutput/bumpparameters_loc_geos_aero.run
                  MPI          6
                )

ecbuild_add_test( TARGET       test_fv3jedi_3denvar_geos_aero
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_var.x
                  ARGS         testinput/3denvar_geos_aero.yaml testoutput/3denvar_geos_aero.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_geos_aero
                 )

endif()


if (NOT SKIP_LARGE_TESTS)

# EDA tests

ecbuild_add_test( TARGET       test_fv3jedi_eda_3dvar_geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_eda.x
                  ARGS         testinput/eda_3dvar_geos.yaml testoutput/eda_3dvar_geos.run
                  MPI          24
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_eda_3dvar_gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_eda.x
                  ARGS         testinput/eda_3dvar_gfs.yaml testoutput/eda_3dvar_gfs.run
                  MPI          24
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_eda_3dvar_gfs_member2
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_eda.x
                  ARGS         testinput/eda_3dvar_gfs_member2.yaml testoutput/eda_3dvar_gfs_member2.run
                  MPI          6
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_eda_4dvar_fv3-geos
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_eda.x
                  ARGS         testinput/eda_4dvar_fv3-geos.yaml testoutput/eda_4dvar_fv3-geos.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

ecbuild_add_test( TARGET       test_fv3jedi_eda_4dvar_fv3-gfs
                  TYPE         EXE
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/fv3jedi_eda.x
                  ARGS         testinput/eda_4dvar_fv3-gfs.yaml testoutput/eda_4dvar_fv3-gfs.run
                  MPI          12
                  TEST_DEPENDS test_fv3jedi_bumpparameters_loc_gfs )

# EnsHofx tests

ecbuild_add_test( TARGET  test_fv3jedi_enshofx_fv3-geos
                  TYPE    EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/fv3jedi_enshofx.x
                  ARGS    testinput/enshofx_fv3-geos.yaml testoutput/enshofx_fv3-geos.run
                  MPI     12 )

ecbuild_add_test( TARGET       test_fv3jedi_enshofx_geos_compare
                  TYPE         SCRIPT
                  COMMAND      ${CMAKE_BINARY_DIR}/bin/compare.py
                  ARGS         testoutput/enshofx_fv3-geos.run testoutput/hofx_fv3-geos.ref 1.0e-2 ${idif}
                  TEST_DEPENDS test_fv3jedi_enshofx_fv3-geos )

endif()

####################################################################
# test time
####################################################################
if ( DEFINED TEST_TIME_REF_FV3 )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/testoutput/${TEST_TIME_REF_FV3}
          ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${TEST_TIME_REF_FV3} )

  ecbuild_add_test( TARGET test_fv3_time
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/run_time_test.py
                  ARGS testoutput/${TEST_TIME_REF_FV3}
                  ${CMAKE_BINARY_DIR}/fv3-jedi/Testing/Temporary/CTestCostData.txt
                  testinput/test_time.yaml)
else()
  message ( STATUS "TEST_TIME not built")
endif()
