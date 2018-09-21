list( APPEND fv3jedi_src_files
ErrorCovarianceFV3JEDI.cc
ErrorCovarianceFV3JEDI.h
ErrorCovarianceFV3JEDIFortran.h
FV3JEDITraits.h
Fortran.h
GeometryFV3JEDI.cc
GeometryFV3JEDI.h
GeometryFV3JEDIFortran.h
GetValuesTrajFV3JEDI.cc
GetValuesTrajFV3JEDI.h
GetValuesTrajFV3JEDIFortran.h
IncrementFV3JEDI.cc
IncrementFV3JEDI.h
IncrementFV3JEDIFortran.h
instantiateFV3JEDIVarChangeFactory.h
LocalizationMatrixFV3JEDI.cc
LocalizationMatrixFV3JEDI.h
LocalizationMatrixFV3JEDIFortran.h
ModelBiasCovarianceFV3JEDI.h
ModelBiasFV3JEDI.h
ModelBiasIncrementFV3JEDI.h
ModelFV3JEDI.cc
ModelFV3JEDI.h
ModelFV3JEDIFortran.h
RunFV3JEDI.cc
RunFV3JEDI.h
RunFV3JEDIFortran.h
StateFV3JEDI.cc
StateFV3JEDI.h
StateFV3JEDIFortran.h
TlmFV3JEDI.cc
TlmFV3JEDI.h
TlmIdFV3JEDI.cc
TlmIdFV3JEDI.h
UtilitiesFV3JEDI.cc
UtilitiesFV3JEDI.h
VarChaC2MFV3JEDI.cc
VarChaC2MFV3JEDI.h
VarChaC2MFV3JEDIFortran.h
fv3jedi_constants.f90
fv3jedi_covariance_interface.F90
fv3jedi_covariance_mod.F90
fv3jedi_geom_mod.F90
fv3jedi_getvaltraj_mod.F90
fv3jedi_increment_interface.F90
fv3jedi_increment_io_mod.F90
fv3jedi_increment_mod.F90
fv3jedi_increment_utils_mod.F90
fv3jedi_kinds_mod.F90
fv3jedi_model_interface.F90
fv3jedi_model_mod.F90
fv3jedi_setup.f90
fv3jedi_state_interface.F90
fv3jedi_state_io_mod.F90
fv3jedi_state_mod.F90
fv3jedi_state_utils_mod.F90
fv3jedi_trajectories.F90
fv3jedi_varcha_c2m_interface.F90
fv3jedi_varcha_c2m_mod.F90
fv3jedi_vars_mod.F90
instantiateLocalizationFactory.h
random_vectors_mod.F90
VariableTransforms/wind_variables.F90
VariableTransforms/moisture_variables.F90
VariableTransforms/pressure_variables.F90
VariableTransforms/surface_variables.F90
VariableTransforms/temperture_variables.F90
VariableTransforms/height_variables.F90
)

set ( FMS_LIBRARIES $ENV{FMS_LIBRARIES} )
set ( FMS_INCLUDES $ENV{FMS_INCLUDES} )

ecbuild_add_library( TARGET   fv3jedi
                     SOURCES  ${fv3jedi_src_files}
                     LIBS     oops ufo ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${FV3_LIBRARIES} ${FMS_LIBRARIES}
                     INCLUDES ${FV3_INCLUDES} ${FMS_INCLUDES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

