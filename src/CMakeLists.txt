list( APPEND fv3jedi_src_files
ErrorCovariance/ErrorCovarianceFV3JEDI.cc
ErrorCovariance/ErrorCovarianceFV3JEDIFortran.h
ErrorCovariance/ErrorCovarianceFV3JEDI.h
ErrorCovariance/fv3jedi_covariance_interface.F90
ErrorCovariance/fv3jedi_covariance_mod.F90
Geometry/fv3jedi_geom_mod.F90
Geometry/GeometryFV3JEDI.cc
Geometry/GeometryFV3JEDIFortran.h
Geometry/GeometryFV3JEDI.h
GetValues/fv3jedi_getvaltraj_mod.F90
GetValues/fv3jedi_getvalues_mod.F90
GetValues/GetValuesTrajFV3JEDI.cc
GetValues/GetValuesTrajFV3JEDIFortran.h
GetValues/GetValuesTrajFV3JEDI.h
Increment/fv3jedi_increment_interface.F90
Increment/fv3jedi_increment_io_mod.F90
Increment/fv3jedi_increment_mod.F90
Increment/fv3jedi_increment_utils_mod.F90
Increment/IncrementFV3JEDI.cc
Increment/IncrementFV3JEDIFortran.h
Increment/IncrementFV3JEDI.h
Localization/instantiateLocalizationFactory.h
Localization/LocalizationMatrixFV3JEDI.cc
Localization/LocalizationMatrixFV3JEDIFortran.h
Localization/LocalizationMatrixFV3JEDI.h
Model/fv3jedi_model_interface_mod.F90
Model/fv3jedi_tlm_interface_mod.F90
Model/fv3jedi_traj_interface_mod.F90
Model/fv3/fv3jedi_model_mod.F90
Model/fv3jedi_tlm_mod.F90
Model/fv3jedi_traj_mod.F90
Model/ModelBiasCovarianceFV3JEDI.h
Model/ModelBiasFV3JEDI.h
Model/ModelBiasIncrementFV3JEDI.h
Model/ModelFV3JEDI.cc
Model/ModelFV3JEDIFortran.h
Model/ModelFV3JEDI.h
Model/TlmFV3JEDI.cc
Model/TlmFV3JEDIFortran.h
Model/TlmFV3JEDI.h
Model/TlmIdFV3JEDI.cc
Model/TlmIdFV3JEDI.h
Run/fv3jedi_setup.f90
Run/RunFV3JEDI.cc
Run/RunFV3JEDIFortran.h
Run/RunFV3JEDI.h
State/fv3jedi_state_interface.F90
State/fv3jedi_state_io_mod.F90
State/fv3jedi_state_mod.F90
State/fv3jedi_state_utils_mod.F90
State/StateFV3JEDI.cc
State/StateFV3JEDIFortran.h
State/StateFV3JEDI.h
Utilities/Fortran.h
Utilities/fv3jedi_constants.f90
Utilities/fv3jedi_kinds_mod.F90
Utilities/FV3JEDITraits.h
Utilities/fv3jedi_vars_mod.F90
Utilities/random_vectors_mod.F90
Utilities/UtilitiesFV3JEDI.cc
Utilities/UtilitiesFV3JEDI.h
VariableChange/fv3jedi_varcha_c2m_interface.F90
VariableChange/fv3jedi_varcha_c2m_mod.F90
VariableChange/height_variables_mod.F90
VariableChange/instantiateFV3JEDIVarChangeFactory.h
VariableChange/moisture_variables_mod.F90
VariableChange/pressure_variables_mod.F90
VariableChange/surface_variables_mod.F90
VariableChange/temperature_variables_mod.F90
VariableChange/VarChaC2MFV3JEDI.cc
VariableChange/VarChaC2MFV3JEDIFortran.h
VariableChange/VarChaC2MFV3JEDI.h
VariableChange/wind_variables_mod.F90
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

