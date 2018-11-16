list( APPEND fv3jedi_src_files
ControlToModel/fv3jedi_varcha_c2m_interface_mod.F90
ControlToModel/fv3jedi_varcha_c2m_mod.f90
ControlToModel/VarChaC2MFV3JEDI.cc
ControlToModel/VarChaC2MFV3JEDIFortran.h
ControlToModel/VarChaC2MFV3JEDI.h
ErrorCovariance/ErrorCovarianceFV3JEDI.cc
ErrorCovariance/ErrorCovarianceFV3JEDI.h
ErrorCovariance/ErrorCovarianceFV3JEDIFortran.h
ErrorCovariance/fv3jedi_covariance_interface_mod.F90
ErrorCovariance/fv3jedi_covariance_mod.f90
Geometry/GeometryFV3JEDI.cc
Geometry/GeometryFV3JEDI.h
Geometry/GeometryFV3JEDIFortran.h
Geometry/fv3jedi_geom_interface_mod.F90
Geometry/fv3jedi_geom_mod.f90
GetValues/GetValuesTrajFV3JEDI.cc
GetValues/GetValuesTrajFV3JEDI.h
GetValues/GetValuesTrajFV3JEDIFortran.h
GetValues/fv3jedi_getvalues_mod.f90
GetValues/fv3jedi_getvalues_traj_mod.F90
Increment/IncrementFV3JEDI.cc
Increment/IncrementFV3JEDI.h
Increment/IncrementFV3JEDIFortran.h
Increment/fv3jedi_increment_interface_mod.f90
Increment/fv3jedi_increment_io_mod.f90
Increment/fv3jedi_increment_mod.f90
Increment/fv3jedi_increment_utils_mod.F90
Localization/LocalizationMatrixFV3JEDI.cc
Localization/LocalizationMatrixFV3JEDI.h
Localization/LocalizationMatrixFV3JEDIFortran.h
Localization/instantiateLocalizationFactory.h
#Model/instantiateFV3JEDIModelFactory.h
Model/fv3/ModelFV3FV3JEDI.h
Model/fv3/fv3jedi_fv3_interface_mod.F90
Model/fv3/fv3jedi_fv3_mod.f90
Model/fv3/ModelFV3FV3JEDI.cc
Model/fv3/ModelFV3FV3JEDIFortran.h
ModelBias/ModelBiasCovarianceFV3JEDI.h
ModelBias/ModelBiasFV3JEDI.h
ModelBias/ModelBiasIncrementFV3JEDI.h
Run/RunFV3JEDI.cc
Run/RunFV3JEDI.h
Run/RunFV3JEDIFortran.h
Run/fv3jedi_setup_mod.f90
State/StateFV3JEDI.cc
State/StateFV3JEDI.h
State/StateFV3JEDIFortran.h
State/fv3jedi_state_interface_mod.f90
State/fv3jedi_state_io_mod.f90
State/fv3jedi_state_mod.f90
State/fv3jedi_state_utils_mod.F90
Utilities/FV3JEDITraits.h
Utilities/Fortran.h
Utilities/UtilitiesFV3JEDI.cc
Utilities/UtilitiesFV3JEDI.h
Utilities/fv3jedi_constants_mod.f90
Utilities/fv3jedi_kinds_mod.f90
Utilities/fv3jedi_vars_mod.f90
Utilities/random_vectors_mod.f90
Utilities/instantiateFV3JEDIVarChangeFactory.h
Utilities/instantiateObsFilterFactory.h
VariableChange/height_variables_mod.f90
VariableChange/moisture_variables_mod.f90
VariableChange/pressure_variables_mod.f90
VariableChange/surface_variables_mod.f90
VariableChange/temperature_variables_mod.f90
VariableChange/wind_variables_mod.f90
Tlm/TlmFV3JEDI.cc
Tlm/TlmFV3JEDIFortran.h
Tlm/TlmFV3JEDI.h
Tlm/TlmIdFV3JEDI.cc
Tlm/TlmIdFV3JEDI.h
Tlm/fv3jedi_tlm_interface_mod.F90
Tlm/fv3jedi_tlm_mod.f90
Tlm/fv3jedi_traj_interface_mod.F90
Tlm/fv3jedi_traj_mod.f90
Tlm/TrajFV3JEDIFortran.h
)

set ( FMS_LIBRARIES $ENV{FMS_LIBRARIES} )
set ( FMS_INCLUDES $ENV{FMS_INCLUDES} )

ecbuild_add_library( TARGET   fv3jedi
                     SOURCES  ${fv3jedi_src_files}
                     LIBS     oops ufo ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${FV3JEDILM_LIBRARIES} ${FMS_LIBRARIES}
                     INCLUDES ${FV3JEDILM_INCLUDES} ${FMS_INCLUDES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

