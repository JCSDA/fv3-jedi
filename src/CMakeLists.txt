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
IO/fv3jedi_io_latlon_mod.f90
Localization/LocalizationMatrixFV3JEDI.cc
Localization/LocalizationMatrixFV3JEDI.h
Localization/LocalizationMatrixFV3JEDIFortran.h
Localization/instantiateLocalizationFactory.h
Model/fv3/fv3jedi_fv3_interface_mod.F90
Model/fv3/fv3jedi_fv3_mod.f90
Model/fv3/ModelFV3FV3JEDI.cc
Model/fv3/ModelFV3FV3JEDI.h
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

list( APPEND fv3jedi_geos_src_files
Model/geos/fv3jedi_geos_interface_mod.F90
Model/geos/fv3jedi_geos_mod.F90
Model/geos/ModelGEOSFV3JEDI.cc
Model/geos/ModelGEOSFV3JEDI.h
Model/geos/ModelGEOSFV3JEDIFortran.h
)

list( APPEND fv3jedi_gfs_src_files
Model/gfs/fv3jedi_gfs_interface_mod.F90
Model/gfs/fv3jedi_gfs_mod.f90
Model/gfs/ModelGFSFV3JEDI.cc
Model/gfs/ModelGFSFV3JEDI.h
Model/gfs/ModelGFSFV3JEDIFortran.h
)

list( APPEND fv3jedi_nuopc_src_files
Model/nuopc/fv3jedi_nuopc_interface_mod.F90
Model/nuopc/fv3jedi_nuopc_mod.F90
Model/nuopc/ModelNUOPCFV3JEDI.cc
Model/nuopc/ModelNUOPCFV3JEDI.h
Model/nuopc/ModelNUOPCFV3JEDIFortran.h
Model/nuopc/esm.F90
Model/nuopc/atm.F90
Model/nuopc/ocn.F90
Model/nuopc/med.F90
)

#Append includes and libraries with particular model 
list (APPEND FV3JEDI_LIB_DEP oops ufo ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${FMS_LIBRARIES} ${FV3JEDILM_LIBRARIES})

if( ${FV3BASEDMODEL_FOUND} )
    list( APPEND FV3JEDI_LIB_DEP ${FV3BASEDMODEL_LIBRARIES})
endif()

#Append src with particular model 
list (APPEND FV3JEDI_SRC_DEP ${fv3jedi_src_files} )

if( ${GEOS_FOUND} )
    list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_geos_src_files} )
    list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_nuopc_src_files} )
elseif( ${GFS_FOUND} )
    list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_gfs_src_files} )
endif()

#Add the library
ecbuild_add_library( TARGET          fv3jedi
                     SOURCES         ${FV3JEDI_SRC_DEP}
                     LIBS            ${FV3JEDI_LIB_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                   )
