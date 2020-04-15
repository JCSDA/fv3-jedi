# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND fv3jedi_src_files
ErrorCovariance/ErrorCovariance.cc
ErrorCovariance/ErrorCovariance.h
ErrorCovariance/ErrorCovariance.interface.h
ErrorCovariance/fv3jedi_covariance_interface_mod.F90
ErrorCovariance/fv3jedi_covariance_mod.f90
FieldMetadata/FieldsMetadata.cc
FieldMetadata/FieldsMetadata.h
FieldMetadata/FieldsMetadata.interface.cc
FieldMetadata/FieldsMetadata.interface.h
FieldMetadata/fields_metadata_mod.f90
Geometry/Geometry.cc
Geometry/Geometry.h
Geometry/Geometry.interface.h
Geometry/fv3jedi_geom_interface_mod.F90
Geometry/fv3jedi_geom_mod.f90
Geometry/fv3geom/fv3geom_constants_mod.f90
Geometry/fv3geom/fv3geom_create_grid.f90
Geometry/fv3geom/fv3geom_gnomonicgrids_mod.f90
Geometry/fv3geom/fv3geom_utils_mod.F90
GeometryIterator/GeometryIterator.cc
GeometryIterator/GeometryIterator.h
GeometryIterator/GeometryIterator.interface.h
GeometryIterator/fv3jedi_geom_iter_interface.F90
GeometryIterator/fv3jedi_geom_iter_mod.F90
GetValues/GetValuesTraj.cc
GetValues/GetValuesTraj.h
GetValues/GetValuesTraj.interface.h
GetValues/fv3jedi_getvalues_mod.f90
GetValues/fv3jedi_getvalues_traj_mod.F90
Increment/Increment.cc
Increment/Increment.h
Increment/Increment.interface.h
Increment/fv3jedi_increment_interface_mod.f90
Increment/fv3jedi_increment_mod.f90
Increment/fv3jedi_increment_utils_mod.F90
IO/fv3jedi_io_geos_mod.f90
IO/fv3jedi_io_gfs_mod.f90
IO/fv3jedi_io_latlon_mod.f90
Localization/LocalizationMatrix.cc
Localization/LocalizationMatrix.h
Localization/LocalizationMatrix.interface.h
Localization/instantiateLocalizationFactory.h
Model/fv3/fv3jedi_fv3_interface_mod.F90
Model/fv3/fv3jedi_fv3_mod.f90
Model/fv3/ModelFV3.cc
Model/fv3/ModelFV3.h
Model/fv3/ModelFV3.interface.h
Model/pseudo/fv3jedi_pseudo_interface_mod.F90
Model/pseudo/fv3jedi_pseudo_mod.f90
Model/pseudo/ModelPseudo.cc
Model/pseudo/ModelPseudo.h
Model/pseudo/ModelPseudo.interface.h
Model/traj/fv3jedi_traj_interface_mod.F90
Model/traj/fv3jedi_traj_mod.f90
Model/traj/ModelTraj.cc
Model/traj/ModelTraj.h
Model/traj/ModelTraj.interface.h
ModelBias/ModelBiasCovariance.h
ModelBias/ModelBias.h
ModelBias/ModelBiasIncrement.h
Run/Run.cc
Run/Run.h
Run/Run.interface.h
Run/fv3jedi_setup_mod.f90
State/State.cc
State/State.h
State/State.interface.h
State/fv3jedi_state_interface_mod.f90
State/fv3jedi_state_mod.f90
State/fv3jedi_state_utils_mod.F90
Utilities/Traits.h
Utilities/interface.h
Utilities/Utilities.cc
Utilities/Utilities.h
Utilities/fv3jedi_bump_mod.f90
Utilities/fv3jedi_communication_mod.f90
Utilities/fv3jedi_constants_mod.f90
Utilities/fv3jedi_field_mod.f90
Utilities/fv3jedi_interpolation_mod.f90
Utilities/fv3jedi_kinds_mod.f90
Utilities/fv3jedi_netcdf_utils_mod.F90
VariableChanges/Analysis2Model/LinVarChaA2M.cc
VariableChanges/Analysis2Model/LinVarChaA2M.h
VariableChanges/Analysis2Model/LinVarChaA2M.interface.h
VariableChanges/Analysis2Model/fv3jedi_linvarcha_a2m.interface.F90
VariableChanges/Analysis2Model/fv3jedi_linvarcha_a2m_mod.f90
VariableChanges/Analysis2Model/VarChaA2M.cc
VariableChanges/Analysis2Model/VarChaA2M.h
VariableChanges/Analysis2Model/VarChaA2M.interface.h
VariableChanges/Analysis2Model/fv3jedi_varcha_a2m.interface.F90
VariableChanges/Analysis2Model/fv3jedi_varcha_a2m_mod.f90
VariableChanges/Control2Analysis/LinVarChaC2A.cc
VariableChanges/Control2Analysis/LinVarChaC2A.h
VariableChanges/Control2Analysis/LinVarChaC2A.interface.h
VariableChanges/Control2Analysis/fv3jedi_linvarcha_c2a.interface.F90
VariableChanges/Control2Analysis/fv3jedi_linvarcha_c2a_mod.f90
VariableChanges/Control2Analysis/VarChaC2A.cc
VariableChanges/Control2Analysis/VarChaC2A.h
VariableChanges/Control2Analysis/VarChaC2A.interface.h
VariableChanges/Control2Analysis/fv3jedi_varcha_c2a.interface.F90
VariableChanges/Control2Analysis/fv3jedi_varcha_c2a_mod.f90
VariableChanges/GEOSRstToBkg/fv3jedi_vc_geosrst2bkg.interface.F90
VariableChanges/GEOSRstToBkg/fv3jedi_vc_geosrst2bkg_mod.f90
VariableChanges/GEOSRstToBkg/VarChaGeosRst2Bkg.cc
VariableChanges/GEOSRstToBkg/VarChaGeosRst2Bkg.h
VariableChanges/GEOSRstToBkg/VarChaGeosRst2Bkg.interface.h
VariableChanges/NMCBalance/LinVarChaNMCBal.cc
VariableChanges/NMCBalance/LinVarChaNMCBal.h
VariableChanges/NMCBalance/LinVarChaNMCBal.interface.h
VariableChanges/NMCBalance/fv3jedi_linvarcha_nmcbal.interface.F90
VariableChanges/NMCBalance/fv3jedi_linvarcha_nmcbal_mod.f90
VariableChanges/Utils/height_variables_mod.f90
VariableChanges/Utils/moisture_variables_mod.f90
VariableChanges/Utils/pressure_variables_mod.f90
VariableChanges/Utils/surface_variables_mod.f90
VariableChanges/Utils/temperature_variables_mod.f90
VariableChanges/Utils/wind_variables_mod.f90
VariableChanges/instantiateVarChangeFactories.h
Tlm/Tlm.cc
Tlm/Tlm.interface.h
Tlm/Tlm.h
Tlm/TlmId.cc
Tlm/TlmId.h
Tlm/fv3jedi_tlm_interface_mod.F90
Tlm/fv3jedi_tlm_mod.f90
)

list( APPEND fv3jedi_geos_src_files
Model/geos/fv3jedi_geos_interface_mod.F90
Model/geos/fv3jedi_geos_mod.F90
Model/geos/ModelGEOS.cc
Model/geos/ModelGEOS.h
Model/geos/ModelGEOS.interface.h
)

list( APPEND fv3jedi_gfs_src_files
Model/gfs/fv3jedi_gfs_interface_mod.F90
Model/gfs/fv3jedi_gfs_mod.f90
Model/gfs/ModelGFS.cc
Model/gfs/ModelGFS.h
Model/gfs/ModelGFS.interface.h
)

list( APPEND fv3jedi_nuopc_src_files
Model/nuopc/fv3jedi_nuopc_interface_mod.F90
Model/nuopc/fv3jedi_nuopc_mod.F90
Model/nuopc/ModelNUOPC.cc
Model/nuopc/ModelNUOPC.h
Model/nuopc/ModelNUOPC.interface.h
Model/nuopc/esm.F90
Model/nuopc/atm.F90
Model/nuopc/ocn.F90
Model/nuopc/med.F90
)

#Append includes and libraries with particular model
list (APPEND FV3JEDI_LIB_DEP oops saber ufo femps atlas atlas_f ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${FMS_LIBRARIES} ${FV3JEDILM_LIBRARIES})

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
                     LINKER_LANGUAGE ${FV3JEDI_LINKER_LANGUAGE}
                   )
