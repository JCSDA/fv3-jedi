# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND fv3jedi_src_files
  ErrorCovariance/ErrorCovariance.cc
  ErrorCovariance/ErrorCovariance.h
  FieldMetadata/FieldsMetadata.cc
  FieldMetadata/FieldsMetadata.h
  FieldMetadata/FieldsMetadata.interface.cc
  FieldMetadata/FieldsMetadata.interface.h
  FieldMetadata/fields_metadata_mod.f90
  Fields/fv3jedi_field_mod.f90
  Fields/fv3jedi_fields_mod.f90
  Geometry/Geometry.cc
  Geometry/Geometry.h
  Geometry/Geometry.interface.h
  Geometry/fv3jedi_geom_interface_mod.F90
  Geometry/fv3jedi_geom_mod.f90
  GeometryIterator/GeometryIterator.cc
  GeometryIterator/GeometryIterator.h
  GeometryIterator/GeometryIterator.interface.h
  GeometryIterator/fv3jedi_geom_iter_interface.F90
  GeometryIterator/fv3jedi_geom_iter_mod.F90
  GetValues/fv3jedi_getvalues_mod.f90
  GetValues/fv3jedi_lineargetvalues_mod.f90
  GetValues/GetValues.cc
  GetValues/GetValues.h
  GetValues/GetValues.interface.F90
  GetValues/GetValues.interface.h
  GetValues/LinearGetValues.cc
  GetValues/LinearGetValues.h
  GetValues/LinearGetValues.interface.F90
  GetValues/LinearGetValues.interface.h
  Increment/Increment.cc
  Increment/Increment.h
  Increment/Increment.interface.h
  Increment/fv3jedi_increment_interface_mod.F90
  Increment/fv3jedi_increment_mod.F90
  IO/fv3jedi_io_geos_mod.f90
  IO/fv3jedi_io_gfs_mod.f90
  IO/fv3jedi_io_latlon_mod.f90
  IO/fv3jedi_io_utils_mod.f90
  Model/fv3lm/fv3jedi_fv3lm_interface_mod.F90
  Model/fv3lm/fv3jedi_fv3lm_mod.f90
  Model/fv3lm/ModelFV3LM.cc
  Model/fv3lm/ModelFV3LM.h
  Model/fv3lm/ModelFV3LM.interface.h
  Model/pseudo/fv3jedi_pseudo_interface_mod.F90
  Model/pseudo/fv3jedi_pseudo_mod.f90
  Model/pseudo/ModelPseudo.cc
  Model/pseudo/ModelPseudo.h
  Model/pseudo/ModelPseudo.interface.h
  ModelBias/ModelBiasCovariance.h
  ModelBias/ModelBias.h
  ModelBias/ModelBiasIncrement.h
  State/State.cc
  State/State.h
  State/State.interface.h
  State/fv3jedi_state_interface_mod.F90
  State/fv3jedi_state_mod.F90
  Utilities/Traits.h
  Utilities/interface.h
  Utilities/Utilities.cc
  Utilities/Utilities.h
  Utilities/fv3jedi_communication_mod.f90
  Utilities/fv3jedi_constants_mod.f90
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
  VariableChanges/ColdStartWinds/fv3jedi_vc_coldstartwinds.interface.F90
  VariableChanges/ColdStartWinds/fv3jedi_vc_coldstartwinds_mod.f90
  VariableChanges/ColdStartWinds/VarChaColdStartWinds.cc
  VariableChanges/ColdStartWinds/VarChaColdStartWinds.h
  VariableChanges/ColdStartWinds/VarChaColdStartWinds.interface.h
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
  VariableChanges/Model2GeoVaLs/fv3jedi_lvc_model2geovals.interface.F90
  VariableChanges/Model2GeoVaLs/fv3jedi_lvc_model2geovals_mod.f90
  VariableChanges/Model2GeoVaLs/fv3jedi_vc_model2geovals.interface.F90
  VariableChanges/Model2GeoVaLs/fv3jedi_vc_model2geovals_mod.f90
  VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.cc
  VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h
  VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.interface.h
  VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.cc
  VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h
  VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.interface.h
  VariableChanges/NMCBalance/LinVarChaNMCBal.cc
  VariableChanges/NMCBalance/LinVarChaNMCBal.h
  VariableChanges/NMCBalance/LinVarChaNMCBal.interface.h
  VariableChanges/NMCBalance/fv3jedi_linvarcha_nmcbal.interface.F90
  VariableChanges/NMCBalance/fv3jedi_linvarcha_nmcbal_mod.f90
  VariableChanges/Utils/fv3jedi_fieldfail_mod.f90
  VariableChanges/Utils/height_variables_mod.f90
  VariableChanges/Utils/moisture_variables_mod.f90
  VariableChanges/Utils/pressure_variables_mod.f90
  VariableChanges/Utils/surface_variables_mod.f90
  VariableChanges/Utils/temperature_variables_mod.f90
  VariableChanges/Utils/wind_variables_mod.f90
  Tlm/fv3jedi_tlm_mod.f90
  Tlm/fv3jedi_traj_mod.f90
  Tlm/Tlm.cc
  Tlm/Tlm.h
  Tlm/TlmId.cc
  Tlm/TlmId.h
  Tlm/Tlm.interface.F90
  Tlm/Tlm.interface.h
  Tlm/Traj.interface.F90
  Tlm/Traj.interface.h
)

# Precision of fv3
IF(FV3_PRECISION MATCHES "DOUBLE" OR NOT FV3_PRECISION)
  list( APPEND fv3jedi_src_files Geometry/fv_prec_r8_mod.f90 )
ELSE()
  list( APPEND fv3jedi_src_files Geometry/fv_prec_r4_mod.f90 )
ENDIF()

# Source files when building with stand alone fv3 (GFDL_atmos_cubed_sphere)
list( APPEND fv3jedi_fv3_src_files
  Geometry/fv_init_control_mod.f90
  VariableChanges/VertRemap/fv3jedi_vc_vertremap.interface.F90
  VariableChanges/VertRemap/fv3jedi_vc_vertremap_mod.f90
  VariableChanges/VertRemap/VarChaVertRemap.cc
  VariableChanges/VertRemap/VarChaVertRemap.h
  VariableChanges/VertRemap/VarChaVertRemap.interface.h
)

# Source files when building with GEOS
list( APPEND fv3jedi_geos_src_files
  Geometry/fv_init_geos_mod.f90
  Model/geos/fv3jedi_geos_interface_mod.F90
  Model/geos/fv3jedi_geos_mod.F90
  Model/geos/ModelGEOS.cc
  Model/geos/ModelGEOS.h
  Model/geos/ModelGEOS.interface.h
)

# Source files when building with GFS
list( APPEND fv3jedi_ufs_src_files
  Geometry/fv_init_mod.f90
  Model/gfs/fv3jedi_gfs_interface_mod.F90
  Model/gfs/fv3jedi_gfs_mod.f90
  Model/gfs/ModelGFS.cc
  Model/gfs/ModelGFS.h
  Model/gfs/ModelGFS.interface.h
)

# Source files when building with UFS
#list( APPEND fv3jedi_ufs_src_files
#)

# Source files when building with UFS prototype
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

#Requirement libraries
set(FV3JEDI_LIB_DEP oops saber ufo femps atlas atlas_f fv3jedilm ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES})

# Optional openMP
if(OpenMP_Fortran_FOUND)
  list( APPEND FV3JEDI_LIB_DEP OpenMP::OpenMP_Fortran)
endif()

#Requirement sources
set(FV3JEDI_SRC_DEP ${fv3jedi_src_files} )

# Optional model - fv3
if( FV3_FORECAST_MODEL MATCHES "FV3CORE" )
  list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_fv3_src_files} )
  list( APPEND FV3JEDI_LIB_DEP fv3 )
endif()

# Optional model - geos
if( FV3_FORECAST_MODEL MATCHES "GEOS" )
  list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_geos_src_files} )
  list( APPEND FV3JEDI_LIB_DEP MPI::MPI_Fortran )
endif()

# Optional model - ufs
if( FV3_FORECAST_MODEL MATCHES "UFS" )
  list( APPEND FV3JEDI_SRC_DEP ${fv3jedi_ufs_src_files} )
endif()


#Add the library
ecbuild_add_library( TARGET          fv3jedi
                     SOURCES         ${FV3JEDI_SRC_DEP}
                     LIBS            ${FV3JEDI_LIB_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${FV3JEDI_LINKER_LANGUAGE}
                   )
