# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDI
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedi C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "fv3-jedi" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( FV3JEDI_LINKER_LANGUAGE CXX )

################################################################################
# Standard compiler flags
################################################################################

include( fv3jedi_extra_macros )

################################################################################
# Set default model to GFDL_atmos_cubed_sphere
################################################################################

if (NOT FV3_FORECAST_MODEL MATCHES GEOS AND NOT FV3_FORECAST_MODEL MATCHES UFS)
  set(FV3_FORECAST_MODEL "FV3CORE")
endif()

################################################################################
# Dependencies
################################################################################

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# atlas
ecbuild_use_package( PROJECT atlas VERSION 0.18.1 REQUIRED )
include_directories( ${ATLAS_INCLUDE_DIRS} )
add_definitions(-DATLASIFIED=1)

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# saber
ecbuild_use_package( PROJECT saber VERSION 0.0.1 REQUIRED )
include_directories( ${SABER_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo REQUIRED )
include_directories( ${UFO_INCLUDE_DIRS} )

# femps
ecbuild_use_package( PROJECT femps REQUIRED )
include_directories( ${FEMPS_INCLUDE_DIRS} )

#OpenMP
find_package( OpenMP COMPONENTS Fortran )

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# Set default forecast model
# --------------------------
if (NOT FV3_FORECAST_MODEL MATCHES GEOS AND NOT FV3_FORECAST_MODEL MATCHES UFS)
  set(FV3_FORECAST_MODEL "FV3CORE")
endif()

# fms
if (FV3_FORECAST_MODEL MATCHES "FV3CORE" OR NOT FV3_FORECAST_MODEL) # OR FV3_FORECAST_MODEL MATCHES "UFS")
  ecbuild_use_package( PROJECT fms REQUIRED )
  include_directories( ${FMS_INCLUDE_DIRS} )
endif()

# fv3
if (FV3_FORECAST_MODEL MATCHES "FV3CORE" OR NOT FV3_FORECAST_MODEL)
  ecbuild_use_package( PROJECT fv3 REQUIRED )
  include_directories( ${FV3_INCLUDE_DIRS} )
endif()

# geos
if (FV3_FORECAST_MODEL MATCHES "GEOS")
  if (NOT DEFINED GEOSgcm_FOUND)
    find_package(MAPL REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake/MAPL )
    find_package(GEOSgcm REQUIRED PATHS ${FV3_FORECAST_MODEL_ROOT}/lib/cmake/GEOSgcm )
  endif()
  find_package(FARGPARSE REQUIRED)
  find_package(FLAP REQUIRED)
  find_package(GFTL REQUIRED)
  find_package(GFTL_SHARED REQUIRED)
  find_package(PFLOGGER REQUIRED)
  find_package(YAFYAML REQUIRED)
  link_libraries(GEOSgcm_GridComp)
endif()


# Linear model
ecbuild_use_package( PROJECT fv3jedilm REQUIRED )
include_directories( ${FV3JEDILM_INCLUDE_DIRS} )

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDI_TPLS atlas atlas_f oops ufo saber )

set( FV3JEDI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_BINARY_DIR}/src
                          ${CMAKE_Fortran_MODULE_DIRECTORY})

set( FV3JEDI_LIBRARIES fv3jedi )

get_directory_property( FV3JEDI_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${FV3JEDI_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND FV3JEDI_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND FV3JEDI_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND FV3JEDI_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################
#set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fdec -fno-range-check -ffree-form ")
include_directories( ${FV3JEDI_INCLUDE_DIRS}
                     ${FV3JEDI_EXTRA_INCLUDE_DIRS} )

# ufs
if (FV3_FORECAST_MODEL MATCHES UFS)
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fdec -fno-range-check -ffree-form ")
  # TODO: fms will be seperated from UFS at some point
  # set CMAKE_PREFIX_PATH to include location of hpc-stack install and ufs install using
  # semi-colon separated lists in quotes e.g.--
  # ecbuild -DCMAKE_PREFIX_PATH="path-to/ufs-weather-model/build/install;/home/mpotts/hpc-stack/install" ...
  find_package(bacio REQUIRED )
  find_package(sigio REQUIRED )
  find_package(sp REQUIRED )
  find_package(w3emc REQUIRED )
  find_package(w3nco REQUIRED )
  find_package(nemsio REQUIRED )
  find_package(fms REQUIRED )
  find_package(fv3atm REQUIRED )
  find_package(ufs REQUIRED )
endif()

add_subdirectory( src )
# ufs
if (FV3_FORECAST_MODEL MATCHES UFS)
  if(POLICY CMP0079)
    cmake_policy(PUSH)
    cmake_policy(SET CMP0079 NEW)
  endif()
  # Link all the libraries
  target_link_libraries( fv3jedi nemsio::nemsio)
  target_link_libraries( fv3jedi w3emc::w3emc_d)
  target_link_libraries( fv3jedi w3nco::w3nco_d)
  target_link_libraries( fv3jedi sp::sp_d)
  target_link_libraries( fv3jedi bacio::bacio_4)
  target_link_libraries( fv3jedi cpl)
  target_link_libraries( fv3jedi fms)
  target_link_libraries( fv3jedi fv3atm)
  target_link_libraries( fv3jedi cpl)
  target_link_libraries( fv3jedi ccppdriver)
  target_link_libraries( fv3jedi ufs)
  target_link_libraries( fv3jedi gfsphysics)
  target_link_libraries( fv3jedi io)
  target_link_libraries( fv3jedi ipd)
  target_link_libraries( fv3jedi stochastic_physics)
  target_link_libraries( fv3jedi NetCDF::NetCDF_Fortran)
  find_package(ESMF MODULE REQUIRED)
  target_include_directories( fv3jedi PRIVATE ${ESMF_F90COMPILEPATHS} )
  add_definitions("-DFRONT_FV3=fv3gfs_cap_mod")
  if(POLICY CMP0079)
    cmake_policy(POP)
  endif()
endif()

add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME fv3jedi )

# print the summary of the configuration
ecbuild_print_summary()
