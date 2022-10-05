# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDI
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedi VERSION 1.2.0 LANGUAGES C CXX Fortran )

option(OPENMP "Build fv3-jedi with OpenMP support" ON)

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

set( FV3JEDI_LINKER_LANGUAGE CXX )

################################################################################
# Standard compiler flags
################################################################################

include( fv3jedi_extra_macros )
option( ENABLE_FV3JEDI_DOC "Build FV3-JEDI documentation" OFF )

################################################################################
# Set default model to GFDL_atmos_cubed_sphere
################################################################################

if (NOT FV3_FORECAST_MODEL MATCHES GEOS AND NOT FV3_FORECAST_MODEL MATCHES UFS)
  set(FV3_FORECAST_MODEL "FV3CORE")
endif()

################################################################################
# Dependencies
################################################################################

# jedi-cmake
find_package( jedicmake QUIET )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )

# atlas
find_package( atlas 0.27.0 REQUIRED )

# oops
find_package( oops 1.3.0 REQUIRED )

# saber
find_package( saber 1.3.0 REQUIRED )

# crtm
find_package( crtm 2.2.3 REQUIRED )

# ropp
find_package( ropp-ufo QUIET )

# geos-aero
find_package( geos-aero QUIET )

# gsibec
find_package( gsibec 1.0.6 QUIET )
if( gsibec_FOUND )
    message( STATUS "FV3-SABER block GSI is enabled" )
    add_definitions(-DENABLE_GSIBEC=1)
endif()

# ncep-sp
find_package( sp QUIET )

# ufo
find_package( ufo 1.3.0 REQUIRED )

# femps
find_package( femps 1.0.0 REQUIRED )

# vader
find_package( vader 1.0.0 REQUIRED )

#OpenMP
if(OPENMP)
  find_package( OpenMP COMPONENTS Fortran )
endif()

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )

# Set default forecast model
# --------------------------
if (NOT FV3_FORECAST_MODEL MATCHES GEOS AND NOT FV3_FORECAST_MODEL MATCHES UFS)
  set(FV3_FORECAST_MODEL "FV3CORE")
endif()

# fms
if (FV3_FORECAST_MODEL MATCHES "FV3CORE" OR NOT FV3_FORECAST_MODEL) # OR FV3_FORECAST_MODEL MATCHES "UFS")
  find_package( fms REQUIRED )
endif()

# fv3
if (FV3_FORECAST_MODEL MATCHES "FV3CORE" OR NOT FV3_FORECAST_MODEL)
  find_package( fv3 REQUIRED )
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
find_package( fv3jedilm 1.0.0 REQUIRED )

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

# ufs
if (FV3_FORECAST_MODEL MATCHES UFS)
  # TODO: fms will be seperated from UFS at some point
  # Requires an installion of ufs-weather-model built with -DJEDI_CONTROL=ON and setting
  # environment variable CMAKE_Platform to linux.custom to ensure INLINE_POST=OFF
  # set CMAKE_PREFIX_PATH to include location of hpc-stack install and ufs install using
  # semi-colon separated lists in quotes e.g.--
  # ecbuild -DCMAKE_PREFIX_PATH="path-to/ufs-weather-model/build/install;path-to/hpc-stack/install" ...
  find_package(bacio REQUIRED )
  find_package(sigio REQUIRED )
  find_package(sp REQUIRED )
  find_package(w3emc REQUIRED )
  find_package(w3nco REQUIRED )
  find_package(nemsio REQUIRED )
  find_package(fms REQUIRED )
  find_package(stochastic_physics REQUIRED )
  find_package(ccpp REQUIRED )
  find_package(ccppphys REQUIRED )
  find_package(fv3atm REQUIRED )
  find_package(ufs REQUIRED )
endif()

add_subdirectory( src )
add_subdirectory( test )

if( ENABLE_FV3JEDI_DOC )
    add_subdirectory( Documents )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME fv3jedi )

# print the summary of the configuration
ecbuild_print_summary()
