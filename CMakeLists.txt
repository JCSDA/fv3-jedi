# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDI
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedi VERSION 1.9.0 LANGUAGES C CXX Fortran )

option(OPENMP "Build fv3-jedi with OpenMP support" ON)

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "fv3-jedi" )

find_package( ecbuild 3.6 REQUIRED )

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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

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
find_package( NetCDF REQUIRED COMPONENTS C Fortran )

# atlas
find_package( atlas 0.35.0 REQUIRED )

# oops
find_package( oops 1.10.0 REQUIRED )

# saber
find_package( saber 1.10.0 REQUIRED )

# crtm
find_package( crtm 2.2.3 REQUIRED )

# ropp
find_package( ropp-ufo QUIET )

# gsibec
find_package( gsibec 1.2.1 QUIET )
if( gsibec_FOUND )
    message( STATUS "FV3-SABER block GSI is enabled" )
    add_definitions(-DENABLE_GSIBEC=1)
endif()

# ncep-sp
find_package( sp QUIET )

# ufo
find_package( ufo 1.10.0 REQUIRED )

# femps
find_package( femps 1.3.0 REQUIRED )

# vader
find_package( vader 1.7.0 REQUIRED )

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
find_package( fv3jedilm 1.5.0 REQUIRED )

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDI_TPLS atlas atlas_f oops ufo saber )

set( FV3JEDI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_BINARY_DIR}/src
                          ${DEPEND_LIB_ROOT}/include
                          ${DEPEND_LIB_ROOT}/mod
                          ${DEPEND_LIB_ROOT}/include/fv3
                          ${CMAKE_Fortran_MODULE_DIRECTORY})
if (FV3_FORECAST_MODEL MATCHES UFS)
  if( UFS_APP MATCHES "^(ATM)$" )
    set( FV3JEDI_LIBRARIES fv3jedi
                           ccpp_framework
                           ccpp_physics )
  elseif( UFS_APP MATCHES "^(ATMAERO)$" )
    find_package (GFTL_SHARED QUIET)
    find_package (FARGPARSE QUIET)
    find_package (PFLOGGER QUIET)
    find_package (MAPL REQUIRED)
    set( FV3JEDI_LIBRARIES fv3jedi
                           ufs_aerosols
                           aerosols_gridcomp
                           gocart_gridcomp
                           ca2g_gridcomp
                           chem_shared
                           gocart_process
                           du2g_gridcomp
                           ni2g_gridcomp
                           ss2g_gridcomp
                           su2g_gridcomp
                           ga_environment
                           MAPL
                           ccpp_framework
                           ccpp_physics
                           GFTL_SHARED::gftl-shared )
  elseif( UFS_APP MATCHES "^(S2S)$" )
    find_package( PIO 2.5.3 REQUIRED COMPONENTS C Fortran)
    set( FV3JEDI_LIBRARIES fv3jedi
                           gsw
                           mom6
                           cice
                           cmeps
                           PIO::PIO_Fortran
                           PIO::PIO_C
                           ccpp_framework
                           ccpp_physics )
  else()
    message(FATAL_ERROR "fv3jedi FV3JEDI_LIBRARIES not configured for UFS_APP ${UFS_APP}")
  endif()
else()
  set( FV3JEDI_LIBRARIES fv3jedi )
endif()

get_directory_property( FV3JEDI_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${FV3JEDI_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND FV3JEDI_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND FV3JEDI_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND FV3JEDI_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Test data
################################################################################

# DH 2023/10/06 - this is unfortunately still true ...
### The fv3-jedi-data package is problematic in bundles because it has some special CMake logic
### at the root bundle location. So, we cannot find the fv3-jedi-data while in fv3-jedi!
### Basically, with the current logic we hack fv3-jedi-data_FOUND to match ENABLE_FV3_JEDI_DATA
#find_package( fv3-jedi-data QUIET )
if (ENABLE_FV3_JEDI_DATA)
  set (fv3-jedi-data_FOUND TRUE)
else()
  # Set values for downloading fv3-jedi test data from tarball
  set (fv3-jedi_data_download_url "https://bin.ssec.wisc.edu/pub/s4/CRTM/file")
  set (fv3-jedi_data_tag "1.9.0")
  set (fv3-jedi_data_tarball_name "fv3-jedi_testinput_tier_1_${fv3-jedi_data_tag}.tar.gz")
  set (fv3-jedi_data_tarball_hash "8ea0904eda4e4b3840536b8c744b3d9e")
endif()

################################################################################
# Sources
################################################################################

# ufs
if (FV3_FORECAST_MODEL MATCHES UFS)
  # Building with the UFS is now accomplished only via the ufs-bundle
  # which will also build the UFS as an external library with a prescribed configuration
  # There is now a unit test which will pull data from an s3 bucket and convert it
  # into a run directory to be used with ctests involving UFS coupling (ATM-only) 
  find_package(bacio REQUIRED )
  find_package(sigio REQUIRED )
  find_package(sp REQUIRED )
  find_package(w3emc REQUIRED )
  find_package(nemsio REQUIRED )
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
