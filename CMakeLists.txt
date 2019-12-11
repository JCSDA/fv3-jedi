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

include( fv3jedi_compiler_flags )
include( fv3jedi_extra_macros )

################################################################################
# Dependencies
################################################################################

# Boost
include_directories( ${Boost_INCLUDE_DIR} )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# fv3-jedi-lm
ecbuild_use_package( PROJECT fv3jedilm REQUIRED )
include_directories( ${FV3JEDILM_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo )
include_directories( ${UFO_INCLUDE_DIRS} )

# femps
ecbuild_use_package( PROJECT femps )
include_directories( ${FEMPS_INCLUDE_DIRS} )

#gfs or geos
if( ${FV3BASEDMODEL_FOUND} )
    if ( ${GFS_FOUND} )
      message("GFS FOUND: Including GFS interface")
    elseif( ${GEOS_FOUND} )
      message("GEOS FOUND: Including GEOS interface")
      include( geos_compiler_flags )
    endif()
    include_directories( ${FV3BASEDMODEL_INCLUDE_DIRS} )
else()
   message("GFS AND GEOS NOT FOUND: Including only FV3-JEDI-LM interface and using local FMS")
   ecbuild_use_package( PROJECT fms REQUIRED )
endif()

include_directories( ${FMS_INCLUDE_DIRS} )

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDI_TPLS ufo oops )

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

include_directories( ${FV3JEDI_INCLUDE_DIRS}
                     ${FV3JEDI_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME fv3jedi )

# print the summary of the configuration
ecbuild_print_summary()
