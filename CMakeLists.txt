# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# FV3JEDI 
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( fv3jedi C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

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
# Dependencies
################################################################################

# Boost
set( Boost_MINIMUM_VERSION "1.47" )
find_package( Boost ${Boost_MINIMUM_VERSION} REQUIRED
              COMPONENTS unit_test_framework )
include_directories( ${Boost_INCLUDE_DIR} )

# fms
ecbuild_use_package( PROJECT fms REQUIRED )
include_directories( ${FMS_INCLUDE_DIRS} )

# fv3
ecbuild_use_package( PROJECT fv3 REQUIRED )
include_directories( ${FV3_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo )
include_directories( ${UFO_INCLUDE_DIRS} )

################################################################################
# Definitions
################################################################################

#Uncomment to compile with TLM/ADM
#set(COMPDEFS "${COMPDEFS};-DTLADPRES")
add_definitions( ${COMPDEFS} )

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDI_TPLS ufo oops )

set( FV3JEDI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
                        ${CMAKE_Fortran_MODULE_DIRECTORY} )
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

include( fv3jedi_compiler_flags )
include_directories( ${FV3JEDI_INCLUDE_DIRS} ${FV3JEDI_EXTRA_INCLUDE_DIRS}
                     ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( mains )
add_subdirectory( src )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME fv3jedi )

# print the summary of the configuration
ecbuild_print_summary()

