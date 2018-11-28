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

# fv3-jedi-lm
ecbuild_use_package( PROJECT fv3-jedi-lm REQUIRED )
include_directories( ${FV3JEDILM_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo )
include_directories( ${UFO_INCLUDE_DIRS} )

# fms
include_directories( ${FMS_INCLUDE_DIRS} )

#gfs or geo
if( ${FV3BASEDMODEL_FOUND} )
    if ( ${GFS_FOUND} )
      message("GFS FOUND: Including GFS interface")
    elseif( ${GEOS_FOUND} )
      message("GEOS FOUND: Including GEOS interface")
    endif()   
    include_directories( ${FV3BASEDMODEL_INCLUDE_DIRS} )
else
   message("GFS AND GEOS NOT FOUND: Including only FV3-JEDI-LM interface")
endif()

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )
ecbuild_include_mpi()
link_libraries(${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})

################################################################################
# Definitions
################################################################################

# Additional user supplied compiler definitions
set( COMPDEFS_ "${COMPDEFS}")
add_definitions( ${COMPDEFS_} )
include( fv3jedi_compiler_flags )

################################################################################
# Export package info
################################################################################

list( APPEND FV3JEDI_TPLS ufo oops fv3-jedi-lm)

set( FV3JEDI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} 
                          ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_BINARY_DIR}/src
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

include_directories( ${FV3JEDI_INCLUDE_DIRS}
                     ${FV3JEDI_EXTRA_INCLUDE_DIRS}
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ErrorCovariance
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ControlToModel
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/GetValues
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Increment
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Localization
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Model
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Model/fv3
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Model/geos
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Model/gfs
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Model/sar
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/ModelBias
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Run
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/State
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Utilities
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/VariableChange
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/Tlm )

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

ecbuild_add_test( TARGET fv3jedi_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

