# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

#######################################################################################
# Fortran
#######################################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -warn unused -init=snan,arrays -traceback -assume realloc_lhs -fPIC -fpe0 -fp-model source -heap-arrays 32 -assume noold_maxminloc -align dcommons" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -xCORE-AVX2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte -traceback -assume realloc_lhs -fPIC -fpe3 -fp-model consistent -g -assume noold_maxminloc -align dcommons" )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()

#######################################################################################
# C
#######################################################################################

# todo

#######################################################################################
# C++
#######################################################################################

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  include( compiler_flags_Clang_CXX )
else()
  message( STATUS "C++ compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()
