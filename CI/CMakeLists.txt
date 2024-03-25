# (C) Copyright 2017-2023 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# ##########
# fv3-bundle
# ##########

cmake_minimum_required( VERSION 3.14 FATAL_ERROR )

project( fv3-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

# Enable MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

# Initialize
ecbuild_bundle_initialize()

# Add git functions from jedi-cmake
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )

# External (optional) observation operators
# ------------------------------
option("BUNDLE_SKIP_ROPP-UFO"  "Don't build ROPP-UFO"  "ON") # Skip ropp-ufo build unless user passes -DBUNDLE_SKIP_ROPP-UFO=OFF

ecbuild_bundle( PROJECT ropp-ufo GIT "https://github.com/jcsda-internal/ropp-test.git" )
ecbuild_bundle( PROJECT crtm GIT "https://github.com/jcsda/CRTMv3.git" )

# Outside JEDI repositories
# -------------------------
option("BUNDLE_SKIP_GSIBEC"  "Don't build GSIBEC"  "ON") # Skip gsibec build unless user passes -DBUNDLE_SKIP_GSIBEC=OFF
ecbuild_bundle( PROJECT gsibec GIT "https://github.com/GEOS-ESM/GSIbec.git" )

# Core JEDI repositories
# ----------------------
option(ENABLE_IODA_DATA "Obtain ufo test data from ioda-data repository (vs tarball)" ON)
option(ENABLE_UFO_DATA "Obtain ufo test data from ufo-data repository (vs tarball)" ON)

ecbuild_bundle( PROJECT oops  GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/jcsda-internal/saber.git" )
# No need to clone, we only run the fv3-jedi ctests
#ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ufo-data.git" )
ecbuild_bundle( PROJECT ioda  GIT "https://github.com/jcsda-internal/ioda.git" )
# No need to clone, we only run the fv3-jedi ctests
#ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" )
ecbuild_bundle( PROJECT ufo   GIT "https://github.com/jcsda-internal/ufo.git" )
ecbuild_bundle( PROJECT vader GIT "https://github.com/jcsda-internal/vader.git" )

# Options for building with certain models
# ----------------------------------------
set(FV3_FORECAST_MODEL "FV3CORE" CACHE STRING "Choose which MODEL to build with")
set_property(CACHE FV3_FORECAST_MODEL PROPERTY STRINGS "FV3CORE" "UFS" "GEOS")
if(NOT FV3_FORECAST_MODEL MATCHES "^(FV3CORE|UFS|GEOS)$")
  ecbuild_error("Valid options for -DFV3_FORECAST_MODEL are FV3CORE|UFS|GEOS")
endif()
ecbuild_info("Building with FV3-based model: ${FV3_FORECAST_MODEL}")

# FMS and FV3 dynamical core model are needed for FV3CORE
# -------------------------------------------------------
if (FV3_FORECAST_MODEL MATCHES "FV3CORE")
  option("BUNDLE_SKIP_FMS" "Don't build FMS" "OFF")
  ecbuild_bundle( PROJECT fms GIT "https://github.com/jcsda/FMS.git" )
  ecbuild_bundle( PROJECT fv3 GIT "https://github.com/jcsda/GFDL_atmos_cubed_sphere.git" )
endif()

# fv3-jedi linear model
# ---------------------
ecbuild_bundle( PROJECT fv3-jedi-lm GIT "https://github.com/jcsda-internal/fv3-jedi-linearmodel.git" )

# fv3-jedi and associated repositories
# ------------------------------------
option(ENABLE_FV3_JEDI_DATA "Obtain fv3-jedi test data from fv3-jedi-data repository (vs tarball)" ON)

ecbuild_bundle( PROJECT femps    GIT "https://github.com/jcsda-internal/femps.git" )
ecbuild_bundle( PROJECT fv3-jedi-data GIT "https://github.com/JCSDA-internal/fv3-jedi-data.git" )
ecbuild_bundle( PROJECT fv3-jedi GIT "https://github.com/jcsda-internal/fv3-jedi.git" )

# Set GIT_BRANCH_FUNC to FV3-JEDI's current branch so that it can be used for fv3-jedi-data
find_branch_name(REPO_DIR_NAME fv3-jedi)

# If FV3-JEDI's current branch is available in fv3-jedi-data repo, that branch will be checked out
branch_checkout (REPO_DIR_NAME fv3-jedi-data
                 BRANCH ${GIT_BRANCH_FUNC} )

# Build Doxygen documentation
# ---------------------------
option(BUILD_FV3_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_FV3_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_FV3_BUNDLE_DOC)

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
