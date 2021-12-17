/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace fv3jedi {

// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Tlm key type
typedef int F90tlm;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// ObOp trajectory
typedef int F90ootrj;
// State key
typedef int F90state;
// Increment key
typedef int F90inc;
// Variable change
typedef int F90varcha;
// GetValues key
typedef int F90getvalues;
typedef int F90lineargetvalues;

}  // namespace fv3jedi
