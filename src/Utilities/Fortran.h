/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_UTILITIES_FORTRAN_H_
#define SRC_UTILITIES_FORTRAN_H_

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
// Localization matrix
typedef int F90lclz;
// ObOp trajectory
typedef int F90ootrj;
// VarChange key
typedef int F90vcc2m;
// State key
typedef int F90state;
// Increment key
typedef int F90inc;

}  // namespace fv3jedi
#endif  // SRC_UTILITIES_FORTRAN_H_
