/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_ERRORCOVARIANCE_ERRORCOVARIANCE_INTERFACE_H_
#define FV3JEDI_ERRORCOVARIANCE_ERRORCOVARIANCE_INTERFACE_H_

#include "fv3jedi/Utilities/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace fv3jedi {

extern "C" {

  void fv3jedi_b_setup_f90(F90bmat &,
                           const eckit::Configuration * const *,
                           const F90geom &);
  void fv3jedi_b_delete_f90(F90bmat &);
  void fv3jedi_b_linearize_f90(const F90bmat &,
                               const eckit::Configuration * const *);
  void fv3jedi_b_mult_f90(const F90bmat &,
                          const F90inc &,
                          const F90inc &);
  void fv3jedi_b_invmult_f90(const F90bmat &,
                             const F90inc &,
                             const F90inc &);
  void fv3jedi_b_randomize_f90(const F90bmat &,
                               const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_ERRORCOVARIANCE_ERRORCOVARIANCE_INTERFACE_H_
