/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_LOCALIZATIONMATRIXFV3JEDIFORTRAN_H_
#define SRC_LOCALIZATIONMATRIXFV3JEDIFORTRAN_H_

#include "Fortran.h"

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
  void fv3jedi_localization_setup_f90(F90lclz &,
                                      const eckit::Configuration * const *,
                                      const F90geom &);
  void fv3jedi_localization_delete_f90(F90lclz &);
  void fv3jedi_localization_mult_f90(const F90lclz &, const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_LOCALIZATIONMATRIXFV3JEDIFORTRAN_H_
