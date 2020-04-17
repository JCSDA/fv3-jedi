/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

  void fv3jedi_lvc_model2geovals_create_f90(const F90lm2g &, const F90geom &, const F90state &,
                                            const F90state &, const eckit::Configuration * const *);
  void fv3jedi_lvc_model2geovals_delete_f90(F90lm2g &);
  void fv3jedi_lvc_model2geovals_multiply_f90(const F90lm2g &, const F90geom &, const F90inc &,
                                              const F90inc &);
  void fv3jedi_lvc_model2geovals_multiplyadjoint_f90(const F90lm2g &, const F90geom &,
                                                     const F90inc &, const F90inc &);
  void fv3jedi_lvc_model2geovals_multiplyinverse_f90(const F90lm2g &, const F90geom &,
                                                     const F90inc &, const F90inc &);
  void fv3jedi_lvc_model2geovals_multiplyinverseadjoint_f90(const F90lm2g &, const F90geom &,
                                                            const F90inc &, const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
