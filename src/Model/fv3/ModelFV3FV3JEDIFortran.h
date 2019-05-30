/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_MODEL_FV3_MODELFV3FV3JEDIFORTRAN_H_
#define SRC_MODEL_FV3_MODELFV3FV3JEDIFORTRAN_H_

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

  void fv3jedi_fv3_create_f90(const eckit::Configuration * const *,
                              const F90geom &,
                              F90model &);
  void fv3jedi_fv3_delete_f90(F90model &);


  void fv3jedi_fv3_initialize_f90(const F90model &,
                                  const F90state &);

  void fv3jedi_fv3_step_f90(const F90model &,
                            const F90state &,
                            const F90geom &,
                            util::DateTime * const *);

  void fv3jedi_fv3_finalize_f90(const F90model &,
                                const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_MODEL_FV3_MODELFV3FV3JEDIFORTRAN_H_
