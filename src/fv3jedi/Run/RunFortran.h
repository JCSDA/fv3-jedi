/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_RUN_RUNFORTRAN_H_
#define FV3JEDI_RUN_RUNFORTRAN_H_

#include "RunFortran.h"

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

  void fv3jedi_setup_f(const eckit::Configuration * const *);
  void fv3jedi_finalize_f();

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_RUN_RUNFORTRAN_H_
