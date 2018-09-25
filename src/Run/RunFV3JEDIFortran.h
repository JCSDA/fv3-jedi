/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_RUN_RUNFV3JEDIFORTRAN_H_
#define SRC_RUN_RUNFV3JEDIFORTRAN_H_

#include "RunFV3JEDIFortran.h"

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
#endif  // SRC_RUN_RUNFV3JEDIFORTRAN_H_
