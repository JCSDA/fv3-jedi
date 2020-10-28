/*
 * (C) Copyright 2017 UCAR
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

  void fv3jedi_traj_set_f90(F90traj &, const F90state &);
  void fv3jedi_traj_wipe_f90(F90traj &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
