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

  void fv3jedi_geos_create_f90(const eckit::Configuration &,
                               const F90geom &,
                               F90model &);
  void fv3jedi_geos_delete_f90(F90model &);


  void fv3jedi_geos_initialize_f90(const F90model &,
                                   const F90state &);

  void fv3jedi_geos_step_f90(const F90model &,
                             const F90state &);

  void fv3jedi_geos_finalize_f90(const F90model &,
                                 const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
