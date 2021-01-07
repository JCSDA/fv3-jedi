/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/Utilities/interface.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace ufo {
  class Locations;
}

namespace fv3jedi {

extern "C" {

  void fv3jedi_getvalues_create_f90(F90getvalues &, const F90geom &, const ufo::Locations &);

  void fv3jedi_getvalues_delete_f90(F90getvalues &);

  void fv3jedi_getvalues_fill_geovals_f90(const F90getvalues &, const F90geom &, const F90state &,
                                          const util::DateTime &, const util::DateTime &,
                                          const ufo::Locations &, const F90goms &);

};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
