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

  void fv3jedi_lineargetvalues_create_f90(F90lineargetvalues &, const F90geom &,
                                          const ufo::Locations &,
                                          const eckit::Configuration * const *);

  void fv3jedi_lineargetvalues_delete_f90(F90lineargetvalues &);

  void fv3jedi_lineargetvalues_set_trajectory_f90(const F90lineargetvalues &, const F90geom &,
                                                  const F90state &,
                                                  const util::DateTime &, const util::DateTime &,
                                                  const ufo::Locations &, const F90goms &);

  void fv3jedi_lineargetvalues_fill_geovals_tl_f90(const F90lineargetvalues &, const F90geom &,
                                                   const F90inc &,
                                                   const util::DateTime &, const util::DateTime &,
                                                   const ufo::Locations &, const F90goms &);

  void fv3jedi_lineargetvalues_fill_geovals_ad_f90(const F90lineargetvalues &, const F90geom &,
                                                   const F90inc &,
                                                   const util::DateTime &, const util::DateTime &,
                                                   const ufo::Locations &, const F90goms &);

};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
