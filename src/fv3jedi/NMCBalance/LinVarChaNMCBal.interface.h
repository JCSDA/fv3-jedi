/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_NMCBALANCE_LINVARCHANMCBAL_INTERFACE_H_
#define FV3JEDI_NMCBALANCE_LINVARCHANMCBAL_INTERFACE_H_

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

  void fv3jedi_linvarcha_nmcbal_create_f90(const F90lvcnmcbal &,
                                           const F90geom &,
                                           const F90state &,
                                           const F90state &,
                                           const eckit::Configuration * const *);
  void fv3jedi_linvarcha_nmcbal_delete_f90(F90lvcnmcbal &);
  void fv3jedi_linvarcha_nmcbal_multiply_f90(const F90lvcnmcbal &,
                                             const F90geom &,
                                             const F90inc &,
                                             const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyadjoint_f90(const F90lvcnmcbal &,
                                                    const F90geom &,
                                                    const F90inc &,
                                                    const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyinverse_f90(const F90lvcnmcbal &,
                                                    const F90geom &,
                                                    const F90inc &,
                                                    const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint_f90(const F90lvcnmcbal &,
                                                           const F90geom &,
                                                           const F90inc &,
                                                           const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_NMCBALANCE_LINVARCHANMCBAL_INTERFACE_H_
