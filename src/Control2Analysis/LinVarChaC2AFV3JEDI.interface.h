/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_INTERFACE_H_
#define SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_INTERFACE_H_

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

  void fv3jedi_linvarcha_c2a_setup_f90(const F90lvcc2a &,
                                       const F90geom &,
                                       const F90state &,
                                       const F90state &,
                                       const eckit::Configuration * const *);
  void fv3jedi_linvarcha_c2a_delete_f90(F90lvcc2a &);
  void fv3jedi_linvarcha_c2a_multiply_f90(const F90lvcc2a &,
                                          const F90geom &,
                                          const F90inc &,
                                          const F90inc &);
  void fv3jedi_linvarcha_c2a_multiplyadjoint_f90(const F90lvcc2a &,
                                                 const F90geom &,
                                                 const F90inc &,
                                                 const F90inc &);
  void fv3jedi_linvarcha_c2a_multiplyinverse_f90(const F90lvcc2a &,
                                                 const F90geom &,
                                                 const F90inc &,
                                                 const F90inc &);
  void fv3jedi_linvarcha_c2a_multiplyinverseadjoint_f90(const F90lvcc2a &,
                                                        const F90geom &,
                                                        const F90inc &,
                                                        const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_INTERFACE_H_
