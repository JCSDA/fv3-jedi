/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_VARCHAC2MFV3JEDIFORTRAN_H_
#define SRC_VARCHAC2MFV3JEDIFORTRAN_H_

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

  void fv3jedi_varcha_c2m_setup_f90(const F90vcc2m &,
                                    const F90geom &,
                                    const F90state &,
                                    const F90state &,
                                    const eckit::Configuration * const *);
  void fv3jedi_varcha_c2m_delete_f90(F90vcc2m &);
  void fv3jedi_varcha_c2m_multiply_f90(const F90vcc2m &,
                                       const F90geom &,
                                       const F90inc &,
                                       const F90inc &);
  void fv3jedi_varcha_c2m_multiplyadjoint_f90(const F90vcc2m &,
                                              const F90geom &,
                                              const F90inc &,
                                              const F90inc &);
  void fv3jedi_varcha_c2m_multiplyinverse_f90(const F90vcc2m &,
                                              const F90geom &,
                                              const F90inc &,
                                              const F90inc &);
  void fv3jedi_varcha_c2m_multiplyinverseadjoint_f90(const F90vcc2m &,
                                                     const F90geom &,
                                                     const F90inc &,
                                                     const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_VARCHAC2MFV3JEDIFORTRAN_H_
