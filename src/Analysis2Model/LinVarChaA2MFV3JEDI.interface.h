/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_ANALYSIS2MODEL_LINVARCHAA2MFV3JEDI_INTERFACE_H_
#define SRC_ANALYSIS2MODEL_LINVARCHAA2MFV3JEDI_INTERFACE_H_

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

  void fv3jedi_linvarcha_a2m_create_f90(const F90lvca2m &,
                                        const F90geom &,
                                        const F90state &,
                                        const F90state &,
                                        const eckit::Configuration * const *);
  void fv3jedi_linvarcha_a2m_delete_f90(F90lvca2m &);
  void fv3jedi_linvarcha_a2m_multiply_f90(const F90lvca2m &,
                                          const F90geom &,
                                          const F90inc &,
                                          const F90inc &);
  void fv3jedi_linvarcha_a2m_multiplyadjoint_f90(const F90lvca2m &,
                                                 const F90geom &,
                                                 const F90inc &,
                                                 const F90inc &);
  void fv3jedi_linvarcha_a2m_multiplyinverse_f90(const F90lvca2m &,
                                                 const F90geom &,
                                                 const F90inc &,
                                                 const F90inc &);
  void fv3jedi_linvarcha_a2m_multiplyinverseadjoint_f90(const F90lvca2m &,
                                                        const F90geom &,
                                                        const F90inc &,
                                                        const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_ANALYSIS2MODEL_LINVARCHAA2MFV3JEDI_INTERFACE_H_
