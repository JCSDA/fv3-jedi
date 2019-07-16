/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_INTERFACE_H_
#define FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_INTERFACE_H_

#include "fv3jedi/Utilities/Fortran.h"

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

  void fv3jedi_varcha_a2m_create_f90(const F90vca2m &,
                                     const F90geom &,
                                     const eckit::Configuration * const *);
  void fv3jedi_varcha_a2m_delete_f90(F90vca2m &);
  void fv3jedi_varcha_a2m_changevar_f90(const F90vca2m &,
                                        const F90geom &,
                                        const F90state &,
                                        const F90state &,
                                        util::DateTime * const *);
  void fv3jedi_varcha_a2m_changevarinverse_f90(const F90vca2m &,
                                               const F90geom &,
                                               const F90state &,
                                               const F90state &,
                                               util::DateTime * const *);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_INTERFACE_H_

