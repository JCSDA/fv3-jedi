/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_INTERFACE_H_
#define FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_INTERFACE_H_

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

  void fv3jedi_varcha_c2a_create_f90(const F90vcc2a &,
                                     const F90geom &,
                                     const eckit::Configuration * const *);
  void fv3jedi_varcha_c2a_delete_f90(F90vcc2a &);
  void fv3jedi_varcha_c2a_changevar_f90(const F90vcc2a &,
                                        const F90geom &,
                                        const F90state &,
                                        const F90state &);
  void fv3jedi_varcha_c2a_changevarinverse_f90(const F90vcc2a &,
                                               const F90geom &,
                                               const F90state &,
                                               const F90state &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_INTERFACE_H_
