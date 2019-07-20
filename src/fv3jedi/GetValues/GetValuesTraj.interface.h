/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GETVALUES_GETVALUESTRAJ_INTERFACE_H_
#define FV3JEDI_GETVALUES_GETVALUESTRAJ_INTERFACE_H_

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

  void fv3jedi_getvalues_traj_setup_f90(const F90ootrj &);
  void fv3jedi_getvalues_traj_delete_f90(const F90ootrj &);

};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_GETVALUES_GETVALUESTRAJ_INTERFACE_H_
