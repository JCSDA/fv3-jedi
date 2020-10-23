/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_UFS_MODELUFS_INTERFACE_H_
#define FV3JEDI_MODEL_UFS_MODELUFS_INTERFACE_H_

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

  void fv3jedi_ufs_create_f90(F90model &, const eckit::Configuration * const *, const F90geom &);
  void fv3jedi_ufs_delete_f90(F90model &);


  void fv3jedi_ufs_initialize_f90(const F90model &, const F90state &, util::DateTime * const *);

  void fv3jedi_ufs_step_f90(const F90model &, const F90state &, util::DateTime * const *,
                            util::DateTime * const *);

  void fv3jedi_ufs_finalize_f90(const F90model &, const F90inc &, util::DateTime * const *);

}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_UFS_MODELUFS_INTERFACE_H_
