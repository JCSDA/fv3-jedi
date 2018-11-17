/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_MODEL_GFS_MODELGFSFV3JEDIFORTRAN_H_
#define SRC_MODEL_GFS_MODELGFSFV3JEDIFORTRAN_H_

#include "Fortran.h"
#include "TrajFV3JEDIFortran.h"

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

  void fv3jedi_gfs_create_f90(const eckit::Configuration * const *,
                               const F90geom &,
                               F90model &);
  void fv3jedi_gfs_delete_f90(F90model &);


  void fv3jedi_gfs_initialize_f90(const F90model &,
                                   const F90state &);

  void fv3jedi_gfs_step_f90(const F90model &,
                             const F90state &,
                             util::DateTime * const *);

  void fv3jedi_gfs_finalize_f90(const F90model &,
                                 const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_MODEL_GFS_MODELGFSFV3JEDIFORTRAN_H_
