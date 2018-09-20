/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_STATEFV3JEDIFORTRAN_H_
#define SRC_STATEFV3JEDIFORTRAN_H_

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

// -----------------------------------------------------------------------------
//  State
// -----------------------------------------------------------------------------
  void fv3jedi_state_create_f90(F90state &,
                                const F90geom &,
                                const eckit::Configuration * const *);

  void fv3jedi_state_delete_f90(F90state &);

  void fv3jedi_state_copy_f90(const F90state &,
                              const F90state &);

  void fv3jedi_state_zero_f90(const F90state &);

  void fv3jedi_state_axpy_f90(const F90state &,
                              const double &,
                              const F90state &);

  void fv3jedi_state_add_incr_f90(const F90geom &,
                                  const F90state &,
                                  const F90inc &);

  void fv3jedi_state_change_resol_f90(const F90state &,
                                      const F90state &);

  void fv3jedi_state_read_file_f90(const F90geom &,
                                   const F90state &,
                                   const eckit::Configuration * const *,
                                   util::DateTime * const *);

  void fv3jedi_state_analytic_init_f90(const F90state &,
                                       const F90geom &,
                                       const eckit::Configuration * const *,
                                       util::DateTime * const *);

  void fv3jedi_state_write_file_f90(const F90geom &,
                                    const F90state &,
                                    const eckit::Configuration * const *,
                                    const util::DateTime * const *);

  void fv3jedi_state_getvalues_notraj_f90(const F90geom &,
                                          const F90state &,
                                          const F90locs &,
                                          const eckit::Configuration * const *,
                                          const F90goms &);

  void fv3jedi_state_getvalues_f90(const F90geom &,
                                   const F90state &,
                                   const F90locs &,
                                   const eckit::Configuration * const *,
                                   const F90goms &,
                                   const F90ootrj &);

  void fv3jedi_state_gpnorm_f90(const F90state &,
                                const int &,
                                double &);

  void fv3jedi_state_sizes_f90(const F90state &,
                               int &,
                               int &,
                               int &);

  void fv3jedi_state_rms_f90(const F90state &,
                             double &);

};  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_STATEFV3JEDIFORTRAN_H_
