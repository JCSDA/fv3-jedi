/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_INCREMENTFV3JEDIFORTRAN_H_
#define SRC_INCREMENTFV3JEDIFORTRAN_H_

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
//  Increment
// -----------------------------------------------------------------------------

  void fv3jedi_increment_create_f90(F90inc &,
                                    const F90geom &,
                                    const eckit::Configuration * const *);

  void fv3jedi_increment_delete_f90(F90inc &);

  void fv3jedi_increment_copy_f90(const F90inc &,
                                  const F90inc &);
 
  void fv3jedi_increment_zero_f90(const F90inc &);

  void fv3jedi_increment_self_add_f90(const F90inc &,
                                      const F90inc &);

  void fv3jedi_increment_self_sub_f90(const F90inc &,
                                      const F90inc &);

  void fv3jedi_increment_self_mul_f90(const F90inc &,
                                      const double &);

  void fv3jedi_increment_axpy_inc_f90(const F90inc &,
                                      const double &,
                                      const F90inc &);

  void fv3jedi_increment_axpy_state_f90(const F90inc &,
                                        const double &,
                                        const F90state &);

  void fv3jedi_increment_dot_prod_f90(const F90inc &,
                                      const F90inc &,
                                      double &);

  void fv3jedi_increment_self_schur_f90(const F90inc &,
                                        const F90inc &);

  void fv3jedi_increment_random_f90(const F90inc &);

  void fv3jedi_increment_diff_incr_f90(const F90inc &,
                                       const F90state &,
                                       const F90state &);

  void fv3jedi_increment_change_resol_f90(const F90inc &,
                                          const F90inc &);

  void fv3jedi_increment_read_file_f90(const F90geom &,
                                       const F90inc &,
                                       const eckit::Configuration * const *,
                                       util::DateTime * const *);

  void fv3jedi_increment_write_file_f90(const F90geom &,
                                        const F90inc &,
                                        const eckit::Configuration * const *,
                                        const util::DateTime * const *);

  void fv3jedi_increment_getvalues_tl_f90(const F90geom &,
                                          const F90inc &,
                                          const F90locs &,
                                          const eckit::Configuration * const *,
                                          const F90goms &,
                                          const F90ootrj &);

  void fv3jedi_increment_getvalues_ad_f90(const F90geom &,
                                          const F90inc &,
                                          const F90locs &,
                                          const eckit::Configuration * const *,
                                          const F90goms &,
                                          const F90ootrj &);

  void fv3jedi_increment_gpnorm_f90(const F90inc &,
                                    const int &,
                                    double &);

  void fv3jedi_increment_sizes_f90(const F90inc &,
                                   int &,
                                   int &,
                                   int &);

  void fv3jedi_increment_rms_f90(const F90inc &,
                                 double &);

  void fv3jedi_increment_ug_coord_f90(const F90inc &,
                                      const int &,
                                      const int &,
                                      const F90geom &);

  void fv3jedi_increment_increment_to_ug_f90(const F90inc &,
                                             const int &,
                                             const int &);

  void fv3jedi_increment_increment_from_ug_f90(const F90inc &,
                                               const int &);

  void fv3jedi_increment_dirac_f90(const F90inc &,
                                   const eckit::Configuration * const *,
                                   const F90geom &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_INCREMENTFV3JEDIFORTRAN_H_
