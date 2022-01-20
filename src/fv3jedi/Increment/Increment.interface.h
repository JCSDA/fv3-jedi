/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"
#include "fv3jedi/Utilities/interface.h"
#include "oops/base/Variables.h"

namespace fv3jedi {
extern "C" {
  void fv3jedi_increment_create_f90(F90inc &, const F90geom &, const oops::Variables &,
                                    const util::DateTime &);
  void fv3jedi_increment_delete_f90(F90inc &);
  void fv3jedi_increment_copy_f90(const F90inc &, const F90inc &);
  void fv3jedi_increment_zero_f90(const F90inc &);
  void fv3jedi_increment_ones_f90(const F90inc &);
  void fv3jedi_increment_self_add_f90(const F90inc &, const F90inc &);
  void fv3jedi_increment_self_sub_f90(const F90inc &, const F90inc &);
  void fv3jedi_increment_self_mul_f90(const F90inc &, const double &);
  void fv3jedi_increment_axpy_inc_f90(const F90inc &, const double &, const F90inc &);
  void fv3jedi_increment_axpy_state_f90(const F90inc &, const double &, const F90state &);
  void fv3jedi_increment_dot_prod_f90(const F90inc &, const F90inc &, double &);
  void fv3jedi_increment_self_schur_f90(const F90inc &, const F90inc &);
  void fv3jedi_increment_random_f90(const F90inc &);
  void fv3jedi_increment_diff_incr_f90(const F90inc &, const F90state &, const F90state &,
                                       const F90geom &);
  void fv3jedi_increment_change_resol_f90(const F90inc &, const F90geom &, const F90inc &,
                                          const F90geom &);
  void fv3jedi_increment_sizes_f90(const F90inc &, int &);
  void fv3jedi_increment_norm_f90(const F90inc &, double &);
  void fv3jedi_increment_update_fields_f90(F90inc &, const F90geom &, const oops::Variables &);
  void fv3jedi_increment_set_atlas_f90(const F90inc &, const F90geom &, const oops::Variables &,
                                       atlas::field::FieldSetImpl *);
  void fv3jedi_increment_to_atlas_f90(const F90inc &, const F90geom &, const oops::Variables &,
                                      atlas::field::FieldSetImpl *);
  void fv3jedi_increment_from_atlas_f90(const F90inc &, const F90geom &, const oops::Variables &,
                                        atlas::field::FieldSetImpl *);
  void fv3jedi_increment_dirac_f90(const F90inc &, const eckit::Configuration &,
                                   const F90geom &);
  void fv3jedi_increment_serialize_f90(const F90inc &, const std::size_t &, double[]);
  void fv3jedi_increment_deserialize_f90(const F90inc &, const std::size_t &, const double[],
                                         const std::size_t &);
  void fv3jedi_increment_getpoint_f90(const F90inc &, const F90iter &, double &, const int &);
  void fv3jedi_increment_setpoint_f90(F90inc &, const F90iter &, const double &, const int &);
  void fv3jedi_increment_getnfieldsncube_f90(const F90state &, int &, int &);
  void fv3jedi_increment_getminmaxrms_f90(const F90state &, int &, const int &, char*, double &);
}  // extern "C"
}  // namespace fv3jedi
