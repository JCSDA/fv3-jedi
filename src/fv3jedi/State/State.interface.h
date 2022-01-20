/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/Utilities/interface.h"
#include "oops/base/Variables.h"

namespace fv3jedi {
extern "C" {
  void fv3jedi_state_create_f90(F90state &, const F90geom &, const oops::Variables &,
                                const util::DateTime &);
  void fv3jedi_state_delete_f90(F90state &);
  void fv3jedi_state_copy_f90(const F90state &, const F90state &);
  void fv3jedi_state_zero_f90(const F90state &);
  void fv3jedi_state_axpy_f90(const F90state &, const double &, const F90state &);
  void fv3jedi_state_add_incr_f90(const F90geom &, const F90state &, const F90inc &);
  void fv3jedi_state_change_resol_f90(const F90state &, const F90geom &, const F90state &,
                                      const F90geom &);
  void fv3jedi_state_update_fields_f90(F90state &, const F90geom &, const oops::Variables &);
  void fv3jedi_state_analytic_init_f90(const F90state &, const F90geom &,
                                       const eckit::Configuration &);
  void fv3jedi_state_fill_geom_orography_f90(const F90state &, F90geom &);
  void fv3jedi_state_sersize_f90(const F90state &, int &);

  void fv3jedi_state_serialize_f90(const F90state &, const std::size_t &, double[]);

  void fv3jedi_state_deserialize_f90(const F90state &, const std::size_t &, const double[],
                                     const std::size_t &);

  void fv3jedi_state_norm_f90(const F90state &, double &);
  void fv3jedi_state_getnfieldsncube_f90(const F90state &, int &, int &);
  void fv3jedi_state_getminmaxrms_f90(const F90state &, int &, const int &, char*, double &);
};  // extern "C"
}  // namespace fv3jedi
