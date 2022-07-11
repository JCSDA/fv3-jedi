/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90vc_M2G;
  extern "C" {
  void fv3jedi_vc_model2geovals_create_f90(const F90vc_M2G &, const F90geom &,
                                           const eckit::LocalConfiguration &);
  void fv3jedi_vc_model2geovals_delete_f90(F90vc_M2G &);
  void fv3jedi_vc_model2geovals_changevar_f90(const F90vc_M2G &, const F90geom &, const F90state &,
                                              const F90state &);
  }  // extern "C"
}  // namespace fv3jedi
