/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90vc_A2M;
  extern "C" {
  void fv3jedi_varcha_a2m_create_f90(const F90vc_A2M &, const F90geom &,
                                     const eckit::LocalConfiguration &);
  void fv3jedi_varcha_a2m_delete_f90(F90vc_A2M &);
  void fv3jedi_varcha_a2m_changevar_f90(const F90vc_A2M &, const F90geom &, const F90state &,
                                        const F90state &);
  void fv3jedi_varcha_a2m_changevarinverse_f90(const F90vc_A2M &, const F90geom &, const F90state &,
                                               const F90state &);
  }  // extern "C"
}  // namespace fv3jedi
