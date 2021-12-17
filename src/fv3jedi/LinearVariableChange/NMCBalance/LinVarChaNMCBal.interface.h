/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/Utilities/interface.h"

namespace fv3jedi {
  typedef int F90lvc_NMC;
  extern "C" {
  void fv3jedi_linvarcha_nmcbal_create_f90(const F90lvc_NMC &, const F90geom &, const F90state &,
                                           const F90state &, const eckit::LocalConfiguration &);
  void fv3jedi_linvarcha_nmcbal_delete_f90(F90lvc_NMC &);
  void fv3jedi_linvarcha_nmcbal_multiply_f90(const F90lvc_NMC &, const F90geom &, const F90inc &,
                                             const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyadjoint_f90(const F90lvc_NMC &, const F90geom &,
                                                    const F90inc &, const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyinverse_f90(const F90lvc_NMC &, const F90geom &,
                                                    const F90inc &, const F90inc &);
  void fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint_f90(const F90lvc_NMC &, const F90geom &,
                                                           const F90inc &, const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
