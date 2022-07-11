/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90vc_CSW;
  extern "C" {
    void fv3jedi_vc_coldstartwinds_create_f90(const F90vc_CSW &, const F90geom &,
                                              const eckit::LocalConfiguration &);
    void fv3jedi_vc_coldstartwinds_delete_f90(F90vc_CSW &);
    void fv3jedi_vc_coldstartwinds_changevar_f90(const F90vc_CSW &, const F90state &,
                                                 const F90state &);
  }  // extern "C"
}  // namespace fv3jedi
