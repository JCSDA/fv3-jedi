/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90iofms;
  extern "C" {
  void fv3jedi_io_fms_create_f90(F90iofms &, const eckit::Configuration &, const F90geom &);
  void fv3jedi_io_fms_delete_f90(F90iofms &);
  void fv3jedi_io_fms_read_state_f90(const F90iofms &, F90state &);
  void fv3jedi_io_fms_read_increment_f90(const F90iofms &, F90inc &);
  void fv3jedi_io_fms_write_state_f90(const F90iofms &, const F90state &);
  void fv3jedi_io_fms_write_increment_f90(const F90iofms &, const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
