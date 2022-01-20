/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90iogfs;
  extern "C" {
  void fv3jedi_io_gfs_create_f90(F90iogfs &, const eckit::Configuration &, const F90geom &);
  void fv3jedi_io_gfs_delete_f90(F90iogfs &);
  void fv3jedi_io_gfs_read_state_f90(const F90iogfs &, F90state &);
  void fv3jedi_io_gfs_read_increment_f90(const F90iogfs &, F90inc &);
  void fv3jedi_io_gfs_write_state_f90(const F90iogfs &, const F90state &);
  void fv3jedi_io_gfs_write_increment_f90(const F90iogfs &, const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
