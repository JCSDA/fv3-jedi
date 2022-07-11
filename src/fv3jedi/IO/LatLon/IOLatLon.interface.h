/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90iolatlon;
  extern "C" {
  void fv3jedi_io_latlon_create_f90(F90iolatlon &, const eckit::Configuration &, const F90geom &);
  void fv3jedi_io_latlon_delete_f90(F90iolatlon &);
  void fv3jedi_io_latlon_write_state_f90(const F90iolatlon &, const F90state &);
  void fv3jedi_io_latlon_write_increment_f90(const F90iolatlon &, const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
