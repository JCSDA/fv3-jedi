/*
 * (C) Copyright 2022 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90ioauxgrid;
  extern "C" {
  void fv3jedi_io_auxgrid_create_f90(F90ioauxgrid &, const eckit::Configuration &, const F90geom &);
  void fv3jedi_io_auxgrid_delete_f90(F90ioauxgrid &);
  void fv3jedi_io_auxgrid_write_state_f90(const F90ioauxgrid &, const F90state &);
  void fv3jedi_io_auxgrid_write_increment_f90(const F90ioauxgrid &, const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
