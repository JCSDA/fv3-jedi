/*
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/IO/CubeSphereHistory/IOCubeSphereHistory.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static IOMaker<IOCubeSphereHistory> makerIOCubeSphereHistory_("cube sphere history");
// -------------------------------------------------------------------------------------------------
IOCubeSphereHistory::IOCubeSphereHistory(const Geometry & geom, const Parameters_ & params)
  : IOBase(geom) {
  util::Timer timer(classname(), "IOCubeSphereHistory");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_io_cube_sphere_history_create_f90(objectKeyForFortran_, params.toConfiguration(),
                                            geom.toFortran());
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOCubeSphereHistory::~IOCubeSphereHistory() {
  util::Timer timer(classname(), "~IOCubeSphereHistory");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_io_cube_sphere_history_delete_f90(objectKeyForFortran_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOCubeSphereHistory::read(State & x) const {
  util::Timer timer(classname(), "read state");
  oops::Log::trace() << classname() << " read state starting" << std::endl;
  fv3jedi_io_cube_sphere_history_read_state_f90(objectKeyForFortran_, x.toFortran());
  oops::Log::trace() << classname() << " read state done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOCubeSphereHistory::read(Increment & dx) const {
  util::Timer timer(classname(), "read increment");
  oops::Log::trace() << classname() << " read increment starting" << std::endl;
  fv3jedi_io_cube_sphere_history_read_increment_f90(objectKeyForFortran_, dx.toFortran());
  oops::Log::trace() << classname() << " read increment done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOCubeSphereHistory::write(const State & x) const {
  util::Timer timer(classname(), "write state");
  oops::Log::trace() << classname() << " write state starting" << std::endl;
  fv3jedi_io_cube_sphere_history_write_state_f90(objectKeyForFortran_, x.toFortran());
  oops::Log::trace() << classname() << " write state done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOCubeSphereHistory::write(const Increment & dx) const {
  util::Timer timer(classname(), "write increment");
  oops::Log::trace() << classname() << " write increment starting" << std::endl;
  fv3jedi_io_cube_sphere_history_write_increment_f90(objectKeyForFortran_, dx.toFortran());
  oops::Log::trace() << classname() << " write increment done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOCubeSphereHistory::print(std::ostream & os) const {
  os << classname() << " IO for Cube Sphere Histories";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
