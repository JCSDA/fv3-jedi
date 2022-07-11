/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/IO/LatLon/IOLatLon.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static IOMaker<IOLatLon> makerIOLatLon_("latlon");
// -------------------------------------------------------------------------------------------------
IOLatLon::IOLatLon(const Geometry & geom, const Parameters_ & params) : IOBase(geom) {
  util::Timer timer(classname(), "IOLatLon");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_io_latlon_create_f90(objectKeyForFortran_, params.toConfiguration(), geom.toFortran());
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOLatLon::~IOLatLon() {
  util::Timer timer(classname(), "~IOLatLon");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_io_latlon_delete_f90(objectKeyForFortran_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOLatLon::read(State & x) const {
  ABORT("IOLatLon::read state not implemented and should not be needed");
}
// -------------------------------------------------------------------------------------------------
void IOLatLon::read(Increment & dx) const {
  ABORT("IOLatLon::read increment not implemented and should not be needed");
}
// -------------------------------------------------------------------------------------------------
void IOLatLon::write(const State & x) const {
  util::Timer timer(classname(), "write state");
  oops::Log::trace() << classname() << " write state starting" << std::endl;
  fv3jedi_io_latlon_write_state_f90(objectKeyForFortran_, x.toFortran());
  oops::Log::trace() << classname() << " write state done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOLatLon::write(const Increment & dx) const {
  util::Timer timer(classname(), "write increment");
  oops::Log::trace() << classname() << " write increment starting" << std::endl;
  fv3jedi_io_latlon_write_increment_f90(objectKeyForFortran_, dx.toFortran());
  oops::Log::trace() << classname() << " write increment done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void IOLatLon::print(std::ostream & os) const {
  os << classname() << " IO to LatLon grid";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
