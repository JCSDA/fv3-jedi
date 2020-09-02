/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/VertRemap/VarChaVertRemap.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
VarChaVertRemap::VarChaVertRemap(const Geometry & resol, const eckit::Configuration & conf):
    geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "VarChaVertRemap");
  oops::Log::trace() << "VarChaVertRemap::VarChaVertRemap start" << std::endl;
  const eckit::Configuration * configc = &conf;

  // Prepare input.nml file
  stageFv3Files(conf, resol.getComm());
  if ( !conf.has("nml_file") ) {
    generateGeomFv3Conf(conf, resol.getComm());
  }
  fv3jedi_vc_vertremap_create_f90(keyFtn_, geom_->toFortran(), &configc);
  // Remove input.nml
  removeFv3Files(resol.getComm());
  oops::Log::trace() << "VarChaVertRemap::VarChaVertRemap done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaVertRemap::~VarChaVertRemap() {
  util::Timer timer(classname(), "~VarChaVertRemap");
  oops::Log::trace() << "VarChaVertRemap::~VarChaVertRemap start" << std::endl;
  fv3jedi_vc_vertremap_delete_f90(keyFtn_);
  oops::Log::trace() << "VarChaVertRemap::~VarChaVertRemap done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaVertRemap::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << "VarChaVertRemap::changeVar starting" << xin << std::endl;
  fv3jedi_vc_vertremap_changevar_f90(keyFtn_, xin.toFortran(), xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << "VarChaVertRemap::changeVar done" << xout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaVertRemap::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << "VarChaVertRemap::changeVarInverse starting" << xin << std::endl;
  xout = xin;  // No inverse required
  xout.validTime() = xin.validTime();
  oops::Log::trace() << "VarChaVertRemap::changeVarInverse done" << xout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaVertRemap::print(std::ostream & os) const {
  os << "VarChaVertRemap";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
