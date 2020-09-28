/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/interface/VariableChange.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/Utilities/Utilities.h"
#include "fv3jedi/VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::VariableChangeMaker<Traits, oops::VariableChange<Traits, VarChaModel2GeoVaLs>>
       makerVarChaModel2GeoVaLs_("Model2GeoVaLs");
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::VarChaModel2GeoVaLs(const Geometry & geom, const eckit::Configuration & conf) :
  geom_(new Geometry(geom))
{
  util::Timer timer(classname(), "VarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  const eckit::Configuration * configc = &conf;
  fv3jedi_vc_model2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::~VarChaModel2GeoVaLs() {
  util::Timer timer(classname(), "~VarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_vc_model2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
  fv3jedi_vc_model2geovals_changevar_f90(keyFtnConfig_, geom_->toFortran(), xin.toFortran(),
                                         xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  xout = xin;
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
