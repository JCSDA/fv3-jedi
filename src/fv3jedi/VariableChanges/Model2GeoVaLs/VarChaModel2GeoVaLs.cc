/*
 * (C) Copyright 2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::VarChaModel2GeoVaLs(const Geometry & geom, const eckit::Configuration & conf) :
  geom_(new Geometry(geom)) {
  oops::Log::trace() << "VarChaModel2GeoVaLs::VarChaModel2GeoVaLs start" << std::endl;
  const eckit::Configuration * configc = &conf;
  fv3jedi_vc_model2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaModel2GeoVaLs::VarChaModel2GeoVaLs done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::~VarChaModel2GeoVaLs() {
  fv3jedi_vc_model2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVar(const State & xm, State & xg) const {
  oops::Log::trace() << "VarChaModel2GeoVaLs::changeVar starting" << xg << std::endl;
  fv3jedi_vc_model2geovals_changevar_f90(keyFtnConfig_, geom_->toFortran(), xm.toFortran(),
                                         xg.toFortran());
  xg.validTime() = xm.validTime();
  oops::Log::trace() << "VarChaModel2GeoVaLs::changeVar done" << xg << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVarInverse(const State & xg, State & xm) const {
  oops::Log::trace() << "VarChaModel2GeoVaLs::changeVarInverse starting" << xg << std::endl;
  fv3jedi_vc_model2geovals_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(), xg.toFortran(),
                                                xm.toFortran());
  xm.validTime() = xg.validTime();
  oops::Log::trace() << "VarChaModel2GeoVaLs::changeVarInverse done" << xg <<
                        std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << "VarChaModel2GeoVaLs";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
