/*
 * (C) Copyright 2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVaLs::LinVarChaModel2GeoVaLs(const State & bg, const State & fg,
                                         const Geometry & resol, const eckit::Configuration & conf):
  geom_(new Geometry(resol))
{
  const eckit::Configuration * configc = &conf;
  fv3jedi_lvc_model2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                       fg.toFortran(), &configc);
  oops::Log::trace() << "LinVarChaModel2GeoVaLs created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVaLs::~LinVarChaModel2GeoVaLs() {
  fv3jedi_lvc_model2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "LinVarChaModel2GeoVaLs destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiply(const Increment & dxm, Increment & dxg) const {
  oops::Log::trace() << "LinVarChaModel2GeoVaLs multiply" << std::endl;
  fv3jedi_lvc_model2geovals_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                         dxm.toFortran(), dxg.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyInverse(const Increment & dxg, Increment & dxm) const {
  oops::Log::trace() << "LinVarChaModel2GeoVaLs multiplyInverse" << std::endl;
  fv3jedi_lvc_model2geovals_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxg.toFortran(), dxm.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyAD(const Increment & dxg, Increment & dxm) const {
  oops::Log::trace() << "LinVarChaModel2GeoVaLs multiplyAD" << std::endl;
  fv3jedi_lvc_model2geovals_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                dxg.toFortran(), dxm.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyInverseAD(const Increment & dxm, Increment & dxg) const {
  oops::Log::trace() << "LinVarChaModel2GeoVaLs multiplyInverseAD" << std::endl;
  fv3jedi_lvc_model2geovals_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                       dxm.toFortran(), dxg.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
