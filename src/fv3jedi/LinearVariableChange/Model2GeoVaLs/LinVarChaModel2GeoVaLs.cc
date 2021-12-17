/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaModel2GeoVaLs> makerLinVarChaMod2Geo_("Model2GeoVaLs");
static LinearVariableChangeMaker<LinVarChaModel2GeoVaLs> makerDefault_("default");
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVaLs::LinVarChaModel2GeoVaLs(const State & bg, const State & fg,
                                    const Geometry & resol, const eckit::LocalConfiguration & conf)
  : LinearVariableChangeBase(), geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "LinVarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_lvc_model2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                       fg.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaModel2GeoVaLs::~LinVarChaModel2GeoVaLs() {
  util::Timer timer(classname(), "~LinVarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_lvc_model2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiply(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  fv3jedi_lvc_model2geovals_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                         dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiply done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyInverse(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " multiplyInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  fv3jedi_lvc_model2geovals_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::multiplyInverseAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  dxout = dxin;
  oops::Log::trace() << classname() << " multiplyInverseAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
