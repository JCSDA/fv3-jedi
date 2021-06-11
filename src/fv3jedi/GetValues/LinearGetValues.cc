/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/GetValues/LinearGetValues.h"

#include "fv3jedi/VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h"
#include "fv3jedi/VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

LinearGetValues::LinearGetValues(const Geometry & geom, const ufo::Locations & locs,
                                 const eckit::Configuration & conf) : locs_(locs),
  geom_(new Geometry(geom)), linearmodel2geovals_(), model2geovals_() {
  oops::Log::trace() << "LinearGetValues::LinearGetValues starting" << std::endl;

  // Create the variable change object
  {
  util::Timer timervc(classname(), "VarChaModel2GeoVaLs");
  model2geovals_.reset(new VarChaModel2GeoVaLs(geom, conf));
  }

  // Call GetValues consructor
  {
  util::Timer timergv(classname(), "LinearGetValues");

  // Pointer to configuration
  const eckit::Configuration * pconf = &conf;

  fv3jedi_lineargetvalues_create_f90(keyLinearGetValues_, geom.toFortran(), locs, &pconf);
  }
  oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

LinearGetValues::~LinearGetValues() {
  oops::Log::trace() << "LinearGetValues::~LinearGetValues starting" << std::endl;

  {
  util::Timer timergv(classname(), "~GetValues");
  fv3jedi_lineargetvalues_delete_f90(keyLinearGetValues_);
  }

  {
  util::Timer timervc(classname(), "~LinVarChaModel2GeoVaLs");
  for (lvcIter jlvc = linearmodel2geovals_.begin(); jlvc != linearmodel2geovals_.end(); ++jlvc) {
    delete jlvc->second;
  }
  }

  oops::Log::trace() << "LinearGetValues::~LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

const LinVarChaModel2GeoVaLs * LinearGetValues::getLinVarCha(const util::DateTime & t1) const {
  lvcIterCnst jlvc = linearmodel2geovals_.find(t1);
  if (jlvc == linearmodel2geovals_.end()) {
    oops::Log::error() << "LinearGetValues.getLinVarCha: linear variable change not available " <<
                          "at time " << t1 << std::endl;
    ABORT("LinearGetValues.getLinVarCha: linear variable change not available");
  }
  return jlvc->second;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::setTrajectory(const State & state, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) {
  oops::Log::trace() << "LinearGetValues::setTrajectory starting" << std::endl;

  // Create state with geovals variables
  State stategeovalvars(*geom_, geovals.getVars(), state.validTime());
  model2geovals_->changeVar(state, stategeovalvars);

  // Create the linear variable change object
  {
  util::Timer timervc(classname(), "LinearVarChaModel2GeoVaLs");
  ASSERT(linearmodel2geovals_.find(t1) == linearmodel2geovals_.end());
  char sep = '.';
  eckit::LocalConfiguration dummyconfig(sep);
  LinVarChaModel2GeoVaLs * linearmodel2geovals = new LinVarChaModel2GeoVaLs(state, state, *geom_,
                                                                            dummyconfig);
  linearmodel2geovals_[t1] = linearmodel2geovals;
  }

  {
  util::Timer timergv(classname(), "SetTrajectory");
  fv3jedi_lineargetvalues_set_trajectory_f90(keyLinearGetValues_, geom_->toFortran(),
                                             stategeovalvars.toFortran(), t1, t2, locs_,
                                             geovals.toFortran());
  }
  oops::Log::trace() << "LinearGetValues::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL starting" << std::endl;

  // Create increment with geovals variables
  Increment incgeovalvars(*geom_, geovals.getVars(), inc.validTime());

  {
  util::Timer timervc(classname(), "multiply");
  const LinVarChaModel2GeoVaLs * linearmodel2geovals = this->getLinVarCha(t1);
  linearmodel2geovals->multiply(inc, incgeovalvars);
  }

  {
  util::Timer timergv(classname(), "fillGeoVaLsTL");
  fv3jedi_lineargetvalues_fill_geovals_tl_f90(keyLinearGetValues_, geom_->toFortran(),
                                              incgeovalvars.toFortran(), t1, t2, locs_,
                                              geovals.toFortran());
  }
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsAD(Increment & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, const ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsAD starting" << std::endl;

  // Create increment with geovals variables
  Increment incgeovalvars(*geom_, geovals.getVars(), inc.validTime());

  {
  util::Timer timergv(classname(), "fillGeoVaLsAD");
  fv3jedi_lineargetvalues_fill_geovals_ad_f90(keyLinearGetValues_, geom_->toFortran(),
                                              incgeovalvars.toFortran(), t1, t2, locs_,
                                              geovals.toFortran());
  }

  // Change variables
  {
  util::Timer timervc(classname(), "multiplyAD");
  const LinVarChaModel2GeoVaLs * linearmodel2geovals = this->getLinVarCha(t1);
  linearmodel2geovals->multiplyAD(incgeovalvars, inc);
  }

  oops::Log::trace() << "LinearGetValues::fillGeovalsAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::print(std::ostream & os) const {
  os << " LinearGetValues for fv3-jedi" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
