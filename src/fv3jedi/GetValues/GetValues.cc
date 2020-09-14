/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/GetValues/GetValues.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

GetValues::GetValues(const Geometry & geom, const ufo::Locations & locs) : locs_(locs),
  geom_(new Geometry(geom)), model2geovals_() {
  oops::Log::trace() << "GetValues::GetValues starting" << std::endl;

  // Create the variable change object
  {
  util::Timer timervc(classname(), "VarChaModel2GeoVaLs");
  char sep = '.';
  eckit::LocalConfiguration dummyconfig(sep);
  model2geovals_.reset(new VarChaModel2GeoVaLs(geom, dummyconfig));
  }

  // Call GetValues consructor
  {
  util::Timer timergv(classname(), "GetValues");
  fv3jedi_getvalues_create_f90(keyGetValues_, geom.toFortran(), locs_.toFortran());
  oops::Log::trace() << "GetValues::GetValues done" << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

GetValues::~GetValues() {
  oops::Log::trace() << "GetValues::~GetValues starting" << std::endl;
  fv3jedi_getvalues_delete_f90(keyGetValues_);
  oops::Log::trace() << "GetValues::~GetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void GetValues::fillGeoVaLs(const State & state, const util::DateTime & t1,
                            const util::DateTime & t2, ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "GetValues::fillGeovals starting" << std::endl;
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;

  // Create state with geovals variables
  State stategeovalvars(*geom_, geovals.getVars(), state.validTime());

  {
  util::Timer timervc(classname(), "changeVar");
  model2geovals_->changeVar(state, stategeovalvars);
  }

  oops::Log::trace() << "GetValues::fillGeovals changeVar done" << stategeovalvars << std::endl;

  // Fill GeoVaLs
  util::Timer timergv(classname(), "fillGeoVaLs");
  fv3jedi_getvalues_fill_geovals_f90(keyGetValues_, geom_->toFortran(),
                                     stategeovalvars.toFortran(), &t1p, &t2p, locs_.toFortran(),
                                     geovals.toFortran());

  oops::Log::trace() << "GetValues::fillGeovals done" << geovals << std::endl;
}

// -------------------------------------------------------------------------------------------------

void GetValues::print(std::ostream & os) const {
  os << " GetValues for fv3-jedi" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
