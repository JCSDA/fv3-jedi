/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/GetValues/LinearGetValues.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

LinearGetValues::LinearGetValues(const Geometry & geom, const ufo::Locations & locs,
                                 const eckit::Configuration & conf) : locs_(locs),
  geom_(new Geometry(geom)) {
  oops::Log::trace() << "LinearGetValues::LinearGetValues starting" << std::endl;

  // Call GetValues consructor
  util::Timer timergv(classname(), "LinearGetValues");

  // Pointer to configuration
  const eckit::Configuration * pconf = &conf;

  fv3jedi_lineargetvalues_create_f90(keyLinearGetValues_, geom.toFortran(), locs, &pconf);
  oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

LinearGetValues::~LinearGetValues() {
  oops::Log::trace() << "LinearGetValues::~LinearGetValues starting" << std::endl;
  util::Timer timergv(classname(), "~GetValues");
  fv3jedi_lineargetvalues_delete_f90(keyLinearGetValues_);
  oops::Log::trace() << "LinearGetValues::~LinearGetValues done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::setTrajectory(const State & state, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) {
  oops::Log::trace() << "LinearGetValues::setTrajectory starting" << std::endl;
  fv3jedi_lineargetvalues_set_trajectory_f90(keyLinearGetValues_, geom_->toFortran(),
                                             state.toFortran(), t1, t2, locs_,
                                             geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL starting" << std::endl;
  util::Timer timergv(classname(), "fillGeoVaLsTL");
  fv3jedi_lineargetvalues_fill_geovals_tl_f90(keyLinearGetValues_, geom_->toFortran(),
                                              inc.toFortran(), t1, t2, locs_,
                                              geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::fillGeoVaLsAD(Increment & inc, const util::DateTime & t1,
                                    const util::DateTime & t2, const ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsAD starting" << std::endl;
  util::Timer timergv(classname(), "fillGeoVaLsAD");
  fv3jedi_lineargetvalues_fill_geovals_ad_f90(keyLinearGetValues_, geom_->toFortran(),
                                              inc.toFortran(), t1, t2, locs_,
                                              geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::fillGeovalsAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearGetValues::print(std::ostream & os) const {
  os << " LinearGetValues for fv3-jedi" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
