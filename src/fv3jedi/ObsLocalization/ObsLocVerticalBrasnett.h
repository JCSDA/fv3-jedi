/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <vector>

#include "fv3jedi/GeometryIterator/GeometryIterator.h"
#include "fv3jedi/Utilities/Traits.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/ObsLocalizationBase.h"
#include "oops/util/missingValues.h"

#include "ufo/ObsTraits.h"

namespace fv3jedi {

/// \brief Options controlling vertical Brasnett observation space localization
/// for snow DA.
class ObsLocBrasnettParameters : public oops::ObsLocalizationParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsLocBrasnettParameters, oops::ObsLocalizationParametersBase)

 public:
  oops::Parameter<double> vertscale{"vertical lengthscale",
                 "lengthscale of vertical localization in meters", 800., this};
};

/// Brasnett 99 observation space localization for snow DA (in vertical).
/// https://doi.org/10.1175/1520-0450(1999)038<0726:AGAOSD>2.0.CO;2
class ObsLocVerticalBrasnett: public oops::ObsLocalizationBase<Traits, ufo::ObsTraits> {
 public:
  typedef ObsLocBrasnettParameters Parameters_;
  ObsLocVerticalBrasnett(const Parameters_ &, const ioda::ObsSpace &);

 protected:
  /// compute localization and update localization values in \p locvector
  /// (missing values for observations outside of localization)
  void computeLocalization(const GeometryIterator &,
                           ioda::ObsVector & locvector) const override;

 private:
  void print(std::ostream &) const override;
  std::vector<float> obsHeight_;  //< height of observations
  double VertScale_;  //< vertical localization scale
};
// -----------------------------------------------------------------------------
ObsLocVerticalBrasnett::ObsLocVerticalBrasnett(const Parameters_ & params,
                                               const ioda::ObsSpace & obsspace):
       obsHeight_(obsspace.nlocs()),
       VertScale_(params.vertscale) {
  oops::Log::trace()<< "VerticalBrasnett localization with: vertical scale=" << VertScale_
                    << std::endl;

  // read height of measurements
  obsspace.get_db("MetaData", "height", obsHeight_);
}
// -----------------------------------------------------------------------------

void ObsLocVerticalBrasnett::computeLocalization(const GeometryIterator & geoiter,
                                                  ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsLocVerticalBrasnett::computeLocalization" << std::endl;

  // retrieve orography for this grid point
  double orog = geoiter.getOrography();
  oops::Log::debug() << "geoiter=" << geoiter << " orog=" << orog << std::endl;

  // compute vertical localization and multiply it by the previously computed localization
  // vloc=exp(- (dz/hfac)^2 )
  const size_t nvars = locvector.nvars();
  const double missing = util::missingValue<double>();
  for (size_t jloc = 0; jloc < locvector.nlocs(); ++jloc) {
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      if (locvector[jvar + jloc * nvars] != missing) {
        double locFactor = std::exp(-std::pow((obsHeight_[jloc] - orog)
                                              /VertScale_, 2));
        locvector[jvar + jloc * nvars] *= locFactor;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ObsLocVerticalBrasnett::print(std::ostream & os) const {
  os << "VerticalBrasnett localization with: vertical scale=" << VertScale_
     << std::endl;
}

}  // namespace fv3jedi
