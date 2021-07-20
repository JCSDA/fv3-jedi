/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_OBSLOCALIZATION_OBSLOCBRASNETT99_H_
#define FV3JEDI_OBSLOCALIZATION_OBSLOCBRASNETT99_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/util/Earth.h"

#include "eckit/config/Configuration.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/UnitSphere.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/soar.h"
#include "oops/util/Printable.h"

#include "ufo/obslocalization/ObsLocParameters.h"
#include "ufo/obslocalization/ObsLocSOAR.h"


namespace fv3jedi {

/// Brasnett 99 observation space localization for snow DA
/// https://doi.org/10.1175/1520-0450(1999)038<0726:AGAOSD>2.0.CO;2
/// Note, Brasnett99 adds vertical localization to the horizontal SOAR function
/// Hence, we inherit from ufo::ObsLocSOAR
template<class MODEL>
class ObsLocBrasnett99: public ufo::ObsLocSOAR<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  ObsLocBrasnett99(const eckit::Configuration &, const ioda::ObsSpace &);

  /// compute localization and save localization values in \p locvector
  /// (missing values for observations outside of localization)
  void computeLocalization(const GeometryIterator_ &,
                           ioda::ObsVector & locvector) const override;

 private:
  void print(std::ostream &) const;
  std::vector<float> obsAltitude_;  //< altitude of observations
  double VertScale_;  //< vertical localization scale
};
// -----------------------------------------------------------------------------
template<typename MODEL>
ObsLocBrasnett99<MODEL>::ObsLocBrasnett99(const eckit::Configuration & config,
                              const ioda::ObsSpace & obsspace):
       ufo::ObsLocSOAR<MODEL>::ObsLocSOAR(config, obsspace),
       obsAltitude_(obsspace.nlocs()),
       VertScale_(config.getDouble("vertical lengthscale", 800)) {
  oops::Log::trace()<< "Brasnett99 localization with: vertical scale=" << VertScale_
                    << std::endl;

  // read altitude of measurements
  obsspace.get_db("MetaData", "altitude", obsAltitude_);
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocBrasnett99<MODEL>::computeLocalization(const GeometryIterator_ & geoiter,
                                                  ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ObsLocBrasnett99::computeLocalization" << std::endl;
  // compute horizontal localization using SOAR
  ufo::ObsLocSOAR<MODEL>::computeLocalization(geoiter, locvector);
  const std::vector<int> & localobs = ufo::ObsLocalization<MODEL>::localobs();

  // retrieve orography for this grid point
  double orog = geoiter.getOrography();
  oops::Log::debug() << "geoiter=" << geoiter << " orog=" << orog << std::endl;

  // compute vertical localization and multiply it by SOAR computed above
  // vloc=exp(- (dz/hfac)^2 )
  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs.size(); ++jlocal) {
    double locFactor = std::exp(-std::pow((obsAltitude_[localobs[jlocal]] - orog)
                                          /VertScale_, 2));
    // obsdist is calculated at each location; need to update R for each variable
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs[jlocal] * nvars] *= locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocBrasnett99<MODEL>::print(std::ostream & os) const {
//  this->print(os);
  os << "Brasnett99 localization with: vertical scale=" << VertScale_
                    << std::endl;
}

}  // namespace fv3jedi

#endif  // FV3JEDI_OBSLOCALIZATION_OBSLOCBRASNETT99_H_
