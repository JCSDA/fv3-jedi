/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_INSTANTIATEOBSFILTERFACTORY_H_
#define SRC_INSTANTIATEOBSFILTERFACTORY_H_

#include "FV3JEDITraits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/BackgroundCheck.h"

namespace fv3jedi {

void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<FV3JEDITraits>();
  static oops::FilterMaker<FV3JEDITraits,
                           oops::ObsFilter<FV3JEDITraits, ufo::BackgroundCheck>
                          >
    makerBkgChk_("Background Check");
}

}  // namespace fv3jedi

#endif  // SRC_INSTANTIATEOBSFILTERFACTORY_H_
