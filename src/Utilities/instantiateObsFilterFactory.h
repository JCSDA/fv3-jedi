/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_
#define SRC_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_

#include "FV3JEDITraits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/BackgroundCheck.h"
#include "ufo/BlackList.h"
#include "ufo/gnssro/QC/ROobserror.h"
#include "ufo/ObsBoundsCheck.h"
#include "ufo/ObsDomainCheck.h"
#include "ufo/ObsPreQC.h"
#include "ufo/Thinning.h"

namespace fv3jedi {

void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<FV3JEDITraits>();
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::ObsPreQC>
                          > makerChk1_("PreQC");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::ObsDomainCheck>
                          > makerChk2_("Domain Check");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::ObsBoundsCheck>
                          > makerChk3_("Bounds Check");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::BlackList>
                          > makerChk4_("BlackList");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::BackgroundCheck>
                          > makerChk5_("Background Check");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::ROobserror>
                          > makerChk6_("ROobserror");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::Thinning>
                          > makerChk7_("Thinning");
}

}  // namespace fv3jedi

#endif  // SRC_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_
