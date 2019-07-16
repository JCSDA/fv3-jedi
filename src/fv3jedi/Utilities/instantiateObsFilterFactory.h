/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_
#define FV3JEDI_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_

#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/filters/BackgroundCheck.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/Thinning.h"
#include "ufo/gnssro/QC/ROobserror.h"

namespace fv3jedi {

void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<FV3JEDITraits>();
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::QCmanager>
                          > makerChk0_("QCmanager");
  static oops::FilterMaker<FV3JEDITraits,
                 oops::ObsFilter<FV3JEDITraits, ufo::PreQC>
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

#endif  // FV3JEDI_UTILITIES_INSTANTIATEOBSFILTERFACTORY_H_
