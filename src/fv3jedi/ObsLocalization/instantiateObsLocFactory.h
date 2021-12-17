/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/ObsLocalization/ObsLocBrasnett99.h"
#include "fv3jedi/Utilities/Traits.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/interface/ObsLocalization.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"

namespace fv3jedi {
template<typename MODEL> void instantiateObsLocFactory() {
  ufo::instantiateObsLocFactory<MODEL>();

  static oops::ObsLocalizationMaker<MODEL, ufo::ObsTraits, fv3jedi::ObsLocBrasnett99<MODEL>>
           makerBrasnett99_("Brasnett99");
}

}  // namespace fv3jedi
