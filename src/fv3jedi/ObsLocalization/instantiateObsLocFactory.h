/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/GeometryIterator/GeometryIterator.h"
#include "fv3jedi/ObsLocalization/ObsLocVerticalBrasnett.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/obslocalization/ObsLocalizationBase.h"

namespace fv3jedi {
void instantiateObsLocFactory() {
  ufo::instantiateObsLocFactory<GeometryIterator>();
  static ufo::ObsLocalizationMaker<GeometryIterator, fv3jedi::ObsLocVerticalBrasnett>
         makerVerticalBrasnett_("Vertical Brasnett");
}

}  // namespace fv3jedi
