/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef FV3JEDI_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_
#define FV3JEDI_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"

#include "fv3jedi/Localization/LocalizationMatrix.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {

void instantiateLocalizationFactory() {
//  static oops::LocalizationMaker<fv3jedi::Traits,
//  LocalizationMatrix> makerRadiosonde_("FV3JEDI");
}

}  // namespace fv3jedi

#endif  // FV3JEDI_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_
