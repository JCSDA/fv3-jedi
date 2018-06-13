/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef FV3_JEDI_SRC_INSTANTIATELOCALIZATIONFACTORY_H_
#define FV3_JEDI_SRC_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"
#include "LocalizationMatrixFV3JEDI.h"
#include "FV3JEDITraits.h"

namespace fv3jedi {

void instantiateLocalizationFactory() {
//  static oops::LocalizationMaker<fv3jedi::FV3JEDITraits,
//  LocalizationMatrix> makerRadiosonde_("FV3JEDI");
}

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_INSTANTIATELOCALIZATIONFACTORY_H_
