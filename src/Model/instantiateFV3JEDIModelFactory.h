/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_MODEL_INSTANTIATEFV3JEDIMODELFACTORY_H_
#define SRC_MODEL_INSTANTIATEFV3JEDIMODELFACTORY_H_

#include "ModelFV3FV3JEDI.h"
/*#include "ModelGEOSFV3JEDI.h"
#include "ModelGFSFV3JEDI.h"
#include "ModelSARFV3JEDI.h"*/
#include "FV3JEDITraits.h"
#include "oops/interface/Model.h"

namespace fv3jedi {

void instantiateFV3JEDIModelFactory() {
  static oops::ModelMaker<fv3jedi::FV3JEDITraits,
         oops::Model<fv3jedi::FV3JEDITraits, fv3jedi::ModelFV3FV3JEDI> >
    makerModelFV3JEDI_("FV3");
}

}  // namespace fv3jedi

#endif  // SRC_MODEL_INSTANTIATEFV3JEDIMODELFACTORY_H_
