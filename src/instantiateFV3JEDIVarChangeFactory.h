  1 /*
  2  * (C) Copyright 2017-2018  UCAR.
  3  * 
  4  * This software is licensed under the terms of the Apache Licence Version 2.0
  5  * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
  6  */

#ifndef FV3JEDI_SRC_INSTANTIATEFV3JEDICHANGEVARFACTORY_H_
#define FV3JEDI_SRC_INSTANTIATEFV3JEDICHANGEVARFACTORY_H_

#include "model/ChangeVar.h"
#include "model/FV3JEDITraits.h"
#include "oops/interface/VariableChange.h"

namespace fv3jedi {

void instantiateFV3JEDIChangeVarFactory() {
  static oops::VariableChangeMaker<fv3jedi::FV3JEDITraits,
         oops::VariableChange<fv3jedi::FV3JEDITraits, fv3jedi::ChangeVar> >
         makerChangeVarFV3JEDI_("ChVarFV3JEDI");
}

}  // namespace fv3jedi

#endif  // FV3JEDI_SRC_INSTANTIATEFV3JEDICHANGEVARFACTORY_H_
