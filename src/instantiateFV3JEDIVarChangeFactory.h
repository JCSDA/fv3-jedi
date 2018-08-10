/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_INSTANTIATEFV3JEDIVARCHANGEFACTORY_H_
#define SRC_INSTANTIATEFV3JEDIVARCHANGEFACTORY_H_

#include "src/VarChangeFV3JEDI.h"
#include "src/FV3JEDITraits.h"
#include "oops/interface/LinearVariableChange.h"

namespace fv3jedi {

void instantiateFV3JEDIVarChangeFactory() {
  static oops::LinearVariableChangeMaker<fv3jedi::FV3JEDITraits,
         oops::LinearVariableChange<fv3jedi::FV3JEDITraits,
                                    fv3jedi::VarChangeFV3JEDI> >
    makerVarChangeFV3JEDI_("VarChangeFV3JEDI");
}

}  // namespace fv3jedi

#endif  // SRC_INSTANTIATEFV3JEDIVARCHANGEFACTORY_H_
