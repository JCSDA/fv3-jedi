/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_UTILITIES_INSTANTIATEFV3JEDIVARCHANGEFACTORIES_H_
#define FV3JEDI_UTILITIES_INSTANTIATEFV3JEDIVARCHANGEFACTORIES_H_

#include "fv3jedi/Analysis2Model/LinVarChaA2MFV3JEDI.h"
#include "fv3jedi/Analysis2Model/VarChaA2MFV3JEDI.h"
#include "fv3jedi/Control2Analysis/LinVarChaC2AFV3JEDI.h"
#include "fv3jedi/Utilities/FV3JEDITraits.h"

#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"

namespace fv3jedi {

void instantiateFV3JEDIVarChangeFactories() {
  static oops::LinearVariableChangeMaker<fv3jedi::FV3JEDITraits,
               oops::LinearVariableChange<fv3jedi::FV3JEDITraits,
               fv3jedi::LinVarChaC2AFV3JEDI> >
                   makerLinVarChaC2AV3JEDI_("Control2Analysis");
  static oops::VariableChangeMaker<fv3jedi::FV3JEDITraits,
               oops::VariableChange<fv3jedi::FV3JEDITraits,
               fv3jedi::VarChaA2MFV3JEDI> >
                   makerVarChaA2MV3JEDI_("Analysis2Model");
  static oops::LinearVariableChangeMaker<fv3jedi::FV3JEDITraits,
               oops::LinearVariableChange<fv3jedi::FV3JEDITraits,
               fv3jedi::LinVarChaA2MFV3JEDI> >
                   makerLinVarChaA2MV3JEDI_("Analysis2Model");
}

}  // namespace fv3jedi

#endif  // FV3JEDI_UTILITIES_INSTANTIATEFV3JEDIVARCHANGEFACTORIES_H_
