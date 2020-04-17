/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_VARIABLECHANGES_INSTANTIATEVARCHANGEFACTORIES_H_
#define FV3JEDI_VARIABLECHANGES_INSTANTIATEVARCHANGEFACTORIES_H_

#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChanges/Analysis2Model/LinVarChaA2M.h"
#include "fv3jedi/VariableChanges/Analysis2Model/VarChaA2M.h"
#include "fv3jedi/VariableChanges/Control2Analysis/LinVarChaC2A.h"
#include "fv3jedi/VariableChanges/Control2Analysis/VarChaC2A.h"
#include "fv3jedi/VariableChanges/GEOSRstToBkg/VarChaGeosRst2Bkg.h"
#include "fv3jedi/VariableChanges/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h"
#include "fv3jedi/VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h"
#include "fv3jedi/VariableChanges/NMCBalance/LinVarChaNMCBal.h"

#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"

namespace fv3jedi {

void instantiateVarChangeFactories() {
  // Control to Analysis
  static oops::VariableChangeMaker<fv3jedi::Traits,
               oops::VariableChange<fv3jedi::Traits,
               fv3jedi::VarChaC2A> >
                   makerVarChaC2A_("Control2Analysis");
  static oops::LinearVariableChangeMaker<fv3jedi::Traits,
               oops::LinearVariableChange<fv3jedi::Traits,
               fv3jedi::LinVarChaC2A> >
                   makerLinVarChaC2A_("Control2Analysis");
  // NMC Balance
  static oops::LinearVariableChangeMaker<fv3jedi::Traits,
               oops::LinearVariableChange<fv3jedi::Traits,
               fv3jedi::LinVarChaNMCBal> >
                   makerLinVarChaNMCBal_("NMCBalance");
  // Analysis to Model
  static oops::VariableChangeMaker<fv3jedi::Traits,
               oops::VariableChange<fv3jedi::Traits,
               fv3jedi::VarChaA2M> >
                   makerVarChaA2M_("Analysis2Model");
  static oops::LinearVariableChangeMaker<fv3jedi::Traits,
               oops::LinearVariableChange<fv3jedi::Traits,
               fv3jedi::LinVarChaA2M> >
                   makerLinVarChaA2M_("Analysis2Model");
  // GEOS Restart to GEOS eta (background)
  static oops::VariableChangeMaker<fv3jedi::Traits,
              oops::VariableChange<fv3jedi::Traits,
              fv3jedi::VarChaGeosRst2Bkg> >
                  makerVarChaGeosRst2Bkg_("GeosRst2Bkg");
  // Model to GeoVaLs
  static oops::VariableChangeMaker<fv3jedi::Traits,
              oops::VariableChange<fv3jedi::Traits,
              fv3jedi::VarChaModel2GeoVaLs> >
                  makerVarChaModel2GeoVaLs_("Model2GeoVaLs");
  static oops::LinearVariableChangeMaker<fv3jedi::Traits,
              oops::LinearVariableChange<fv3jedi::Traits,
              fv3jedi::LinVarChaModel2GeoVaLs> >
                  makerLinVarChaModel2GeoVaLs_("Model2GeoVaLs");
}

}  // namespace fv3jedi

#endif  // FV3JEDI_VARIABLECHANGES_INSTANTIATEVARCHANGEFACTORIES_H_
