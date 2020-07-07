/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChanges/instantiateVarChangeFactories.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"
#include "saber/oops/instantiateVariableChangeFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

#include "oops/runs/EDA.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  fv3jedi::instantiateVarChangeFactories();
  saber::instantiateCovarFactory<fv3jedi::Traits>();
  saber::instantiateLocalizationFactory<fv3jedi::Traits>();
  saber::instantiateVariableChangeFactory<fv3jedi::Traits>();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::EDA<fv3jedi::Traits, ufo::ObsTraits> eda;
  return run.execute(eda);
}
