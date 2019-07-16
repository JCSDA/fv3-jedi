/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Localization/instantiateLocalizationFactory.h"
#include "fv3jedi/Run/RunFV3JEDI.h"
#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "fv3jedi/Utilities/instantiateFV3JEDIVarChangeFactories.h"
#include "fv3jedi/Utilities/instantiateObsFilterFactory.h"

#include "oops/runs/Variational.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  fv3jedi::instantiateLocalizationFactory();
  fv3jedi::instantiateFV3JEDIVarChangeFactories();
  fv3jedi::instantiateObsFilterFactory();
  oops::Variational<fv3jedi::FV3JEDITraits> var;
  run.execute(var);
  return 0;
}
