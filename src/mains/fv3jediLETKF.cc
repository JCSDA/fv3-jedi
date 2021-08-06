/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/ObsLocalization/instantiateObsLocFactory.h"
#include "fv3jedi/Utilities/Traits.h"
#include "oops/runs/LocalEnsembleDA.h"
#include "oops/runs/Run.h"
#include "saber/oops/instantiateVariableChangeFactory.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsLocFactory<fv3jedi::Traits>();
  ufo::instantiateObsErrorFactory<ufo::ObsTraits>();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  saber::instantiateVariableChangeFactory<fv3jedi::Traits>();
  fv3jedi::instantiateObsLocFactory<fv3jedi::Traits>();
  oops::LocalEnsembleDA<fv3jedi::Traits, ufo::ObsTraits> letkf;
  return run.execute(letkf);
}
