/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/Traits.h"
#include "oops/base/State.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/EnsembleGETKFApplication.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<fv3jedi::Traits>();
  oops::EnsembleForecastApplication<oops::Forecast<fv3jedi::Traits>, fv3jedi::Traits > ensfc;
  return run.execute(ensfc);
}
