/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChanges/instantiateVarChangeFactories.h"
#include "oops/runs/Dirac.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"
#include "saber/oops/instantiateVariableChangeFactory.h"

int main(int argc,  char ** argv) {
  fv3jedi::Run run(argc, argv);
  fv3jedi::instantiateVarChangeFactories();
  saber::instantiateCovarFactory<fv3jedi::Traits>();
  saber::instantiateLocalizationFactory<fv3jedi::Traits>();
  saber::instantiateVariableChangeFactory<fv3jedi::Traits>();
  oops::Dirac<fv3jedi::Traits> dir;
  run.execute(dir);
  return 0;
}
