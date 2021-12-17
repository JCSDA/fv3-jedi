/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "fv3jedi/ErrorCovariance/ErrorCovariance.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.h"
#include "fv3jedi/GetValues/GetValues.h"
#include "fv3jedi/GetValues/LinearGetValues.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/LinearVariableChange/LinearVariableChange.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/VariableChange/VariableChange.h"

#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/ModelBias/ModelBiasIncrement.h"

#include "fv3jedi/ModelBias/ModelBiasCovariance.h"

namespace fv3jedi {

struct Traits {
  static std::string name() {return "FV3JEDI";}
  static std::string nameCovar() {return "FV3JEDI-ID";}

  typedef fv3jedi::ErrorCovariance      Covariance;
  typedef fv3jedi::Increment            Increment;
  typedef fv3jedi::Geometry             Geometry;
  typedef fv3jedi::GeometryIterator     GeometryIterator;
  typedef fv3jedi::GetValues            GetValues;
  typedef fv3jedi::LinearGetValues      LinearGetValues;
  typedef fv3jedi::LinearVariableChange LinearVariableChange;
  typedef fv3jedi::ModelBias            ModelAuxControl;
  typedef fv3jedi::ModelBiasIncrement   ModelAuxIncrement;
  typedef fv3jedi::ModelBiasCovariance  ModelAuxCovariance;
  typedef fv3jedi::State                State;
  typedef fv3jedi::VariableChange       VariableChange;
};

}  // namespace fv3jedi
