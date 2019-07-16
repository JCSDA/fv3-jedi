/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_UTILITIES_FV3JEDITRAITS_H_
#define FV3JEDI_UTILITIES_FV3JEDITRAITS_H_

#include <string>

#include "fv3jedi/ErrorCovariance/ErrorCovarianceFV3JEDI.h"
#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/GetValues/GetValuesTrajFV3JEDI.h"
#include "fv3jedi/Increment/IncrementFV3JEDI.h"
#include "fv3jedi/Localization/LocalizationMatrixFV3JEDI.h"
#include "fv3jedi/State/StateFV3JEDI.h"

#include "fv3jedi/ModelBias/ModelBiasFV3JEDI.h"
#include "fv3jedi/ModelBias/ModelBiasIncrementFV3JEDI.h"

#include "fv3jedi/ModelBias/ModelBiasCovarianceFV3JEDI.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsOperator.h"

namespace fv3jedi {

struct FV3JEDITraits {
  static std::string name() {return "FV3JEDI";}
  static std::string nameCovar() {return "FV3JEDIstatic";}
  static std::string nameCovar4D() {return "FV3JEDIstatic";}

  typedef fv3jedi::ErrorCovarianceFV3JEDI      Covariance;
  typedef fv3jedi::IncrementFV3JEDI            Increment;
  typedef fv3jedi::GeometryFV3JEDI             Geometry;
  typedef fv3jedi::GetValuesTrajFV3JEDI        InterpolatorTraj;
  typedef fv3jedi::LocalizationMatrixFV3JEDI   LocalizationMatrix;
  typedef fv3jedi::ModelBiasFV3JEDI            ModelAuxControl;
  typedef fv3jedi::ModelBiasIncrementFV3JEDI   ModelAuxIncrement;
  typedef fv3jedi::ModelBiasCovarianceFV3JEDI  ModelAuxCovariance;
  typedef fv3jedi::StateFV3JEDI                State;

  typedef ufo::GeoVaLs                         GeoVaLs;
  typedef ufo::LinearObsOperator               LinearObsOperator;
  typedef ufo::Locations                       Locations;
  typedef ufo::ObsBias                         ObsAuxControl;
  typedef ufo::ObsBiasCovariance               ObsAuxCovariance;
  typedef ufo::ObsBiasIncrement                ObsAuxIncrement;
  typedef ufo::ObsOperator                     ObsOperator;

  typedef ioda::ObsSpace                       ObsSpace;
  typedef ioda::ObsVector                      ObsVector;
  template <typename DATA> using ObsDataVector = ioda::ObsDataVector<DATA>;
};

}  // namespace fv3jedi

#endif  // FV3JEDI_UTILITIES_FV3JEDITRAITS_H_
