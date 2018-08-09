/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_FV3JEDITRAITS_H_
#define FV3_JEDI_SRC_FV3JEDITRAITS_H_

#include <string>

#include "ChangeVarFV3JEDI.h"
#include "ErrorCovarianceFV3JEDI.h"
#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "LocalizationMatrixFV3JEDI.h"
#include "ModelFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "ModelBiasIncrementFV3JEDI.h"
#include "ModelBiasCovarianceFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsCheck.h"
#include "ufo/ObsOperator.h"
#include "ufo/LinearObsOperator.h"
#include "GetValuesTrajFV3JEDI.h"

namespace fv3jedi {

struct FV3JEDITraits {
  static std::string name() {return "FV3JEDI";}
  static std::string nameCovar() {return "FV3JEDIstatic";}

  typedef fv3jedi::GeometryFV3JEDI             Geometry;

  typedef fv3jedi::StateFV3JEDI                State;
  typedef fv3jedi::ModelFV3JEDI                Model;
  typedef fv3jedi::IncrementFV3JEDI            Increment;
  typedef fv3jedi::ErrorCovarianceFV3JEDI      Covariance;
  typedef fv3jedi::ChangeVarFV3JEDI            ChangeVar;

  typedef fv3jedi::ModelBiasFV3JEDI            ModelAuxControl;
  typedef fv3jedi::ModelBiasIncrementFV3JEDI   ModelAuxIncrement;
  typedef fv3jedi::ModelBiasCovarianceFV3JEDI  ModelAuxCovariance;

  typedef fv3jedi::LocalizationMatrixFV3JEDI   LocalizationMatrix;

  typedef fv3jedi::GetValuesTrajFV3JEDI        InterpolatorTraj;

  typedef ufo::ObsOperator                     ObsOperator;
  typedef ufo::LinearObsOperator               LinearObsOperator;
  typedef ioda::ObsSpace                       ObsSpace;
  typedef ioda::ObsVector                      ObsVector;

  typedef ufo::ObsBias                         ObsAuxControl;
  typedef ufo::ObsBiasIncrement                ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance               ObsAuxCovariance;
  typedef ufo::ObsCheck                        ObsCheck;

  typedef ufo::GeoVaLs                         GeoVaLs;
  typedef ioda::Locations                      Locations;
};

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_FV3JEDITRAITS_H_
