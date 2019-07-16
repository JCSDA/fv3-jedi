/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <unistd.h>

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Model/geos/ModelGEOSFV3JEDI.h"
#include "fv3jedi/ModelBias/ModelBiasFV3JEDI.h"
#include "fv3jedi/State/StateFV3JEDI.h"
#include "fv3jedi/Utilities/UtilitiesFV3JEDI.h"
#include "ModelGEOSFV3JEDIFortran.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelGEOSFV3JEDI> makermodel_("GEOS");
// -----------------------------------------------------------------------------
ModelGEOSFV3JEDI::ModelGEOSFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf)
{
  oops::Log::trace() << "ModelGEOSFV3JEDI::ModelGEOSFV3JEDI" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;

  // JEDI to GEOS directory
  getcwd(jedidir_, 10000);

  std::string sGEOSSCRDIR = mconf.getString("GEOSSCRDIR");
  strcpy(geosscrdir_, sGEOSSCRDIR.c_str());
  chdir(geosscrdir_);

  // Create the model
  fv3jedi_geos_create_f90(&configc, geom_.toFortran(), keyConfig_);

  // GEOS to JEDI directory
  chdir(jedidir_);

  oops::Log::trace() << "ModelGEOSFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelGEOSFV3JEDI::~ModelGEOSFV3JEDI() {
  chdir(geosscrdir_);
  fv3jedi_geos_delete_f90(keyConfig_);
  chdir(jedidir_);
  oops::Log::trace() << "ModelGEOSFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::initialize(StateFV3JEDI & xx) const {
  chdir(geosscrdir_);
  fv3jedi_geos_initialize_f90(keyConfig_, xx.toFortran());
  chdir(jedidir_);
  oops::Log::debug() << "ModelGEOSFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  chdir(geosscrdir_);
  fv3jedi_geos_step_f90(keyConfig_, xx.toFortran(), &dtp);
  chdir(jedidir_);
  oops::Log::debug() << "ModelGEOSFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::finalize(StateFV3JEDI & xx) const {
  chdir(geosscrdir_);
  fv3jedi_geos_finalize_f90(keyConfig_, xx.toFortran());
  chdir(jedidir_);
  oops::Log::debug() << "ModelGEOSFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelGEOSFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  ABORT("Model:GEOS should not be used for the trajecotry");
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::print(std::ostream & os) const {
  os << "ModelGEOSFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
