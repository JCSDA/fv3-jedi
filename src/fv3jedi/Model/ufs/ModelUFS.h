/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_UFS_MODELUFS_H_
#define FV3JEDI_MODEL_UFS_MODELUFS_H_

#include <ostream>
#include <string>

#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Utilities/Traits.h"
#include "ModelUFS.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -----------------------------------------------------------------------------

class ModelUFS: public oops::ModelBase<Traits>,
                        private util::ObjectCounter<ModelUFS> {
 public:
  static const std::string classname() {return "fv3jedi::ModelUFS";}

  ModelUFS(const Geometry &, const eckit::Configuration &);
  ~ModelUFS();

  void initialize(State &) const;
  void step(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;
  int saveTrajectory(State &, const ModelBias &) const;

  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
  char jedidir_[10000];
  char ufsdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_UFS_MODELUFS_H_
