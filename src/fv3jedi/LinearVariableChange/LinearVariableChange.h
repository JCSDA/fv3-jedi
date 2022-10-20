/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "vader/vader.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/LinearVariableChange/Base/LinearVariableChangeBase.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::LinearVariableChange";}

  typedef LinearVariableChangeParametersWrapper Parameters_;

  explicit LinearVariableChange(const Geometry &, const Parameters_ &);
  ~LinearVariableChange();

  void changeVarTraj(const State &, const oops::Variables &);

  void changeVarTL(Increment &, const oops::Variables &) const;
  void changeVarInverseTL(Increment &, const oops::Variables &) const;
  void changeVarAD(Increment &, const oops::Variables &) const;
  void changeVarInverseAD(Increment &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  Parameters_ params_;
  const Geometry & geom_;
  std::unique_ptr<LinearVariableChangeBase> linearVariableChange_;
  FieldsMetadata fieldsMetadata_;
  vader::Vader vader_;
  oops::Variables varsVaderPopulates_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
