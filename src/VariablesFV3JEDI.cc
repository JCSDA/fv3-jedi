/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "VariablesFV3JEDI.h"

#include<vector>

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------

VariablesFV3JEDI::VariablesFV3JEDI(const oops::Variables & oopsvars) {
  oops::Log::debug() << "VariablesFV3JEDI oopsvar:" << oopsvars.variables() << std::endl;
  this->setF90(oopsvars.variables());
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

VariablesFV3JEDI::VariablesFV3JEDI(const eckit::Configuration & config) {
  oops::Log::debug() << "VariablesFV3JEDI config:" << config << std::endl;
  std::vector<std::string> vars;
  config.get("variables", vars);
  this->setF90(vars);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

void VariablesFV3JEDI::setF90(const std::vector<std::string> vars) {
  size_t nv = vars.size();
  oops::Log::debug() << "setF90 " << nv << " vars = " << vars << std::endl;
  fvars_.resize(nv + 2);
  fvars_[0] = nv;
  for (size_t jj = 0; jj < nv; ++jj) {
     int ii = 0;
     if (vars[jj]=="u") ii = 1;
     if (vars[jj]=="v") ii = 2;
     if (vars[jj]=="pt") ii = 3;
     if (vars[jj]=="delp") ii = 4;
     if (vars[jj]=="q") ii = 5;
     ASSERT(ii > 0);
     fvars_[jj+1] = ii;
  }
  fvars_[nv+1] = 999;  // just for checking
  oops::Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

VariablesFV3JEDI::~VariablesFV3JEDI() {}

// -----------------------------------------------------------------------------

VariablesFV3JEDI::VariablesFV3JEDI(const VariablesFV3JEDI & other): fvars_(other.fvars_) {}

// -----------------------------------------------------------------------------

void VariablesFV3JEDI::print(std::ostream & os) const {
  os << "fv3jedi::VariablesFV3JEDI: vars = " << fvars_;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
