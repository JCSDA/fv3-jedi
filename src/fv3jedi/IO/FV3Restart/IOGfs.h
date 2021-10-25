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

#include "eckit/config/Configuration.h"
#include "fv3jedi/IO/Utils/IOBase.h"
#include "IOGfs.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOGfs: public IOBase {
 public:
  static const std::string classname() {return "fv3jedi::IOGfs";}
  IOGfs(const eckit::Configuration &, const Geometry &);
  ~IOGfs();
  void read(State &) const override;
  void read(Increment &) const override;
  void write(const State &) const override;
  void write(const Increment &) const override;

 private:
  F90iogfs objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
