/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Printable.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOBase : public util::Printable, private boost::noncopyable {
 public:
  explicit IOBase(const eckit::Configuration &, const Geometry &) {}
  virtual ~IOBase() {}
  virtual void read(State &) const = 0;
  virtual void read(Increment &) const = 0;
  virtual void write(const State &) const = 0;
  virtual void write(const Increment &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -------------------------------------------------------------------------------------------------

class IOFactory;

// -------------------------------------------------------------------------------------------------

class IOFactory {
 public:
  static IOBase * create(const eckit::Configuration &, const Geometry &);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~IOFactory() = default;

 protected:
  explicit IOFactory(const std::string &name);

 private:
  virtual IOBase * make(const eckit::Configuration &, const Geometry &) = 0;

  static std::map < std::string, IOFactory * > & getMakers() {
    static std::map < std::string, IOFactory * > makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template<class T>
class IOMaker : public IOFactory {
  IOBase * make(const eckit::Configuration & conf, const Geometry & geom) override {
    return new T(conf, geom);
  }

 public:
  explicit IOMaker(const std::string & name) : IOFactory(name) {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
