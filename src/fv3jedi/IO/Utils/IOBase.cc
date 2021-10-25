/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/IO/Utils/IOBase.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

IOFactory::IOFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in fv3jedi::IOFactory." << std::endl;
    ABORT("Element already registered in fv3jedi::IOFactory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

IOBase * IOFactory::create(const eckit::Configuration & conf, const Geometry & geom) {
  oops::Log::trace() << "fv3jedi::IOBase constructor starting" << std::endl;
  const std::string &id = conf.getString("filetype");
  typename std::map<std::string, IOFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in fv3jedi::IOFactory." << std::endl;
    ABORT("Element does not exist in fv3jedi::IOFactory.");
  }
  IOBase * ptr = jloc->second->make(conf, geom);
  oops::Log::trace() << "fv3jedi::IOBase constructor done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
