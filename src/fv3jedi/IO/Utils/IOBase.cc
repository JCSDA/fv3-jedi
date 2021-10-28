/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/IO/Utils/IOBase.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

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

IOBase * IOFactory::create(const Geometry & geom, const IOParametersBase & params) {
  oops::Log::trace() << "IOBase::create starting" << std::endl;
  const std::string &id = params.filetype.value().value();
  typename std::map<std::string, IOFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in fv3jedi::IOFactory." << std::endl;
    ABORT("Element does not exist in fv3jedi::IOFactory.");
  }
  IOBase * ptr = jloc->second->make(geom, params);
  oops::Log::trace() << "IOBase::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

std::unique_ptr<IOParametersBase>
IOFactory::createParameters(const std::string &name) {
  typename std::map<std::string, IOFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in fv3jedi::IOFactory");
  }
  return it->second->makeParameters();
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
