/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_UTILITIES_UTILITIESFV3JEDI_H_
#define SRC_UTILITIES_UTILITIESFV3JEDI_H_

#include "eckit/config/Configuration.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

  void stageFv3Files(const eckit::Configuration &);
  void removeFv3Files();

  void delete_file(const char *);

}  // namespace fv3jedi

#endif  // SRC_UTILITIES_UTILITIESFV3JEDI_H_

