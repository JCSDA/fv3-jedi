/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_UTILITIESFV3JEDI_H_
#define FV3_JEDI_SRC_UTILITIESFV3JEDI_H_

#include "eckit/config/Configuration.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

  void stageFv3Files(const eckit::Configuration &conf);
  void removeFv3Files();

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_UTILITIESFV3JEDI_H_

