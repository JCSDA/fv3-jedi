/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_LOCALIZATIONMATRIXFV3JEDI_H_
#define FV3_JEDI_SRC_LOCALIZATIONMATRIXFV3JEDI_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "GeometryFV3JEDI.h"
#include "eckit/config/Configuration.h"
#include "oops/interface/LocalizationBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "Fortran.h"

// Forward declarations
namespace fv3jedi {
  class GeometryFV3JEDI;
  class IncrementFV3JEDI;

/// Localization matrix for FV3JEDI model.

// -----------------------------------------------------------------------------
class LocalizationMatrixFV3JEDI: public util::Printable,
                        private boost::noncopyable,
                        private util::ObjectCounter<LocalizationMatrixFV3JEDI> {
 public:
  static const std::string classname()
                             {return "fv3jedi::LocalizationMatrixFV3JEDI";}

  LocalizationMatrixFV3JEDI(const GeometryFV3JEDI &,
                            const eckit::Configuration &);
  ~LocalizationMatrixFV3JEDI();
  void multiply(IncrementFV3JEDI &) const;

 private:
  void print(std::ostream &) const;
  F90lclz keyFtnConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_LOCALIZATIONMATRIXFV3JEDI_H_
