/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_LOCALIZATION_LOCALIZATIONMATRIX_H_
#define FV3JEDI_LOCALIZATION_LOCALIZATIONMATRIX_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Localization/LocalizationMatrixFortran.h"
#include "oops/interface/LocalizationBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"


// Forward declarations
namespace fv3jedi {
  class Geometry;
  class Increment;

/// Localization matrix for FV3JEDI model.

// -----------------------------------------------------------------------------
class LocalizationMatrix: public util::Printable,
                        private boost::noncopyable,
                        private util::ObjectCounter<LocalizationMatrix> {
 public:
  static const std::string classname()
                             {return "fv3jedi::LocalizationMatrix";}

  LocalizationMatrix(const Geometry &,
                            const eckit::Configuration &);
  ~LocalizationMatrix();
  void multiply(Increment &) const;

 private:
  void print(std::ostream &) const;
  F90lclz keyFtnConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_LOCALIZATION_LOCALIZATIONMATRIX_H_
