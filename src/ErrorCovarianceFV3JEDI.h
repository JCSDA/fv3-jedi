/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_ERRORCOVARIANCEFV3JEDI_H_
#define FV3_JEDI_SRC_ERRORCOVARIANCEFV3JEDI_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "Fortran.h"
#include "GeometryFV3JEDI.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace fv3jedi {
  class IncrementFV3JEDI;
  class StateFV3JEDI;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for FV3JEDI

class ErrorCovarianceFV3JEDI : public util::Printable,
                           private boost::noncopyable,
                           private util::ObjectCounter<ErrorCovarianceFV3JEDI> {
 public:
  static const std::string classname()
                                  {return "fv3jedi::ErrorCovarianceFV3JEDI";}

  ErrorCovarianceFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &,
                       const eckit::Configuration &, const StateFV3JEDI &);
  ~ErrorCovarianceFV3JEDI();

  void linearize(const StateFV3JEDI &, const GeometryFV3JEDI &);
  void multiply(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void inverseMultiply(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void randomize(IncrementFV3JEDI &) const;

 private:
  void print(std::ostream &) const;
  F90bmat keyFtnConfig_;
  boost::scoped_ptr<const GeometryFV3JEDI> geom_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3_JEDI_SRC_ERRORCOVARIANCEFV3JEDI_H_
