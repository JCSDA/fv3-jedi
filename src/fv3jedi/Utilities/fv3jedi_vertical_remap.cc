/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <optional>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/ModelData/ModelData.h"
#include "fv3jedi/Utilities/Constants.h"
#include "fv3jedi/Utilities/fv3jedi_vertical_remap.h"
#include "fv3jedi/VariableChange/VaderCookbook.h"
#include "fv3jedi/VariableChange/VariableChange.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "vader/recipes/LnAirPressure.h"

namespace fv3jedi {

// 1D interpolation inline helper functions
inline double interpZ0(double & p, double & tv,
                       double & zd, double & pd,
                       double & gor) {
  return zd + tv/gor*log(pd/p);
}
inline double interpZ1(double & p, double & tv,
                       double & zd, double & pd,
                       double & gamma, double & gor) {
  return zd - tv/gamma*( pow(pd/p, -gamma/gor) - 1. );
}
inline double interpP0(const double & z, double & zu, \
                       double & pu, double & tvu,
                       double & gor) {
  return pu*exp( -gor/tvu*( z - zu ) );
}
inline double interpP1(const double & z, double & zu, \
                       double & pu, double & tvu,
                       double & gamma, double & gor) {
  return pu*pow( 1. + gamma/tvu*( z - zu ), -gor/gamma );
}
inline double computeGamma(double & pu, double & tvu, \
                           double & pd, double & tvd,
                           double & gor) {
  return -gor*log(tvd/tvu)/log(pd/pu);
}

// Adjusts surface pressure based on differences between interpolated and grid terrain
atlas::Field newPs(atlas::FieldSet & fset, atlas::Field & zsOrog,
                   const atlas::FieldSet & fsetThermo,
                   const double epsilon, const double beta) {
  // This method computes a new surface pressure given a new orography.
  // The new pressure is computed assuming a hydrostatic balance
  // and a constant temperature lapse rate. Below ground, the
  // lapse rate is assumed to be -6.5 k/km.

  bool found_surface;
  double pd;
  double pu;
  double tvd;
  double tvu;
  double gamma;
  double zu;

  auto viewZs = atlas::array::make_view<double, 2>(fset["geopotential_height_at_surface"]);
  auto viewPs = atlas::array::make_view<double, 2>(fset["air_pressure_at_surface"]);

  auto viewP = atlas::array::make_view<double, 2>(fsetThermo["air_pressure"]);
  auto viewTv = atlas::array::make_view<double, 2>(fsetThermo["virtual_temperature"]);

  auto viewZsNew = atlas::array::make_view<double, 2>(zsOrog);

  size_t nGrid = viewP.shape(0);
  int npz = viewP.shape(1);

  double gor = fv3jedi::getConstant("grav")/fv3jedi::getConstant("rdry");

  atlas::Field fieldPsNew = zsOrog.functionspace().createField<double>(
                              atlas::option::name("air_pressure_at_surface") |
                              atlas::option::levels(1));

  auto viewPsNew = atlas::array::make_view<double, 2>(fieldPsNew);

  for ( size_t iGrid = 0; iGrid < nGrid ; iGrid++ ) {
    found_surface = false;

    // Compute surface pressure below the original ground
    // --------------------------------------------------

    // Assume a specified lapse rate below the surface
    gamma = beta;

    // Compute values at lowest grid cell
    pu = viewP(iGrid, npz-1);
    tvu = viewTv(iGrid, npz-1);
    zu = interpZ1(pu, tvu, viewZs(iGrid, 0), viewPs(iGrid, 0), gamma, gor);

    if ( viewZsNew(iGrid, 0) <= zu ) {
      found_surface = true;
      if ( fabs(gamma) > epsilon ) {
        viewPsNew(iGrid, 0) = interpP1(viewZsNew(iGrid, 0), zu, pu, tvu, gamma, gor);
      } else {
        viewPsNew(iGrid, 0) = interpP0(viewZsNew(iGrid, 0), zu, pu, tvu, gor);
      }
    }

    // Compute surface pressure above the original ground
    // --------------------------------------------------

    if ( !found_surface ) {
      int iLevel = npz-2;
      do {
        // Replace lower grid cell values with previous upper grid cell values
        pd = pu;
        tvd = tvu;

        // Compute upper grid cell values
        pu = viewP(iGrid, iLevel);
        tvu = viewTv(iGrid, iLevel);

        // Compute lapse rate of virtual temperature
        gamma = computeGamma(pu, tvu, pd, tvd, gor);

        // Interpolate height at upper grid cell
        if ( fabs(gamma) > epsilon ) {
          zu = interpZ1(pu, tvu, zu, pd, gamma, gor);
        } else {
          zu = interpZ0(pu, tvu, zu, pd, gor);
        }

        // Interpolate surface pressure
        if ( viewZsNew(iGrid, 0) <= zu ) {
          found_surface = true;
          if ( fabs(gamma) > epsilon ) {
            viewPsNew(iGrid, 0) = interpP1(viewZsNew(iGrid, 0), zu, pu, tvu, gamma, gor);
          } else {
            viewPsNew(iGrid, 0) = interpP0(viewZsNew(iGrid, 0), zu, pu, tvu, gor);
          }
        }

        iLevel--;
      } while ( !found_surface && iLevel > 0 );
    }

    // Compute surface pressure above the top
    if ( !found_surface ) {
      // Interpolate surface pressure
      viewPsNew(iGrid, 0) = interpP0(viewZsNew(iGrid, 0), zu, pu, tvu, gor);
    }
  }

  // Return recalculated surface pressure
  return fieldPsNew;
}

// Performs vertical interpolation from one grid to another
atlas::FieldSet interpolate(atlas::FieldSet & fset,
                            const atlas::FieldSet & fsetThermo,
                            const atlas::FieldSet & fsetThermoAdjusted,
                            atlas::Field & fieldZsNew, const atlas::Field & fieldPsNew,
                            const double beta) {
  // This method vertically interpolates upper-air fields.
  // Wind, temperature, humidity and other tracers are interpolated.
  // The interpolation is cubic lagrangian in log pressure
  // with a monotonic constraint in the center of the domain.
  // In the outer intervals it is linear in log pressure.
  // outside the domain, fields are generally held constant,
  // except for temperature and humidity below the input domain,
  // where the temperature lapse rate is held fixed at -6.5 k/km and
  // the relative humidity is held constant. Pressure fields are
  // computed using a variable change.

  int k;
  double dz;
  double z2s, q2s;
  double z1a, z1b, z1c, z1d;
  double q1a, q1b, q1c, q1d;

  auto viewLogP = atlas::array::make_view<double, 2>(fsetThermo["ln_air_pressure"]);
  auto viewLogPRemap = atlas::array::make_view<double, 2>(fsetThermoAdjusted["ln_air_pressure"]);

  auto viewZsNew = atlas::array::make_view<double, 2>(fieldZsNew);
  auto viewPsNew = atlas::array::make_view<double, 2>(fieldPsNew);

  size_t nGrid = viewLogP.shape(0);
  int npz = viewLogP.shape(1);

  std::vector<int> k1s(nGrid*npz);
  std::vector<double> ffa(nGrid*npz), ffb(nGrid*npz), ffc(nGrid*npz), ffd(nGrid*npz);

  // Initialize output field set
  atlas::FieldSet fsetRemap = fset.clone();

  for ( size_t iGrid = 0; iGrid < nGrid ; iGrid++ ) {
    // Match old and remapped vertical coordinates
    for ( int iLevel = npz-1; iLevel >= 0; iLevel-- ) {
      k = npz;
      do {
        if ( -viewLogPRemap(iGrid, iLevel) < -viewLogP(iGrid, k-1) ) break;
        k--;
        if ( k == 0 ) break;
      } while ( true );
      k1s[iGrid*npz+iLevel] = k;
    }

    // Compute interpolation coefficients
    for ( int iLevel = npz-1; iLevel >= 0; iLevel-- ) {
      if ( k1s[iGrid*npz+iLevel] == 1 || k1s[iGrid*npz+iLevel] == npz-1 ) {
        z2s = -viewLogPRemap(iGrid, iLevel);

        z1a = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]);
        z1b = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]-1);

        ffa[iGrid*npz+iLevel] = (z2s - z1b)/(z1a - z1b);
        ffb[iGrid*npz+iLevel] = (z2s - z1a)/(z1b - z1a);
      } else if ( k1s[iGrid*npz+iLevel] > 1 && k1s[iGrid*npz+iLevel] < npz-1 )  {
        z2s = -viewLogPRemap(iGrid, iLevel);

        z1a = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]+1);
        z1b = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]);
        z1c = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]-1);
        z1d = -viewLogP(iGrid, k1s[iGrid*npz+iLevel]-2);

        ffa[iGrid*npz+iLevel] = (z2s - z1b)/(z1a - z1b)* \
                                (z2s - z1c)/(z1a - z1c)* \
                                (z2s - z1d)/(z1a - z1d);
        ffb[iGrid*npz+iLevel] = (z2s - z1a)/(z1b - z1a)* \
                                (z2s - z1c)/(z1b - z1c)* \
                                (z2s - z1d)/(z1b - z1d);
        ffc[iGrid*npz+iLevel] = (z2s - z1a)/(z1c - z1a)* \
                                (z2s - z1b)/(z1c - z1b)* \
                                (z2s - z1d)/(z1c - z1d);
        ffd[iGrid*npz+iLevel] = (z2s - z1a)/(z1d - z1a)* \
                                (z2s - z1b)/(z1d - z1b)* \
                                (z2s - z1c)/(z1d - z1c);
      }
    }
  }

  // Loop through fields and interpolate or substitute
  for ( auto & fieldRemap : fsetRemap ) {
    if ( fieldRemap.shape(1) == 1 ) {
      if (fieldRemap.name() == "geopotential_height_at_surface") {
        // Set surface height to new grid's orogography
        fieldRemap = fieldZsNew.clone();
      } else if (fieldRemap.name() == "air_pressure_at_surface") {
        // Set surface pressure to recalculated pressure
        fieldRemap = fieldPsNew.clone();
      } else {
        // Other surface fields not allowed
        oops::Log::error()
            << fieldRemap.name()
            << ": vertical remapping of surface fields other than pressure "
               "and geopotential height not defined";
      }

    } else if ( fieldRemap.name() == "air_pressure" ||
                fieldRemap.name() == "air_pressure_levels" ||
                fieldRemap.name() == "air_pressure_thickness" ||
                fieldRemap.name() == "ln_air_pressure" ||
                fieldRemap.name() == "ln_air_pressure_at_interface" ||
                fieldRemap.name() == "air_pressure_to_kappa" ) {
      // Substitute 3D pressure from thermodynamic fieldset into output fieldset
      fieldRemap = fsetThermoAdjusted[fieldRemap.name()].clone();

    } else {
      // Perform interpolation
      auto viewField = atlas::array::make_view<double, 2>(fset[fieldRemap.name()]);
      auto viewFieldRemap = atlas::array::make_view<double, 2>(fieldRemap);

      for ( size_t iGrid = 0; iGrid < nGrid ; iGrid++ ) {
        for ( int iLevel = npz-1; iLevel >= 0; iLevel-- ) {
          if ( k1s[iGrid*npz+iLevel] == npz ) {
            // Constant below domain
            viewFieldRemap(iGrid, iLevel) = viewField(iGrid, npz-1);
          } else if ( k1s[iGrid*npz+iLevel] == 0 ) {
            // Constant above domain
            viewFieldRemap(iGrid, iLevel) = viewField(iGrid, 0);
          } else if ( k1s[iGrid*npz+iLevel] == npz-1 || k1s[iGrid*npz+iLevel] == 1 ) {
            // Linear in log pressure in outside intervals of domain
            q1a = viewField(iGrid, k1s[iGrid*npz+iLevel]);
            q1b = viewField(iGrid, k1s[iGrid*npz+iLevel]-1);
            viewFieldRemap(iGrid, iLevel) = ffa[iGrid*npz+iLevel]*q1a + ffb[iGrid*npz+iLevel]*q1b;
          } else {
            // Cubic Lagrangian in log pressure with monotonoic constraint in center of domain
            q1a = viewField(iGrid, k1s[iGrid*npz+iLevel]+1);
            q1b = viewField(iGrid, k1s[iGrid*npz+iLevel]);
            q1c = viewField(iGrid, k1s[iGrid*npz+iLevel]-1);
            q1d = viewField(iGrid, k1s[iGrid*npz+iLevel]-2);
            q2s = ffa[iGrid*npz+iLevel]*q1a + ffb[iGrid*npz+iLevel]*q1b +
                  ffc[iGrid*npz+iLevel]*q1c + ffd[iGrid*npz+iLevel]*q1d;
            if ( q2s < std::min(q1b, q1c) ) {
              viewFieldRemap(iGrid, iLevel) = std::min(q1b, q1c);
            } else if ( q2s > std::max(q1b, q1c) ) {
              viewFieldRemap(iGrid, iLevel) = std::max(q1b, q1c);
            } else {
              viewFieldRemap(iGrid, iLevel) = q2s;
            }
          }
        }
      }
    }
  }

  // Compute temperature and humidity before the input domain

  auto viewTRemap = atlas::array::make_view<double, 2>(fsetRemap["air_temperature"]);
  auto viewQRemap = atlas::array::make_view<double, 2>(
                      fsetRemap["water_vapor_mixing_ratio_wrt_moist_air"]);

  //
  const double dltdz   = beta*fv3jedi::getConstant("rdry")/fv3jedi::getConstant("grav");
  const double dlpvdrt = -2.5e6/fv3jedi::getConstant("rvap");

  for ( size_t iGrid = 0; iGrid < nGrid ; iGrid++ ) {
    double tSurfRemap = viewTRemap(iGrid, npz-1);
    double qvSurfRemap = viewQRemap(iGrid, npz-1);
    for ( int iLevel = npz-1; iLevel >= 0; iLevel-- ) {
      dz = -viewLogPRemap(iGrid, iLevel) - (-viewLogP(iGrid, npz-1));

      if ( dz < 0. ) {
        viewTRemap(iGrid, iLevel) = tSurfRemap*exp(dltdz*dz);
        viewQRemap(iGrid, iLevel) = qvSurfRemap*exp(dlpvdrt*(1./viewTRemap(iGrid, iLevel) - \
                                                             1./tSurfRemap) - dz);
      }
    }
  }

  // Return interpolated fieldset
  return fsetRemap;
}

// Computes necessary fields for vertical remapping using Vader
atlas::FieldSet computeThermoFields(atlas::FieldSet & fset, const atlas::Field & ps,
                                    const ModelData & modelData, bool needNewPs) {
  // This method computes variables that are required by other methods using a
  // Vader variable change. If needNewPs is true, then log air pressure and virtual
  // temperature are computed. Otherwise, log air pressure and any other pressure or
  // hydrostatic variables in the original input fieldset are computed.

  // Initialize Vader cookbook config
  eckit::LocalConfiguration vaderCookbookConfig, vaderConfig;
  vader::VaderParameters vaderParams;

  // Configure Vader cookbook with FV3-JEDI custom cookbook
  std::map<std::string, std::vector<std::string>> fv3Cookbook = vaderFV3CustomCookbook();

  // Initialize output fieldset
  atlas::FieldSet fsetThermo;

  // In all cases we need ln_air_pressure
  fsetThermo.add(ps.clone());
  oops::Variables varsIn(std::vector<std::string>{"ln_air_pressure"});
  vaderCookbookConfig.set("air_pressure_levels", fv3Cookbook["air_pressure_levels"]);
  vaderCookbookConfig.set("air_pressure", fv3Cookbook["air_pressure"]);
  vaderCookbookConfig.set("ln_air_pressure", fv3Cookbook["ln_air_pressure"]);

  // Prepare output fieldset for variable change
  if ( needNewPs ) {
    // If trying to obtain new surface pressure, we need virtual temperature
    if ( fset.has("virtual_temperature") ) {
      fsetThermo.add(fset["virtual_temperature"].clone());
    } else {
      fsetThermo.add(fset["air_temperature"].clone());
      fsetThermo.add(fset["water_vapor_mixing_ratio_wrt_moist_air"].clone());
      varsIn += oops::Variables(std::vector<std::string>{"virtual_temperature"});
      vaderCookbookConfig.set("virtual_temperature", fv3Cookbook["virtual_temperature"]);
    }
  } else {
    // Now that we have new surface pressure, we need to get any additional pressure variables
    // since they can be obtain precisely and don't need to be interpolated
    for ( auto name : std::vector<std::string>{"air_pressure_thickness",
                                               "ln_air_pressure_at_interface",
                                               "air_pressure_to_kappa"} ) {
      if ( fset.has(name) ) {
        varsIn += oops::Variables(std::vector<std::string>{name});
        vaderCookbookConfig.set(name, fv3Cookbook[name]);
      }
    }
  }

  // Initialize Vader
  vaderConfig.set(vader::configCookbookKey, vaderCookbookConfig);
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());
  vader::Vader vader(vaderParams, vaderConfig);

  // Performed variable change
  oops::Variables varsOut = vader.changeVar(fsetThermo, varsIn);

  // Return 3D thermodyamic state fields
  return fsetThermo;
}

// VertRemap constructor
VertRemap::VertRemap(const fv3jedi::Geometry & geom,
                     const atlas::FieldSet & fsetOrogNew)
  : modelData_(geom) {
  oops::Log::trace() << "Entering VertRemap constructor" << std::endl;

  // Save new surface heights
  fieldZsNew_ = fsetOrogNew["geopotential_height_at_surface"].clone();

  oops::Log::trace() << "Leaving VertRemap constructor" << std::endl;
}


// Print member variables;
void VertRemap::print(std::ostream & os) const {
  os << "VertRemap: modelData_=" << modelData_ << ", fieldZsNew_=" << fieldZsNew_;
}

// Remapping method
atlas::FieldSet VertRemap::remap(atlas::FieldSet & fset) {
  oops::Log::trace() << "Entering VertRemap::remap method" << std::endl;

  // A minimum set of fields is required for the remapping
  ASSERT(fset.has("air_pressure_at_surface"));
  ASSERT(fset.has("geopotential_height_at_surface"));
  ASSERT(fset.has("air_temperature"));
  ASSERT(fset.has("water_vapor_mixing_ratio_wrt_moist_air"));

  // Compute thermodynamic profiles using the interpolated surface pressures
  const atlas::FieldSet fsetThermo = computeThermoFields(fset, fset["air_pressure_at_surface"],
                                                         modelData_, true);

  // Adjust surface pressure based on differences between interpolated and grid terrain
  const atlas::Field fieldPsNew = newPs(fset, fieldZsNew_, fsetThermo, epsilon_, beta_);

  // Compute thermdynamic profiles using the adjusted surface pressure
  const atlas::FieldSet fsetThermoAdjusted = computeThermoFields(fset, fieldPsNew,
                                                                 modelData_, false);

  // Perform interpolation of input fieldset to remapped vertical coordinates
  atlas::FieldSet fsetRemap = interpolate(fset, fsetThermo, fsetThermoAdjusted,
                                          fieldZsNew_, fieldPsNew, beta_);

  oops::Log::trace() << "Leaving VertRemap::remap method" << std::endl;

  // Return remapped
  return fsetRemap;
}

}  // namespace fv3jedi
