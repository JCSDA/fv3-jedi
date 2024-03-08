/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Geometry/TimeInvariantFieldsHelpers.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

namespace fv3jedi {

extern "C" {
  void fv3jedi_nominal_surface_pressure_f90(atlas::field::FieldSetImpl *);
}

void insertDerivedTimeInvariantFields(atlas::FieldSet & fset,
                                      const oops::Variables & varsToAdd) {
  for (const auto & var : varsToAdd.variables()) {
    if (var == "nominal_surface_pressure") {
      // Check requested field doesn't already exist
      ASSERT(!fset.has(var));
      // Check necessary inputs are present
      ASSERT(fset.has("filtered_orography"));

      // Add new atlas::Field for nominal surface pressure
      const auto & orog = fset.field("filtered_orography");
      atlas::Field nsp =
          orog.functionspace().template createField<double>(atlas::option::name(var)
                                                            | atlas::option::levels(1));
      fset.add(nsp);

      // Call fortran to set values in the new NSP field
      fv3jedi_nominal_surface_pressure_f90(fset.get());

    } else if (var == "water_mask"
               || var == "land_mask"
               || var == "land_uncovered_by_snow_mask"
               || var == "ice_uncovered_by_snow_mask"
               || var == "snow_mask") {
      // Set up several surface-type based masking fields:
      // water_mask
      // - where 1, marks ocean and lakes; this is where slmsk == 0
      // - used for interpolating CRTM surface fields
      //
      // land_mask
      // - where 1, marks land; this is where slmsk == 1 AND land type is not "glacial land ice"
      //   excludes land marked as glacier, but includes land covered by snow
      // - used for land DA when snow covered/uncovered is to be treated the same
      // land_uncovered_by_snow_mask
      // - as above, but excludes land with sufficient snow
      // - used for interpolating CRTM surface fields
      //
      // ice_uncovered_by_snow_mask
      // - where 1, marks ice; this is where slmsk == 2 OR (slmsk == 1 and land type is glacial)
      //   includes sea ice, frozen lakes, and glaciers over land, but excludes ice with sufficient
      //   snow cover
      // - used for interpolating CRTM surface fields
      //
      // snow_mask
      // - where 1, marks snow; these are land/ice points with sufficient snow as determined by
      //   the snow water equivalent ("sheleg" in background files)
      // - used for interpolating CRTM surface fields
      ASSERT(!fset.has(var));
      ASSERT(fset.has("slmsk"));
      ASSERT(fset.has("sheleg"));
      ASSERT(fset.has("stype"));
      ASSERT(fset.has("vtype"));

      // Wraps std::round with a safety check we're not TOO far from an integer, which may be a bug
      const auto & nearestInt = [](const double x) -> int {
        const double r = std::round(x);
        // But for now we have to use a very loose tolerance, because fv3-jedi-data test data has
        // an interpolation bug such that the slmsk integer field is not an integer. See:
        // https://github.com/JCSDA-internal/fv3-jedi/issues/1144
        ASSERT(std::abs(x - r) < 0.499);
        return static_cast<int>(r);
      };

      const atlas::Field slmsk = fset.field("slmsk");
      const auto slmsk_view = atlas::array::make_view<double, 2>(slmsk);
      const auto swe_view = atlas::array::make_view<double, 2>(fset.field("sheleg"));
      const auto stype_view = atlas::array::make_view<double, 2>(fset.field("stype"));
      const auto vtype_view = atlas::array::make_view<double, 2>(fset.field("vtype"));

      atlas::Field mask = slmsk.functionspace().template createField<double>(
          atlas::option::name(var) | atlas::option::levels(1));
      auto mask_view = atlas::array::make_view<double, 2>(mask);

      // From the fv3-jedi CRTM variable change code
      constexpr double minswe = 0.1;
      constexpr int vtype_glacial = 15;
      constexpr int stype_glacial = 16;

      for (size_t i = 0; i < mask.shape(0); ++i) {
        // The GFS slmsk has values {0,1,2} denoting {sea,land,ice}.
        int surfType = nearestInt(slmsk_view(i, 0));
        ASSERT(surfType >= 0 && surfType <= 2);

        // If land with vtype or stype matching "glacial land ice", recategorize as ice
        if (surfType == 1 && (nearestInt(vtype_view(i, 0)) == vtype_glacial
                              || nearestInt(stype_view(i, 0)) == stype_glacial)) {
          surfType = 2;
        }

        // If land/ice are covered in snow (as determined from SWE > min), flag for snow
        const bool snowCover = (surfType > 0 && swe_view(i, 0) > minswe);

        mask_view(i, 0) = 0.0;
        if ((var == "water_mask" && surfType == 0)
            || (var == "land_mask" && surfType == 1)
            || (var == "land_uncovered_by_snow_mask" && surfType == 1 && !snowCover)
            || (var == "ice_uncovered_by_snow_mask" && surfType == 2 && !snowCover)
            || (var == "snow_mask" && snowCover)) {
          mask_view(i, 0) = 1.0;
        }
      }
      fset.add(mask);

    } else {
      throw eckit::UserError("Unsupported time-invariant derived field: " + var);
    }
  }
}

}  // namespace fv3jedi
