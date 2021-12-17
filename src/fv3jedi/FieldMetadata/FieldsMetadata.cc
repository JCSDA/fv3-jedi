/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// -------------------------------------------------------------------------------------------------

/*!
  Configuration files should be formatted as e.g.:

  Fields:
    - FieldName: ud
      FieldIONames: [u, ud, U]
      InterpType: nearest
      Kind: double
      Levels: full
      LongName: u_component_of_native_D_grid_wind
      Space: vector
      StaggerLoc: northsouth
      Tracer: false
      Units: ms-1
      IOFile: core

    - FieldName: t
      FieldIONames: [t, T, air_temperature]
      Kind: double
      InterpType: nearest
      Levels: full
      LongName: air_temperature
      Space: magnitude
      StaggerLoc: center
      Tracer: false
      Units: K
      IOFile: core

  FieldName:    The name that the interface recognizes a field as. E.g. if a field needs to be
                accessed in a variable change this is the name that is used to do so.
  FieldIONames: Name used by the model and in the files. There can be multiple options per
                FieldName.
  InterpType:   Type of interpolation used.
  Kind:         What kind of data, 'double', 'integer' (default: double)
  Levels:       Can be 'full' (nz), 'half' (nz+1) or a number. (Default: full)
  LongName:     Long name written to output files (default: FieldName)
  Space:        Where the field represents a 'vector', a 'magnitude' or a 'direction'
                (default: magnitude)
  StaggerLoc:   Horizontal staggering location, options are 'center', 'northsouth', 'eastwest' or
                'corner' (default: center)
  Tracer:       Boolean of whether or nor the field is a tracer (default: false)
  Units:        Units of the field, can be 1 for dimensionless
  IOFile:       Optional; can be core, tracer, sfcd, sfcw; determines which restart file to
                read the field from. Default will be defined in IO/fv3jedi_io_<gfs|geos>_mod.f90

  If no default the metadata must be provided below.
*/

// -------------------------------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "fv3jedi/FieldMetadata/FieldsMetadata.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  FieldsMetadata::FieldsMetadata(const Parameters_ & params, int & nlev) {
    // Local versions
    std::string fieldIOName;
    std::string fieldName;
    std::string kind;
    std::string levels_str;
    std::string longName;
    std::string space;
    std::string staggerLoc;
    bool tracer;
    std::string units;
    std::string interpType;
    std::string io_file;

    //  List of field sets
    std::vector<FieldsetsParameters> paramsFieldSets = params.fieldSets.value();
    int paramsFieldSetsSize = paramsFieldSets.size();

    //  Loop over sets of fields, i.e. each yaml file
    for (int km = 0; km < paramsFieldSetsSize; ++km) {
      // Create new config for each field set
      eckit::PathName pathNameFieldSet(paramsFieldSets[km].fieldSet.value());
      const eckit::YAMLConfiguration paramsFieldSet(pathNameFieldSet);

      FieldsParameters fieldsparams;
      fieldsparams.validateAndDeserialize(paramsFieldSet);

      //  List of fields
      std::vector<FieldParameters> paramsFields = fieldsparams.fields.value();
      int paramsFieldsSize = paramsFields.size();

      // Loop over fields
      for (int jm = 0; jm < paramsFieldsSize; ++jm) {
        //  List of potential names, default is FieldName
        fieldName = paramsFields[jm].fieldName.value();
        std::vector<std::string> fieldIONames;

        // Always add the Field name to the map keys
        fieldIONames.push_back(fieldName);

        // Add any user provided IO names to the map keys
        for (auto &IOName : paramsFields[jm].fieldIONames.value()) {
            fieldIONames.push_back(IOName);
        }

        // Add any user provided long name to the map keys
        if (paramsFields[jm].longName.value() != boost::none) {
          fieldIONames.push_back(*paramsFields[jm].longName.value());
        }

        int confFieldIONamesSize = fieldIONames.size();

        // Loop over potential field names
        for (int im = 0; im < confFieldIONamesSize; ++im) {
          // Insert new field into the map
          fieldIOName = fieldIONames[im];

          // Only add if not already in the map
          if ( fields_.find(fieldIOName) == fields_.end() ) {
            fields_.insert(std::pair<std::string, FieldMetadata>(fieldIOName,
                                                                 FieldMetadata(fieldIOName)));

            // Pointer to current
            FieldMetadata& fieldmetadata = fields_.find(fieldIOName)->second;

            // Push other metadata
            // -------------------

            kind = paramsFields[jm].precisionKind.value();
            levels_str = paramsFields[jm].levels.value();
            if (paramsFields[jm].longName.value() == boost::none) {
              longName = fieldName;
            } else {
              longName = *paramsFields[jm].longName.value();
            }
            space = paramsFields[jm].space.value();
            staggerLoc = paramsFields[jm].staggerLoc.value();
            tracer = paramsFields[jm].tracer.value();
            units = paramsFields[jm].units.value();
            interpType = paramsFields[jm].interpType.value();
            io_file = paramsFields[jm].ioFile.value();

            // Check for valid choices
            fieldmetadata.checkKindValid(fieldIOName, kind);
            fieldmetadata.checkSpaceValid(fieldIOName, space);
            fieldmetadata.checkStaggerLocValid(fieldIOName, staggerLoc);
            fieldmetadata.checkLevelValid(fieldIOName, levels_str);
            fieldmetadata.checkInterpTypeValid(fieldIOName, interpType);

            // Convert string levels to integer levels
            int levels;

            if (levels_str == "full") {
              levels = nlev;
            } else if (levels_str == "half") {
              levels = nlev + 1;
            } else {
              try {
                levels = std::stoi(levels_str);
              } catch (std::invalid_argument& e) {
                ABORT("FieldMetadata::FieldsMetadata levels neither full, half or an integer");
              }
            }

            fieldmetadata.setFieldName(fieldName);
            fieldmetadata.setKind(kind);
            fieldmetadata.setLevels(levels);
            fieldmetadata.setLongName(longName);
            fieldmetadata.setSpace(space);
            fieldmetadata.setStaggerLoc(staggerLoc);
            fieldmetadata.setTracer(tracer);
            fieldmetadata.setUnits(units);
            fieldmetadata.setInterpType(interpType);
            fieldmetadata.setIOFile(io_file);
          }  // End if not already in the map
        }  // End loop over potential field names
      }  // End loop over fields
    }  // End loop over field sets
  }

  // -----------------------------------------------------------------------------------------------

  FieldMetadata FieldsMetadata::getField(const std::string & fieldIOName) const {
    if ( fields_.find(fieldIOName) == fields_.end() ) {
      ABORT("FieldMetadata::getField: Field "+fieldIOName+" not found in FieldsMetadata");
    }
    return fields_.find(fieldIOName)->second;
  }

  // -----------------------------------------------------------------------------------------------

  size_t FieldsMetadata::getLevelsFromLongName(const std::string & longName) const {
    // Iterate over map and search for longname. When found return levels
    size_t levels = 0;
    for ( auto it = fields_.begin(); it != fields_.end(); ++it ) {
       if (longName == it->second.getLongName()) {
         levels = it->second.getLevels();
       }
    }
    if (levels == 0) {
      oops::Log::error() << "getLevelsFromLongName did not find " << longName << std::endl;
      throw eckit::BadParameter("Long name not found");
    }
    return levels;
  }

  // -----------------------------------------------------------------------------------------------

  oops::Variables FieldsMetadata::LongNameFromIONameLongNameOrFieldName(
                                                               const oops::Variables & vars) const {
    // Method takes as input an oops Variables object containing either IO names, field names or
    // long names. It returns the corresponding long name.

    // Convert variables to vector of strings
    const std::vector<std::string>& vecStrInName = vars.variables();

    // Vector of long name strings
    std::vector<std::string> vecStrLongName;

    // Make sure a field is found
    bool found;

    // Iterate over vars and find equivalent long name
    for (auto &StrInName : vecStrInName) {
      // First check if the incoming name is an IO name:
      auto it = fields_.find(StrInName);

      // Make sure something is found
      found = false;

      if (it != fields_.end()) {
        // If IO name extract from element (it)
        vecStrLongName.push_back(it->second.getLongName());
        found = true;
      } else {
        // If not IO name iterate through elements and check long names or field names
        for ( auto ite = fields_.begin(); ite != fields_.end(); ++ite ) {
          if (StrInName == ite->second.getLongName() || StrInName == ite->second.getFieldName()) {
            vecStrLongName.push_back(ite->second.getLongName());
            found = true;
          }
        }
      }

      // Make sure element was found
      ASSERT(found);
    }

    // Return long name vars
    return oops::Variables(vecStrLongName);
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi

// -------------------------------------------------------------------------------------------------
