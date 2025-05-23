/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"

namespace fv3jedi {

extern "C" {
  void get_field_metadata_f(const FieldsMetadata* fieldsMetadata,
                            const char longNameC[], char varUnitsC[], char dataKindC[],
                            bool& tracer, int & levels, char mathSpacC[]);
}

}  // namespace fv3jedi
