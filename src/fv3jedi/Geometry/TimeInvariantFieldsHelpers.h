/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace atlas {
class FieldSet;
}  // namespace atlas

namespace oops {
class Variables;
}  // namespace oops

namespace fv3jedi {
/// Computes new time-invariant fields, and inserts them into the FieldSet
///
/// varsToAdd lists the new time-invariant fields to be computed and inserted
/// fset on input should contain all the inputs needed to compute varsToAdd; on output also
/// contains the varsToAdd.
void insertDerivedTimeInvariantFields(atlas::FieldSet & fset,
                                      const oops::Variables & varsToAdd);
}  // namespace fv3jedi

