! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_time_invariant_helpers_interface_mod

use atlas_module, only: atlas_fieldset
use iso_c_binding

use fv3jedi_time_invariant_helpers_mod
use fv3jedi_geom_interface_mod
use fv3jedi_geom_mod

implicit none

private

contains

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_nominal_surface_pressure(c_afieldset) &
                                              bind(c, name='fv3jedi_nominal_surface_pressure_f90')

!Arguments
type(c_ptr), intent(in), value :: c_afieldset

type(atlas_fieldset) :: afieldset

! Fortran APIs
! ------------
afieldset = atlas_fieldset(c_afieldset)

! Call implementation
! -------------------
call calculate_nominal_surface_pressure(afieldset)

end subroutine c_fv3jedi_nominal_surface_pressure

! --------------------------------------------------------------------------------------------------

end module fv3jedi_time_invariant_helpers_interface_mod
