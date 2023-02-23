! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_time_invariant_helpers_mod

use atlas_module, only: atlas_field, atlas_fieldset

use fv3jedi_constants_mod, only: ps, lapse_rate, lapse_exponent
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_geom_mod

implicit none

private
public :: calculate_nominal_surface_pressure

contains

! --------------------------------------------------------------------------------------------------

subroutine calculate_nominal_surface_pressure(afieldset)
  type(atlas_fieldset), intent(in) :: afieldset

  integer :: n, nmax
  type(atlas_field) :: orog_field
  type(atlas_field) :: nsp_field
  real(kind_real), pointer :: orog_ptr(:,:)
  real(kind_real), pointer :: nsp_ptr(:,:)

  real(kind=kind_real), parameter :: ps_at_msl = ps
  real(kind=kind_real) :: orog_factor

  orog_field = afieldset%field("filtered_orography")
  call orog_field%data(orog_ptr)
  if (afieldset%has_field("nominal_surface_pressure")) then
    nsp_field = afieldset%field("nominal_surface_pressure")
    call nsp_field%data(nsp_ptr)
  else
    call abor1_ftn('Need to include nominal_surface_pressure in yaml extra_fields')
  end if

  nmax = size(orog_ptr)

  n = 0
  do n=1, nmax
    if(orog_ptr(1,n) <= 0.0) then
      nsp_ptr(1,n) = ps_at_msl
    else
      orog_factor = 1.0 + lapse_rate/288.15_kind_real*orog_ptr(1,n)
      nsp_ptr(1,n) = ps_at_msl*orog_factor**lapse_exponent
    end if
  enddo

end subroutine calculate_nominal_surface_pressure

! --------------------------------------------------------------------------------------------------

end module fv3jedi_time_invariant_helpers_mod
