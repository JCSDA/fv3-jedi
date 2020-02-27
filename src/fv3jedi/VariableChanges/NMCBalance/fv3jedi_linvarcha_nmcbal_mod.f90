! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_linvarcha_nmcbal_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fv3jedi_kinds_mod

use fv3jedi_field_mod, only: copy_subset, has_field, pointer_field_array

use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use wind_vt_mod

use type_bump, only: bump_type
!upe type_gauss, only: gauss_type

implicit none
private

public :: fv3jedi_linvarcha_nmcbal
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_linvarcha_nmcbal
 integer :: nx,ny,nz
 ! regression coeffs for GSI unbalanced variables
 real(kind=kind_real), allocatable :: agvi(:,:,:) ! coeffs for psi and t
 real(kind=kind_real), allocatable :: wgvi(:,:)   ! coeffs for psi and ln(ps)
 real(kind=kind_real), allocatable :: bvi(:,:)    ! coeffs for psi and chi
 type(bump_type) :: bump
! type(gauss_type) :: gauss
end type fv3jedi_linvarcha_nmcbal

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, conf)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
type(fv3jedi_geom), target,     intent(in)    :: geom
type(fv3jedi_state), target,    intent(in)    :: bg
type(fv3jedi_state), target,    intent(in)    :: fg
type(fckit_configuration),      intent(in)    :: conf

integer :: nx,ny,nz
character(len=256) :: path_to_gsi_coeffs
character(len=:),allocatable :: str

! 1. Get nx/ny/nz from config
nx = 192
if (conf%has("nx")) then
  call conf%get_or_die("nx",nx)
endif
!self%gauss%nx = nx

ny = 96
if (conf%has("ny")) then
  call conf%get_or_die("ny",ny)
endif
!self%gauss%nx = ny

nz = 64
if (conf%has("nz")) then
  call conf%get_or_die("nz",nz)
endif
!self%gauss%nx = nz

! 2. Create Gaussian grid in saber
!call self%gauss%create()

! 3. Create interplation weights (bump)
call setup_grid_interp(self)

! 4. Read balance coeffs from file
path_to_gsi_coeffs = 'Data/nmcbal/nmcbal_stats.dat'
if (conf%has("path_to_gsi_coeffs")) then
  call conf%get_or_die("path_to_gsi_coeffs",str)
  path_to_gsi_coeffs = str
  deallocate(str)
endif
allocate(self%agvi(ny+1,nz,nz))
allocate(self%wgvi(ny+1,nz))
allocate(self%bvi(ny+1,nz))
call read_balance_coeffs(self,path_to_gsi_coeffs)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

! Deallacote everthing
deallocate(self%agvi)
deallocate(self%wgvi)
deallocate(self%bvi)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xuba,xbal)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xuba
type(fv3jedi_increment),        intent(inout) :: xbal

integer :: n

real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: psi
real(kind=kind_real), pointer, dimension(:,:,:) :: chi
real(kind=kind_real), pointer, dimension(:,:,:) :: tv
real(kind=kind_real), pointer, dimension(:,:,:) :: rh

! Unbalanced to balanced

! 1. Intep xuba -> Gaussian (self%bump)
! 2. Balance operator (Gaussian)
! 3. Intep back to CS = xbal

! Copy variables that are the same in both increments
call copy_subset(xuba%fields, xbal%fields)

! Pointers to increments
call pointer_field_array(xbal%fields, 'ua' , ua)
call pointer_field_array(xbal%fields, 'va' , va)
call pointer_field_array(xbal%fields, 't'  , t )
call pointer_field_array(xbal%fields, 'q'  , q )
call pointer_field_array(xuba%fields, 'psi', psi)
call pointer_field_array(xuba%fields, 'chi', chi)
call pointer_field_array(xuba%fields, 'tv' , tv )
call pointer_field_array(xuba%fields, 'rh' , rh )

! Identity for now
psi = ua
chi = va
tv  = t
rh  = q

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xbal,xuba)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xbal
type(fv3jedi_increment),        intent(inout) :: xuba

real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: psi
real(kind=kind_real), pointer, dimension(:,:,:) :: chi
real(kind=kind_real), pointer, dimension(:,:,:) :: tv
real(kind=kind_real), pointer, dimension(:,:,:) :: rh

! Copy variables that are the same in both increments
call copy_subset(xbal%fields, xuba%fields)

! Pointers to increments
call pointer_field_array(xbal%fields, 'ua' , ua)
call pointer_field_array(xbal%fields, 'va' , va)
call pointer_field_array(xbal%fields, 't'  , t )
call pointer_field_array(xbal%fields, 'q'  , q )
call pointer_field_array(xuba%fields, 'psi', psi)
call pointer_field_array(xuba%fields, 'chi', chi)
call pointer_field_array(xuba%fields, 'tv' , tv )
call pointer_field_array(xuba%fields, 'rh' , rh )

! Identity for now
ua = psi
va = chi
t  = tv
q  = rh

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xbal,xuba)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xbal
type(fv3jedi_increment),        intent(inout) :: xuba

real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: psi
real(kind=kind_real), pointer, dimension(:,:,:) :: chi
real(kind=kind_real), pointer, dimension(:,:,:) :: tv
real(kind=kind_real), pointer, dimension(:,:,:) :: rh

! Copy variables that are the same in both increments
call copy_subset(xbal%fields, xuba%fields)

! Pointers to increments
call pointer_field_array(xbal%fields, 'ua' , ua)
call pointer_field_array(xbal%fields, 'va' , va)
call pointer_field_array(xbal%fields, 't'  , t )
call pointer_field_array(xbal%fields, 'q'  , q )
call pointer_field_array(xuba%fields, 'psi', psi)
call pointer_field_array(xuba%fields, 'chi', chi)
call pointer_field_array(xuba%fields, 'tv' , tv )
call pointer_field_array(xuba%fields, 'rh' , rh )

! Identity for now
ua = psi
va = chi
t  = tv
q  = rh

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xuba,xbal)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xuba
type(fv3jedi_increment),        intent(inout) :: xbal

real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: psi
real(kind=kind_real), pointer, dimension(:,:,:) :: chi
real(kind=kind_real), pointer, dimension(:,:,:) :: tv
real(kind=kind_real), pointer, dimension(:,:,:) :: rh

! Copy variables that are the same in both increments
call copy_subset(xuba%fields, xbal%fields)

! Pointers to increments
call pointer_field_array(xbal%fields, 'ua' , ua)
call pointer_field_array(xbal%fields, 'va' , va)
call pointer_field_array(xbal%fields, 't'  , t )
call pointer_field_array(xbal%fields, 'q'  , q )
call pointer_field_array(xuba%fields, 'psi', psi)
call pointer_field_array(xuba%fields, 'chi', chi)
call pointer_field_array(xuba%fields, 'tv' , tv )
call pointer_field_array(xuba%fields, 'rh' , rh )

! Identity for now
psi = ua
chi = va
tv  = t
rh  = q

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

subroutine setup_grid_interp(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout)    :: self


end subroutine setup_grid_interp

! ------------------------------------------------------------------------------

subroutine read_balance_coeffs(self,path_to_gsi_coeffs)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout)    :: self
character(len=256), intent(in) :: path_to_gsi_coeffs


end subroutine read_balance_coeffs

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_nmcbal_mod
