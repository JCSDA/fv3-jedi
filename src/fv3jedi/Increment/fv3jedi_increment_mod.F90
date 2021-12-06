! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_increment_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use oops_variables_mod, only: oops_variables

use random_mod
use fckit_mpi_module

use fields_metadata_mod, only: field_metadata

use fv3jedi_field_mod,           only: fv3jedi_field, checksame, hasfield, get_field
use fv3jedi_fields_mod,          only: fv3jedi_fields
use fv3jedi_constants_mod,       only: constoz, cp, alhl, rgas
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_iter_mod,       only: fv3jedi_geom_iter
use fv3jedi_kinds_mod,           only: kind_real

use wind_vt_mod, only: d_to_a

use mpp_domains_mod, only: mpp_global_sum, bitwise_efp_sum, center, east, north, center

implicit none
private
public :: fv3jedi_increment, fv3jedi_increment_registry

type, extends(fv3jedi_fields) :: fv3jedi_increment
contains
  procedure, public :: ones
  procedure, public :: random
  procedure, public :: self_add
  procedure, public :: self_schur
  procedure, public :: self_sub
  procedure, public :: self_mul
  procedure, public :: dot_prod
  procedure, public :: diff_incr
  procedure, public :: dirac
  procedure, public :: getpoint
  procedure, public :: setpoint
end type fv3jedi_increment

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_increment

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_increment_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine ones(self)

implicit none
class(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 1.0_kind_real
enddo

end subroutine ones

! --------------------------------------------------------------------------------------------------

subroutine random(self)

implicit none
class(fv3jedi_increment), intent(inout) :: self

integer :: var
integer, parameter :: rseed = 7

do var = 1,self%nf
  call normal_distribution(self%fields(var)%array, 0.0_kind_real, 1.0_kind_real, rseed)
enddo

end subroutine random

! --------------------------------------------------------------------------------------------------

subroutine self_add(self,rhs)

implicit none
class(fv3jedi_increment), intent(inout) :: self
class(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_add")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + rhs%fields(var)%array
enddo

end subroutine self_add

! --------------------------------------------------------------------------------------------------

subroutine self_schur(self,rhs)

implicit none
class(fv3jedi_increment), intent(inout) :: self
class(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_schur")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array * rhs%fields(var)%array
enddo

end subroutine self_schur

! --------------------------------------------------------------------------------------------------

subroutine self_sub(self,rhs)

implicit none
class(fv3jedi_increment), intent(inout) :: self
class(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_sub")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array - rhs%fields(var)%array
enddo

end subroutine self_sub

! --------------------------------------------------------------------------------------------------

subroutine self_mul(self,zz)

implicit none
class(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),     intent(in)    :: zz

integer :: var

do var = 1,self%nf
  self%fields(var)%array = zz * self%fields(var)%array
enddo

end subroutine self_mul

! --------------------------------------------------------------------------------------------------

subroutine dot_prod(self,other,zprod)

implicit none
class(fv3jedi_increment), intent(in)    :: self
class(fv3jedi_increment), intent(in)    :: other
real(kind=kind_real),     intent(inout) :: zprod

real(kind=kind_real) :: zp
integer :: i,j,k
integer :: var

call checksame(self%fields,other%fields,"fv3jedi_increment_mod.dot_prod")

zp=0.0_kind_real
do var = 1,self%nf
  do k = 1,self%fields(var)%npz
    do j = self%jsc,self%jec
      do i = self%isc,self%iec
        zp = zp + self%fields(var)%array(i,j,k) * other%fields(var)%array(i,j,k)
      enddo
    enddo
  enddo
enddo

!Get global dot product
call self%f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
!if (self%f_comm%rank() == 0) print*, "Dot product test result: ", zprod

end subroutine dot_prod

! --------------------------------------------------------------------------------------------------

subroutine diff_incr(self, x1_fields, x2_fields, geom)

implicit none
class(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_field),      intent(in)    :: x1_fields(:)
type(fv3jedi_field),      intent(in)    :: x2_fields(:)
type(fv3jedi_geom),       intent(inout) :: geom

integer :: var, isc, iec, jsc, jec, npz
type(fv3jedi_field), pointer :: x1p, x2p
real(kind=kind_real), allocatable :: x1_ua(:,:,:), x1_va(:,:,:)
real(kind=kind_real), allocatable :: x2_ua(:,:,:), x2_va(:,:,:)
real(kind=kind_real), pointer :: x1_ud(:,:,:), x1_vd(:,:,:)
real(kind=kind_real), pointer :: x2_ud(:,:,:), x2_vd(:,:,:)
real(kind=kind_real), pointer :: x1_delp(:,:,:)
real(kind=kind_real), pointer :: x2_delp(:,:,:)
real(kind=kind_real), allocatable :: x1_ps(:,:,:), x2_ps(:,:,:)

! Make sure two states have same fields and in same position
call checksame(x1_fields, x2_fields, "fv3jedi_increment_mod.diff_incr x1 vs x2")

! Convenience
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


!D-Grid to A-Grid (if needed)
if (self%has_field('ua')) then
  allocate(x1_ua(isc:iec,jsc:jec,1:npz))
  allocate(x1_va(isc:iec,jsc:jec,1:npz))
  allocate(x2_ua(isc:iec,jsc:jec,1:npz))
  allocate(x2_va(isc:iec,jsc:jec,1:npz))
  if (hasfield(x1_fields, 'ua')) then
    call get_field(x1_fields, 'ua', x1_ua)
    call get_field(x1_fields, 'va', x1_va)
    call get_field(x2_fields, 'ua', x2_ua)
    call get_field(x2_fields, 'va', x2_va)
  elseif (hasfield(x1_fields, 'ud')) then
    call get_field(x1_fields, 'ud', x1_ud)
    call get_field(x1_fields, 'vd', x1_vd)
    call get_field(x2_fields, 'ud', x2_ud)
    call get_field(x2_fields, 'vd', x2_vd)
    call d_to_a(geom, x1_ud, x1_vd, x1_ua, x1_va)
    call d_to_a(geom, x2_ud, x2_vd, x2_ua, x2_va)
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: no way to determine A grid winds")
  endif
endif

!delp to ps
if (self%has_field('ps')) then

  allocate(x1_ps(isc:iec,jsc:jec,1))
  allocate(x2_ps(isc:iec,jsc:jec,1))

  if (hasfield(x1_fields, 'delp')) then
    call get_field(x1_fields, 'delp', x1_delp)
    x1_ps(:,:,1) = sum(x1_delp,3)
  elseif (hasfield(x1_fields, 'ps')) then
    call get_field(x1_fields, 'ps', x1_ps)
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x1")
  endif

  if (hasfield(x2_fields, 'delp')) then
    call get_field(x2_fields, 'delp', x2_delp)
    x2_ps(:,:,1) = sum(x2_delp,3)
  elseif (hasfield(x2_fields, 'ps')) then
    call get_field(x2_fields, 'ps', x2_ps)
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x2")
  endif

endif

do var = 1,self%nf

  !A-Grid winds can be a special case
  if (self%fields(var)%fv3jedi_name == 'ua') then

    self%fields(var)%array = x1_ua - x2_ua

  elseif (self%fields(var)%fv3jedi_name == 'va') then

    self%fields(var)%array = x1_va - x2_va

  !Ps can be a special case
  elseif (self%fields(var)%fv3jedi_name == 'ps') then

    self%fields(var)%array = x1_ps - x2_ps

  else

    !Get pointer to states
    call get_field(x1_fields, self%fields(var)%fv3jedi_name, x1p)
    call get_field(x2_fields, self%fields(var)%fv3jedi_name, x2p)

    !inc = state - state
    self%fields(var)%array = x1p%array - x2p%array

    !Nullify pointers
    nullify(x1p,x2p)

  endif

enddo

if (allocated(x1_ua)) deallocate(x1_ua)
if (allocated(x1_va)) deallocate(x1_va)
if (allocated(x2_ua)) deallocate(x2_ua)
if (allocated(x2_va)) deallocate(x2_va)
if (allocated(x1_ps)) deallocate(x1_ps)
if (allocated(x2_ps)) deallocate(x2_ps)

end subroutine diff_incr

! --------------------------------------------------------------------------------------------------

subroutine dirac(self, conf, geom)

implicit none
class(fv3jedi_increment),  intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fv3jedi_geom),        intent(in)    :: geom

integer :: ndir,idir,var

integer, allocatable :: ixdir(:),iydir(:),ildir(:),itdir(:)
character(len=32), allocatable :: ifdir(:)

logical :: found

character(len=:), allocatable :: str
character(len=:), allocatable :: str_array(:)

! Get Diracs positions
call conf%get_or_die("ndir",ndir)

allocate(ixdir(ndir))
allocate(iydir(ndir))
allocate(ildir(ndir))
allocate(itdir(ndir))

if ((conf%get_size("ixdir")/=ndir) .or. &
    (conf%get_size("iydir")/=ndir) .or. &
    (conf%get_size("ildir")/=ndir) .or. &
    (conf%get_size("itdir")/=ndir) .or. &
    (conf%get_size("ifdir")/=ndir)) &
  call abor1_ftn("fv3jedi_increment_mod.diracL=: dimension inconsistency")

call conf%get_or_die("ixdir",ixdir)
call conf%get_or_die("iydir",iydir)
call conf%get_or_die("ildir",ildir)
call conf%get_or_die("itdir",itdir)

call conf%get_or_die("ifdir",str_array)
ifdir = str_array
deallocate(str_array)

print*, ndir, ixdir, iydir, ildir, itdir, ifdir

! Setup Diracs
call self%zero()

! only u, v, T, ps and tracers allowed
do idir=1,ndir

  found = .false.

  ! is specified grid point, tile number on this processor
  if (geom%ntile == itdir(idir) .and. &
    ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
    iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then

    ! If so, perturb desired increment and level
    do var = 1,self%nf
      if (trim(self%fields(var)%fv3jedi_name) == trim(ifdir(idir))) then
        found = .true.
        self%fields(var)%array(ixdir(idir),iydir(idir),ildir(idir)) = 1.0
      endif
    enddo

    if (.not.found) call abor1_ftn("fv3jedi_increment_mod.dirac: dirac not found")

  endif

enddo

end subroutine dirac

! --------------------------------------------------------------------------------------------------

subroutine getpoint(self, geoiter, values)

implicit none

class(fv3jedi_increment), intent(in)    :: self
type(fv3jedi_geom_iter),  intent(in)    :: geoiter
real(kind=kind_real),     intent(inout) :: values(:)

integer :: var, nz, ii

ii = 0
do var = 1,self%nf
  nz = self%fields(var)%npz
  values(ii+1:ii+nz) = self%fields(var)%array(geoiter%iindex, geoiter%jindex,:)
  ii = ii + nz
enddo

end subroutine getpoint

! --------------------------------------------------------------------------------------------------

subroutine setpoint(self, geoiter, values)

implicit none

! Passed variables
class(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom_iter),  intent(in)    :: geoiter
real(kind=kind_real),     intent(in)    :: values(:)

integer :: var, nz, ii

ii = 0
do var = 1,self%nf
  nz = self%fields(var)%npz
  self%fields(var)%array(geoiter%iindex, geoiter%jindex,:) = values(ii+1:ii+nz)
  ii = ii + nz
enddo

end subroutine setpoint

! --------------------------------------------------------------------------------------------------

end module fv3jedi_increment_mod
