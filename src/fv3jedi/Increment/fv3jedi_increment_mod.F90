! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_increment_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module,  only: fckit_configuration
use fckit_mpi_module,            only: fckit_mpi_sum

! oops
use random_mod,                  only: normal_distribution

! fv3jedi
use fv3jedi_field_mod,           only: fv3jedi_field, checksame, get_field
use fv3jedi_fields_mod,          only: fv3jedi_fields
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_iter_mod,       only: fv3jedi_geom_iter
use fv3jedi_kinds_mod,           only: kind_real

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
  procedure, public :: diff_states
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

class(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 1.0_kind_real
enddo

end subroutine ones

! --------------------------------------------------------------------------------------------------

subroutine random(self)

class(fv3jedi_increment), intent(inout) :: self

integer :: var
integer, parameter :: rseed = 7

do var = 1,self%nf
  call normal_distribution(self%fields(var)%array, 0.0_kind_real, 1.0_kind_real, rseed)
enddo

end subroutine random

! --------------------------------------------------------------------------------------------------

subroutine self_add(self,rhs)

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

class(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),     intent(in)    :: zz

integer :: var

do var = 1,self%nf
  self%fields(var)%array = zz * self%fields(var)%array
enddo

end subroutine self_mul

! --------------------------------------------------------------------------------------------------

subroutine dot_prod(self,other,zprod)

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

subroutine diff_states(self, state1_fields, state2_fields, geom)

! Arguments
class(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_field),      intent(in)    :: state1_fields(:)
type(fv3jedi_field),      intent(in)    :: state2_fields(:)
type(fv3jedi_geom),       intent(inout) :: geom

! Locals
integer :: f
type(fv3jedi_field), pointer :: state1p, state2p

! Loop over increment vars and set data to diff of two states
do f = 1, self%nf

  !Get pointers to states
  call get_field(state1_fields, self%fields(f)%short_name, state1p)
  call get_field(state2_fields, self%fields(f)%short_name, state2p)

  !inc = state - state
  self%fields(f)%array = state1p%array - state2p%array

  !Nullify pointers
  nullify(state1p, state2p)

enddo

end subroutine diff_states

! --------------------------------------------------------------------------------------------------

subroutine dirac(self, conf, geom)

! Arguments
class(fv3jedi_increment),  intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fv3jedi_geom),        intent(in)    :: geom

! Locals
integer :: ndir,idir
integer, allocatable :: ixdir(:),iydir(:),ildir(:),itdir(:)
character(len=32), allocatable :: ifdir(:)
character(len=:), allocatable :: str_array(:)
type(fv3jedi_field), pointer :: dirac_field

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

! Setup Diracs
call self%zero()

! only u, v, T, ps and tracers allowed
do idir=1,ndir

  ! Get the field
  call get_field(self%fields, trim(ifdir(idir)), dirac_field)

  ! is specified grid point, tile number on this processor
  if (geom%ntile == itdir(idir) .and. &
    ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
    iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then

    ! Make perturbation
    dirac_field%array(ixdir(idir),iydir(idir),ildir(idir)) = 1.0_kind_real

  endif

enddo

end subroutine dirac

! --------------------------------------------------------------------------------------------------

subroutine getpoint(self, geoiter, values)

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
