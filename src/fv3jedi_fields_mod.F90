! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Handle fields for the FV3JEDI model

module fv3jedi_fields_mod

use iso_c_binding
use config_mod
use datetime_mod
use fv3jedi_geom_mod
use ufo_vars_mod
use fv3jedi_kinds
use ioda_locs_mod
use ufo_geovals_mod

use mpp_domains_mod,   only: EAST, NORTH
use mpp_mod,           only: mpp_pe, mpp_npes, mpp_root_pe
use fms_io_mod,        only: restart_file_type, register_restart_field, &
                             free_restart_type, restore_state, save_restart

use fv3jedi_mod,       only: fv_atmos_type, allocate_fv_atmos_type
use fv3jedi_constants, only: deg2rad

implicit none
private

public :: fv3jedi_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, interp, interp_tl, interp_ad, &
        & convert_to_ug, convert_from_ug, dirac, &
        & analytic_IC
public :: fv3jedi_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold FV3JEDI fields
type :: fv3jedi_field
  type(fv_atmos_type) :: Atm
  type(fv3jedi_geom), pointer :: geom
  integer :: nf
  integer :: isc, iec, jsc, jec
  integer :: root_pe
end type fv3jedi_field

#define LISTED_TYPE fv3jedi_field

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_field_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(ufo_vars),  intent(in)    :: vars

! Allocate main fields
call allocate_fv_atmos_type(self%Atm, &
                            geom%bd%isd, geom%bd%ied, geom%bd%jsd, geom%bd%jed, &
                            geom%bd%isc, geom%bd%iec, geom%bd%jsc, geom%bd%jec, &
                            geom%npz, geom%ntracers, &
                            geom%hydrostatic, geom%wind_type)

! Pointer to geometry
self%geom => geom

! Initialize all domain arrays to zero
call zeros(self)
self%Atm%phis   = 0.0_kind_real

! Number of fields
self%nf = 5

! For convenience
self%isc = geom%bd%isc
self%iec = geom%bd%iec
self%jsc = geom%bd%jsc
self%jec = geom%bd%jec

self%root_pe = 0
if (mpp_pe() == mpp_root_pe()) self%root_pe = 1

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(fv3jedi_field), intent(inout) :: self

if (allocated(   self%Atm%u)) deallocate (    self%Atm%u )
if (allocated(   self%Atm%v)) deallocate (    self%Atm%v )
if (allocated(  self%Atm%pt)) deallocate (   self%Atm%pt )
if (allocated(self%Atm%delp)) deallocate ( self%Atm%delp )
if (allocated(   self%Atm%q)) deallocate (    self%Atm%q )
if (allocated(self%Atm%phis)) deallocate ( self%Atm%phis )
if (allocated(   self%Atm%w)) deallocate (    self%Atm%w )
if (allocated(self%Atm%delz)) deallocate ( self%Atm%delz )

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(fv3jedi_field), intent(inout) :: self

!Zero out the entire domain

self%Atm%u = 0.0_kind_real
self%Atm%v = 0.0_kind_real
self%Atm%pt = 0.0_kind_real
self%Atm%delp = 0.0_kind_real
self%Atm%q = 0.0_kind_real
if (.not. self%Atm%hydrostatic) then
   self%Atm%w = 0.0_kind_real
   self%Atm%delz = 0.0_kind_real
endif

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine ones(self)
implicit none
type(fv3jedi_field), intent(inout) :: self

call zeros(self)

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = 1.0_kind_real
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = 1.0_kind_real
endif

end subroutine ones

! ------------------------------------------------------------------------------

subroutine random(self)
use random_vectors_mod
implicit none
type(fv3jedi_field), intent(inout) :: self
integer :: nq

call random_vector(self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:))
call random_vector(self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:))
call random_vector(self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:))
call random_vector(self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:))
do nq = 1,self%geom%ntracers
   call random_vector(self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,nq))
enddo
if (.not. self%Atm%hydrostatic) then
   call random_vector(self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:))
   call random_vector(self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:))
endif

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

self%Atm%calendar_type = rhs%Atm%calendar_type
self%Atm%date = rhs%Atm%date
self%Atm%date_init = rhs%Atm%date_init

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif
self%Atm%phis(self%isc:self%iec,self%jsc:self%jec) = rhs%Atm%phis(self%isc:self%iec,self%jsc:self%jec)

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) + rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif

return
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) * rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) * rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif

return
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) - rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) - rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif

return
end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(fv3jedi_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = zz * self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = zz * self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif

return
end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_field), intent(in)    :: rhs

self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) + zz * rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
if (.not. self%Atm%hydrostatic) then
   self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) + zz * rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
endif

return
end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)

use mpi,             only: mpi_real8,mpi_comm_world,mpi_sum

implicit none
type(fv3jedi_field), intent(in) :: fld1, fld2
real(kind=kind_real), intent(inout) :: zprod
real(kind=kind_real) :: zp
integer :: i,j,k,l
integer :: ierr

zp=0.0_kind_real

!u
do i = fld1%geom%bd%isc,fld1%geom%bd%iec
   do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do k = 1,fld1%geom%npz
         zp = zp + fld1%Atm%u(i,j,k) * fld2%Atm%u(i,j,k)
      enddo
   enddo
enddo

!v
do i = fld1%geom%bd%isc,fld1%geom%bd%iec
   do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do k = 1,fld1%geom%npz
         zp = zp + fld1%Atm%v(i,j,k) * fld2%Atm%v(i,j,k)
      enddo
   enddo
enddo

!pt
do i = fld1%geom%bd%isc,fld1%geom%bd%iec
   do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do k = 1,fld1%geom%npz
         zp = zp + fld1%Atm%pt(i,j,k) * fld2%Atm%pt(i,j,k)
      enddo
   enddo
enddo

!delp
do i = fld1%geom%bd%isc,fld1%geom%bd%iec
   do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do k = 1,fld1%geom%npz
         zp = zp + fld1%Atm%delp(i,j,k) * fld2%Atm%delp(i,j,k)
      enddo
   enddo
enddo

!tracers
do i = fld1%geom%bd%isc,fld1%geom%bd%iec
   do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do k = 1,fld1%geom%npz
         do l = 1,fld1%geom%ntracers
            zp = zp + fld1%Atm%q(i,j,k,l) * fld2%Atm%q(i,j,k,l)
         enddo
      enddo
   enddo
enddo

!Get global dot product
call mpi_allreduce(zp,zprod,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

!For debugging print result:
if (fld1%root_pe == 1) print*, "Dot product test result: ", zprod

return
end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

if (self%geom%size_cubic_grid==rhs%geom%size_cubic_grid) then
   self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%u(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%u(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%v(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%v(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%pt(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%delp(self%isc:self%iec,self%jsc:self%jec,:)
   self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) = self%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:) + rhs%Atm%q(self%isc:self%iec,self%jsc:self%jec,:,:)
   if (.not. self%Atm%hydrostatic) then
      self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%w(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%w(self%isc:self%iec,self%jsc:self%jec,:)
      self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) = self%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:) + rhs%Atm%delz(self%isc:self%iec,self%jsc:self%jec,:)
   endif
else
   call abor1_ftn("fv3jedi fields:  add_incr not implemented for low res increment yet")
endif

return
end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(fv3jedi_field), intent(inout) :: lhs
type(fv3jedi_field), intent(in)    :: x1
type(fv3jedi_field), intent(in)    :: x2

call zeros(lhs)
if (x1%geom%size_cubic_grid==x2%geom%size_cubic_grid) then
   lhs%Atm%u(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%u(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%u(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
   lhs%Atm%v(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%v(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%v(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
   lhs%Atm%pt(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%pt(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%pt(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
   lhs%Atm%delp(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%delp(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%delp(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
   lhs%Atm%q(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:,:) = x1%Atm%q(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:,:) - x2%Atm%q(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:,:)
   if (.not. lhs%Atm%hydrostatic) then
      lhs%Atm%w(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%w(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%w(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
      lhs%Atm%delz(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) = x1%Atm%delz(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:) - x2%Atm%delz(lhs%isc:lhs%iec,lhs%jsc:lhs%jec,:)
   endif
else
   call abor1_ftn("fv3jedi fields:  diff_incr not implemented for low res increment yet")
endif

return
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: fld
type(fv3jedi_field), intent(in)    :: rhs

if (fld%geom%size_cubic_grid==rhs%geom%size_cubic_grid) then
   call copy(fld, rhs)
else
   call abor1_ftn("fv3jedi_fields: change_resol not implmeneted yet")
endif

return
end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)

  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log
  use variable_transforms, only: d2a2c_vect
  use mpp_domains_mod, only: mpp_update_domains, DGRID_NE

  use field_manager_mod,       only: MODEL_ATMOS
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, set_tracer_profile

  implicit none

  type(fv3jedi_field), intent(inout) :: fld  !< Fields
  type(c_ptr), intent(in)       :: c_conf   !< Configuration
  type(datetime), intent(inout) :: vdate    !< DateTime

  type(restart_file_type) :: Fv_restart
  type(restart_file_type) :: Tr_restart

  integer :: id_restart, iounit, io_status, layout(2)
  integer :: date_init(6), date(6), calendar_type
  integer(kind=c_int) :: idate, isecs

  character(len=255) :: datapath_in, datapath_ti
  character(len=255) :: filename_core
  character(len=255) :: filename_trcr
  character(len=255) :: filename_cplr

  character(len=20) :: sdate,validitydate
  character(len=1024)  :: buf

  integer :: k

  character(len=64):: tracer_name
  integer :: ntracers, ntprog, nt, ierr

  integer :: print_read_info = 0


  !Set filenames
  !--------------
  filename_core = 'fv_core.res.nc'
  filename_trcr = 'fv_tracer.res.nc'
  filename_cplr = 'coupler.res'

  if (config_element_exists(c_conf,"filename_core")) then
     filename_core = config_get_string(c_conf,len(filename_core),"filename_core")
  endif
  if (config_element_exists(c_conf,"filename_trcr")) then
     filename_trcr = config_get_string(c_conf,len(filename_trcr),"filename_trcr")
  endif
  if (config_element_exists(c_conf,"filename_cplr")) then
     filename_cplr = config_get_string(c_conf,len(filename_cplr),"filename_cplr")
  endif

  datapath_in = config_get_string(c_conf,len(datapath_in),"datapath_read")
  datapath_ti = config_get_string(c_conf,len(datapath_ti),"datapath_tile")

  ! Register the variables that should be read
  ! ------------------------------------------
  !Winds
  if (trim(fld%geom%wind_type) == 'D-grid') then
     id_restart = register_restart_field(Fv_restart, filename_core, 'u', fld%Atm%u, &
                  domain=fld%geom%domain, position=NORTH)
     id_restart = register_restart_field(Fv_restart, filename_core, 'v', fld%Atm%v, &
                  domain=fld%geom%domain, position=EAST)
  elseif (trim(fld%geom%wind_type) == 'A-grid') then
     id_restart = register_restart_field(Fv_restart, filename_core, 'ua', fld%Atm%u, &
                                         domain=fld%geom%domain)
     id_restart = register_restart_field(Fv_restart, filename_core, 'va', fld%Atm%v, &
                                         domain=fld%geom%domain)
  endif

  !phis
  id_restart = register_restart_field(Fv_restart, filename_core, 'phis', fld%Atm%phis, &
               domain=fld%geom%domain)

  !Temperature
  id_restart = register_restart_field(Fv_restart, filename_core, 'T', fld%Atm%pt, &
               domain=fld%geom%domain)

  !Pressure thickness
  id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', fld%Atm%delp, &
               domain=fld%geom%domain)

  !Nonhydrostatic variables
  if (.not. fld%Atm%hydrostatic) then
     id_restart =  register_restart_field(Fv_restart, filename_core, 'W', fld%Atm%w, &
                   domain=fld%geom%domain)
     id_restart =  register_restart_field(Fv_restart, filename_core, 'DZ', fld%Atm%delz, &
                   domain=fld%geom%domain)
  endif


  ! Read file and fill variables
  ! ----------------------------
  call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
  call free_restart_type(Fv_restart)
 

! !Register and read tracers
! !-------------------------
! call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)
!
! do nt = 1, ntprog
!    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
!    call set_tracer_profile (MODEL_ATMOS, nt, fld%Atm%q(fld%geom%bd%isc:fld%geom%bd%iec,fld%geom%bd%isc:fld%geom%bd%iec,:,nt) )
!    id_restart = register_restart_field(Tr_restart, filename_trcr, tracer_name, fld%Atm%q(:,:,:,nt), &
!                                        domain=fld%geom%domain)
! enddo
! call restore_state(Tr_restart, directory=trim(adjustl(datapath_ti)))
! call free_restart_type(Tr_restart)

  id_restart =  register_restart_field(Tr_restart, filename_trcr, 'sphum', fld%Atm%q(:,:,:,1), &
                                       domain=fld%geom%domain)

  call restore_state(Tr_restart, directory=trim(adjustl(datapath_ti)))
  call free_restart_type(Tr_restart)


  ! Get dates from file
  !--------------------

  ! read date from coupler.res text file.
  sdate = config_get_string(c_conf,len(sdate),"date")
  iounit = 101
  open(iounit, file=trim(adjustl(datapath_in))//filename_cplr, form='formatted')
  read(iounit, '(i6)')  calendar_type
  read(iounit, '(6i6)') date_init
  read(iounit, '(6i6)') date
  close(iounit)
  fld%Atm%date = date
  fld%Atm%date_init = date_init
  fld%Atm%calendar_type = calendar_type
  idate=date(1)*10000+date(2)*100+date(3)
  isecs=date(4)*3600+date(5)*60+date(6)

  if (fld%root_pe == 1 .and. print_read_info == 1 ) then
     print *,'read_file: integer time from coupler.res: ',date,idate,isecs
  endif

  call datetime_from_ifs(vdate, idate, isecs)
  call datetime_to_string(vdate, validitydate)

  if (fld%root_pe == 1 .and. print_read_info == 1 ) then
     print *,'read_file: validity date: ',trim(validitydate)
     print *,'read_file: expected validity date: ',trim(sdate)
  endif

  return

end subroutine read_file

! ------------------------------------------------------------------------------
!> Analytic Initialization for the FV3 Model
!!
!! \details **analytic_IC()** initializes the FV3JEDI Field and State objects using one of
!! several alternative idealized analytic models.  This is intended to facilitate testing by
!! eliminating the need to read in the initial state from a file and by providing exact expressions
!! to test interpolations.  This function is activated by setting the "analytic_init" field in the
!! "initial" or "StateFile" section of the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * invent-state: Backward compatibility with original analytic init option
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!!
!! \author M. Miesch (adapted from a pre-existing call to invent_state)
!! \date March 13, 2018: Created
!!
!! \warning This routine initializes the fv3jedi_field object.  However, since the fv_atmos_type
!! component of fv3jedi_field is a subset of the corresponding object in the fv3 model,
!! this initialization routine is not sufficient to comprehensively define the full fv3 state.
!! So, this intitialization can be used for interpolation and other tests within JEDI but it is
!! cannot currently be used to initiate a forecast with fv3gfs.
!!
!! \warning This routine does not initialize the fv3jedi_interp member of the fv3jedi_field object
!!
!! \warning Though an input state file is not required for these analytic initialization routines,
!! some grid information (in particular the hybrid vertical grid coefficients ak and bk)
!! is still read in from an input file when creating the geometry object that is a required
!! member of fv3jedi_field; see c_fv3jedi_geo_setup() in fv3jedi_geom_mod.F90.
!!
!! \warning It's unclear whether the pt member of the fv_atmos_type structure is potential temperature
!! or temperature.  This routine assumes the latter.  If this is not correct, then we will need to
!! implement a conversion
!!
subroutine analytic_IC(fld, geom, c_conf, vdate)

  use kinds
  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log
  use constants_mod, only: pi=>pi_8
!  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, test1_advection_hadley
  
  !FV3 Test Cases
  use fv_arrays_mod,  only: fv_atmos_type, deallocate_fv_atmos_type
  use test_cases_mod, only: init_case, test_case
  use fv_control_mod, only: fv_init, pelist_all

  implicit none

  type(fv3jedi_field), intent(inout)     :: fld !< Fields
  type(fv3jedi_geom), target, intent(in) :: geom    !< Geometry 
  type(c_ptr), intent(in)                :: c_conf   !< Configuration
  type(datetime), intent(inout)          :: vdate    !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: i,j,k
  real(kind=kind_real) :: deg_to_rad = pi/180.0_kind_real
  real(kind=kind_real) :: rlat, rlon, z
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  type(fv_atmos_type), allocatable :: FV_AtmIC(:)
  real(kind=kind_real)             :: DTdummy = 900.0
  logical, allocatable             :: grids_on_this_pe(:)
  integer                          :: p_split = 1

  ! Poitner to geometry component of field object
  fld%geom => geom

  If (config_element_exists(c_conf,"analytic_init")) Then
     IC = Trim(config_get_string(c_conf,len(IC),"analytic_init"))
  Else
     ! This default value is for backward compatibility
     IC = "invent-state"
  EndIf

  call log%warning("fv3jedi_fields:analytic_init: "//IC)
  sdate = config_get_string(c_conf,len(sdate),"date")
  WRITE(buf,*) 'validity date is: '//sdate
  call log%info(buf)
  call datetime_set(sdate, vdate)

  !===================================================================
  int_option: Select Case (IC)

     Case("invent-state")

        call invent_state(fld,c_conf)

     Case("fv3_init_case")

        !Initialize temporary FV_Atm fv3 construct
        call fv_init(FV_AtmIC, DTdummy, grids_on_this_pe, p_split)
        deallocate(pelist_all)

        !Test case to run, see fv3: /tools/test_cases.F90 for possibilities
        test_case = config_get_int(c_conf,"fv3_test_case")

        call init_case( FV_AtmIC(1)%u,FV_AtmIC(1)%v,FV_AtmIC(1)%w,FV_AtmIC(1)%pt,FV_AtmIC(1)%delp,FV_AtmIC(1)%q, &
                        FV_AtmIC(1)%phis, FV_AtmIC(1)%ps,FV_AtmIC(1)%pe, FV_AtmIC(1)%peln,FV_AtmIC(1)%pk,FV_AtmIC(1)%pkz, &
                        FV_AtmIC(1)%uc,FV_AtmIC(1)%vc, FV_AtmIC(1)%ua,FV_AtmIC(1)%va,        & 
                        FV_AtmIC(1)%ak, FV_AtmIC(1)%bk, FV_AtmIC(1)%gridstruct, FV_AtmIC(1)%flagstruct,&
                        FV_AtmIC(1)%npx, FV_AtmIC(1)%npy, FV_AtmIC(1)%npz, FV_AtmIC(1)%ng, &
                        FV_AtmIC(1)%flagstruct%ncnst, FV_AtmIC(1)%flagstruct%nwat,  &
                        FV_AtmIC(1)%flagstruct%ndims, FV_AtmIC(1)%flagstruct%ntiles, &
                        FV_AtmIC(1)%flagstruct%dry_mass, &
                        FV_AtmIC(1)%flagstruct%mountain,       &
                        FV_AtmIC(1)%flagstruct%moist_phys, FV_AtmIC(1)%flagstruct%hydrostatic, &
                        FV_AtmIC(1)%flagstruct%hybrid_z, FV_AtmIC(1)%delz, FV_AtmIC(1)%ze0, &
                        FV_AtmIC(1)%flagstruct%adiabatic, FV_AtmIC(1)%ks, FV_AtmIC(1)%neststruct%npx_global, &
                        FV_AtmIC(1)%ptop, FV_AtmIC(1)%domain, FV_AtmIC(1)%tile, FV_AtmIC(1)%bd )

        !Copy from temporary structure into fields
        fld%Atm%u = FV_AtmIC(1)%u
        fld%Atm%v = FV_AtmIC(1)%v
        fld%Atm%pt = FV_AtmIC(1)%pt
        fld%Atm%delp = FV_AtmIC(1)%delp
        fld%Atm%q = FV_AtmIC(1)%q
        fld%Atm%phis = FV_AtmIC(1)%phis
        fld%geom%ak = FV_AtmIC(1)%ak
        fld%geom%ak = FV_AtmIC(1)%ak
        fld%geom%ptop = FV_AtmIC(1)%ptop
        if (.not. fld%Atm%hydrostatic) then
           fld%Atm%w = FV_AtmIC(1)%w
           fld%Atm%delz = FV_AtmIC(1)%delz
        endif

        !Deallocate temporary FV_Atm fv3 structure
        call deallocate_fv_atmos_type(FV_AtmIC(1))
        deallocate(FV_AtmIC)
        deallocate(grids_on_this_pe)

     Case ("dcmip-test-1-1")

        do i = geom%bd%isc,geom%bd%iec
           do j = geom%bd%jsc,geom%bd%jec
              rlat = deg_to_rad*geom%grid_lat(i,j)
              rlon = deg_to_rad*geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              !Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
              !                                 phis0,ps,rho0,hum0,q1,q2,q3,q4)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 !Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                 !                                 phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 !fld%Atm%u(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 !fld%Atm%v(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
                 If (.not.fld%Atm%hydrostatic) fld%Atm%w(i,j,k) = w0
                 fld%Atm%pt(i,j,k) = t0
                 fld%Atm%delp(i,j,k) = pe2-pe1
                 If (geom%ntracers >= 1) fld%Atm%q(i,j,k,1) = hum0
                 If (geom%ntracers >= 2) fld%Atm%q(i,j,k,2) = q1
                 If (geom%ntracers >= 3) fld%Atm%q(i,j,k,3) = q2
                 If (geom%ntracers >= 4) fld%Atm%q(i,j,k,4) = q3
                 If (geom%ntracers >= 5) fld%Atm%q(i,j,k,5) = q4
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do i = geom%bd%isc,geom%bd%iec
           do j = geom%bd%jsc,geom%bd%jec
              rlat = deg_to_rad*geom%grid_lat(i,j)
              rlon = deg_to_rad*geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              !Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
              !                            t0,phis0,ps,rho0,hum0,q1)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 !Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                 !                            t0,phis0,ps,rho0,hum0,q1)

                 !fld%Atm%u(i,j,k) = u0 !ATTN comment above
                 !fld%Atm%v(i,j,k) = v0
                 If (.not.fld%Atm%hydrostatic) fld%Atm%w(i,j,k) = w0
                 fld%Atm%pt(i,j,k) = t0
                 fld%Atm%delp(i,j,k) = pe2-pe1
                 If (geom%ntracers >= 1) fld%Atm%q(i,j,k,1) = hum0
                 If (geom%ntracers >= 2) fld%Atm%q(i,j,k,2) = q1
                 
              enddo
           enddo
        enddo

        WRITE(*,*) "DCMIP TEST 1-2"

     Case Default

        call invent_state(fld,c_conf)

     End Select int_option
        
end subroutine analytic_IC
  
! ------------------------------------------------------------------------------
subroutine invent_state(flds,config)

use kinds

implicit none

type(fv3jedi_field), intent(inout) :: flds    !< Model fields
type(c_ptr), intent(in)           :: config  !< Configuration structure

integer :: i,j,k

!u
do i = flds%geom%bd%isc,flds%geom%bd%iec
   do j = flds%geom%bd%jsc,flds%geom%bd%jec
      do k = 1,flds%geom%npz
         flds%Atm%u(i,j,k) = cos(0.25*flds%geom%grid_lon(i,j)) + cos(0.25*flds%geom%grid_lat(i,j))
      enddo
   enddo
enddo

!v
do i = flds%geom%bd%isc,flds%geom%bd%iec
   do j = flds%geom%bd%jsc,flds%geom%bd%jec
      do k = 1,flds%geom%npz
         flds%Atm%v(i,j,k) = 1.0_kind_real
      enddo
   enddo
enddo

!pt
do i = flds%geom%bd%isc,flds%geom%bd%iec
   do j = flds%geom%bd%jsc,flds%geom%bd%jec
      do k = 1,flds%geom%npz
         flds%Atm%pt(i,j,k) = cos(0.25*flds%geom%grid_lon(i,j)) + cos(0.25*flds%geom%grid_lat(i,j))
      enddo
   enddo
enddo

!delp
do i = flds%geom%bd%isc,flds%geom%bd%iec
   do j = flds%geom%bd%jsc,flds%geom%bd%jec
      do k = 1,flds%geom%npz
         flds%Atm%delp(i,j,k) = k
      enddo
   enddo
enddo

!q
do i = flds%geom%bd%isc,flds%geom%bd%iec
   do j = flds%geom%bd%jsc,flds%geom%bd%jec
      do k = 1,flds%geom%npz
         flds%Atm%q(i,j,k,1) = 0.0
      enddo
   enddo
enddo

return
end subroutine invent_state

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)
  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log

  implicit none

  type(fv3jedi_field), intent(in) :: fld  !< Fields
  type(c_ptr), intent(in)       :: c_conf   !< Configuration
  type(datetime), intent(in) :: vdate    !< DateTime

  type(restart_file_type) :: Fv_restart
  type(restart_file_type) :: Tr_restart

  integer :: id_restart, iounit, io_status
  integer :: date_init(6), date(6), calendar_type, layout(2)
  integer(kind=c_int) :: idate, isecs

  character(len=64) :: filename_core
  character(len=64) :: filename_trcr
  character(len=64) :: filename_cplr
  character(len=64) :: datefile

  character(len=20) :: sdate,validitydate

  character(len=255) :: datapath_out


  ! Place to save restarts
  ! ----------------------
  datapath_out = "Data/"
  if (config_element_exists(c_conf,"datapath_write")) then
     datapath_out = config_get_string(c_conf,len(datapath_out),"datapath_write")
  endif


  ! Current date
  ! ------------
  call datetime_to_ifs(vdate, idate, isecs)
  date(1) = idate/10000
  date(2) = idate/100 - date(1)*100
  date(3) = idate - (date(1)*10000 + date(2)*100)
  date(4) = isecs/3600
  date(5) = (isecs - date(4)*3600)/60
  date(6) = isecs - (date(4)*3600 + date(5)*60)
 

  ! Naming convection for the file
  ! ------------------------------
  filename_core = 'fv_core.res.nc'
  if (config_element_exists(c_conf,"filename_core")) then
    filename_core = config_get_string(c_conf,len(filename_core),"filename_core")
  endif

  filename_trcr = 'fv_tracer.res.nc'
  if (config_element_exists(c_conf,"filename_trcr")) then
    filename_trcr = config_get_string(c_conf,len(filename_trcr),"filename_trcr")
  endif

  filename_cplr = 'coupler.res'
  if (config_element_exists(c_conf,"filename_cplr")) then
    filename_cplr = config_get_string(c_conf,len(filename_cplr),"filename_cplr")
  endif

  !Append with the date
  write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2,A1)') date(1),date(2),date(3),".",date(4),date(5),date(6),"."
  filename_core = trim(datefile)//trim(filename_core)
  filename_trcr = trim(datefile)//trim(filename_trcr)
  filename_cplr = trim(datefile)//trim(filename_cplr)


  ! Register the variables that should be written
  ! ---------------------------------------------
  !Winds
  if (trim(fld%geom%wind_type) == 'D-grid') then
     id_restart = register_restart_field( Fv_restart, filename_core, 'u', fld%Atm%u, &
                                          domain=fld%geom%domain,position=NORTH )
     id_restart = register_restart_field( Fv_restart, filename_core, 'v', fld%Atm%v, &
                                          domain=fld%geom%domain,position=EAST )
  elseif (trim(fld%geom%wind_type) == 'A-grid') then
      id_restart =  register_restart_field(Fv_restart, filename_core, 'ua', fld%Atm%u, &
                                            domain=fld%geom%domain )
      id_restart =  register_restart_field(Fv_restart, filename_core, 'va', fld%Atm%v, &
                                            domain=fld%geom%domain )
  endif

  !phis
  id_restart = register_restart_field( Fv_restart, filename_core, 'phis', fld%Atm%phis, &
                                       domain=fld%geom%domain )

  !Temperature
  id_restart = register_restart_field( Fv_restart, filename_core, 'T', fld%Atm%pt, &
                                       domain=fld%geom%domain )

  !Pressure thickness
  id_restart = register_restart_field( Fv_restart, filename_core, 'DELP', fld%Atm%delp, &
                                       domain=fld%geom%domain )

  !Nonhydrostatic fields
  if (.not. fld%Atm%hydrostatic) then
      id_restart =  register_restart_field( Fv_restart, filename_core, 'W', fld%Atm%w, &
                                            domain=fld%geom%domain )
      id_restart =  register_restart_field( Fv_restart, filename_core, 'DZ', fld%Atm%delz, &
                                            domain=fld%geom%domain )
  endif

  !Cell center lat/lon
  id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lat', fld%geom%grid_lat, &
                                       domain=fld%geom%domain )
  id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lon', fld%geom%grid_lon, &
                                       domain=fld%geom%domain )

  id_restart =  register_restart_field( Fv_restart, filename_core, 'ua', fld%Atm%ua, &
                                        domain=fld%geom%domain )
  id_restart =  register_restart_field( Fv_restart, filename_core, 'va', fld%Atm%va, &
                                        domain=fld%geom%domain )

  ! Write variables to file
  ! -----------------------
  call save_restart(Fv_restart, directory=trim(adjustl(datapath_out))//'RESTART')
  call free_restart_type(Fv_restart)


  !Write tracers to file
  !---------------------
  id_restart = register_restart_field( Tr_restart, filename_trcr, 'sphum', fld%Atm%q(:,:,:,1), &
                                       domain=fld%geom%domain )

  call save_restart(Tr_restart, directory=trim(adjustl(datapath_out))//'RESTART')
  call free_restart_type(Tr_restart)


  !Write date/time info in coupler.res
  !-----------------------------------
  iounit = 101
  if (fld%root_pe == 1) then
     print *,'write_file: date model init = ',fld%Atm%date_init
     print *,'write_file: date model now  = ',fld%Atm%date
     print *,'write_file: date vdate      = ',date
     open(iounit, file=trim(adjustl(datapath_out))//'RESTART/'//trim(adjustl(filename_cplr)), form='formatted')
     write( iounit, '(i6,8x,a)' ) fld%Atm%calendar_type, &
          '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
     write( iounit, '(6i6,8x,a)' )date, &
           'Model start time:   year, month, day, hour, minute, second'
     write( iounit, '(6i6,8x,a)' )date, &
           'Current model time: year, month, day, hour, minute, second'
     close(iounit)
  endif

return

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)
implicit none
type(fv3jedi_field), intent(in) :: fld
integer, intent(in) :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)

integer :: isc, iec, jsc, jec, gs

!1. Min
!2. Max
!3. RMS

isc = fld%geom%bd%isc
iec = fld%geom%bd%iec
jsc = fld%geom%bd%jsc
jec = fld%geom%bd%jec

gs = (iec-isc+1)*(jec-jsc+1)*fld%geom%npz

!u
pstat(1,1) = minval(fld%Atm%u(isc:iec,jsc:jec,:))
pstat(2,1) = maxval(fld%Atm%u(isc:iec,jsc:jec,:))
pstat(3,1) = sqrt((sum(fld%Atm%u(isc:iec,jsc:jec,:))/gs)**2)

!v
pstat(1,2) = minval(fld%Atm%v(isc:iec,jsc:jec,:))
pstat(2,2) = maxval(fld%Atm%v(isc:iec,jsc:jec,:))
pstat(3,2) = sqrt((sum(fld%Atm%v(isc:iec,jsc:jec,:))/gs)**2)

!pt
pstat(1,3) = minval(fld%Atm%pt(isc:iec,jsc:jec,:))
pstat(2,3) = maxval(fld%Atm%pt(isc:iec,jsc:jec,:))
pstat(3,3) = sqrt((sum(fld%Atm%pt(isc:iec,jsc:jec,:))/gs)**2)

!delp
pstat(1,4) = minval(fld%Atm%delp(isc:iec,jsc:jec,:))
pstat(2,4) = maxval(fld%Atm%delp(isc:iec,jsc:jec,:))
pstat(3,4) = sqrt((sum(fld%Atm%delp(isc:iec,jsc:jec,:))/gs)**2)

!q
pstat(1,5) = minval(fld%Atm%q(isc:iec,jsc:jec,:,1))
pstat(2,5) = maxval(fld%Atm%q(isc:iec,jsc:jec,:,1))
pstat(3,5) = sqrt((sum(fld%Atm%q(isc:iec,jsc:jec,:,1))/gs)**2)

return

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
use mpi,             only: mpi_real8,mpi_comm_world,mpi_sum,mpi_int
implicit none
type(fv3jedi_field), intent(in) :: fld
real(kind=kind_real), intent(out) :: prms
real(kind=kind_real) :: zz
integer jx,jy,jz,ii,nt,ierr,npes,iisum

zz = 0.0_kind_real
prms = 0.0_kind_real
ii = 0

do jy=fld%geom%bd%jsc,fld%geom%bd%jec
   do jx=fld%geom%bd%isc,fld%geom%bd%iec
      do jz=1,fld%geom%npz
         zz = zz + fld%Atm%pt(jx,jy,jz)**2
         zz = zz + fld%Atm%delp(jx,jy,jz)**2
         if (.not. fld%Atm%hydrostatic) then
            zz = zz + fld%Atm%w(jx,jy,jz)**2
            zz = zz + fld%Atm%delz(jx,jy,jz)**2
         endif
         do nt=1,fld%geom%ntracers
            zz = zz + fld%Atm%q(jx,jy,jz,nt)**2
         enddo
         ii = ii + 1
      enddo
   enddo
enddo

!Get global values
call mpi_allreduce(zz,prms,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
call mpi_allreduce(ii,iisum,1,mpi_int,mpi_sum,mpi_comm_world,ierr)

if (ierr .ne. 0) then
   print *,'error in fldrms/mpi_allreduce, error code=',ierr
endif
prms = sqrt(prms/real(iisum,kind_real))

!Print for debugging
if (fld%root_pe == 1) print *,'fldrms: prms = ', prms
if (fld%root_pe == 1) print *,'fldrms: iisum = ', iisum

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)

use iso_c_binding
use mpi

implicit none
type(fv3jedi_field), intent(inout) :: self
type(c_ptr), intent(in)       :: c_conf   !< Configuration
integer :: ndir,idir,ildir,ifdir,itiledir
integer,allocatable :: ixdir(:),iydir(:)
character(len=3) :: idirchar

! Get Diracs positions
ndir = config_get_int(c_conf,"ndir")
allocate(ixdir(ndir))
allocate(iydir(ndir))

do idir=1,ndir
   write(idirchar,'(i3)') idir
   ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
   iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
end do
ildir = config_get_int(c_conf,"ildir")
ifdir = config_get_int(c_conf,"ifdir")
itiledir = config_get_int(c_conf,"itiledir")


! Check
if (ndir<1) call abor1_ftn("fv3jedi_fields:dirac non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid ixdir")
endif
if (any(iydir<1).or.any(iydir>self%geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid iydir")
endif
if ((ildir<1).or.(ildir>self%geom%npz)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid ildir")
endif
if ((ifdir<1).or.(ifdir>5)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid ifdir")
endif
if ((itiledir<1).or.(itiledir>6)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid itiledir")
endif

! Setup Diracs
call zeros(self)

! only u,v,theta,delp and humidity allowed
do idir=1,ndir

   ! is specified grid point, tile number on this processor
   if (self%geom%ntile == itiledir .and. &
       ixdir(idir) >= self%geom%bd%isc .and. ixdir(idir) <= self%geom%bd%iec .and. &
       iydir(idir) >= self%geom%bd%jsc .and. iydir(idir) <= self%geom%bd%jec) then
       ! If so, perturb desired field and level
       if (ifdir == 1) then
          self%Atm%u(ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 2) then
          self%Atm%v(ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 3) then
          self%Atm%pt(ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 4) then
          self%Atm%delp(ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 5) then
          self%Atm%q(ixdir(idir),iydir(idir),ildir,1) = 1.0
       endif
   endif
end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine convert_to_ug(self, ug)
use unstructured_grid_mod
implicit none
type(fv3jedi_field), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug

integer :: nmga,ic0a,jx,jy,jl,nf
integer,allocatable :: imask(:,:)
real(kind=kind_real),allocatable :: lon(:),lat(:),area(:),vunit(:,:)
real(kind=kind_real) :: sigmaup,sigmadn

! Define local number of gridpoints
nmga = (self%geom%bd%iec - self%geom%bd%isc + 1) * (self%geom%bd%jec - self%geom%bd%jsc + 1)

! Allocation
allocate(lon(nmga))
allocate(lat(nmga))
allocate(area(nmga))
allocate(vunit(nmga,self%geom%npz))
allocate(imask(nmga,self%geom%npz))

! Copy coordinates
ic0a = 0
do jy=self%geom%bd%jsc,self%geom%bd%jec
  do jx=self%geom%bd%isc,self%geom%bd%iec
    ic0a = ic0a+1
    lon(ic0a) = deg2rad*self%geom%grid_lon(jx,jy)
    lat(ic0a) = deg2rad*self%geom%grid_lat(jx,jy)
    area(ic0a) = self%geom%area(jx,jy)
  enddo
enddo
imask = 1

! Define vertical unit
do jl=1,self%geom%npz
  sigmaup = self%geom%ak(jl+1)/101300.0+self%geom%bk(jl+1) ! si are now sigmas
  sigmadn = self%geom%ak(jl  )/101300.0+self%geom%bk(jl  )
  vunit(:,jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
enddo

! Should this come from self/vars?
nf = 5
if (.not. self%geom%hydrostatic) nf = 7

! Create unstructured grid
call create_unstructured_grid(ug, nmga, self%geom%npz, nf, 1, lon, lat, area, vunit, imask)

! Copy field
ic0a = 0
do jy=self%geom%bd%jsc,self%geom%bd%jec
  do jx=self%geom%bd%isc,self%geom%bd%iec
    ic0a = ic0a+1
    do jl=1,self%geom%npz
        ug%fld(ic0a,jl,1,1) = self%Atm%u(jx,jy,jl)
        ug%fld(ic0a,jl,2,1) = self%Atm%v(jx,jy,jl)
        ug%fld(ic0a,jl,3,1) = self%Atm%pt(jx,jy,jl)
        ug%fld(ic0a,jl,4,1) = self%Atm%q(jx,jy,jl,1)
        ug%fld(ic0a,jl,5,1) = self%Atm%delp(jx,jy,jl)
        if (.not. self%geom%hydrostatic) ug%fld(ic0a,jl,6,1) = self%Atm%w(jx,jy,jl)
        if (.not. self%geom%hydrostatic) ug%fld(ic0a,jl,7,1) = self%Atm%delz(jx,jy,jl)
    enddo
  enddo
enddo

end subroutine convert_to_ug

! -----------------------------------------------------------------------------

subroutine convert_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(fv3jedi_field), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug

integer :: ic0a,jx,jy,jl

! Copy field
ic0a = 0
do jy=self%geom%bd%jsc,self%geom%bd%jec
  do jx=self%geom%bd%isc,self%geom%bd%iec
    ic0a = ic0a+1
    do jl=1,self%geom%npz
        self%Atm%u(jx,jy,jl)    = ug%fld(ic0a,jl,1,1)
        self%Atm%v(jx,jy,jl)    = ug%fld(ic0a,jl,2,1)
        self%Atm%pt(jx,jy,jl)   = ug%fld(ic0a,jl,3,1)
        self%Atm%q(jx,jy,jl,1)  = ug%fld(ic0a,jl,4,1)
        self%Atm%delp(jx,jy,jl) = ug%fld(ic0a,jl,5,1)
        if (.not. self%geom%hydrostatic) self%Atm%w(jx,jy,jl) = ug%fld(ic0a,jl,6,1)
        if (.not. self%geom%hydrostatic) self%Atm%delz(jx,jy,jl) = ug%fld(ic0a,jl,7,1)
    enddo
  enddo
enddo

end subroutine convert_from_ug

! ------------------------------------------------------------------------------

subroutine interp(fld, locs, vars, gom)

use type_bump, only: bump_type
use variable_transforms

implicit none
type(fv3jedi_field), intent(inout) :: fld 
type(ioda_locs),     intent(in)     :: locs 
type(ufo_vars),  intent(in)        :: vars
type(ufo_geovals),  intent(inout)  :: gom

character(len=*), parameter :: myname = 'interp'

type(bump_type), pointer :: pbump

integer :: ii, jj, ji, jvar, jlev, ngrid, nobs
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)
real(kind=kind_real), allocatable :: geoval(:,:,:)

! Initialize the interpolation
! ----------------------------
call initialize_interp( fld, locs, vars, gom, myname, &
                        pbump, ngrid, nobs )

! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(ngrid,1))
allocate(obs_field(nobs,1))

! Interpolate fields to obs locations using pre-calculated weights
! ----------------------------------------------------------------

!Temporary variable in case of variable transform
allocate(geoval(fld%geom%bd%isd:fld%geom%bd%ied,fld%geom%bd%jsd:fld%geom%bd%jed,fld%geom%npz))

!write(*,*)'interp model    t min, max= ',minval(fld%Atm%pt),maxval(fld%Atm%pt)
!write(*,*)'interp model delp min, max= ',minval(fld%Atm%delp),maxval(fld%Atm%delp)

do jvar = 1, vars%nv

  select case (trim(vars%fldnames(jvar)))

  case ("wind_u")

  case ("wind_v")
                
  case ("virtual_temperature")

    !Convert temperature to virtual temperature
    call T_to_Tv(fld%geom,fld%Atm%pt,fld%Atm%q,geoval)

    do jlev = 1, fld%geom%npz
      ii = 0
      do jj = fld%geom%bd%jsc, fld%geom%bd%jec
        do ji = fld%geom%bd%isc, fld%geom%bd%iec
          ii = ii + 1
          mod_field(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call pbump%apply_obsop(mod_field,obs_field)
      gom%geovals(jvar)%vals(jlev,:) = obs_field(:,1)
    enddo

  case ("humidity_mixing_ratio")

  case ("atmosphere_ln_pressure_coordinate")

    !Convert pressure thickness to log pressure level mid point
    call delp_to_logP(fld%geom,fld%Atm%delp,geoval)

    geoval = log(0.001) + geoval !Convert to kPa

    do jlev = 1, fld%geom%npz
      ii = 0
      do jj = fld%geom%bd%jsc, fld%geom%bd%jec
        do ji = fld%geom%bd%isc, fld%geom%bd%iec
          ii = ii + 1
          mod_field(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call pbump%apply_obsop(mod_field,obs_field)
      gom%geovals(jvar)%vals(jlev,:) = obs_field(:,1)
    enddo

  case ("air_pressure")

  case ("air_pressure_levels")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

enddo

deallocate(geoval)
deallocate(mod_field)
deallocate(obs_field)

!write(*,*)'interp geovals t min, max= ',minval(gom%geovals(1)%vals(:,:)),maxval(gom%geovals(1)%vals(:,:))
!write(*,*)'interp geovals p min, max= ',minval(gom%geovals(2)%vals(:,:)),maxval(gom%geovals(2)%vals(:,:))

end subroutine interp

! ------------------------------------------------------------------------------

subroutine interp_tl(fld, locs, vars, gom)!, traj)

use type_bump, only: bump_type
use variable_transforms
use fv3jedi_trajectories, only: fv3jedi_trajectory

implicit none
type(fv3jedi_field),      intent(inout) :: fld 
type(ioda_locs),           intent(in)    :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_trajectory), pointer       :: traj

character(len=*), parameter :: myname = 'interp_tl'

type(bump_type), pointer :: pbump

integer :: ii, jj, ji, jvar, jlev, ngrid, nobs
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)

! Initialize the interpolation
! ----------------------------
call initialize_interp( fld, locs, vars, gom, myname, &
                        pbump, ngrid, nobs, traj )

! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(ngrid,1))
allocate(obs_field(nobs,1))

!write(*,*)'interp_tl model    t min, max= ',minval(fld%Atm%pt),maxval(fld%Atm%pt)
!write(*,*)'interp_tl model delp min, max= ',minval(fld%Atm%delp),maxval(fld%Atm%delp)

! Interpolate fields to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv
 
  gom%geovals(jvar)%vals(:,:) = 0.0

  select case (trim(vars%fldnames(jvar)))
   
  case ("wind_v")
  
  case ("virtual_temperature")

    !Tangent linear of temperature to virtual temperature
    call T_to_Tv_tl(fld%geom, traj%pt, fld%Atm%pt, traj%q(:,:,:,1), fld%Atm%q (:,:,:,1) )

    do jlev = 1, fld%geom%npz
      ii = 0
      do jj = fld%geom%bd%jsc, fld%geom%bd%jec
        do ji = fld%geom%bd%isc, fld%geom%bd%iec
          ii = ii + 1
          mod_field(ii, 1) = fld%Atm%pt(ji, jj, jlev)
        enddo
      enddo
      call pbump%apply_obsop(mod_field,obs_field)
      gom%geovals(jvar)%vals(jlev,:) = obs_field(:,1)
    enddo
  
  case ("humidity_mixing_ratio")
  
  case ("atmosphere_ln_pressure_coordinate")

  case ("air_pressure")

  case ("air_pressure_levels")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

enddo


deallocate(mod_field)
deallocate(obs_field)

!write(*,*)'interp_tl geovals t min, max= ',minval(gom%geovals(1)%vals(:,:)),maxval(gom%geovals(1)%vals(:,:))
!write(*,*)'interp_tl geovals p min, max= ',minval(gom%geovals(2)%vals(:,:)),maxval(gom%geovals(2)%vals(:,:))

end subroutine interp_tl

! ------------------------------------------------------------------------------

subroutine interp_ad(fld, locs, vars, gom)!, traj)

use type_bump, only: bump_type
use variable_transforms
use fv3jedi_trajectories, only: fv3jedi_trajectory

implicit none
type(fv3jedi_field),      intent(inout) :: fld 
type(ioda_locs),           intent(in)    :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_trajectory), pointer       :: traj

character(len=*), parameter :: myname = 'interp_ad'

type(bump_type), pointer :: pbump

integer :: ii, jj, ji, jvar, jlev, ngrid, nobs
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)

! Initialize the interpolation
! ----------------------------
call initialize_interp( fld, locs, vars, gom, myname, &
                        pbump,  ngrid, nobs, traj )

! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(ngrid,1))
allocate(obs_field(nobs,1))

!write(*,*)'interp_ad geovals t min, max= ',minval(gom%geovals(1)%vals(:,:)),maxval(gom%geovals(1)%vals(:,:))
!write(*,*)'interp_ad geovals p min, max= ',minval(gom%geovals(2)%vals(:,:)),maxval(gom%geovals(2)%vals(:,:))

! Interpolate fields to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv
 
  select case (trim(vars%fldnames(jvar)))
 
  case ("wind_u")
  
  case ("wind_v")
  
  case ("virtual_temperature")

    do jlev = 1, fld%geom%npz
      obs_field(:,1) = gom%geovals(jvar)%vals(jlev,:)
      call pbump%apply_obsop_ad(obs_field,mod_field)
      ii = 0
      do jj = fld%geom%bd%jsc, fld%geom%bd%jec
        do ji = fld%geom%bd%isc, fld%geom%bd%iec
          ii = ii + 1
          fld%Atm%pt(ji, jj, jlev) = fld%Atm%pt(ji, jj, jlev) + mod_field(ii, 1)
        enddo
      enddo
    enddo 
  
    !Adjoint of temperature to virtual temperature
    call T_to_Tv_ad(fld%geom, traj%pt, fld%Atm%pt, traj%q(:,:,:,1), fld%Atm%q (:,:,:,1) )
  
  case ("humidity_mixing_ratio")
  
  case ("atmosphere_ln_pressure_coordinate")

  case ("air_pressure")

  case ("air_pressure_levels")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

enddo

deallocate(mod_field)
deallocate(obs_field)

!write(*,*)'interp_ad model    t min, max= ',minval(fld%Atm%pt),maxval(fld%Atm%pt)
!write(*,*)'interp_ad model delp min, max= ',minval(fld%Atm%delp),maxval(fld%Atm%delp)

end subroutine interp_ad

! ------------------------------------------------------------------------------

subroutine initialize_interp( fld, locs, vars, gom, myname, &
                                   pbump, ngrid, nobs, trajp )

use type_bump, only: bump_type
use variable_transforms
use fv3jedi_trajectories, only: fv3jedi_trajectory

implicit none
type(fv3jedi_field),                intent(inout) :: fld 
type(ioda_locs),                     intent(in   ) :: locs 
type(ufo_vars),                     intent(in   ) :: vars
type(ufo_geovals),                  intent(inout) :: gom
type(fv3jedi_trajectory), optional, pointer, intent(inout) :: trajp
character(len=*),                   intent(in   ) :: myname
type(bump_type), pointer,           intent(inout) :: pbump
integer,                            intent(  out) :: ngrid, nobs
integer :: jvar

!HACK: read a trajectory from file
!---------------------------------

  type(fv3jedi_trajectory), save, target :: traj
  type(restart_file_type) :: Fv_restart
  type(restart_file_type) :: Tr_restart
  integer :: id_restart
  character(len=255) :: datapath_in, datapath_ti
  character(len=255) :: filename_core
  character(len=255) :: filename_trcr
  integer, save :: gottraj = 0

  if (present(trajp)) then

    if (gottraj == 0) then

      if (fld%geom%npx == 97) then

         if (fld%root_pe == 1) print*, 'HACK TO PROVIDE TRAJECOTORY FOR INTERPOLATION (C96)'

         datapath_in = 'Data/C96_RESTART_2016-01-01-06/'
         datapath_ti = 'Data/C96_RESTART_2016-01-01-06/INPUT/'

         filename_core = 'fv_core.res.nc'
         filename_trcr = 'fv_tracer.res.nc'

      elseif (fld%geom%npx == 49) then

         if (fld%root_pe == 1) print*, 'HACK TO PROVIDE TRAJECOTORY FOR INTERPOLATION (C48)'

         datapath_in = 'Data/C48_RESTART_2017-08-01-00/ENSEMBLE/mem001/RESTART/'
         datapath_ti = 'Data/C48_RESTART_2017-08-01-00/ENSEMBLE/mem001/RESTART/'

         filename_core = '20170801.000000.fv_core.res.nc'
         filename_trcr = '20170801.000000.fv_tracer.res.nc'

      else

         call abor1_ftn("Resolution not supported for trajectory hack")

      endif

      allocate(traj%pt  (fld%geom%bd%isd:fld%geom%bd%ied,fld%geom%bd%jsd:fld%geom%bd%jed,fld%geom%npz))
      allocate(traj%delp(fld%geom%bd%isd:fld%geom%bd%ied,fld%geom%bd%jsd:fld%geom%bd%jed,fld%geom%npz))
      allocate(traj%q   (fld%geom%bd%isd:fld%geom%bd%ied,fld%geom%bd%jsd:fld%geom%bd%jed,fld%geom%npz,fld%geom%ntracers))

      id_restart = register_restart_field(Fv_restart, filename_core, 'T',    traj%pt,   domain=fld%geom%domain)
      id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', traj%delp, domain=fld%geom%domain)

      call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
      call free_restart_type(Fv_restart)

      id_restart =  register_restart_field(Tr_restart, filename_trcr, 'sphum', traj%q(:,:,:,1), domain=fld%geom%domain)

      call restore_state(Tr_restart, directory=trim(adjustl(datapath_ti)))
      call free_restart_type(Tr_restart) 

      trajp => traj

      gottraj = 1

    else

      trajp => traj

    endif

  endif
!END HACK
!--------

! Get grid dimensions and checks
! ------------------------------
ngrid = (fld%geom%bd%iec - fld%geom%bd%isc + 1)*(fld%geom%bd%jec - fld%geom%bd%jsc + 1)
nobs = locs%nlocs 

! Make sure the return values are allocated and set
! -------------------------------------------------
if (trim(myname)/="interp_ad") then
   do jvar=1,vars%nv
      if (.not.allocated(gom%geovals(jvar)%vals)) then
         if (fld%root_pe == 1) print*, trim(myname), ': allocating ', trim(vars%fldnames(jvar)),' GeoVaLs'
         gom%geovals(jvar)%nval = fld%geom%npz
         gom%geovals(jvar)%nobs = nobs
         allocate( gom%geovals(jvar)%vals(fld%geom%npz,nobs) )
         gom%geovals(jvar)%vals = 0.0
      endif
   enddo
   gom%linit = .true.
   gom%lalloc = .true.
endif

!Run some basic checks on the interpolation
!------------------------------------------
call interp_checks(myname, fld, locs, vars, gom)

! Calculate interpolation weight using BUMP
! -----------------------------------------
call initialize_bump(fld, fld%geom, locs, pbump)

end subroutine initialize_interp

! ------------------------------------------------------------------------------

subroutine initialize_bump(fld, grid, locs, pbump)

use fv3jedi_geom_mod, only: fv3jedi_geom
use type_bump, only: bump_type
use mpi, only: mpi_comm_world

implicit none
type(fv3jedi_field), intent(in) :: fld
type(fv3jedi_geom), intent(in) :: grid
type(ioda_locs), intent(in)    :: locs
type(bump_type), pointer, intent(out) :: pbump

logical, save :: interp_initialized = .FALSE.
type(bump_type), save, target :: bump

integer :: mod_nx,mod_ny,mod_nz,mod_num,obs_num
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 

real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

integer :: ii, jj, ji, jvar, jlev

!Get the Solution dimensions
!---------------------------
mod_nx  = grid%bd%iec - grid%bd%isc + 1  
mod_ny  = grid%bd%jec - grid%bd%jsc + 1
mod_num = mod_nx * mod_ny
mod_nz  = grid%npz
obs_num = locs%nlocs 

!Calculate interpolation weight using BUMP
!-----------------------------------------
if (.NOT.interp_initialized) then

   write(*,*)'initialize_bump mod_num,obs_num = ',mod_num,obs_num

   allocate( mod_lat(mod_num), mod_lon(mod_num) )
   mod_lat = reshape( grid%grid_lat(grid%bd%isc:grid%bd%iec,      &
                                    grid%bd%jsc:grid%bd%jec),     &
                                   [mod_num] )  
   mod_lon = reshape( grid%grid_lon(grid%bd%isc:grid%bd%iec,      &
                                    grid%bd%jsc:grid%bd%jec),     &
                                   [mod_num] )

   !Important namelist options
   bump%nam%prefix = 'oops_data'   ! Prefix for files output
   bump%nam%nobs = obs_num         ! Number of observations
   bump%nam%obsop_interp = 'bilin' ! Interpolation type (bilinear)
   bump%nam%obsdis = 'local'       ! Observation distribution parameter ('random','local' or 'adjusted')
   bump%nam%diag_interp = 'bilin'

   !Less important namelist options (should not be changed)
   bump%nam%default_seed = .true.
   bump%nam%new_hdiag = .false.
   bump%nam%new_param = .false.
   bump%nam%check_adjoints = .false.
   bump%nam%check_pos_def = .false.
   bump%nam%check_sqrt = .false.
   bump%nam%check_dirac = .false.
   bump%nam%check_randomization = .false.
   bump%nam%check_consistency = .false.
   bump%nam%check_optimality = .false.
   bump%nam%new_lct = .false.
   bump%nam%new_obsop = .true.

   !Initialize geometry
   allocate(area(mod_num))
   allocate(vunit(mod_num,1))
   allocate(lmask(mod_num,1))
   area = 1.0           ! Dummy area
   vunit = 1.0          ! Dummy vertical unit
   lmask = .true.       ! Mask

   !Initialize BUMP
   call bump%setup_online( mpi_comm_world,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                           nobs=obs_num,lonobs=locs%lon(:),latobs=locs%lat(:) )

   !Release memory
   deallocate(area)
   deallocate(vunit)
   deallocate(lmask)
   deallocate( mod_lat, mod_lon )

   interp_initialized = .TRUE. 

endif

pbump => bump

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)
implicit none
character(len=*), intent(in) :: cop
type(fv3jedi_field), intent(in) :: fld
type(ioda_locs), intent(in)    :: locs
type(ufo_vars), intent(in)    :: vars
type(ufo_geovals), intent(in) :: gom

integer :: jvar
character(len=255) :: cinfo

cinfo="fv3jedi_fields:checks "//trim(cop)//" : "

!Check things are the sizes we expect
!------------------------------------
if (gom%nobs /= locs%nlocs ) then
   call abor1_ftn(trim(cinfo)//"geovals wrong size")
endif
if( gom%nvar .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//"nvar wrong size")
endif
if( .not. allocated(gom%geovals) )then
   call abor1_ftn(trim(cinfo)//"geovals unallocated")
endif
if( size(gom%geovals) .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//"geovals wrong size")
endif
if (trim(cop)/="interp_tl" .and. trim(cop)/="interp_ad") then
if (.not.gom%linit) then
   call abor1_ftn(trim(cinfo)//"geovals not initialized")
endif
do jvar=1,vars%nv
   if (allocated(gom%geovals(jvar)%vals)) then  
      if( gom%geovals(jvar)%nval .ne. fld%geom%npz )then
         call abor1_ftn(trim(cinfo)//"nval wrong size")
      endif
      if( gom%geovals(jvar)%nobs .ne. locs%nlocs )then
         call abor1_ftn(trim(cinfo)//"nobs wrong size")
      endif
      if( size(gom%geovals(jvar)%vals, 1) .ne. fld%geom%npz )then
         call abor1_ftn(trim(cinfo)//"vals wrong size 1")
      endif
      if( size(gom%geovals(jvar)%vals, 2) .ne. locs%nlocs )then
         call abor1_ftn(trim(cinfo)//"vals wrong size 2")
      endif       
   else
     call abor1_ftn(trim(cinfo)//"vals not allocated")
   endif 
enddo
endif

write(*,*)'interp_checks ',trim(cinfo),' done'

end subroutine interp_checks

! ------------------------------------------------------------------------------

subroutine check_resolution(x1, x2)

implicit none
type(fv3jedi_field), intent(in) :: x1, x2

if (x1%geom%size_cubic_grid /= x2%geom%size_cubic_grid .or.  x1%geom%npz /= x2%geom%npz) then
  call abor1_ftn ("fv3jedi_fields: resolution error")
endif
call check(x1)
call check(x2)

end subroutine check_resolution

! ------------------------------------------------------------------------------

subroutine check(self)
implicit none
type(fv3jedi_field), intent(in) :: self
logical :: bad

bad = .false.

if (bad) then
   call abor1_ftn ("fv3jedi_fields: field not consistent")
endif

end subroutine check

! ------------------------------------------------------------------------------

end module fv3jedi_fields_mod
