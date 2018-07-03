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

use mpp_domains_mod,    only: EAST, NORTH
use mpp_mod,            only: mpp_pe, mpp_npes, mpp_root_pe
use fms_io_mod,         only: restart_file_type, register_restart_field, &
                              free_restart_type, restore_state, save_restart
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index, get_tracer_names

use fv3jedi_fields_utils_mod
use fv3jedi_fields_io_mod
use fv3jedi_constants, only: deg2rad, constoz

use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj

implicit none
private

public :: create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, fldrms, &
          change_resol, getvalues, getvalues_tl, getvalues_ad, &
          convert_to_ug, convert_from_ug, dirac, &
          analytic_IC

public :: fv3jedi_field

public :: fv3jedi_field_registry

! ------------------------------------------------------------------------------

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

! Index in q for the required tracers
self%ti_q  = get_tracer_index (MODEL_ATMOS, 'sphum')
self%ti_ql = get_tracer_index (MODEL_ATMOS, 'liq_wat')
self%ti_qi = get_tracer_index (MODEL_ATMOS, 'ice_wat')
self%ti_o3 = get_tracer_index (MODEL_ATMOS, 'o3mr')

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

if (allocated(self%Atm%slmsk )) deallocate(self%Atm%slmsk )
if (allocated(self%Atm%sheleg)) deallocate(self%Atm%sheleg)
if (allocated(self%Atm%tsea  )) deallocate(self%Atm%tsea  )
if (allocated(self%Atm%vtype )) deallocate(self%Atm%vtype )
if (allocated(self%Atm%stype )) deallocate(self%Atm%stype )
if (allocated(self%Atm%vfrac )) deallocate(self%Atm%vfrac )
if (allocated(self%Atm%stc   )) deallocate(self%Atm%stc   )
if (allocated(self%Atm%smc   )) deallocate(self%Atm%smc   )
if (allocated(self%Atm%u_srf )) deallocate(self%Atm%snwdph)
if (allocated(self%Atm%u_srf )) deallocate(self%Atm%u_srf )
if (allocated(self%Atm%v_srf )) deallocate(self%Atm%v_srf )
if (allocated(self%Atm%f10m  )) deallocate(self%Atm%f10m  )

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

self%Atm%hydrostatic = rhs%Atm%hydrostatic

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

self%Atm%ua(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%ua(self%isc:self%iec,self%jsc:self%jec,:)
self%Atm%va(self%isc:self%iec,self%jsc:self%jec,:) = rhs%Atm%va(self%isc:self%iec,self%jsc:self%jec,:)

self%Atm%calendar_type = rhs%Atm%calendar_type
self%Atm%date = rhs%Atm%date
self%Atm%date_init = rhs%Atm%date_init

self%Atm%slmsk  = rhs%Atm%slmsk
self%Atm%sheleg = rhs%Atm%sheleg
self%Atm%tsea   = rhs%Atm%tsea
self%Atm%vtype  = rhs%Atm%vtype
self%Atm%stype  = rhs%Atm%stype
self%Atm%vfrac  = rhs%Atm%vfrac
self%Atm%stc    = rhs%Atm%stc
self%Atm%smc    = rhs%Atm%smc
self%Atm%snwdph = rhs%Atm%snwdph
self%Atm%u_srf  = rhs%Atm%u_srf
self%Atm%v_srf  = rhs%Atm%v_srf
self%Atm%f10m   = rhs%Atm%f10m

self%geom => rhs%geom

self%nf = rhs%nf

self%isc = rhs%isc
self%iec = rhs%iec
self%jsc = rhs%jsc
self%jec = rhs%jec

self%root_pe = rhs%root_pe

self%havecrtmfields = rhs%havecrtmfields

self%ti_q  = rhs%ti_q 
self%ti_ql = rhs%ti_ql
self%ti_qi = rhs%ti_qi
self%ti_o3 = rhs%ti_o3


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

!u,v,T,delp
do k = 1,fld1%geom%npz
  do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
    do i = fld1%geom%bd%isc,fld1%geom%bd%iec
      zp = zp + fld1%Atm%u   (i,j,k) * fld2%Atm%u   (i,j,k)
      zp = zp + fld1%Atm%v   (i,j,k) * fld2%Atm%v   (i,j,k)
      zp = zp + fld1%Atm%pt  (i,j,k) * fld2%Atm%pt  (i,j,k)
      zp = zp + fld1%Atm%delp(i,j,k) * fld2%Atm%delp(i,j,k)
    enddo
  enddo
enddo

!Tracers
do l = 1,fld1%geom%ntracers
  do k = 1,fld1%geom%npz
    do j = fld1%geom%bd%jsc,fld1%geom%bd%jec
      do i = fld1%geom%bd%isc,fld1%geom%bd%iec
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
!! * dcmip-test-3-1: Non-hydrostatic gravity wave
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \author M. Miesch (adapted from a pre-existing call to invent_state)
!! \date March, 2018: Created
!! \date May, 2018: Added dcmip-test-3-1
!! \date June, 2018: Added dcmip-test-4-0
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
  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

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

  ! Pointer to geometry component of field object
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
              Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
                                               phis0,ps,rho0,hum0,q1,q2,q3,q4)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                                                  phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 fld%Atm%u(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 fld%Atm%v(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
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
              Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                          t0,phis0,ps,rho0,hum0,q1)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                             t0,phis0,ps,rho0,hum0,q1)

                 fld%Atm%u(i,j,k) = u0 !ATTN comment above
                 fld%Atm%v(i,j,k) = v0
                 If (.not.fld%Atm%hydrostatic) fld%Atm%w(i,j,k) = w0
                 fld%Atm%pt(i,j,k) = t0
                 fld%Atm%delp(i,j,k) = pe2-pe1
                 If (geom%ntracers >= 1) fld%Atm%q(i,j,k,1) = hum0
                 If (geom%ntracers >= 2) fld%Atm%q(i,j,k,2) = q1
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do i = geom%bd%isc,geom%bd%iec
           do j = geom%bd%jsc,geom%bd%jec
              rlat = deg_to_rad*geom%grid_lat(i,j)
              rlon = deg_to_rad*geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                      t0,phis0,ps,rho0,hum0)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0)

                 fld%Atm%u(i,j,k) = u0 !ATTN comment above
                 fld%Atm%v(i,j,k) = v0
                 If (.not.fld%Atm%hydrostatic) fld%Atm%w(i,j,k) = w0
                 fld%Atm%pt(i,j,k) = t0
                 fld%Atm%delp(i,j,k) = pe2-pe1
                 If (geom%ntracers >= 1) fld%Atm%q(i,j,k,1) = hum0
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do i = geom%bd%isc,geom%bd%iec
           do j = geom%bd%jsc,geom%bd%jec
              rlat = deg_to_rad*geom%grid_lat(i,j)
              rlon = deg_to_rad*geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

              fld%Atm%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

                 fld%Atm%u(i,j,k) = u0 !ATTN comment above
                 fld%Atm%v(i,j,k) = v0
                 If (.not.fld%Atm%hydrostatic) fld%Atm%w(i,j,k) = w0
                 fld%Atm%pt(i,j,k) = t0
                 fld%Atm%delp(i,j,k) = pe2-pe1
                 If (geom%ntracers >= 1) fld%Atm%q(i,j,k,1) = hum0
                 
              enddo
           enddo
        enddo

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
do k = 1,flds%geom%npz
  do j = flds%geom%bd%jsc,flds%geom%bd%jec
    do i = flds%geom%bd%isc,flds%geom%bd%iec
      flds%Atm%u(i,j,k) = cos(0.25*flds%geom%grid_lon(i,j)) + cos(0.25*flds%geom%grid_lat(i,j))
    enddo
  enddo
enddo

!v
do k = 1,flds%geom%npz
  do j = flds%geom%bd%jsc,flds%geom%bd%jec
    do i = flds%geom%bd%isc,flds%geom%bd%iec
      flds%Atm%v(i,j,k) = 1.0_kind_real
    enddo
  enddo
enddo

!pt
do k = 1,flds%geom%npz
  do j = flds%geom%bd%jsc,flds%geom%bd%jec
    do i = flds%geom%bd%isc,flds%geom%bd%iec
      flds%Atm%pt(i,j,k) = cos(0.25*flds%geom%grid_lon(i,j)) + cos(0.25*flds%geom%grid_lat(i,j))
    enddo
  enddo
enddo

!delp
do k = 1,flds%geom%npz
  do j = flds%geom%bd%jsc,flds%geom%bd%jec
    do i = flds%geom%bd%isc,flds%geom%bd%iec
      flds%Atm%delp(i,j,k) = k
    enddo
  enddo
enddo

!q
do k = 1,flds%geom%npz
  do j = flds%geom%bd%jsc,flds%geom%bd%jec
    do i = flds%geom%bd%isc,flds%geom%bd%iec
      flds%Atm%q(i,j,k,1) = 0.0
    enddo
  enddo
enddo

return
end subroutine invent_state

! ------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)

  implicit none

  type(fv3jedi_field), intent(inout) :: fld      !< Fields
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'fv3gfs') then
     call read_fms_restart(fld, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call read_geos_restart(fld, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_fields read: restart type not supported")
  endif

  return

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)

  implicit none

  type(fv3jedi_field), intent(in)    :: fld      !< Fields
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'fv3gfs') then
     call write_fms_restart(fld, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call write_geos_restart(fld, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_fields write: restart type not supported")
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

subroutine getvalues(fld, locs, vars, gom, traj)

use surface_vt_mod
use pressure_vt_mod
use tmprture_vt_mod
use moisture_vt_mod, only: crtm_ade_efr, crtm_mixratio
use wind_vt_mod
use field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod,only: get_tracer_index, get_tracer_names
use type_bump, only: bump_type

implicit none
type(fv3jedi_field),                        intent(inout) :: fld 
type(ioda_locs),                            intent(in)    :: locs 
type(ufo_vars),                             intent(in)    :: vars
type(ufo_geovals),                          intent(inout) :: gom
type(fv3jedi_getvaltraj), optional, target, intent(inout) :: traj

character(len=*), parameter :: myname = 'interp'

type(bump_type), target  :: bump
type(bump_type), pointer :: pbump
logical, target :: bump_alloc
logical, pointer :: pbumpa

integer :: ii, jj, ji, jvar, jlev, ngrid, nobs
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
integer :: nvl
logical :: do_interp

integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,npz,i,j

integer :: nt, trcount
character(len=20) :: trname

!Local pressure variables
real(kind=kind_real), allocatable :: prsi(:,:,:) !Pressure Pa, interfaces
real(kind=kind_real), allocatable :: prs (:,:,:) !Pressure Pa, midpoint
real(kind=kind_real), allocatable :: logp(:,:,:) !Log(pressue), (Pa) midpoint

!Local CRTM moisture variables
real(kind=kind_real), allocatable :: ql_ade(:,:,:) !Cloud liq water kgm^2
real(kind=kind_real), allocatable :: qi_ade(:,:,:) !Cloud ice water kgm^2
real(kind=kind_real), allocatable :: ql_efr(:,:,:) !Cloud effective radius microns
real(kind=kind_real), allocatable :: qi_efr(:,:,:) !Cloud effective radium microns
real(kind=kind_real), allocatable :: qmr(:,:,:)    !Moisture mixing ratio
real(kind=kind_real), allocatable :: water_coverage_m(:,:) !Water coverage, model grid

!Local CRTM surface variables
integer             , allocatable :: vegetation_type(:)          !Index of vege type              | surface(1)%Vegetation_Type
integer             , allocatable :: land_type(:)                !Index of land type              | surface(1)%Land_Type
integer             , allocatable :: soil_type(:)                !Index of soil type              | surface(1)%Soil_Type
real(kind=kind_real), allocatable :: wind_speed(:)               !10 meter wind speed m/s         | surface(1)%wind_speed
real(kind=kind_real), allocatable :: wind_direction(:)           !10 meter wind direction degrees | surface(1)%wind_direction
real(kind=kind_real), allocatable :: water_coverage(:)           !Fraction of water coverage      | surface(1)%water_coverage
real(kind=kind_real), allocatable :: land_coverage(:)            !Fraction of land coverage       | surface(1)%land_coverage
real(kind=kind_real), allocatable :: ice_coverage(:)             !Fraction of ice coverage        | surface(1)%ice_coverage
real(kind=kind_real), allocatable :: snow_coverage(:)            !Fraction of snow coverage       | surface(1)%snow_coverage
real(kind=kind_real), allocatable :: lai(:)                      !Leaf area index                 ! surface(1)%lai
real(kind=kind_real), allocatable :: water_temperature(:)        !Water temp (K)                  | surface(1)%water_temperature    
real(kind=kind_real), allocatable :: land_temperature(:)         !Land temp (K)                   | surface(1)%land_temperature     
real(kind=kind_real), allocatable :: ice_temperature(:)          !Ice temp (K)                    | surface(1)%ice_temperature      
real(kind=kind_real), allocatable :: snow_temperature(:)         !Snow temp (K)                   | surface(1)%snow_temperature     
real(kind=kind_real), allocatable :: soil_moisture_content(:)    !Soil moisture content           | surface(1)%soil_moisture_content
real(kind=kind_real), allocatable :: vegetation_fraction(:)      !Vegetation fraction             | surface(1)%vegetation_fraction  
real(kind=kind_real), allocatable :: soil_temperature(:)         !Soil temperature                | surface(1)%soil_temperature     
real(kind=kind_real), allocatable :: snow_depth(:)               !Snow depth                      | surface(1)%snow_depth           


! Grid convenience
! ----------------
isc = fld%geom%bd%isc
iec = fld%geom%bd%iec
jsc = fld%geom%bd%jsc
jec = fld%geom%bd%jec
isd = fld%geom%bd%isd
ied = fld%geom%bd%ied
jsd = fld%geom%bd%jsd
jed = fld%geom%bd%jed
npz = fld%geom%npz

ngrid = (iec-isc+1)*(jec-jsc+1)
nobs = locs%nlocs 

! Initialize the interpolation
! ----------------------------
if (present(traj)) then

  pbump => traj%bump(1)

  if (.not. traj%lalloc) then
  
     traj%ngrid = ngrid
     traj%nobs = nobs
   
     if (.not.allocated(traj%pt)) allocate(traj%pt(isd:ied,jsd:jed,1:npz))
     if (.not.allocated(traj%q )) allocate(traj%q (isd:ied,jsd:jed,1:npz,1:fld%geom%ntracers))
  
     traj%pt = fld%Atm%pt
     traj%q  = fld%Atm%q
 
     pbumpa => traj%lalloc

  endif

else

  pbump => bump
  bump_alloc = .false.
  pbumpa => bump_alloc

endif

if (.not. pbumpa) then
   call initialize_bump(fld%geom, locs, pbump)
   pbumpa = .true.
endif


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(ngrid,1))
allocate(obs_field(nobs,1))

! Local GeoVals
! -------------
allocate(geovale(isd:ied,jsd:jed,npz+1))
allocate(geovalm(isd:ied,jsd:jed,npz))

! Get pressures at edge, center & log center
! ------------------------------------------
allocate(prsi(isd:ied,jsd:jed,npz+1))
allocate(prs (isd:ied,jsd:jed,npz  ))
allocate(logp(isd:ied,jsd:jed,npz  ))

call delp_to_pe_p_logp(fld%geom,fld%Atm%delp,prsi,prs,logp)

! Get CRTM surface variables
! ----------------------
allocate(wind_speed(nobs))
allocate(wind_direction(nobs))
allocate(land_type(nobs))
allocate(vegetation_type(nobs))
allocate(soil_type(nobs))
allocate(water_coverage(nobs))
allocate(land_coverage(nobs))
allocate(ice_coverage(nobs))
allocate(snow_coverage(nobs))
allocate(lai(nobs))
allocate(water_temperature(nobs))
allocate(land_temperature(nobs))
allocate(ice_temperature(nobs))
allocate(snow_temperature(nobs))
allocate(soil_moisture_content(nobs))
allocate(vegetation_fraction(nobs))
allocate(soil_temperature(nobs))
allocate(snow_depth(nobs))

wind_speed = 0.0_kind_real
wind_direction = 0.0_kind_real
land_type = 0
vegetation_type = 0
soil_type = 0
water_coverage = 0.0_kind_real
land_coverage = 0.0_kind_real
ice_coverage = 0.0_kind_real
snow_coverage = 0.0_kind_real
lai = 0.0_kind_real
water_temperature = 0.0_kind_real
land_temperature = 0.0_kind_real
ice_temperature = 0.0_kind_real
snow_temperature = 0.0_kind_real
soil_moisture_content = 0.0_kind_real
vegetation_fraction = 0.0_kind_real
soil_temperature = 0.0_kind_real
snow_depth = 0.0_kind_real

if (fld%havecrtmfields) then
  !TODO only if a radiance
  call crtm_surface( fld%geom, nobs, ngrid, locs%lat(:), locs%lon(:), &
                     fld%Atm%slmsk, fld%Atm%sheleg, fld%Atm%tsea, fld%Atm%vtype, &
                     fld%Atm%stype, fld%Atm%vfrac, fld%Atm%stc, fld%Atm%smc, fld%Atm%snwdph, &
                     fld%Atm%u_srf,fld%Atm%v_srf,fld%Atm%f10m, &
                     land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                     snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                     snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                     wind_speed, wind_direction )
endif


! Get CRTM moisture variables
! ---------------------------
allocate(ql_ade(isd:ied,jsd:jed,npz))
allocate(qi_ade(isd:ied,jsd:jed,npz))
allocate(ql_efr(isd:ied,jsd:jed,npz))
allocate(qi_efr(isd:ied,jsd:jed,npz))
allocate(qmr(isd:ied,jsd:jed,npz))
allocate(water_coverage_m(isd:ied,jsd:jed))

ql_ade = 0.0_kind_real
qi_ade = 0.0_kind_real
ql_efr = 0.0_kind_real
qi_efr = 0.0_kind_real

if (fld%havecrtmfields) then

  !TODO Is it water_coverage or sea_coverage fed in here?
  water_coverage_m = 0.0_kind_real
  do j = jsc,jec
    do i = isc,iec
      if (fld%Atm%slmsk(i,j) == 0) water_coverage_m(i,j) = 1.0_kind_real
    enddo
  enddo
  
  call crtm_ade_efr( fld%geom,prsi,fld%Atm%pt,fld%Atm%delp,water_coverage_m,fld%Atm%q(:,:,:,fld%ti_q), &
                     fld%Atm%q(:,:,:,fld%ti_ql),fld%Atm%q(:,:,:,fld%ti_qi),ql_ade,qi_ade,ql_efr,qi_efr )
  
  call crtm_mixratio(fld%geom,fld%Atm%q(:,:,:,fld%ti_q),qmr)

endif


!write(*,*)'interp model    t min, max= ',minval(fld%Atm%pt),maxval(fld%Atm%pt)
!write(*,*)'interp model delp min, max= ',minval(fld%Atm%delp),maxval(fld%Atm%delp)

! Variable transforms and interpolate to obs locations
! ----------------------------------------------------

do jvar = 1, vars%nv

  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  ! Convert to observation variables/units
  ! --------------------------------------
  select case (trim(vars%fldnames(jvar)))

  case ("temperature")

    nvl = npz
    do_interp = .true.
    geovalm = fld%Atm%pt
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    call T_to_Tv(fld%geom,fld%Atm%pt,fld%Atm%q,geovalm)
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

    nvl = npz
    do_interp = .true.
    geovalm = log(0.001_kind_real) + logp !to kPa
    geoval => geovalm

  case ("humidity_mixing_ratio")

    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm  

  case ("air_pressure")

    nvl = npz
    do_interp = .true.
    geovalm = prs / 100.0_kind_real !to hPa
    geoval => geovalm

  case ("air_pressure_levels")

    nvl = npz + 1
    do_interp = .true.
    geovale = prsi / 100.0_kind_real !to hPa
    geoval => geovale

  case ("mass_concentration_of_ozone_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = fld%Atm%q(:,:,:,fld%ti_o3) * constoz
   geoval => geovalm

  case ("mass_concentration_of_carbon_dioxide_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = 407.0_kind_real !Just a constant for now
   geoval => geovalm

  case ("atmosphere_mass_content_of_cloud_liquid_water")

   nvl = npz
   do_interp = .true.
   geovalm = ql_ade
   geoval => geovalm

  case ("atmosphere_mass_content_of_cloud_ice")

   nvl = npz
   do_interp = .true.
   geovalm = qi_ade
   geoval => geovalm

  case ("effective_radius_of_cloud_liquid_water_particle")

   nvl = npz
   do_interp = .true.
   geovalm = ql_efr
   geoval => geovalm

  case ("effective_radius_of_cloud_ice_particle")

   nvl = npz
   do_interp = .true.
   geovalm = qi_efr
   geoval => geovalm

  case ("Water_Fraction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = water_coverage

  case ("Land_Fraction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = land_coverage

  case ("Ice_Fraction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = ice_coverage

  case ("Snow_Fraction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = snow_coverage

  case ("Water_Temperature")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = water_temperature

  case ("Land_Temperature")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = land_temperature

  case ("Ice_Temperature")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = ice_temperature

  case ("Snow_Temperature")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = snow_temperature

  case ("Snow_Depth")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = snow_depth

  case ("Vegetation_Fraction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = vegetation_fraction

  case ("Sfc_Wind_Speed")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = wind_speed

  case ("Sfc_Wind_Direction")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = wind_direction

  case ("Lai")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = lai

  case ("Soil_Moisture")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = soil_moisture_content

  case ("Soil_Temperature")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = soil_temperature

  case ("Land_Type_Index")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = real(land_type,kind_real)

  case ("Vegetation_Type")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = real(vegetation_type,kind_real)

  case ("Soil_Type")

   nvl = 1
   do_interp = .false.
   obs_field(:,1) = real(soil_type,kind_real)

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select


  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nobs,nvl)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call interp_checks(myname, fld, locs, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = 1, nvl
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          mod_field(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call pbump%apply_obsop(mod_field,obs_field)
      gom%geovals(jvar)%vals(jlev,:) = obs_field(:,1)
    enddo
  else
    gom%geovals(jvar)%vals(nvl,:) = obs_field(:,1)
  endif

  nullify(geoval)

enddo

nullify(pbump)

deallocate(mod_field)
deallocate(obs_field)
deallocate(geovale)
deallocate(geovalm)
deallocate(prsi)
deallocate(prs )
deallocate(logp)
deallocate(wind_speed)
deallocate(wind_direction)
deallocate(land_type)
deallocate(vegetation_type)
deallocate(soil_type)
deallocate(water_coverage)
deallocate(land_coverage)
deallocate(ice_coverage)
deallocate(snow_coverage)
deallocate(lai)
deallocate(water_temperature)
deallocate(land_temperature)
deallocate(ice_temperature)
deallocate(snow_temperature)
deallocate(soil_moisture_content)
deallocate(vegetation_fraction)
deallocate(soil_temperature)
deallocate(snow_depth)
deallocate(ql_ade)
deallocate(qi_ade)
deallocate(ql_efr)
deallocate(qi_efr)
deallocate(qmr)
deallocate(water_coverage_m)

!write(*,*)'interp geovals t min, max= ',minval(gom%geovals(1)%vals(:,:)),maxval(gom%geovals(1)%vals(:,:))
!write(*,*)'interp geovals p min, max= ',minval(gom%geovals(2)%vals(:,:)),maxval(gom%geovals(2)%vals(:,:))

end subroutine getvalues

! ------------------------------------------------------------------------------

subroutine getvalues_tl(fld, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_tl

implicit none
type(fv3jedi_field),      intent(inout) :: fld 
type(ioda_locs),          intent(in)    :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_tl'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


! Grid convenience
! ----------------
isc = fld%geom%bd%isc
iec = fld%geom%bd%iec
jsc = fld%geom%bd%jsc
jec = fld%geom%bd%jec
isd = fld%geom%bd%isd
ied = fld%geom%bd%ied
jsd = fld%geom%bd%jsd
jed = fld%geom%bd%jed
npz = fld%geom%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(traj%ngrid,1))
allocate(obs_field(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isd:ied,jsd:jed,npz+1))
allocate(geovalm(isd:ied,jsd:jed,npz))


! Interpolate fields to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv
 
  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  select case (trim(vars%fldnames(jvar)))
   
  case ("temperature")
  
    nvl = npz
    do_interp = .true.
    geovalm = fld%Atm%pt
    geoval => geovalm

  case ("virtual_temperature")

    nvl = fld%geom%npz
    do_interp = .true.
    call T_to_Tv_tl(fld%geom, traj%pt, fld%Atm%pt, traj%q(:,:,:,1), fld%Atm%q (:,:,:,1) )
    geovalm = fld%Atm%pt
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = fld%geom%npz
    do_interp = .true.
    call crtm_mixratio_tl(fld%geom, traj%q(:,:,:,1), fld%Atm%q(:,:,:,1), geovalm)
    geoval => geovalm  

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,traj%nobs,nvl)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call interp_checks(myname, fld, locs, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = 1, nvl
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          mod_field(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call traj%bump(1)%apply_obsop(mod_field,obs_field)
      gom%geovals(jvar)%vals(jlev,:) = obs_field(:,1)
    enddo
  else
    gom%geovals(jvar)%vals(nvl,:) = obs_field(:,1)
  endif

  nullify(geoval)

enddo

deallocate(geovalm,geovale)

deallocate(mod_field)
deallocate(obs_field)

end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(fld, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_ad

implicit none
type(fv3jedi_field),      intent(inout) :: fld 
type(ioda_locs),           intent(in)   :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_ad'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_field(:,:)
real(kind=kind_real), allocatable :: obs_field(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


! Grid convenience
! ----------------
isc = fld%geom%bd%isc
iec = fld%geom%bd%iec
jsc = fld%geom%bd%jsc
jec = fld%geom%bd%jec
isd = fld%geom%bd%isd
ied = fld%geom%bd%ied
jsd = fld%geom%bd%jsd
jed = fld%geom%bd%jed
npz = fld%geom%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_field(traj%ngrid,1))
allocate(obs_field(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isd:ied,jsd:jed,npz+1))
allocate(geovalm(isd:ied,jsd:jed,npz))

geovale = 0.0_kind_real
geovalm = 0.0_kind_real

! Interpolate fields to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv

  ! PART 1, do_interp flag
  ! ----------------------
  do_interp = .false.

  select case (trim(vars%fldnames(jvar)))
   
  case ("temperature")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  !Part 2, apply adjoint of interp
  !-------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = nvl, 1, -1
      obs_field(:,1) = gom%geovals(jvar)%vals(jlev,:)
      gom%geovals(jvar)%vals(jlev,:) = 0.0_kind_real
      call traj%bump(1)%apply_obsop_ad(obs_field,mod_field)
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          geoval(ji, jj, jlev) = mod_field(ii, 1)
        enddo
      enddo
    enddo
  else
    obs_field(:,1) = gom%geovals(jvar)%vals(nvl,:)
  endif

  !Part 3, back to state variables
  !-------------------------------
 
  select case (trim(vars%fldnames(jvar)))
 
  case ("temperature")
  
    fld%Atm%pt = geovalm 

  case ("virtual_temperature")
    
    fld%Atm%pt = geovalm
    call T_to_Tv_ad(fld%geom, traj%pt, fld%Atm%pt, traj%q(:,:,:,1), fld%Atm%q (:,:,:,1) )

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    call crtm_mixratio_ad(fld%geom, traj%q(:,:,:,1), fld%Atm%q(:,:,:,1), geovalm)

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  geovale = 0.0_kind_real
  geovalm = 0.0_kind_real

enddo

deallocate(mod_field)
deallocate(obs_field)

end subroutine getvalues_ad

! ------------------------------------------------------------------------------

subroutine allocate_geovals_vals(gom,jvar,nobs,gvlev)

implicit none
integer, intent(in) :: jvar, nobs, gvlev
type(ufo_geovals), intent(inout) :: gom

! Allocate geovals for this jvar
if (allocated(gom%geovals(jvar)%vals)) deallocate(gom%geovals(jvar)%vals)

allocate(gom%geovals(jvar)%vals(gvlev,nobs))

gom%geovals(jvar)%nval = gvlev
gom%geovals(jvar)%nobs = nobs
gom%geovals(jvar)%vals = 0.0_kind_real

gom%linit  = .true.

end subroutine allocate_geovals_vals

! ------------------------------------------------------------------------------

subroutine initialize_bump(geom, locs, bump)

use fv3jedi_geom_mod, only: fv3jedi_geom
use type_bump, only: bump_type
use mpi, only: mpi_comm_world

implicit none

!Arguments
type(fv3jedi_geom), intent(in)    :: geom
type(ioda_locs),    intent(in)    :: locs
type(bump_type),    intent(inout) :: bump

!Locals
integer :: mod_num
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

integer, save :: bumpcount = 0
character(len=5) :: cbumpcount
character(len=16) :: bump_nam_prefix


! Each bump%nam%prefix must be distinct
! -------------------------------------
bumpcount = bumpcount + 1
write(cbumpcount,"(I0.5)") bumpcount
bump_nam_prefix = 'fv3jedi_bump_data_'//cbumpcount


!Get the Solution dimensions
!---------------------------
mod_num = (geom%bd%iec - geom%bd%isc + 1) * (geom%bd%jec - geom%bd%jsc + 1)


!Calculate interpolation weight using BUMP
!-----------------------------------------
allocate( mod_lat(mod_num), mod_lon(mod_num) )
mod_lat = reshape( geom%grid_lat(geom%bd%isc:geom%bd%iec,      &
                                 geom%bd%jsc:geom%bd%jec),     &
                                [mod_num] )  
mod_lon = reshape( geom%grid_lon(geom%bd%isc:geom%bd%iec,      &
                                 geom%bd%jsc:geom%bd%jec),     &
                                [mod_num] ) - 180.0_kind_real

!Important namelist options
call bump%nam%init

bump%nam%prefix = bump_nam_prefix   ! Prefix for files output
bump%nam%nobs = locs%nlocs          ! Number of observations
bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)
bump%nam%obsdis = 'random'          ! Observation distribution parameter ('random','local' or 'adjusted')
bump%nam%diag_interp = 'bilin'

bump%nam%local_diag = .false.
bump%nam%nldwv = 10

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
                                nobs=locs%nlocs,lonobs=locs%lon(:)-180.0_kind_real,latobs=locs%lat(:) )

!Release memory
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate( mod_lat, mod_lon )

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom, jvar)
implicit none
character(len=*), intent(in) :: cop
type(fv3jedi_field), intent(in) :: fld
type(ioda_locs), intent(in)     :: locs
type(ufo_vars), intent(in)      :: vars
type(ufo_geovals), intent(in)   :: gom
integer, intent(in)             :: jvar

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
   call abor1_ftn(trim(cinfo)//"geovals not allocated")
endif
if( size(gom%geovals) .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//"geovals wrong size")
endif
if (.not.gom%linit) then
   call abor1_ftn(trim(cinfo)//"geovals not initialized")
endif
if (allocated(gom%geovals(jvar)%vals)) then  
   if( gom%geovals(jvar)%nobs .ne. locs%nlocs )then
      call abor1_ftn(trim(cinfo)//"nobs wrong size")
   endif
   if( size(gom%geovals(jvar)%vals, 2) .ne. locs%nlocs )then
      call abor1_ftn(trim(cinfo)//"vals wrong size 2")
   endif       
else
  call abor1_ftn(trim(cinfo)//"vals not allocated")
endif 

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
