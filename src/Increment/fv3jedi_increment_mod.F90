! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Handle increment for the FV3JEDI odel

module fv3jedi_increment_mod

use iso_c_binding
use config_mod
use datetime_mod

use ioda_locs_mod
use ufo_vars_mod
use ufo_geovals_mod

use fv3jedi_constants, only: rad2deg, constoz
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_increment_io_mod
use fv3jedi_kinds, only: kind_real
use fv3jedi_state_utils_mod, only: fv3jedi_state
use fv3jedi_vars_mod, only: fv3jedi_vars
use fv3jedi_getvalues_mod, only: getvalues_tl, getvalues_ad

implicit none
private

public :: create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, incrms, &
          change_resol, getvalues_tl, getvalues_ad, &
          ug_coord, increment_to_ug, increment_from_ug, dirac, svnorm
public :: fv3jedi_increment
public :: fv3jedi_increment_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_increment

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_increment_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(fv3jedi_vars),  intent(in)    :: vars

integer :: var

! Copy the variable names
self%vars%nv = vars%nv
allocate(self%vars%fldnames(self%vars%nv))
self%vars%fldnames = vars%fldnames

! Allocate variables based on names
do var = 1, self%vars%nv

   select case (trim(self%vars%fldnames(var)))

     case("ud")
      !if (.not.allocated(  self%ud)) allocate (  self%ud(geom%isc:geom%iec,  geom%jsc:geom%jec+1, geom%npz))
     case("vd")
      !if (.not.allocated(  self%vd)) allocate (  self%vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  , geom%npz))
     case("ua")
       if (.not.allocated(  self%ua)) allocate (  self%ua(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("va")
       if (.not.allocated(  self%va)) allocate (  self%va(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("t")
       if (.not.allocated(   self%t)) allocate (   self%t(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("delp")
       if (.not.allocated(self%delp)) allocate (self%delp(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("q")
       if (.not.allocated(   self%q)) allocate (   self%q(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("qi")
       if (.not.allocated(  self%qi)) allocate (  self%qi(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("ql")
       if (.not.allocated(  self%ql)) allocate (  self%ql(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("o3")
       if (.not.allocated(  self%o3)) allocate (  self%o3(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("psi")
       if (.not.allocated( self%psi)) allocate ( self%psi(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("chi")
       if (.not.allocated( self%chi)) allocate ( self%chi(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("tv")
       if (.not.allocated(  self%tv)) allocate (  self%tv(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("ps")
       if (.not.allocated(  self%ps)) allocate (  self%ps(geom%isc:geom%iec,  geom%jsc:geom%jec            ))
     case("qc")
       if (.not.allocated(  self%qc)) allocate (  self%qc(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("qic")
       if (.not.allocated( self%qic)) allocate ( self%qic(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("qlc")
       if (.not.allocated( self%qlc)) allocate ( self%qlc(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("o3c")
       if (.not.allocated( self%o3c)) allocate ( self%o3c(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("w")
       if (.not.allocated(   self%w)) allocate (   self%w(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("delz")
       if (.not.allocated(self%delz)) allocate (self%delz(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case default 
       call abor1_ftn("Create: unknown variable "//trim(self%vars%fldnames(var)))

   end select

enddo

self%hydrostatic = .true.
if (allocated(self%w).and.allocated(self%delz)) self%hydrostatic = .false.

! Initialize all domain arrays to zero
call zeros(self)

! For convenience
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%isd = geom%isd
self%ied = geom%ied
self%jsd = geom%jsd
self%jed = geom%jed
self%npx = geom%npx
self%npy = geom%npy
self%npz = geom%npz

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

if (allocated(self%ua  )) deallocate (self%ua  )
if (allocated(self%va  )) deallocate (self%va  )
if (allocated(self%t   )) deallocate (self%t   )
if (allocated(self%delp)) deallocate (self%delp)
if (allocated(self%q   )) deallocate (self%q   )
if (allocated(self%qi  )) deallocate (self%qi  )
if (allocated(self%ql  )) deallocate (self%ql  )
if (allocated(self%o3  )) deallocate (self%o3  )

if (allocated(self%psi )) deallocate(self%psi )
if (allocated(self%chi )) deallocate(self%chi )
if (allocated(self%tv  )) deallocate(self%tv  )
if (allocated(self%ps  )) deallocate(self%ps  )
if (allocated(self%qc  )) deallocate(self%qc  )
if (allocated(self%qic )) deallocate(self%qic )
if (allocated(self%qlc )) deallocate(self%qlc )
if (allocated(self%o3c )) deallocate(self%o3c )

if (allocated(self%w   )) deallocate (self%w   )
if (allocated(self%delz)) deallocate (self%delz)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

!Zero out the entire domain

if(allocated(self%ua  )) self%ua   = 0.0_kind_real
if(allocated(self%va  )) self%va   = 0.0_kind_real
if(allocated(self%t   )) self%t    = 0.0_kind_real
if(allocated(self%delp)) self%delp = 0.0_kind_real
if(allocated(self%q   )) self%q    = 0.0_kind_real
if(allocated(self%qi  )) self%qi   = 0.0_kind_real
if(allocated(self%ql  )) self%ql   = 0.0_kind_real
if(allocated(self%o3  )) self%o3   = 0.0_kind_real

if(allocated(self%psi)) self%psi   = 0.0_kind_real
if(allocated(self%chi)) self%chi   = 0.0_kind_real
if(allocated(self%tv )) self%tv    = 0.0_kind_real
if(allocated(self%ps )) self%ps    = 0.0_kind_real
if(allocated(self%qc )) self%qc    = 0.0_kind_real
if(allocated(self%qic)) self%qic   = 0.0_kind_real
if(allocated(self%qlc)) self%qlc   = 0.0_kind_real
if(allocated(self%o3c)) self%o3c   = 0.0_kind_real

if(allocated(self%w   )) self%w    = 0.0_kind_real
if(allocated(self%delz)) self%delz = 0.0_kind_real

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine ones(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

call zeros(self)

if(allocated(self%ua  )) self%ua   = 1.0_kind_real
if(allocated(self%va  )) self%va   = 1.0_kind_real
if(allocated(self%t   )) self%t    = 1.0_kind_real
if(allocated(self%delp)) self%delp = 1.0_kind_real
if(allocated(self%q   )) self%q    = 1.0_kind_real
if(allocated(self%qi  )) self%qi   = 1.0_kind_real
if(allocated(self%ql  )) self%ql   = 1.0_kind_real
if(allocated(self%o3  )) self%o3   = 1.0_kind_real

if(allocated(self%psi)) self%psi   = 1.0_kind_real
if(allocated(self%chi)) self%chi   = 1.0_kind_real
if(allocated(self%tv )) self%tv    = 1.0_kind_real
if(allocated(self%ps )) self%ps    = 1.0_kind_real
if(allocated(self%qc )) self%qc    = 1.0_kind_real
if(allocated(self%qic)) self%qic   = 1.0_kind_real
if(allocated(self%qlc)) self%qlc   = 1.0_kind_real
if(allocated(self%o3c)) self%o3c   = 1.0_kind_real

if(allocated(self%w   )) self%w    = 1.0_kind_real
if(allocated(self%delz)) self%delz = 1.0_kind_real

end subroutine ones

! ------------------------------------------------------------------------------

subroutine random(self)
use random_vectors_mod
implicit none
type(fv3jedi_increment), intent(inout) :: self
integer :: nq

if(allocated(self%ua  )) call random_vector(self%ua  )
if(allocated(self%va  )) call random_vector(self%va  )
if(allocated(self%t   )) call random_vector(self%t   )
if(allocated(self%delp)) call random_vector(self%delp)
if(allocated(self%q   )) call random_vector(self%q   )
if(allocated(self%qi  )) call random_vector(self%qi  )
if(allocated(self%ql  )) call random_vector(self%ql  )
if(allocated(self%o3  )) call random_vector(self%o3  )

if(allocated(self%psi )) call random_vector(self%psi )
if(allocated(self%chi )) call random_vector(self%chi )
if(allocated(self%tv  )) call random_vector(self%tv  )
if(allocated(self%ps  )) call random_vector(self%ps  )
if(allocated(self%qc  )) call random_vector(self%qc  )
if(allocated(self%qic )) call random_vector(self%qic )
if(allocated(self%qlc )) call random_vector(self%qlc )
if(allocated(self%o3c )) call random_vector(self%o3c )

if(allocated(self%w   )) call random_vector(self%w   )
if(allocated(self%delz)) call random_vector(self%delz)

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

self%isc            = rhs%isc           
self%iec            = rhs%iec           
self%jsc            = rhs%jsc           
self%jec            = rhs%jec           
self%isd            = rhs%isd           
self%ied            = rhs%ied           
self%jsd            = rhs%jsd           
self%jed            = rhs%jed           
self%npx            = rhs%npx           
self%npy            = rhs%npy           
self%npz            = rhs%npz           
self%hydrostatic    = rhs%hydrostatic   
self%calendar_type  = rhs%calendar_type 
self%date           = rhs%date          
self%date_init      = rhs%date_init     

if(allocated(self%ua  )) self%ua   = rhs%ua
if(allocated(self%va  )) self%va   = rhs%va
if(allocated(self%t   )) self%t    = rhs%t
if(allocated(self%delp)) self%delp = rhs%delp
if(allocated(self%q   )) self%q    = rhs%q
if(allocated(self%qi  )) self%qi   = rhs%qi
if(allocated(self%ql  )) self%ql   = rhs%ql
if(allocated(self%o3  )) self%o3   = rhs%o3

if(allocated(self%psi )) self%psi  = rhs%psi
if(allocated(self%chi )) self%chi  = rhs%chi
if(allocated(self%tv  )) self%tv   = rhs%tv
if(allocated(self%ps  )) self%ps   = rhs%ps
if(allocated(self%qc  )) self%qc   = rhs%qc
if(allocated(self%qic )) self%qic  = rhs%qic
if(allocated(self%qlc )) self%qlc  = rhs%qlc
if(allocated(self%o3c )) self%o3c  = rhs%o3c

if(allocated(self%w   )) self%w    = rhs%w
if(allocated(self%delz)) self%delz = rhs%delz

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ua  )) self%ua   = self%ua   + rhs%ua  
if(allocated(self%va  )) self%va   = self%va   + rhs%va  
if(allocated(self%t   )) self%t    = self%t    + rhs%t   
if(allocated(self%delp)) self%delp = self%delp + rhs%delp
if(allocated(self%q   )) self%q    = self%q    + rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   + rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   + rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   + rhs%o3  

if(allocated(self%psi )) self%psi  = self%psi  + rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  + rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   + rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   + rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   + rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  + rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  + rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  + rhs%o3c

if(allocated(self%w   )) self%w    = self%w    + rhs%w   
if(allocated(self%delz)) self%delz = self%delz + rhs%delz

return
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ua  )) self%ua   = self%ua   * rhs%ua  
if(allocated(self%va  )) self%va   = self%va   * rhs%va  
if(allocated(self%t   )) self%t    = self%t    * rhs%t   
if(allocated(self%delp)) self%delp = self%delp * rhs%delp
if(allocated(self%q   )) self%q    = self%q    * rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   * rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   * rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   * rhs%o3  

if(allocated(self%psi )) self%psi  = self%psi  * rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  * rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   * rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   * rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   * rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  * rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  * rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  * rhs%o3c

if(allocated(self%w   )) self%w    = self%w    * rhs%w   
if(allocated(self%delz)) self%delz = self%delz * rhs%delz

return
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ua  )) self%ua   = self%ua   - rhs%ua  
if(allocated(self%va  )) self%va   = self%va   - rhs%va  
if(allocated(self%t   )) self%t    = self%t    - rhs%t   
if(allocated(self%delp)) self%delp = self%delp - rhs%delp
if(allocated(self%q   )) self%q    = self%q    - rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   - rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   - rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   - rhs%o3  

if(allocated(self%psi )) self%psi  = self%psi  - rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  - rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   - rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   - rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   - rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  - rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  - rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  - rhs%o3c

if(allocated(self%w   )) self%w    = self%w    - rhs%w   
if(allocated(self%delz)) self%delz = self%delz - rhs%delz

return
end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz

if(allocated(self%ua  )) self%ua   = zz * self%ua  
if(allocated(self%va  )) self%va   = zz * self%va  
if(allocated(self%t   )) self%t    = zz * self%t   
if(allocated(self%delp)) self%delp = zz * self%delp
if(allocated(self%q   )) self%q    = zz * self%q   
if(allocated(self%qi  )) self%qi   = zz * self%qi  
if(allocated(self%ql  )) self%ql   = zz * self%ql  
if(allocated(self%o3  )) self%o3   = zz * self%o3  

if(allocated(self%psi )) self%psi  = zz * self%psi
if(allocated(self%chi )) self%chi  = zz * self%chi
if(allocated(self%tv  )) self%tv   = zz * self%tv 
if(allocated(self%ps  )) self%ps   = zz * self%ps 
if(allocated(self%qc  )) self%qc   = zz * self%qc 
if(allocated(self%qic )) self%qic  = zz * self%qic
if(allocated(self%qlc )) self%qlc  = zz * self%qlc
if(allocated(self%o3c )) self%o3c  = zz * self%o3c

if(allocated(self%w   )) self%w    = zz * self%w   
if(allocated(self%delz)) self%delz = zz * self%delz

return
end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy_inc(self,zz,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ua  )) self%ua   = self%ua   + zz * rhs%ua  
if(allocated(self%va  )) self%va   = self%va   + zz * rhs%va  
if(allocated(self%t   )) self%t    = self%t    + zz * rhs%t   
if(allocated(self%delp)) self%delp = self%delp + zz * rhs%delp
if(allocated(self%q   )) self%q    = self%q    + zz * rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   + zz * rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   + zz * rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   + zz * rhs%o3  

if(allocated(self%psi )) self%psi  = self%psi  + zz * rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  + zz * rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   + zz * rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   + zz * rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   + zz * rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  + zz * rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  + zz * rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  + zz * rhs%o3c

if(allocated(self%w   )) self%w    = self%w    + zz * rhs%w   
if(allocated(self%delz)) self%delz = self%delz + zz * rhs%delz

return
end subroutine axpy_inc

! ------------------------------------------------------------------------------

subroutine axpy_state(self,zz,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_state), intent(in)    :: rhs

if(allocated(self%ua  )) self%ua   = self%ua   + zz * rhs%ua  
if(allocated(self%va  )) self%va   = self%va   + zz * rhs%va  
if(allocated(self%t   )) self%t    = self%t    + zz * rhs%t   
if(allocated(self%delp)) self%delp = self%delp + zz * rhs%delp
if(allocated(self%q   )) self%q    = self%q    + zz * rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   + zz * rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   + zz * rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   + zz * rhs%o3  

if(allocated(self%w   )) self%w    = self%w    + zz * rhs%w   
if(allocated(self%delz)) self%delz = self%delz + zz * rhs%delz

return
end subroutine axpy_state

! ------------------------------------------------------------------------------

subroutine dot_prod(inc1,inc2,zprod)

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

implicit none
type(fv3jedi_increment), intent(in) :: inc1, inc2
real(kind=kind_real), intent(inout) :: zprod
real(kind=kind_real) :: zp
integer :: i,j,k
integer :: ierr
type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

zp=0.0_kind_real

!ua
if (allocated(inc1%ua)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%ua(i,j,k) * inc2%ua(i,j,k)
      enddo
    enddo
  enddo
endif

!va
if (allocated(inc1%va)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%va(i,j,k) * inc2%va(i,j,k)
      enddo
    enddo
  enddo
endif

!t
if (allocated(inc1%t)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%t(i,j,k) * inc2%t(i,j,k)
      enddo
    enddo
  enddo
endif

!delp
if (allocated(inc1%delp)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%delp(i,j,k) * inc2%delp(i,j,k)
      enddo
    enddo
  enddo
endif

!q
if (allocated(inc1%q)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%q(i,j,k) * inc2%q(i,j,k)
      enddo
    enddo
  enddo
endif

!qi
if (allocated(inc1%qi)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%qi(i,j,k) * inc2%qi(i,j,k)
      enddo
    enddo
  enddo
endif

!ql
if (allocated(inc1%ql)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%ql(i,j,k) * inc2%ql(i,j,k)
      enddo
    enddo
  enddo
endif

!o3
if (allocated(inc1%o3)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%o3(i,j,k) * inc2%o3(i,j,k)
      enddo
    enddo
  enddo
endif

!psi
if (allocated(inc1%psi)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%psi(i,j,k) * inc2%psi(i,j,k)
      enddo
    enddo
  enddo
endif

!chi
if (allocated(inc1%chi)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%chi(i,j,k) * inc2%chi(i,j,k)
      enddo
    enddo
  enddo
endif

!tv
if (allocated(inc1%tv)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%tv(i,j,k) * inc2%tv(i,j,k)
      enddo
    enddo
  enddo
endif

!ps
if (allocated(inc1%ps)) then
  do j = inc1%jsc,inc1%jec
    do i = inc1%isc,inc1%iec
      zp = zp + inc1%ps(i,j) * inc2%ps(i,j)
    enddo
  enddo
endif

!qc
if (allocated(inc1%qc)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%qc(i,j,k) * inc2%qc(i,j,k)
      enddo
    enddo
  enddo
endif

!qic
if (allocated(inc1%qic)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%qic(i,j,k) * inc2%qic(i,j,k)
      enddo
    enddo
  enddo
endif

!qlc
if (allocated(inc1%qlc)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%qlc(i,j,k) * inc2%qlc(i,j,k)
      enddo
    enddo
  enddo
endif

!o3c
if (allocated(inc1%o3c)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%o3c(i,j,k) * inc2%o3c(i,j,k)
      enddo
    enddo
  enddo
endif

!delz
if (allocated(inc1%delz)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%delz(i,j,k) * inc2%delz(i,j,k)
      enddo
    enddo
  enddo
endif

!w
if (allocated(inc1%w)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%w(i,j,k) * inc2%w(i,j,k)
      enddo
    enddo
  enddo
endif

!Get global dot product
call f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
if (f_comm%rank() == 0) print*, "Dot product test result: ", zprod

return
end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: check
check = (rhs%iec-rhs%isc+1) - (self%iec-self%isc+1)

if (check==0) then
  if(allocated(rhs%ua  )) self%ua   = self%ua   + rhs%ua  
  if(allocated(rhs%va  )) self%va   = self%va   + rhs%va  
  if(allocated(rhs%t   )) self%t    = self%t    + rhs%t   
  if(allocated(rhs%delp)) self%delp = self%delp + rhs%delp
  if(allocated(rhs%q   )) self%q    = self%q    + rhs%q   
  if(allocated(rhs%qi  )) self%qi   = self%qi   + rhs%qi  
  if(allocated(rhs%ql  )) self%ql   = self%ql   + rhs%ql  
  if(allocated(rhs%o3  )) self%o3   = self%o3   + rhs%o3  

  if(allocated(rhs%psi )) self%psi  = self%psi  + rhs%psi
  if(allocated(rhs%chi )) self%chi  = self%chi  + rhs%chi
  if(allocated(rhs%tv  )) self%tv   = self%tv   + rhs%tv 
  if(allocated(rhs%ps  )) self%ps   = self%ps   + rhs%ps 
  if(allocated(rhs%qc  )) self%qc   = self%qc   + rhs%qc 
  if(allocated(rhs%qic )) self%qic  = self%qic  + rhs%qic
  if(allocated(rhs%qlc )) self%qlc  = self%qlc  + rhs%qlc
  if(allocated(rhs%o3c )) self%o3c  = self%o3c  + rhs%o3c

  if(allocated(rhs%w   )) self%w    = self%w    + rhs%w   
  if(allocated(rhs%delz)) self%delz = self%delz + rhs%delz 
else
   call abor1_ftn("fv3jedi increment:  add_incr not implemented for low res increment yet")
endif

return
end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(fv3jedi_increment), intent(inout) :: lhs
type(fv3jedi_state), intent(in)    :: x1
type(fv3jedi_state), intent(in)    :: x2

integer :: check
check = (x1%iec-x1%isc+1) - (x2%iec-x2%isc+1)

call zeros(lhs)
if (check==0) then

  if(allocated(lhs%ua  )) lhs%ua   = x1%ua   - x2%ua  
  if(allocated(lhs%va  )) lhs%va   = x1%va   - x2%va  
  if(allocated(lhs%t   )) lhs%t    = x1%t    - x2%t   
  if(allocated(lhs%delp)) lhs%delp = x1%delp - x2%delp
  if(allocated(lhs%q   )) lhs%q    = x1%q    - x2%q   
  if(allocated(lhs%qi  )) lhs%qi   = x1%qi   - x2%qi  
  if(allocated(lhs%ql  )) lhs%ql   = x1%ql   - x2%ql  
  if(allocated(lhs%o3  )) lhs%o3   = x1%o3   - x2%o3  

  if(allocated(lhs%w   )) lhs%w    = x1%w    - x2%w   
  if(allocated(lhs%delz)) lhs%delz = x1%delz - x2%delz

else

   call abor1_ftn("fv3jedi increment:  diff_incr not implemented for low res increment yet")

endif

return
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(inc,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_increment), intent(in)    :: rhs

integer :: check
check = (rhs%iec-rhs%isc+1) - (inc%iec-inc%isc+1)

if (check==0) then
   call copy(inc, rhs)
else
   call abor1_ftn("fv3jedi_increment: change_resol not implmeneted yet")
endif

return
end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(geom, inc, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_increment), intent(inout) :: inc      !< Increment
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call read_fms_restart(geom, inc, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call read_geos_restart(geom, inc, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_increment read: restart type not supported")
  endif

  return

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, inc, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_increment), intent(in)    :: inc      !< Increment
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call write_fms_restart(geom, inc, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call write_geos_restart(geom, inc, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_increment write: restart type not supported")
  endif

  return

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(inc, nf, pstat)
implicit none
type(fv3jedi_increment), intent(in) :: inc
integer, intent(in) :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)

integer :: isc, iec, jsc, jec, gs

!1. Min
!2. Max
!3. RMS

isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec

gs = (iec-isc+1)*(jec-jsc+1)*inc%npz

!ua
if (allocated(inc%ua)) then
  pstat(1,1) = minval(inc%ua(isc:iec,jsc:jec,:))
  pstat(2,1) = maxval(inc%ua(isc:iec,jsc:jec,:))
  pstat(3,1) = sqrt((sum(inc%ua(isc:iec,jsc:jec,:))/gs)**2)
endif

!va
if (allocated(inc%va)) then
  pstat(1,2) = minval(inc%va(isc:iec,jsc:jec,:))
  pstat(2,2) = maxval(inc%va(isc:iec,jsc:jec,:))
  pstat(3,2) = sqrt((sum(inc%va(isc:iec,jsc:jec,:))/gs)**2)
endif
  
!t
if (allocated(inc%t)) then
  pstat(1,3) = minval(inc%t(isc:iec,jsc:jec,:))
  pstat(2,3) = maxval(inc%t(isc:iec,jsc:jec,:))
  pstat(3,3) = sqrt((sum(inc%t(isc:iec,jsc:jec,:))/gs)**2)
endif
  
!delp
if (allocated(inc%delp)) then
  pstat(1,4) = minval(inc%delp(isc:iec,jsc:jec,:))
  pstat(2,4) = maxval(inc%delp(isc:iec,jsc:jec,:))
  pstat(3,4) = sqrt((sum(inc%delp(isc:iec,jsc:jec,:))/gs)**2)
endif
  
!q
if (allocated(inc%q)) then
  pstat(1,5) = minval(inc%q(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(inc%q(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(inc%q(isc:iec,jsc:jec,:))/gs)**2)
endif
  
!qi
if (allocated(inc%qi)) then
  pstat(1,5) = minval(inc%qi(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(inc%qi(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(inc%qi(isc:iec,jsc:jec,:))/gs)**2)
endif

!ql
if (allocated(inc%ql)) then
  pstat(1,5) = minval(inc%ql(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(inc%ql(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(inc%ql(isc:iec,jsc:jec,:))/gs)**2)
endif
  
!o3
if (allocated(inc%o3)) then
  pstat(1,5) = minval(inc%o3(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(inc%o3(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(inc%o3(isc:iec,jsc:jec,:))/gs)**2)
endif

return

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine incrms(inc, prms)
use fckit_mpi_module, only : fckit_mpi_comm, fckit_mpi_sum
implicit none
type(fv3jedi_increment), intent(in) :: inc
real(kind=kind_real), intent(out) :: prms

real(kind=kind_real) :: zz
integer i,j,k,ii,nt,ierr,npes,iisum
integer :: isc,iec,jsc,jec,npz
type(fckit_mpi_comm) :: f_comm

isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
npz = inc%npz

f_comm = fckit_mpi_comm()

zz = 0.0_kind_real
prms = 0.0_kind_real
ii = 0

!ua
if (allocated(inc%ua)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%ua(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!va
if (allocated(inc%va)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%va(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!t
if (allocated(inc%t)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%t(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!delp
if (allocated(inc%delp)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%delp(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!q
if (allocated(inc%q)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%q(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qi
if (allocated(inc%qi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%qi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ql
if (allocated(inc%ql)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%ql(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!o3
if (allocated(inc%o3)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%o3(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!psi
if (allocated(inc%psi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%psi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!chi
if (allocated(inc%chi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%chi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!tv
if (allocated(inc%tv)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%tv(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ps
if (allocated(inc%ps)) then
  do j = jsc,jec
    do i = isc,iec
      zz = zz + inc%ps(i,j)**2
      ii = ii + 1
    enddo
  enddo
endif

!qc
if (allocated(inc%qc)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%qc(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qic
if (allocated(inc%qic)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%qic(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qlc
if (allocated(inc%qlc)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%qlc(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!o3c
if (allocated(inc%o3c)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%o3c(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!w
if (allocated(inc%w)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%w(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!delz
if (allocated(inc%delz)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%delz(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!Get global values
call f_comm%allreduce(zz,prms,fckit_mpi_sum())
call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

!if (ierr .ne. 0) then
!   print *,'error in incrms/mpi_allreduce, error code=',ierr
!endif
prms = sqrt(prms/real(iisum,kind_real))

end subroutine incrms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf, geom)

use iso_c_binding

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(c_ptr), intent(in)       :: c_conf   !< Configuration
type(fv3jedi_geom), intent(in) :: geom
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
if (ndir<1) call abor1_ftn("fv3jedi_increment:dirac non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%npx)) then
   call abor1_ftn("fv3jedi_increment:dirac invalid ixdir")
endif
if (any(iydir<1).or.any(iydir>geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_increment:dirac invalid iydir")
endif
if ((ildir<1).or.(ildir>self%npz)) then
   call abor1_ftn("fv3jedi_increment:dirac invalid ildir")
endif
if ((ifdir<1).or.(ifdir>5)) then
   call abor1_ftn("fv3jedi_increment:dirac invalid ifdir")
endif
if ((itiledir<1).or.(itiledir>6)) then
   call abor1_ftn("fv3jedi_increment:dirac invalid itiledir")
endif

! Setup Diracs
call zeros(self)

! only u,v,theta,delp and humidity allowed
do idir=1,ndir

   ! is specified grid point, tile number on this processor
   if (geom%ntile == itiledir .and. &
       ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
       iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then
       ! If so, perturb desired increment and level
       if (ifdir == 1) then
          self%ua  (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 2) then
          self%va  (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 3) then
          self%t   (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 4) then
          self%delp(ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 5) then
          self%q   (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 6) then
          self%qi  (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 7) then
          self%ql  (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 8) then
          self%o3  (ixdir(idir),iydir(idir),ildir) = 1.0
       endif
   endif
end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine ug_size(self, ug)
use unstructured_grid_mod
implicit none
type(fv3jedi_increment), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug

integer :: igrid

! Set number of grids
if (ug%colocated==1) then
   ! Colocatd
   ug%ngrid = 1
else
   ! Not colocatedd
   ug%ngrid = 1
end if

! Allocate grid instances
if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

if (ug%colocated==1) then
  ! colocatedd

  ! Set local number of points
  ug%grid(1)%nmga = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1) 

  ! Set number of levels
  ug%grid(1)%nl0 = self%npz

  ! Set number of variables
  ug%grid(1)%nv = 8

  ! Set number of timeslots
  ug%grid(1)%nts = 1
else
  ! Not colocatedd
  do igrid=1,ug%ngrid
     ! Set local number of points
     ug%grid(igrid)%nmga = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1)

     ! Set number of levels
     ug%grid(igrid)%nl0 = self%npz

     ! Set number of variables
     ug%grid(igrid)%nv = 8

     ! Set number of timeslots
     ug%grid(igrid)%nts = 1
  enddo
end if

end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine ug_coord(self, ug, colocated, geom)
use unstructured_grid_mod
implicit none
type(fv3jedi_increment), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
integer, intent(in) :: colocated
type(fv3jedi_geom), intent(in) :: geom

integer :: imga,jx,jy,jl
real(kind=kind_real),allocatable :: lon(:),lat(:),area(:),vunit(:,:)
real(kind=kind_real) :: sigmaup,sigmadn
integer :: igrid

! Copy colocated
ug%colocated = colocated

! Define size
call ug_size(self, ug)

! Allocate unstructured grid coordinates
call allocate_unstructured_grid_coord(ug)

if (ug%colocated==1) then
  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      ug%grid(1)%lon(imga) = rad2deg*geom%grid_lon(jx,jy)
      ug%grid(1)%lat(imga) = rad2deg*geom%grid_lat(jx,jy)
      ug%grid(1)%area(imga) = geom%area(jx,jy)
      do jl=1,self%npz
        sigmaup = geom%ak(jl+1)/101300.0+geom%bk(jl+1) ! si are now sigmas
        sigmadn = geom%ak(jl  )/101300.0+geom%bk(jl  )
        ug%grid(1)%vunit(imga,jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
        ug%grid(1)%lmask(imga,jl) = .true.
      enddo
    enddo
  enddo 
else
  !Not colocated
  do igrid=1,ug%ngrid
    imga = 0
    do jy=self%jsc,self%jec
      do jx=self%isc,self%iec
        imga = imga+1
        ug%grid(igrid)%lon(imga) = rad2deg*geom%grid_lon(jx,jy)
        ug%grid(igrid)%lat(imga) = rad2deg*geom%grid_lat(jx,jy)
        ug%grid(igrid)%area(imga) = geom%area(jx,jy)
        do jl=1,self%npz
          sigmaup = geom%ak(jl+1)/101300.0+geom%bk(jl+1) ! si are now sigmas
          sigmadn = geom%ak(jl  )/101300.0+geom%bk(jl  )
          ug%grid(igrid)%vunit(imga,jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
          ug%grid(igrid)%lmask(imga,jl) = .true.
        enddo
      enddo
    enddo
  enddo
endif

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine increment_to_ug(self, ug, colocated)
use unstructured_grid_mod
implicit none
type(fv3jedi_increment), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
integer, intent(in) :: colocated

integer :: imga,jx,jy,jl
real(kind=kind_real),allocatable :: ptmp(:,:,:)
integer :: igrid

! Copy colocated
ug%colocated = colocated

! Define size
call ug_size(self, ug)

! Allocate unstructured grid increment
call allocate_unstructured_grid_field(ug)

! Copy increment

if (ug%colocated==1) then

if (allocated(self%ua)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(1)%fld(imga,jl,1,1) = self%ua  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,2,1) = self%va  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,3,1) = self%t   (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,4,1) = self%delp(jx,jy,jl)
          ug%grid(1)%fld(imga,jl,5,1) = self%q   (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,6,1) = self%qi  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,7,1) = self%ql  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,8,1) = self%o3  (jx,jy,jl)
      enddo
    enddo
  enddo

elseif (allocated(self%psi)) then

  allocate(ptmp(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  ptmp = 0.0
  ptmp(:,:,self%npz) = self%ps(self%isc:self%iec,self%jsc:self%jec)

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(1)%fld(imga,jl,1,1) = self%psi (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,2,1) = self%chi (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,3,1) = self%tv  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,4,1) = ptmp     (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,5,1) = self%qc  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,6,1) = self%qic (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,7,1) = self%qlc (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,8,1) = self%o3c (jx,jy,jl)
      enddo
    enddo
  enddo

  deallocate(ptmp)

endif

else

do igrid=1,ug%ngrid

if (allocated(self%ua)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(igrid)%fld(imga,jl,1,1) = self%ua  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,2,1) = self%va  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,3,1) = self%t   (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,4,1) = self%delp(jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,5,1) = self%q   (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,6,1) = self%qi  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,7,1) = self%ql  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,8,1) = self%o3  (jx,jy,jl)
      enddo
    enddo
  enddo

elseif (allocated(self%psi)) then

  allocate(ptmp(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  ptmp = 0.0
  ptmp(:,:,self%npz) = self%ps(self%isc:self%iec,self%jsc:self%jec)

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(igrid)%fld(imga,jl,1,1) = self%psi (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,2,1) = self%chi (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,3,1) = self%tv  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,4,1) = ptmp     (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,5,1) = self%qc  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,6,1) = self%qic (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,7,1) = self%qlc (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,8,1) = self%o3c (jx,jy,jl)
      enddo
    enddo
  enddo

  deallocate(ptmp)

endif

enddo

endif

end subroutine increment_to_ug

! -----------------------------------------------------------------------------

subroutine increment_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug

integer :: imga,jx,jy,jl
real(kind=kind_real),allocatable :: ptmp(:,:,:)
integer :: igrid

! Copy increment

if (ug%colocated==1) then

if (allocated(self%ua)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%ua  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,1,1)
          self%va  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,2,1)
          self%t   (jx,jy,jl) = ug%grid(1)%fld(imga,jl,3,1)
          self%delp(jx,jy,jl) = ug%grid(1)%fld(imga,jl,4,1)
          self%q   (jx,jy,jl) = ug%grid(1)%fld(imga,jl,5,1)
          self%qi  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,6,1)
          self%ql  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,7,1)
          self%o3  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,8,1)
      enddo
    enddo
  enddo

elseif (allocated(self%psi)) then

  allocate(ptmp(self%isc:self%iec,self%jsc:self%jec,1:self%npz))

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%psi (jx,jy,jl) = ug%grid(1)%fld(imga,jl,1,1)
          self%chi (jx,jy,jl) = ug%grid(1)%fld(imga,jl,2,1)
          self%tv  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,3,1)
          ptmp     (jx,jy,jl) = ug%grid(1)%fld(imga,jl,4,1)
          self%qc  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,5,1)
          self%qic (jx,jy,jl) = ug%grid(1)%fld(imga,jl,6,1)
          self%qlc (jx,jy,jl) = ug%grid(1)%fld(imga,jl,7,1)
          self%o3c (jx,jy,jl) = ug%grid(1)%fld(imga,jl,8,1)
      enddo
    enddo
  enddo

  self%ps(self%isc:self%iec,self%jsc:self%jec) = ptmp(:,:,self%npz)

  deallocate(ptmp)

endif

else

do igrid=1,ug%ngrid

if (allocated(self%ua)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%ua  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,1,1)
          self%va  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,2,1)
          self%t   (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,3,1)
          self%delp(jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,4,1)
          self%q   (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,5,1)
          self%qi  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,6,1)
          self%ql  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,7,1)
          self%o3  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,8,1)
      enddo
    enddo
  enddo

elseif (allocated(self%psi)) then

  allocate(ptmp(self%isc:self%iec,self%jsc:self%jec,1:self%npz))

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%psi (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,1,1)
          self%chi (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,2,1)
          self%tv  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,3,1)
          ptmp     (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,4,1)
          self%qc  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,5,1)
          self%qic (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,6,1)
          self%qlc (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,7,1)
          self%o3c (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,8,1)
      enddo
    enddo
  enddo

  self%ps(self%isc:self%iec,self%jsc:self%jec) = ptmp(:,:,self%npz)

  deallocate(ptmp)

endif

enddo

endif


end subroutine increment_from_ug

! ------------------------------------------------------------------------------

subroutine check_resolution(x1, x2)

implicit none
type(fv3jedi_increment), intent(in) :: x1, x2

if (x1%npx /= x2%npx .or.  x1%npz /= x2%npz) then
  call abor1_ftn ("fv3jedi_increment: resolution error")
endif
call check(x1)
call check(x2)

end subroutine check_resolution

! ------------------------------------------------------------------------------

subroutine check(self)
implicit none
type(fv3jedi_increment), intent(in) :: self
logical :: bad

bad = .false.

if (bad) then
   call abor1_ftn ("fv3jedi_increment: increment not consistent")
endif

end subroutine check

! ------------------------------------------------------------------------------

subroutine svnorm(self,geom,ref,c_conf)

implicit none
type(fv3jedi_increment) :: self
type(fv3jedi_geom)      :: geom
type(fv3jedi_state)     :: ref !To linearize around if nl
type(c_ptr)             :: c_conf

integer :: i,j,k
integer :: isc,iec,jsc,jec,npz
real(kind=kind_real) :: psref

!Code to compute a vector norm for an increment, e.g. the energy norm for FSOI

isc = self%isc
iec = self%iec
jsc = self%isc
jec = self%iec
npz = self%npz

allocate(ps_ref(isc:iec,jsc:jec))

psref = sum(ref%delp,3)



deallocate(psref)

end subroutine svnorm

! ------------------------------------------------------------------------------

end module fv3jedi_increment_mod
