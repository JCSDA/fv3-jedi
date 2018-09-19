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

implicit none
private

public :: create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, incrms, &
          change_resol, getvalues_tl, getvalues_ad, &
          ug_coord, increment_to_ug, increment_from_ug, dirac
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

integer :: isd,ied,jsd,jed,npz,hydroi
integer :: var

! Grid convenience
isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed
npz = geom%npz

! Copy the variable names
self%vars%nv = vars%nv
allocate(self%vars%fldnames(self%vars%nv))
self%vars%fldnames = vars%fldnames

! Allocate variables based on names
do var = 1, self%vars%nv

   select case (trim(self%vars%fldnames(var)))

     case("ud")
       if (.not.allocated(  self%ud)) allocate (  self%ud(isd:ied,  jsd:jed+1, npz))
     case("vd")
       if (.not.allocated(  self%vd)) allocate (  self%vd(isd:ied+1,jsd:jed  , npz))
     case("ua")
       if (.not.allocated(  self%ua)) allocate (  self%ua(isd:ied,  jsd:jed  , npz))
     case("va")
       if (.not.allocated(  self%va)) allocate (  self%va(isd:ied,  jsd:jed  , npz))
     case("t")
       if (.not.allocated(   self%t)) allocate (   self%t(isd:ied,  jsd:jed  , npz))
     case("delp")
       if (.not.allocated(self%delp)) allocate (self%delp(isd:ied,  jsd:jed  , npz))
     case("q")
       if (.not.allocated(   self%q)) allocate (   self%q(isd:ied,  jsd:jed  , npz))
     case("qi")
       if (.not.allocated(  self%qi)) allocate (  self%qi(isd:ied,  jsd:jed  , npz))
     case("ql")
       if (.not.allocated(  self%ql)) allocate (  self%ql(isd:ied,  jsd:jed  , npz))
     case("o3")
       if (.not.allocated(  self%o3)) allocate (  self%o3(isd:ied,  jsd:jed  , npz))
     case("psi")
       if (.not.allocated( self%psi)) allocate ( self%psi(isd:ied,  jsd:jed  , npz))
     case("chi")
       if (.not.allocated( self%chi)) allocate ( self%chi(isd:ied,  jsd:jed  , npz))
     case("tv")
       if (.not.allocated(  self%tv)) allocate (  self%tv(isd:ied,  jsd:jed  , npz))
     case("ps")
       if (.not.allocated(  self%ps)) allocate (  self%ps(isd:ied,  jsd:jed       ))
     case("qc")
       if (.not.allocated(  self%qc)) allocate (  self%qc(isd:ied,  jsd:jed  , npz))
     case("qic")
       if (.not.allocated( self%qic)) allocate ( self%qic(isd:ied,  jsd:jed  , npz))
     case("qlc")
       if (.not.allocated( self%qlc)) allocate ( self%qlc(isd:ied,  jsd:jed  , npz))
     case("o3c")
       if (.not.allocated( self%o3c)) allocate ( self%o3c(isd:ied,  jsd:jed  , npz))
     case("w")
       if (.not.allocated(   self%w)) allocate (   self%w(isd:ied,  jsd:jed  , npz))
     case("delz")
       if (.not.allocated(self%delz)) allocate (self%delz(isd:ied,  jsd:jed  , npz))
     case default 
       call abor1_ftn("Create: unknown variable "//trim(self%vars%fldnames(var)))

   end select

enddo

self%hydrostatic = .true.
if (allocated(self%w).and.allocated(self%delz)) self%hydrostatic = .false.

!Extra dynamical variables
if (.not.allocated(self%phis)) allocate(self%phis(isd:ied,jsd:jed    ))

if (.not.allocated(self%vort)) allocate(self%vort(isd:ied,jsd:jed,npz))
if (.not.allocated(self%divg)) allocate(self%divg(isd:ied,jsd:jed,npz))

!CRTM surface variables
if (.not.allocated(self%slmsk )) allocate(self%slmsk (isd:ied,jsd:jed))
if (.not.allocated(self%sheleg)) allocate(self%sheleg(isd:ied,jsd:jed))
if (.not.allocated(self%tsea  )) allocate(self%tsea  (isd:ied,jsd:jed))
if (.not.allocated(self%vtype )) allocate(self%vtype (isd:ied,jsd:jed))
if (.not.allocated(self%stype )) allocate(self%stype (isd:ied,jsd:jed))
if (.not.allocated(self%vfrac )) allocate(self%vfrac (isd:ied,jsd:jed))
if (.not.allocated(self%stc   )) allocate(self%stc   (isd:ied,jsd:jed,4))
if (.not.allocated(self%smc   )) allocate(self%smc   (isd:ied,jsd:jed,4))
if (.not.allocated(self%snwdph)) allocate(self%snwdph(isd:ied,jsd:jed))
if (.not.allocated(self%u_srf )) allocate(self%u_srf (isd:ied,jsd:jed))
if (.not.allocated(self%v_srf )) allocate(self%v_srf (isd:ied,jsd:jed))
if (.not.allocated(self%f10m  )) allocate(self%f10m  (isd:ied,jsd:jed))

! Initialize all domain arrays to zero
call zeros(self)
self%phis   = 0.0_kind_real

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
self%am_i_root_pe = geom%am_i_root_pe

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

if (allocated(self%ud  )) deallocate (self%ud  )
if (allocated(self%vd  )) deallocate (self%vd  )
if (allocated(self%ua  )) deallocate (self%ua  )
if (allocated(self%va  )) deallocate (self%va  )
if (allocated(self%t   )) deallocate (self%t   )
if (allocated(self%delp)) deallocate (self%delp)
if (allocated(self%q   )) deallocate (self%q   )
if (allocated(self%qi  )) deallocate (self%qi  )
if (allocated(self%ql  )) deallocate (self%ql  )
if (allocated(self%o3  )) deallocate (self%o3  )
if (allocated(self%phis)) deallocate (self%phis)
if (allocated(self%w   )) deallocate (self%w   )
if (allocated(self%delz)) deallocate (self%delz)

if (allocated(self%slmsk )) deallocate(self%slmsk )
if (allocated(self%sheleg)) deallocate(self%sheleg)
if (allocated(self%tsea  )) deallocate(self%tsea  )
if (allocated(self%vtype )) deallocate(self%vtype )
if (allocated(self%stype )) deallocate(self%stype )
if (allocated(self%vfrac )) deallocate(self%vfrac )
if (allocated(self%stc   )) deallocate(self%stc   )
if (allocated(self%smc   )) deallocate(self%smc   )
if (allocated(self%snwdph)) deallocate(self%snwdph)
if (allocated(self%u_srf )) deallocate(self%u_srf )
if (allocated(self%v_srf )) deallocate(self%v_srf )
if (allocated(self%f10m  )) deallocate(self%f10m  )

if (allocated(self%vort)) deallocate(self%vort)
if (allocated(self%divg)) deallocate(self%divg)
if (allocated(self%psi )) deallocate(self%psi )
if (allocated(self%chi )) deallocate(self%chi )
if (allocated(self%tv  )) deallocate(self%tv  )
if (allocated(self%ps  )) deallocate(self%ps  )
if (allocated(self%qc  )) deallocate(self%qc  )
if (allocated(self%qic )) deallocate(self%qic )
if (allocated(self%qlc )) deallocate(self%qlc )
if (allocated(self%o3c )) deallocate(self%o3c )

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

!Zero out the entire domain

!Model
if(allocated(self%ud  )) self%ud   = 0.0_kind_real
if(allocated(self%vd  )) self%vd   = 0.0_kind_real
if(allocated(self%ua  )) self%ua   = 0.0_kind_real
if(allocated(self%va  )) self%va   = 0.0_kind_real
if(allocated(self%t   )) self%t    = 0.0_kind_real
if(allocated(self%delp)) self%delp = 0.0_kind_real
if(allocated(self%q   )) self%q    = 0.0_kind_real
if(allocated(self%qi  )) self%qi   = 0.0_kind_real
if(allocated(self%ql  )) self%ql   = 0.0_kind_real
if(allocated(self%o3  )) self%o3   = 0.0_kind_real
if(allocated(self%w   )) self%w    = 0.0_kind_real
if(allocated(self%delz)) self%delz = 0.0_kind_real

!Control
if(allocated(self%psi)) self%psi = 0.0_kind_real
if(allocated(self%chi)) self%chi = 0.0_kind_real
if(allocated(self%tv )) self%tv  = 0.0_kind_real
if(allocated(self%ps )) self%ps  = 0.0_kind_real
if(allocated(self%qc )) self%qc  = 0.0_kind_real
if(allocated(self%qic)) self%qic = 0.0_kind_real
if(allocated(self%qlc)) self%qlc = 0.0_kind_real
if(allocated(self%o3c)) self%o3c = 0.0_kind_real

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine ones(self)
implicit none
type(fv3jedi_increment), intent(inout) :: self

call zeros(self)

!Model
if(allocated(self%ud  )) self%ud   = 1.0_kind_real
if(allocated(self%vd  )) self%vd   = 1.0_kind_real
if(allocated(self%ua  )) self%ua   = 1.0_kind_real
if(allocated(self%va  )) self%va   = 1.0_kind_real
if(allocated(self%t   )) self%t    = 1.0_kind_real
if(allocated(self%delp)) self%delp = 1.0_kind_real
if(allocated(self%q   )) self%q    = 1.0_kind_real
if(allocated(self%qi  )) self%qi   = 1.0_kind_real
if(allocated(self%ql  )) self%ql   = 1.0_kind_real
if(allocated(self%o3  )) self%o3   = 1.0_kind_real
if(allocated(self%w   )) self%w    = 1.0_kind_real
if(allocated(self%delz)) self%delz = 1.0_kind_real

!Control
if(allocated(self%psi)) self%psi = 1.0_kind_real
if(allocated(self%chi)) self%chi = 1.0_kind_real
if(allocated(self%tv )) self%tv  = 1.0_kind_real
if(allocated(self%ps )) self%ps  = 1.0_kind_real
if(allocated(self%qc )) self%qc  = 1.0_kind_real
if(allocated(self%qic)) self%qic = 1.0_kind_real
if(allocated(self%qlc)) self%qlc = 1.0_kind_real
if(allocated(self%o3c)) self%o3c = 1.0_kind_real


end subroutine ones

! ------------------------------------------------------------------------------

subroutine random(self)
use random_vectors_mod
implicit none
type(fv3jedi_increment), intent(inout) :: self
integer :: nq

if(allocated(self%ud  )) call random_vector(self%ud  )
if(allocated(self%vd  )) call random_vector(self%vd  )
if(allocated(self%ua  )) call random_vector(self%ua  )
if(allocated(self%va  )) call random_vector(self%va  )
if(allocated(self%t   )) call random_vector(self%t   )
if(allocated(self%delp)) call random_vector(self%delp)
if(allocated(self%q   )) call random_vector(self%q   )
if(allocated(self%qi  )) call random_vector(self%qi  )
if(allocated(self%ql  )) call random_vector(self%ql  )
if(allocated(self%o3  )) call random_vector(self%o3  )

if(allocated(self%w   )) call random_vector(self%w   )
if(allocated(self%delz)) call random_vector(self%delz)

if(allocated(self%psi )) call random_vector(self%psi )
if(allocated(self%chi )) call random_vector(self%chi )
if(allocated(self%tv  )) call random_vector(self%tv  )
if(allocated(self%ps  )) call random_vector(self%ps  )
if(allocated(self%qc  )) call random_vector(self%qc  )
if(allocated(self%qic )) call random_vector(self%qic )
if(allocated(self%qlc )) call random_vector(self%qlc )
if(allocated(self%o3c )) call random_vector(self%o3c )

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
self%havecrtmfields = rhs%havecrtmfields
self%hydrostatic    = rhs%hydrostatic   
self%calendar_type  = rhs%calendar_type 
self%date           = rhs%date          
self%date_init      = rhs%date_init     

if(allocated(self%ud  )) self%ud   = rhs%ud  
if(allocated(self%vd  )) self%vd   = rhs%vd  
if(allocated(self%ua  )) self%ua   = rhs%ua  
if(allocated(self%va  )) self%va   = rhs%va  
if(allocated(self%t   )) self%t    = rhs%t   
if(allocated(self%delp)) self%delp = rhs%delp
if(allocated(self%q   )) self%q    = rhs%q   
if(allocated(self%qi  )) self%qi   = rhs%qi  
if(allocated(self%ql  )) self%ql   = rhs%ql  
if(allocated(self%o3  )) self%o3   = rhs%o3  
if(allocated(self%w   )) self%w    = rhs%w   
if(allocated(self%delz)) self%delz = rhs%delz

if(allocated(self%psi )) self%psi  = rhs%psi
if(allocated(self%chi )) self%chi  = rhs%chi
if(allocated(self%tv  )) self%tv   = rhs%tv 
if(allocated(self%ps  )) self%ps   = rhs%ps 
if(allocated(self%qc  )) self%qc   = rhs%qc 
if(allocated(self%qic )) self%qic  = rhs%qic
if(allocated(self%qlc )) self%qlc  = rhs%qlc
if(allocated(self%o3c )) self%o3c  = rhs%o3c

self%phis   = rhs%phis
self%slmsk  = rhs%slmsk
self%sheleg = rhs%sheleg
self%tsea   = rhs%tsea
self%vtype  = rhs%vtype
self%stype  = rhs%stype
self%vfrac  = rhs%vfrac
self%stc    = rhs%stc
self%smc    = rhs%smc
self%snwdph = rhs%snwdph
self%u_srf  = rhs%u_srf
self%v_srf  = rhs%v_srf
self%f10m   = rhs%f10m

self%vort = rhs%vort
self%divg = rhs%divg

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ud  )) self%ud   = self%ud   + rhs%ud  
if(allocated(self%vd  )) self%vd   = self%vd   + rhs%vd  
if(allocated(self%ua  )) self%ua   = self%ua   + rhs%ua  
if(allocated(self%va  )) self%va   = self%va   + rhs%va  
if(allocated(self%t   )) self%t    = self%t    + rhs%t   
if(allocated(self%delp)) self%delp = self%delp + rhs%delp
if(allocated(self%q   )) self%q    = self%q    + rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   + rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   + rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   + rhs%o3  
if(allocated(self%w   )) self%w    = self%w    + rhs%w   
if(allocated(self%delz)) self%delz = self%delz + rhs%delz

if(allocated(self%psi )) self%psi  = self%psi  + rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  + rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   + rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   + rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   + rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  + rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  + rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  + rhs%o3c

return
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ud  )) self%ud   = self%ud   * rhs%ud  
if(allocated(self%vd  )) self%vd   = self%vd   * rhs%vd  
if(allocated(self%ua  )) self%ua   = self%ua   * rhs%ua  
if(allocated(self%va  )) self%va   = self%va   * rhs%va  
if(allocated(self%t   )) self%t    = self%t    * rhs%t   
if(allocated(self%delp)) self%delp = self%delp * rhs%delp
if(allocated(self%q   )) self%q    = self%q    * rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   * rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   * rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   * rhs%o3  
if(allocated(self%w   )) self%w    = self%w    * rhs%w   
if(allocated(self%delz)) self%delz = self%delz * rhs%delz

if(allocated(self%psi )) self%psi  = self%psi  * rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  * rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   * rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   * rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   * rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  * rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  * rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  * rhs%o3c

return
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ud  )) self%ud   = self%ud   - rhs%ud  
if(allocated(self%vd  )) self%vd   = self%vd   - rhs%vd  
if(allocated(self%ua  )) self%ua   = self%ua   - rhs%ua  
if(allocated(self%va  )) self%va   = self%va   - rhs%va  
if(allocated(self%t   )) self%t    = self%t    - rhs%t   
if(allocated(self%delp)) self%delp = self%delp - rhs%delp
if(allocated(self%q   )) self%q    = self%q    - rhs%q   
if(allocated(self%qi  )) self%qi   = self%qi   - rhs%qi  
if(allocated(self%ql  )) self%ql   = self%ql   - rhs%ql  
if(allocated(self%o3  )) self%o3   = self%o3   - rhs%o3  
if(allocated(self%w   )) self%w    = self%w    - rhs%w   
if(allocated(self%delz)) self%delz = self%delz - rhs%delz

if(allocated(self%psi )) self%psi  = self%psi  - rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  - rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   - rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   - rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   - rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  - rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  - rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  - rhs%o3c

return
end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz

if(allocated(self%ud  )) self%ud   = zz * self%ud  
if(allocated(self%vd  )) self%vd   = zz * self%vd  
if(allocated(self%ua  )) self%ua   = zz * self%ua  
if(allocated(self%va  )) self%va   = zz * self%va  
if(allocated(self%t   )) self%t    = zz * self%t   
if(allocated(self%delp)) self%delp = zz * self%delp
if(allocated(self%q   )) self%q    = zz * self%q   
if(allocated(self%qi  )) self%qi   = zz * self%qi  
if(allocated(self%ql  )) self%ql   = zz * self%ql  
if(allocated(self%o3  )) self%o3   = zz * self%o3  
if(allocated(self%w   )) self%w    = zz * self%w   
if(allocated(self%delz)) self%delz = zz * self%delz

if(allocated(self%psi )) self%psi  = zz * self%psi
if(allocated(self%chi )) self%chi  = zz * self%chi
if(allocated(self%tv  )) self%tv   = zz * self%tv 
if(allocated(self%ps  )) self%ps   = zz * self%ps 
if(allocated(self%qc  )) self%qc   = zz * self%qc 
if(allocated(self%qic )) self%qic  = zz * self%qic
if(allocated(self%qlc )) self%qlc  = zz * self%qlc
if(allocated(self%o3c )) self%o3c  = zz * self%o3c

return
end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy_inc(self,zz,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_increment), intent(in)    :: rhs

if(allocated(self%ud  )) self%ud   = self%ud   + zz * rhs%ud  
if(allocated(self%vd  )) self%vd   = self%vd   + zz * rhs%vd  
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

if(allocated(self%psi )) self%psi  = self%psi  + zz * rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  + zz * rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   + zz * rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   + zz * rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   + zz * rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  + zz * rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  + zz * rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  + zz * rhs%o3c

return
end subroutine axpy_inc

! ------------------------------------------------------------------------------

subroutine axpy_state(self,zz,rhs)
implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_state), intent(in)    :: rhs

if(allocated(self%ud  )) self%ud   = self%ud   + zz * rhs%ud  
if(allocated(self%vd  )) self%vd   = self%vd   + zz * rhs%vd  
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

if(allocated(self%psi )) self%psi  = self%psi  + zz * rhs%psi
if(allocated(self%chi )) self%chi  = self%chi  + zz * rhs%chi
if(allocated(self%tv  )) self%tv   = self%tv   + zz * rhs%tv 
if(allocated(self%ps  )) self%ps   = self%ps   + zz * rhs%ps 
if(allocated(self%qc  )) self%qc   = self%qc   + zz * rhs%qc 
if(allocated(self%qic )) self%qic  = self%qic  + zz * rhs%qic
if(allocated(self%qlc )) self%qlc  = self%qlc  + zz * rhs%qlc
if(allocated(self%o3c )) self%o3c  = self%o3c  + zz * rhs%o3c

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

!ud
if (allocated(inc1%ud)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%ud(i,j,k) * inc2%ud(i,j,k)
      enddo
    enddo
  enddo
endif

!vd
if (allocated(inc1%vd)) then
  do k = 1,inc1%npz
    do j = inc1%jsc,inc1%jec
      do i = inc1%isc,inc1%iec
        zp = zp + inc1%vd(i,j,k) * inc2%vd(i,j,k)
      enddo
    enddo
  enddo
endif

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

!Get global dot product
call f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
if (inc1%am_i_root_pe) print*, "Dot product test result: ", zprod

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
  if(allocated(rhs%ud  )) self%ud   = self%ud   + rhs%ud  
  if(allocated(rhs%vd  )) self%vd   = self%vd   + rhs%vd  
  if(allocated(rhs%ua  )) self%ua   = self%ua   + rhs%ua  
  if(allocated(rhs%va  )) self%va   = self%va   + rhs%va  
  if(allocated(rhs%t   )) self%t    = self%t    + rhs%t   
  if(allocated(rhs%delp)) self%delp = self%delp + rhs%delp
  if(allocated(rhs%q   )) self%q    = self%q    + rhs%q   
  if(allocated(rhs%qi  )) self%qi   = self%qi   + rhs%qi  
  if(allocated(rhs%ql  )) self%ql   = self%ql   + rhs%ql  
  if(allocated(rhs%o3  )) self%o3   = self%o3   + rhs%o3  
  if(allocated(rhs%w   )) self%w    = self%w    + rhs%w   
  if(allocated(rhs%delz)) self%delz = self%delz + rhs%delz 
  if(allocated(rhs%psi )) self%psi  = self%psi  + rhs%psi
  if(allocated(rhs%chi )) self%chi  = self%chi  + rhs%chi
  if(allocated(rhs%tv  )) self%tv   = self%tv   + rhs%tv 
  if(allocated(rhs%ps  )) self%ps   = self%ps   + rhs%ps 
  if(allocated(rhs%qc  )) self%qc   = self%qc   + rhs%qc 
  if(allocated(rhs%qic )) self%qic  = self%qic  + rhs%qic
  if(allocated(rhs%qlc )) self%qlc  = self%qlc  + rhs%qlc
  if(allocated(rhs%o3c )) self%o3c  = self%o3c  + rhs%o3c
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

  if(allocated(lhs%ud  )) lhs%ud   = x1%ud   - x2%ud  
  if(allocated(lhs%vd  )) lhs%vd   = x1%vd   - x2%vd  
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

  if(allocated(lhs%psi )) lhs%psi  = x1%psi  - x2%psi
  if(allocated(lhs%chi )) lhs%chi  = x1%chi  - x2%chi
  if(allocated(lhs%tv  )) lhs%tv   = x1%tv   - x2%tv 
  if(allocated(lhs%ps  )) lhs%ps   = x1%ps   - x2%ps 
  if(allocated(lhs%qc  )) lhs%qc   = x1%qc   - x2%qc 
  if(allocated(lhs%qic )) lhs%qic  = x1%qic  - x2%qic
  if(allocated(lhs%qlc )) lhs%qlc  = x1%qlc  - x2%qlc
  if(allocated(lhs%o3c )) lhs%o3c  = x1%o3c  - x2%o3c

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

!u
pstat(1,1) = minval(inc%ud(isc:iec,jsc:jec,:))
pstat(2,1) = maxval(inc%ud(isc:iec,jsc:jec,:))
pstat(3,1) = sqrt((sum(inc%ud(isc:iec,jsc:jec,:))/gs)**2)

!v
pstat(1,2) = minval(inc%vd(isc:iec,jsc:jec,:))
pstat(2,2) = maxval(inc%vd(isc:iec,jsc:jec,:))
pstat(3,2) = sqrt((sum(inc%vd(isc:iec,jsc:jec,:))/gs)**2)

!pt
pstat(1,3) = minval(inc%t(isc:iec,jsc:jec,:))
pstat(2,3) = maxval(inc%t(isc:iec,jsc:jec,:))
pstat(3,3) = sqrt((sum(inc%t(isc:iec,jsc:jec,:))/gs)**2)

!delp
pstat(1,4) = minval(inc%delp(isc:iec,jsc:jec,:))
pstat(2,4) = maxval(inc%delp(isc:iec,jsc:jec,:))
pstat(3,4) = sqrt((sum(inc%delp(isc:iec,jsc:jec,:))/gs)**2)

!q
pstat(1,5) = minval(inc%q(isc:iec,jsc:jec,:))
pstat(2,5) = maxval(inc%q(isc:iec,jsc:jec,:))
pstat(3,5) = sqrt((sum(inc%q(isc:iec,jsc:jec,:))/gs)**2)

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

!ud
if (allocated(inc%ud)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%ud(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!vd
if (allocated(inc%vd)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + inc%vd(i,j,k)**2
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

!Get global values
call f_comm%allreduce(zz,prms,fckit_mpi_sum())
call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

!if (ierr .ne. 0) then
!   print *,'error in incrms/mpi_allreduce, error code=',ierr
!endif
prms = sqrt(prms/real(iisum,kind_real))

!Print for debugging
if (inc%am_i_root_pe) print *,'incrms: prms = ', prms
if (inc%am_i_root_pe) print *,'incrms: iisum = ', iisum

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
          self%ud  (ixdir(idir),iydir(idir),ildir) = 1.0
       else if (ifdir == 2) then
          self%vd  (ixdir(idir),iydir(idir),ildir) = 1.0
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

if (allocated(self%ud)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(1)%fld(imga,jl,1,1) = self%ud  (jx,jy,jl)
          ug%grid(1)%fld(imga,jl,2,1) = self%vd  (jx,jy,jl)
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

if (allocated(self%ud)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          ug%grid(igrid)%fld(imga,jl,1,1) = self%ud  (jx,jy,jl)
          ug%grid(igrid)%fld(imga,jl,2,1) = self%vd  (jx,jy,jl)
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

if (allocated(self%ud)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%ud  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,1,1)
          self%vd  (jx,jy,jl) = ug%grid(1)%fld(imga,jl,2,1)
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

if (allocated(self%ud)) then

  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      do jl=1,self%npz
          self%ud  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,1,1)
          self%vd  (jx,jy,jl) = ug%grid(igrid)%fld(imga,jl,2,1)
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

subroutine getvalues_tl(geom, inc, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_tl
use pressure_vt_mod
implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
type(fv3jedi_increment),      intent(inout) :: inc 
type(ioda_locs),          intent(in)    :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_tl'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

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
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
isd = inc%isd
ied = inc%ied
jsd = inc%jsd
jed = inc%jed
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isd:ied,jsd:jed,npz+1))
allocate(geovalm(isd:ied,jsd:jed,npz))


! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv
 
  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  select case (trim(vars%fldnames(jvar)))
   
  case ("temperature")
  
    nvl = npz
    do_interp = .true.
    geovalm = inc%t
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geovalm = inc%q
    geoval => geovalm

  case ("virtual_temperature")

    nvl = inc%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, inc%t, traj%q, inc%q )
    geovalm = inc%t
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = inc%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, inc%q, geovalm)
    geoval => geovalm  

  case ("air_pressure")

    nvl = geom%npz
    do_interp = .true.
    call delp_to_p_tl(geom,inc%delp,geovalm)
    geoval => geovalm

  case ("air_pressure_levels")
 
  case ("geopotential_height")

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
  call interp_checks(myname, inc, locs, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = 1, nvl
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          mod_increment(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call traj%bump%apply_obsop(mod_increment,obs_increment)
      gom%geovals(jvar)%vals(jlev,:) = obs_increment(:,1)
    enddo
  else
    gom%geovals(jvar)%vals(nvl,:) = obs_increment(:,1)
  endif

  nullify(geoval)

enddo

deallocate(geovalm,geovale)

deallocate(mod_increment)
deallocate(obs_increment)

end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(geom, inc, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_ad
use pressure_vt_mod
implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
type(fv3jedi_increment),      intent(inout) :: inc 
type(ioda_locs),           intent(in)   :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_ad'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

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
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
isd = inc%isd
ied = inc%ied
jsd = inc%jsd
jed = inc%jed
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isd:ied,jsd:jed,npz+1))
allocate(geovalm(isd:ied,jsd:jed,npz))

geovale = 0.0_kind_real
geovalm = 0.0_kind_real

! Interpolate increment to obs locations using pre-calculated weights
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

  case ("specific_humidity")

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
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_pressure_levels")

  case ("geopotential_height")

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
      obs_increment(:,1) = gom%geovals(jvar)%vals(jlev,:)
      gom%geovals(jvar)%vals(jlev,:) = 0.0_kind_real
      call traj%bump%apply_obsop_ad(obs_increment,mod_increment)
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          geoval(ji, jj, jlev) = mod_increment(ii, 1)
        enddo
      enddo
    enddo
  else
    obs_increment(:,1) = gom%geovals(jvar)%vals(nvl,:)
  endif

  !Part 3, back to state variables
  !-------------------------------
 
  select case (trim(vars%fldnames(jvar)))
 
  case ("temperature")
    inc%t = geovalm

  case ("specific_humidity")
    inc%q = geovalm

  case ("virtual_temperature")
    
    inc%t = geovalm
    call T_to_Tv_ad(geom, traj%t, inc%t, traj%q, inc%q )

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    call crtm_mixratio_ad(geom, traj%q, inc%q, geovalm)

  case ("air_pressure")

    call delp_to_p_ad(geom,inc%delp,geovalm)

  case ("air_pressure_levels")
 
  case ("geopotential_height")

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

deallocate(mod_increment)
deallocate(obs_increment)

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

use fckit_mpi_module, only: fckit_mpi_comm
use fv3jedi_geom_mod, only: fv3jedi_geom
use type_bump, only: bump_type

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

type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

! Each bump%nam%prefix must be distinct
! -------------------------------------
bumpcount = bumpcount + 1
write(cbumpcount,"(I0.5)") bumpcount
bump_nam_prefix = 'fv3jedi_bump_data_'//cbumpcount


!Get the Solution dimensions
!---------------------------
mod_num = (geom%iec - geom%isc + 1) * (geom%jec - geom%jsc + 1)


!Calculate interpolation weight using BUMP
!-----------------------------------------
allocate( mod_lat(mod_num), mod_lon(mod_num) )
mod_lat = reshape( rad2deg*geom%grid_lat(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] )  
mod_lon = reshape( rad2deg*geom%grid_lon(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] ) - 180.0_kind_real

!Important namelist options
call bump%nam%init

bump%nam%prefix = bump_nam_prefix   ! Prefix for files output
bump%nam%nobs = locs%nlocs          ! Number of observations
bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)
bump%nam%obsdis = 'local'           ! Local or BUMP may try to redistribute obs
bump%nam%diag_interp = 'bilin'
bump%nam%local_diag = .false.

!Less important namelist options (should not be changed)
bump%nam%default_seed = .true.
bump%nam%new_hdiag = .false.
bump%nam%new_nicas = .false.
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
call bump%setup_online( f_comm%communicator(),mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                                nobs=locs%nlocs,lonobs=locs%lon(:)-180.0_kind_real,latobs=locs%lat(:) )

!Release memory
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate( mod_lat, mod_lon )

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, inc, locs, vars, gom, jvar)
implicit none
character(len=*), intent(in) :: cop
type(fv3jedi_increment), intent(in) :: inc
type(ioda_locs), intent(in)     :: locs
type(ufo_vars), intent(in)      :: vars
type(ufo_geovals), intent(in)   :: gom
integer, intent(in)             :: jvar

character(len=255) :: cinfo

cinfo="fv3jedi_increment:checks "//trim(cop)//" : "

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

end module fv3jedi_increment_mod
