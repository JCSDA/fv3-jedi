! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Handle fields for the FV3JEDI odel

module fv3jedi_fields_mod

use iso_c_binding
use config_mod
use datetime_mod
use fv3jedi_geom_mod
use ufo_vars_mod
use fv3jedi_kinds
use ioda_locs_mod
use ufo_geovals_mod

use fv3jedi_fields_utils_mod
use fv3jedi_fields_io_mod
use fv3jedi_constants, only: rad2deg, constoz

use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj

use fv3jedi_vars_mod, only: fv3jedi_vars

implicit none
private

public :: create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, fldrms, &
          change_resol, getvalues, getvalues_tl, getvalues_ad, &
          ug_coord, field_to_ug, field_from_ug, dirac, &
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
type(fv3jedi_field), intent(inout) :: self

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
type(fv3jedi_field), intent(inout) :: self

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
type(fv3jedi_field), intent(inout) :: self

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
type(fv3jedi_field), intent(inout) :: self
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
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

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
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

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
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

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
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

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
type(fv3jedi_field), intent(inout) :: self
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

subroutine axpy(self,zz,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_field), intent(in)    :: rhs

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
end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

implicit none
type(fv3jedi_field), intent(in) :: fld1, fld2
real(kind=kind_real), intent(inout) :: zprod
real(kind=kind_real) :: zp
integer :: i,j,k
integer :: ierr
type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

zp=0.0_kind_real

!ud
if (allocated(fld1%ud)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%ud(i,j,k) * fld2%ud(i,j,k)
      enddo
    enddo
  enddo
endif

!vd
if (allocated(fld1%vd)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%vd(i,j,k) * fld2%vd(i,j,k)
      enddo
    enddo
  enddo
endif

!ua
if (allocated(fld1%ua)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%ua(i,j,k) * fld2%ua(i,j,k)
      enddo
    enddo
  enddo
endif

!va
if (allocated(fld1%va)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%va(i,j,k) * fld2%va(i,j,k)
      enddo
    enddo
  enddo
endif

!t
if (allocated(fld1%t)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%t(i,j,k) * fld2%t(i,j,k)
      enddo
    enddo
  enddo
endif

!delp
if (allocated(fld1%delp)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%delp(i,j,k) * fld2%delp(i,j,k)
      enddo
    enddo
  enddo
endif

!q
if (allocated(fld1%q)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%q(i,j,k) * fld2%q(i,j,k)
      enddo
    enddo
  enddo
endif

!qi
if (allocated(fld1%qi)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%qi(i,j,k) * fld2%qi(i,j,k)
      enddo
    enddo
  enddo
endif

!ql
if (allocated(fld1%ql)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%ql(i,j,k) * fld2%ql(i,j,k)
      enddo
    enddo
  enddo
endif

!o3
if (allocated(fld1%o3)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%o3(i,j,k) * fld2%o3(i,j,k)
      enddo
    enddo
  enddo
endif

!psi
if (allocated(fld1%psi)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%psi(i,j,k) * fld2%psi(i,j,k)
      enddo
    enddo
  enddo
endif

!chi
if (allocated(fld1%chi)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%chi(i,j,k) * fld2%chi(i,j,k)
      enddo
    enddo
  enddo
endif

!tv
if (allocated(fld1%tv)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%tv(i,j,k) * fld2%tv(i,j,k)
      enddo
    enddo
  enddo
endif

!ps
if (allocated(fld1%ps)) then
  do j = fld1%jsc,fld1%jec
    do i = fld1%isc,fld1%iec
      zp = zp + fld1%ps(i,j) * fld2%ps(i,j)
    enddo
  enddo
endif

!qc
if (allocated(fld1%qc)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%qc(i,j,k) * fld2%qc(i,j,k)
      enddo
    enddo
  enddo
endif

!qic
if (allocated(fld1%qic)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%qic(i,j,k) * fld2%qic(i,j,k)
      enddo
    enddo
  enddo
endif

!qlc
if (allocated(fld1%qlc)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%qlc(i,j,k) * fld2%qlc(i,j,k)
      enddo
    enddo
  enddo
endif

!o3c
if (allocated(fld1%o3c)) then
  do k = 1,fld1%npz
    do j = fld1%jsc,fld1%jec
      do i = fld1%isc,fld1%iec
        zp = zp + fld1%o3c(i,j,k) * fld2%o3c(i,j,k)
      enddo
    enddo
  enddo
endif

!Get global dot product
call f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
if (fld1%am_i_root_pe) print*, "Dot product test result: ", zprod

return
end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: self
type(fv3jedi_field), intent(in)    :: rhs

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

   call abor1_ftn("fv3jedi fields:  diff_incr not implemented for low res increment yet")

endif

return
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(fv3jedi_field), intent(inout) :: fld
type(fv3jedi_field), intent(in)    :: rhs

integer :: check
check = (rhs%iec-rhs%isc+1) - (fld%iec-fld%isc+1)

if (check==0) then
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
!! cannot currently be used to initiate a forecast with gfs.
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

  type(fv3jedi_field), intent(inout) :: fld !< Fields
  type(fv3jedi_geom),  intent(inout) :: geom    !< Geometry 
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: i,j,k
  real(kind=kind_real) :: rlat, rlon, z
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  type(fv_atmos_type), allocatable :: FV_AtmIC(:)
  real(kind=kind_real)             :: DTdummy = 900.0
  logical, allocatable             :: grids_on_this_pe(:)
  integer                          :: p_split = 1

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

        call invent_state(fld,c_conf,geom)

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
        fld%ud = FV_AtmIC(1)%u
        fld%vd = FV_AtmIC(1)%v
        fld%t = FV_AtmIC(1)%pt
        fld%delp = FV_AtmIC(1)%delp
        fld%q = FV_AtmIC(1)%q(:,:,:,1)
        fld%phis = FV_AtmIC(1)%phis
        geom%ak = FV_AtmIC(1)%ak
        geom%ak = FV_AtmIC(1)%ak
        geom%ptop = FV_AtmIC(1)%ptop
        if (.not. fld%hydrostatic) then
           fld%w = FV_AtmIC(1)%w
           fld%delz = FV_AtmIC(1)%delz
        endif

        !Deallocate temporary FV_Atm fv3 structure
        call deallocate_fv_atmos_type(FV_AtmIC(1))
        deallocate(FV_AtmIC)
        deallocate(grids_on_this_pe)

     Case ("dcmip-test-1-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
                                               phis0,ps,rho0,hum0,q1,q2,q3,q4)

              fld%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                                                  phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 fld%ud(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 fld%vd(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
                 If (.not.fld%hydrostatic) fld%w(i,j,k) = w0
                 fld%t(i,j,k) = t0
                 fld%delp(i,j,k) = pe2-pe1
                 fld%q (i,j,k) = hum0
                 fld%qi(i,j,k) = q1
                 fld%ql(i,j,k) = q2
                 fld%o3(i,j,k) = q3
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                          t0,phis0,ps,rho0,hum0,q1)

              fld%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                             t0,phis0,ps,rho0,hum0,q1)

                 fld%ud(i,j,k) = u0 !ATTN comment above
                 fld%vd(i,j,k) = v0
                 If (.not.fld%hydrostatic) fld%w(i,j,k) = w0
                 fld%t(i,j,k) = t0
                 fld%delp(i,j,k) = pe2-pe1
                 fld%q (i,j,k) = hum0
                 fld%qi(i,j,k) = q1
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                      t0,phis0,ps,rho0,hum0)

              fld%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0)

                 fld%ud(i,j,k) = u0 !ATTN comment above
                 fld%vd(i,j,k) = v0
                 If (.not.fld%hydrostatic) fld%w(i,j,k) = w0
                 fld%t(i,j,k) = t0
                 fld%delp(i,j,k) = pe2-pe1
                 fld%q(i,j,k) = hum0
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

              fld%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

                 fld%ud(i,j,k) = u0 !ATTN comment above
                 fld%vd(i,j,k) = v0
                 If (.not.fld%hydrostatic) fld%w(i,j,k) = w0
                 fld%t(i,j,k) = t0
                 fld%delp(i,j,k) = pe2-pe1
                 fld%q(i,j,k) = hum0
                 
              enddo
           enddo
        enddo

     Case Default

        call invent_state(fld,c_conf,geom)

     End Select int_option
        
end subroutine analytic_IC
  
! ------------------------------------------------------------------------------
subroutine invent_state(flds,config,geom)

use kinds

implicit none

type(fv3jedi_field), intent(inout) :: flds    !< Model fields
type(c_ptr), intent(in)            :: config  !< Configuration structure
type(fv3jedi_geom),  intent(in)    :: geom

integer :: i,j,k

!ua
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      flds%ua(i,j,k) = cos(0.25*geom%grid_lon(i,j)) + cos(0.25*geom%grid_lat(i,j))
    enddo
  enddo
enddo

!va
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      flds%va(i,j,k) = 1.0_kind_real
    enddo
  enddo
enddo

!t
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      flds%t(i,j,k) = cos(0.25*geom%grid_lon(i,j)) + cos(0.25*geom%grid_lat(i,j))
    enddo
  enddo
enddo

!delp
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      flds%delp(i,j,k) = k
    enddo
  enddo
enddo

!q
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      flds%q(i,j,k) = 0.0
    enddo
  enddo
enddo

return
end subroutine invent_state

! ------------------------------------------------------------------------------

subroutine read_file(geom, fld, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_field), intent(inout) :: fld      !< Fields
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call read_fms_restart(geom, fld, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call read_geos_restart(geom, fld, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_fields read: restart type not supported")
  endif

  return

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, fld, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_field), intent(in)    :: fld      !< Fields
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call write_fms_restart(geom, fld, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call write_geos_restart(geom, fld, c_conf, vdate)
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

isc = fld%isc
iec = fld%iec
jsc = fld%jsc
jec = fld%jec

gs = (iec-isc+1)*(jec-jsc+1)*fld%npz

!u
pstat(1,1) = minval(fld%ud(isc:iec,jsc:jec,:))
pstat(2,1) = maxval(fld%ud(isc:iec,jsc:jec,:))
pstat(3,1) = sqrt((sum(fld%ud(isc:iec,jsc:jec,:))/gs)**2)

!v
pstat(1,2) = minval(fld%vd(isc:iec,jsc:jec,:))
pstat(2,2) = maxval(fld%vd(isc:iec,jsc:jec,:))
pstat(3,2) = sqrt((sum(fld%vd(isc:iec,jsc:jec,:))/gs)**2)

!pt
pstat(1,3) = minval(fld%t(isc:iec,jsc:jec,:))
pstat(2,3) = maxval(fld%t(isc:iec,jsc:jec,:))
pstat(3,3) = sqrt((sum(fld%t(isc:iec,jsc:jec,:))/gs)**2)

!delp
pstat(1,4) = minval(fld%delp(isc:iec,jsc:jec,:))
pstat(2,4) = maxval(fld%delp(isc:iec,jsc:jec,:))
pstat(3,4) = sqrt((sum(fld%delp(isc:iec,jsc:jec,:))/gs)**2)

!q
pstat(1,5) = minval(fld%q(isc:iec,jsc:jec,:))
pstat(2,5) = maxval(fld%q(isc:iec,jsc:jec,:))
pstat(3,5) = sqrt((sum(fld%q(isc:iec,jsc:jec,:))/gs)**2)

return

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
use fckit_mpi_module, only : fckit_mpi_comm, fckit_mpi_sum
implicit none
type(fv3jedi_field), intent(in) :: fld
real(kind=kind_real), intent(out) :: prms

real(kind=kind_real) :: zz
integer i,j,k,ii,nt,ierr,npes,iisum
integer :: isc,iec,jsc,jec,npz
type(fckit_mpi_comm) :: f_comm

isc = fld%isc
iec = fld%iec
jsc = fld%jsc
jec = fld%jec
npz = fld%npz

f_comm = fckit_mpi_comm()

zz = 0.0_kind_real
prms = 0.0_kind_real
ii = 0

!ud
if (allocated(fld%ud)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%ud(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!vd
if (allocated(fld%vd)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%vd(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!t
if (allocated(fld%t)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%t(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!delp
if (allocated(fld%delp)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%delp(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!q
if (allocated(fld%q)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%q(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qi
if (allocated(fld%qi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%qi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ql
if (allocated(fld%ql)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%ql(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!o3
if (allocated(fld%o3)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%o3(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!psi
if (allocated(fld%psi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%psi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!chi
if (allocated(fld%chi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%chi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!tv
if (allocated(fld%tv)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%tv(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ps
if (allocated(fld%ps)) then
  do j = jsc,jec
    do i = isc,iec
      zz = zz + fld%ps(i,j)**2
      ii = ii + 1
    enddo
  enddo
endif

!qc
if (allocated(fld%qc)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%qc(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qic
if (allocated(fld%qic)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%qic(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qlc
if (allocated(fld%qlc)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%qlc(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!o3c
if (allocated(fld%o3c)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + fld%o3c(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!Get global values
call f_comm%allreduce(zz,prms,fckit_mpi_sum())
call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

!if (ierr .ne. 0) then
!   print *,'error in fldrms/mpi_allreduce, error code=',ierr
!endif
prms = sqrt(prms/real(iisum,kind_real))

!Print for debugging
if (fld%am_i_root_pe) print *,'fldrms: prms = ', prms
if (fld%am_i_root_pe) print *,'fldrms: iisum = ', iisum

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf, geom)

use iso_c_binding

implicit none
type(fv3jedi_field), intent(inout) :: self
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
if (ndir<1) call abor1_ftn("fv3jedi_fields:dirac non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%npx)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid ixdir")
endif
if (any(iydir<1).or.any(iydir>geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_fields:dirac invalid iydir")
endif
if ((ildir<1).or.(ildir>self%npz)) then
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
   if (geom%ntile == itiledir .and. &
       ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
       iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then
       ! If so, perturb desired field and level
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
type(fv3jedi_field), intent(in) :: self
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
type(fv3jedi_field), intent(in) :: self
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

subroutine field_to_ug(self, ug, colocated)
use unstructured_grid_mod
implicit none
type(fv3jedi_field), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
integer, intent(in) :: colocated

integer :: imga,jx,jy,jl
real(kind=kind_real),allocatable :: ptmp(:,:,:)
integer :: igrid

! Copy colocated
ug%colocated = colocated

! Define size
call ug_size(self, ug)

! Allocate unstructured grid field
call allocate_unstructured_grid_field(ug)

! Copy field

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

end subroutine field_to_ug

! -----------------------------------------------------------------------------

subroutine field_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(fv3jedi_field), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug

integer :: imga,jx,jy,jl
real(kind=kind_real),allocatable :: ptmp(:,:,:)
integer :: igrid

! Copy field

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


end subroutine field_from_ug

! ------------------------------------------------------------------------------

subroutine getvalues(geom, fld, locs, vars, gom, traj)

use surface_vt_mod
use pressure_vt_mod
use tmprture_vt_mod
use moisture_vt_mod, only: crtm_ade_efr, crtm_mixratio
use wind_vt_mod
use height_vt_mod,   only: geop_height
use type_bump, only: bump_type

implicit none
type(fv3jedi_geom),                         intent(inout) :: geom 
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
logical,  parameter                ::use_compress = .true.  !!could be a fv3 namelist option?


! Grid convenience
! ----------------
isc = fld%isc
iec = fld%iec
jsc = fld%jsc
jec = fld%jec
isd = fld%isd
ied = fld%ied
jsd = fld%jsd
jed = fld%jed
npz = fld%npz

ngrid = (iec-isc+1)*(jec-jsc+1)
nobs = locs%nlocs 

! Initialize the interpolation
! ----------------------------
if (present(traj)) then

  pbump => traj%bump

  if (.not. traj%lalloc) then
  
     traj%ngrid = ngrid
     traj%nobs = nobs
   
     if (.not.allocated(traj%t)) allocate(traj%t(isd:ied,jsd:jed,1:npz))
     if (.not.allocated(traj%q)) allocate(traj%q(isd:ied,jsd:jed,1:npz))
  
     traj%t = fld%t
     traj%q = fld%q
 
     pbumpa => traj%lalloc

  endif

else

  pbump => bump
  bump_alloc = .false.
  pbumpa => bump_alloc

endif

if (.not. pbumpa) then
   call initialize_bump(geom, locs, pbump)
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

call delp_to_pe_p_logp(geom,fld%delp,prsi,prs,logp)

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
  call crtm_surface( geom, nobs, ngrid, locs%lat(:), locs%lon(:), &
                     fld%slmsk, fld%sheleg, fld%tsea, fld%vtype, &
                     fld%stype, fld%vfrac, fld%stc, fld%smc, fld%snwdph, &
                     fld%u_srf,fld%v_srf,fld%f10m, &
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
      if (fld%slmsk(i,j) == 0) water_coverage_m(i,j) = 1.0_kind_real
    enddo
  enddo
  
  call crtm_ade_efr( geom,prsi,fld%t,fld%delp,water_coverage_m,fld%q, &
                     fld%ql,fld%qi,ql_ade,qi_ade,ql_efr,qi_efr )
  
  call crtm_mixratio(geom,fld%q,qmr)

endif


!write(*,*)'interp model    t min, max= ',minval(fld%t),maxval(fld%t)
!write(*,*)'interp model delp min, max= ',minval(fld%delp),maxval(fld%delp)

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
    geovalm = fld%t
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geovalm = fld%q
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    call T_to_Tv(geom,fld%t,fld%q,geovalm)
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

  case ("geopotential_height")
    call geop_height(geom,prs,prsi,fld%t,fld%q,fld%phis,use_compress,geovalm)
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("mass_concentration_of_ozone_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = fld%o3 * constoz
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

if (.not. present(traj)) then
  call pbump%dealloc
endif

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

subroutine getvalues_tl(geom, fld, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_tl
use pressure_vt_mod
implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
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
isc = fld%isc
iec = fld%iec
jsc = fld%jsc
jec = fld%jec
isd = fld%isd
ied = fld%ied
jsd = fld%jsd
jed = fld%jed
npz = fld%npz


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
    geovalm = fld%t
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geovalm = fld%q
    geoval => geovalm

  case ("virtual_temperature")

    nvl = fld%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, fld%t, traj%q, fld%q )
    geovalm = fld%t
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = fld%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, fld%q, geovalm)
    geoval => geovalm  

  case ("air_pressure")

    nvl = geom%npz
    do_interp = .true.
    call delp_to_p_tl(geom,fld%delp,geovalm)
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
      call traj%bump%apply_obsop(mod_field,obs_field)
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

subroutine getvalues_ad(geom, fld, locs, vars, gom, traj)

use tmprture_vt_mod
use moisture_vt_mod, only: crtm_mixratio_ad
use pressure_vt_mod
implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
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
isc = fld%isc
iec = fld%iec
jsc = fld%jsc
jec = fld%jec
isd = fld%isd
ied = fld%ied
jsd = fld%jsd
jed = fld%jed
npz = fld%npz


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
      obs_field(:,1) = gom%geovals(jvar)%vals(jlev,:)
      gom%geovals(jvar)%vals(jlev,:) = 0.0_kind_real
      call traj%bump%apply_obsop_ad(obs_field,mod_field)
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
    fld%t = geovalm

  case ("specific_humidity")
    fld%q = geovalm

  case ("virtual_temperature")
    
    fld%t = geovalm
    call T_to_Tv_ad(geom, traj%t, fld%t, traj%q, fld%q )

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    call crtm_mixratio_ad(geom, traj%q, fld%q, geovalm)

  case ("air_pressure")

    call delp_to_p_ad(geom,fld%delp,geovalm)

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

if (x1%npx /= x2%npx .or.  x1%npz /= x2%npz) then
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
