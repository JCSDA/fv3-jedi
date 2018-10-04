! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_tlm_mod

use iso_c_binding
use config_mod
use duration_mod

use kinds
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 
use fv3jedi_traj_mod, only: fv3jedi_traj

use fv3jedi_lm_mod, only: fv3jedi_lm_type

implicit none
private

public :: fv3jedi_tlm 
public :: tlm_create
public :: tlm_delete
public :: tlm_initialize_tl
public :: tlm_initialize_ad
public :: tlm_step_tl
public :: tlm_step_ad
public :: tlm_finalize_tl
public :: tlm_finalize_ad

! ------------------------------------------------------------------------------

!> Fortran derived type to hold tlm definition
type:: fv3jedi_tlm
  type(fv3jedi_lm_type) :: fv3jedi_lm  !<Linearized model object
end type fv3jedi_tlm

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine tlm_create(self, geom, c_conf)

implicit none
type(fv3jedi_tlm),  intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom
type(c_ptr),        intent(in)    :: c_conf

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt


! Model time step
! ---------------
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)


! Model configuration and creation
! --------------------------------
self%fv3jedi_lm%conf%do_dyn     = config_get_int(c_conf,"lm_do_dyn")
self%fv3jedi_lm%conf%do_phy_trb = config_get_int(c_conf,"lm_do_trb")
self%fv3jedi_lm%conf%do_phy_mst = config_get_int(c_conf,"lm_do_mst")

call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

end subroutine tlm_create

! ------------------------------------------------------------------------------

subroutine tlm_delete(self)

implicit none
type(fv3jedi_tlm), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine tlm_delete

! ------------------------------------------------------------------------------

subroutine tlm_initialize_ad(self, inc)

implicit none
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%init_ad()

end subroutine tlm_initialize_ad

! ------------------------------------------------------------------------------

subroutine tlm_initialize_tl(self, inc)

implicit none
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%init_tl()

end subroutine tlm_initialize_tl

! ------------------------------------------------------------------------------

subroutine tlm_step_ad(geom, self, incr, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: incr
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call lm_to_incr_ad(geom,self%fv3jedi_lm,incr)
call self%fv3jedi_lm%step_ad()
call incr_to_lm_ad(geom,incr,self%fv3jedi_lm)

end subroutine tlm_step_ad

! ------------------------------------------------------------------------------

subroutine tlm_step_tl(geom, self, incr, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: incr
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call incr_to_lm_tl(geom,incr,self%fv3jedi_lm)
call self%fv3jedi_lm%step_tl()
call lm_to_incr_tl(geom,self%fv3jedi_lm,incr)

end subroutine tlm_step_tl

! ------------------------------------------------------------------------------

subroutine tlm_finalize_ad(self, inc)

implicit none
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%final_ad()

end subroutine tlm_finalize_ad

! ------------------------------------------------------------------------------

subroutine tlm_finalize_tl(self, inc)

implicit none
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%final_tl()

end subroutine tlm_finalize_tl

! ------------------------------------------------------------------------------

subroutine incr_to_lm_tl(geom, incr, lm)

use wind_vt_mod, only: a2d

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: incr
type(fv3jedi_lm_type),   intent(inout) :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

allocate(ud(incr%isc:incr%iec  ,incr%jsc:incr%jec+1,1:incr%npz))
allocate(vd(incr%isc:incr%iec+1,incr%jsc:incr%jec  ,1:incr%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

call a2d(geom, incr%ua(incr%isc:incr%iec,incr%jsc:incr%jec,:), &
               incr%va(incr%isc:incr%iec,incr%jsc:incr%jec,:), &
                    ud(incr%isc:incr%iec  ,incr%jsc:incr%jec+1,:), &
                    vd(incr%isc:incr%iec+1,incr%jsc:incr%jec  ,:))

lm%pert%u    = ud(incr%isc:incr%iec,incr%jsc:incr%jec,:)
lm%pert%v    = vd(incr%isc:incr%iec,incr%jsc:incr%jec,:)
lm%pert%ua   = incr%ua
lm%pert%va   = incr%va
lm%pert%t    = incr%t
lm%pert%delp = incr%delp
lm%pert%qv   = incr%q
lm%pert%qi   = incr%qi
lm%pert%ql   = incr%ql
lm%pert%o3   = incr%o3

if (.not. incr%hydrostatic) then
   lm%pert%delz = incr%delz
   lm%pert%w    = incr%w
endif

deallocate(ud,vd)

end subroutine incr_to_lm_tl

! ------------------------------------------------------------------------------

subroutine lm_to_incr_tl(geom,lm,incr)

use wind_vt_mod, only: d2a

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(inout) :: incr
type(fv3jedi_lm_type),   intent(in)    :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ua,va,ud,vd

allocate(ud(incr%isd:incr%ied  ,incr%jsd:incr%jed+1,1:incr%npz))
allocate(vd(incr%isd:incr%ied+1,incr%jsd:incr%jed  ,1:incr%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

allocate(ua(incr%isd:incr%ied  ,incr%jsd:incr%jed  ,1:incr%npz))
allocate(va(incr%isd:incr%ied  ,incr%jsd:incr%jed  ,1:incr%npz))
ua = 0.0_kind_real
va = 0.0_kind_real

ud(incr%isc:incr%iec,incr%jsc:incr%jec,:) = lm%pert%u
vd(incr%isc:incr%iec,incr%jsc:incr%jec,:) = lm%pert%v

call d2a(geom, ud, vd, ua, va)

incr%ua   = ua(incr%isc:incr%iec,incr%jsc:incr%jec,:)
incr%va   = va(incr%isc:incr%iec,incr%jsc:incr%jec,:)
incr%t    = lm%pert%t
incr%delp = lm%pert%delp
incr%q    = lm%pert%qv
incr%qi   = lm%pert%qi
incr%ql   = lm%pert%ql
incr%o3   = lm%pert%o3
if (.not. incr%hydrostatic) then
  incr%delz = lm%pert%delz
  incr%w    = lm%pert%w
endif

deallocate(ud,vd,ua,va)

end subroutine lm_to_incr_tl

! ------------------------------------------------------------------------------

subroutine lm_to_incr_ad(geom,lm,incr)

use wind_vt_mod, only: d2a_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: incr
type(fv3jedi_lm_type),   intent(inout) :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ua,va,ud,vd

allocate(ud(incr%isd:incr%ied  ,incr%jsd:incr%jed+1,1:incr%npz))
allocate(vd(incr%isd:incr%ied+1,incr%jsd:incr%jed  ,1:incr%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

allocate(ua(incr%isd:incr%ied  ,incr%jsd:incr%jed  ,1:incr%npz))
allocate(va(incr%isd:incr%ied  ,incr%jsd:incr%jed  ,1:incr%npz))
ua = 0.0_kind_real
va = 0.0_kind_real

ua(incr%isc:incr%iec,incr%jsc:incr%jec,:) = incr%ua
va(incr%isc:incr%iec,incr%jsc:incr%jec,:) = incr%va

lm%pert%t    = incr%t
lm%pert%delp = incr%delp
lm%pert%qv   = incr%q
lm%pert%qi   = incr%qi
lm%pert%ql   = incr%ql
lm%pert%o3   = incr%o3
if (.not. incr%hydrostatic) then
   lm%pert%delz = incr%delz
   lm%pert%w    = incr%w
endif

call d2a_ad(geom, ud, vd, ua, va)

lm%pert%u = ud(incr%isc:incr%iec,incr%jsc:incr%jec,:)
lm%pert%v = vd(incr%isc:incr%iec,incr%jsc:incr%jec,:)

lm%pert%ua(incr%isc:incr%iec,incr%jsc:incr%jec,:) = 0.0
lm%pert%va(incr%isc:incr%iec,incr%jsc:incr%jec,:) = 0.0

deallocate(ud,vd,ua,va)

end subroutine lm_to_incr_ad

! ------------------------------------------------------------------------------

subroutine incr_to_lm_ad(geom,incr,lm)

use wind_vt_mod, only: a2d_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(inout) :: incr
type(fv3jedi_lm_type),   intent(in)    :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

allocate(ud(incr%isc:incr%iec  ,incr%jsc:incr%jec+1,1:incr%npz))
allocate(vd(incr%isc:incr%iec+1,incr%jsc:incr%jec  ,1:incr%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

incr%ua   = 0.0
incr%va   = 0.0
incr%t    = 0.0
incr%delp = 0.0
incr%q    = 0.0
incr%qi   = 0.0
incr%ql   = 0.0
incr%o3   = 0.0
if (.not. incr%hydrostatic) then
   incr%delz = 0.0
   incr%w    = 0.0
endif

ud(incr%isc:incr%iec,incr%jsc:incr%jec,:) = lm%pert%u
vd(incr%isc:incr%iec,incr%jsc:incr%jec,:) = lm%pert%v
incr%t    = lm%pert%t
incr%delp = lm%pert%delp
incr%q    = lm%pert%qv
incr%qi   = lm%pert%qi
incr%ql   = lm%pert%ql
incr%o3   = lm%pert%o3
if (.not. incr%hydrostatic) then
  incr%delz = lm%pert%delz
  incr%w    = lm%pert%w
endif

!Convert A to D
call a2d_ad(geom, incr%ua, incr%va, ud, vd)

deallocate(ud,vd)

end subroutine incr_to_lm_ad

! ------------------------------------------------------------------------------

subroutine traj_to_traj( traj, lm )

implicit none
type(fv3jedi_traj),    intent(in)    :: traj
type(fv3jedi_lm_type), intent(inout) :: lm

lm%traj%u    = traj%u   
lm%traj%v    = traj%v   
lm%traj%ua   = traj%ua  
lm%traj%va   = traj%va  
lm%traj%t    = traj%t   
lm%traj%delp = traj%delp
lm%traj%qv   = traj%qv  
lm%traj%ql   = traj%ql  
lm%traj%qi   = traj%qi  
lm%traj%o3   = traj%o3  
                                        
if (.not. lm%conf%hydrostatic) then
  lm%traj%w    = traj%w   
  lm%traj%delz = traj%delz
endif
                                        
if (lm%conf%do_phy_mst .ne. 0) then
  lm%traj%qls  = traj%qls 
  lm%traj%qcn  = traj%qcn 
  lm%traj%cfcn = traj%cfcn
endif

!> Rank two
lm%traj%phis    = traj%phis   
lm%traj%frocean = traj%frocean
lm%traj%frland  = traj%frland 
lm%traj%varflt  = traj%varflt 
lm%traj%ustar   = traj%ustar  
lm%traj%bstar   = traj%bstar  
lm%traj%zpbl    = traj%zpbl   
lm%traj%cm      = traj%cm     
lm%traj%ct      = traj%ct     
lm%traj%cq      = traj%cq     
lm%traj%kcbl    = traj%kcbl   
lm%traj%ts      = traj%ts     
lm%traj%khl     = traj%khl    
lm%traj%khu     = traj%khu    

end subroutine traj_to_traj

! ------------------------------------------------------------------------------

end module fv3jedi_tlm_mod
