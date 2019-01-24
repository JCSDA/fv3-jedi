! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_tlm_mod

use iso_c_binding
use config_mod
use duration_mod

use fv3jedi_kinds_mod
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
  logical :: fsoi_mode
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
integer :: tmp

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

self%fsoi_mode = .false.
if (config_element_exists(c_conf,"fsoi_mode")) then
   tmp = config_get_int(c_conf,"fsoi_mode")
   if (tmp==1) self%fsoi_mode = .true.
endif

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

subroutine tlm_initialize_ad(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%init_ad()

end subroutine tlm_initialize_ad

! ------------------------------------------------------------------------------

subroutine tlm_initialize_tl(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

if (self%fsoi_mode) call inc_to_lm_tl(geom,inc,self%fv3jedi_lm)

call self%fv3jedi_lm%init_tl()

end subroutine tlm_initialize_tl

! ------------------------------------------------------------------------------

subroutine tlm_step_ad(geom, self, inc, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

if (.not. self%fsoi_mode) call lm_to_inc_ad(geom,self%fv3jedi_lm,inc)
call self%fv3jedi_lm%step_ad()
call inc_to_lm_ad(geom,inc,self%fv3jedi_lm)

end subroutine tlm_step_ad

! ------------------------------------------------------------------------------

subroutine tlm_step_tl(geom, self, inc, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

if (.not. self%fsoi_mode) call inc_to_lm_tl(geom,inc,self%fv3jedi_lm)
call self%fv3jedi_lm%step_tl()
call lm_to_inc_tl(geom,self%fv3jedi_lm,inc)

end subroutine tlm_step_tl

! ------------------------------------------------------------------------------

subroutine tlm_finalize_ad(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%final_ad()
if (.not. self%fsoi_mode) call lm_to_inc_ad(geom,self%fv3jedi_lm,inc)

end subroutine tlm_finalize_ad

! ------------------------------------------------------------------------------

subroutine tlm_finalize_tl(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: inc

call self%fv3jedi_lm%final_tl()

end subroutine tlm_finalize_tl

! ------------------------------------------------------------------------------

subroutine inc_to_lm_tl(geom, inc, lm)

use wind_vt_mod, only: a2d

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: k
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

allocate(ud(inc%isc:inc%iec  ,inc%jsc:inc%jec+1,1:inc%npz))
allocate(vd(inc%isc:inc%iec+1,inc%jsc:inc%jec  ,1:inc%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

call a2d(geom, inc%fields(inc%ua)%field(inc%isc:inc%iec,inc%jsc:inc%jec,:), &
               inc%fields(inc%va)%field(inc%isc:inc%iec,inc%jsc:inc%jec,:), &
                    ud(inc%isc:inc%iec  ,inc%jsc:inc%jec+1,:), &
                    vd(inc%isc:inc%iec+1,inc%jsc:inc%jec  ,:))

lm%pert%u    = ud(inc%isc:inc%iec,inc%jsc:inc%jec,:)
lm%pert%v    = vd(inc%isc:inc%iec,inc%jsc:inc%jec,:)
lm%pert%ua   = inc%fields(inc%ua)%field
lm%pert%va   = inc%fields(inc%va)%field
lm%pert%t    = inc%fields(inc%t)%field
do k = 1,geom%npz
  lm%pert%delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*inc%fields(inc%ps)%field(:,:,1)
enddo
lm%pert%qv   = inc%fields(inc%q)%field
lm%pert%qi   = inc%fields(inc%qi)%field
lm%pert%ql   = inc%fields(inc%ql)%field
lm%pert%o3   = inc%fields(inc%o3)%field

if (.not. inc%hydrostatic) then
   lm%pert%delz = inc%fields(inc%delz)%field
   lm%pert%w    = inc%fields(inc%w)%field
endif

deallocate(ud,vd)

end subroutine inc_to_lm_tl

! ------------------------------------------------------------------------------

subroutine lm_to_inc_tl(geom,lm,inc)

use wind_vt_mod, only: d2a

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_lm_type),   intent(in)    :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ua,va,ud,vd

allocate(ud(inc%isd:inc%ied  ,inc%jsd:inc%jed+1,1:inc%npz))
allocate(vd(inc%isd:inc%ied+1,inc%jsd:inc%jed  ,1:inc%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

allocate(ua(inc%isd:inc%ied  ,inc%jsd:inc%jed  ,1:inc%npz))
allocate(va(inc%isd:inc%ied  ,inc%jsd:inc%jed  ,1:inc%npz))
ua = 0.0_kind_real
va = 0.0_kind_real

ud(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%u
vd(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%v

call d2a(geom, ud, vd, ua, va)

inc%fields(inc%ua)%field   = ua(inc%isc:inc%iec,inc%jsc:inc%jec,:)
inc%fields(inc%va)%field   = va(inc%isc:inc%iec,inc%jsc:inc%jec,:)
inc%fields(inc%t)%field    = lm%pert%t
inc%fields(inc%ps)%field(:,:,1)   = sum(lm%pert%delp,3)
inc%fields(inc%q)%field    = lm%pert%qv
inc%fields(inc%qi)%field   = lm%pert%qi
inc%fields(inc%ql)%field   = lm%pert%ql
inc%fields(inc%o3)%field   = lm%pert%o3
if (.not. inc%hydrostatic) then
  inc%fields(inc%delz)%field = lm%pert%delz
  inc%fields(inc%w)%field    = lm%pert%w
endif

deallocate(ud,vd,ua,va)

end subroutine lm_to_inc_tl

! ------------------------------------------------------------------------------

subroutine lm_to_inc_ad(geom,lm,inc)

use wind_vt_mod, only: d2a_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: k
real(kind=kind_real), allocatable, dimension(:,:,:) :: ua,va,ud,vd

allocate(ud(inc%isd:inc%ied  ,inc%jsd:inc%jed+1,1:inc%npz))
allocate(vd(inc%isd:inc%ied+1,inc%jsd:inc%jed  ,1:inc%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

allocate(ua(inc%isd:inc%ied  ,inc%jsd:inc%jed  ,1:inc%npz))
allocate(va(inc%isd:inc%ied  ,inc%jsd:inc%jed  ,1:inc%npz))
ua = 0.0_kind_real
va = 0.0_kind_real

ua(inc%isc:inc%iec,inc%jsc:inc%jec,:) = inc%fields(inc%ua)%field
va(inc%isc:inc%iec,inc%jsc:inc%jec,:) = inc%fields(inc%va)%field

lm%pert%t    = inc%fields(inc%t)%field
do k = 1,geom%npz
  lm%pert%delp(:,:,k) = inc%fields(inc%ps)%field(:,:,1)
enddo
lm%pert%qv   = inc%fields(inc%q)%field
lm%pert%qi   = inc%fields(inc%qi)%field
lm%pert%ql   = inc%fields(inc%ql)%field
lm%pert%o3   = inc%fields(inc%o3)%field
if (.not. inc%hydrostatic) then
   lm%pert%delz = inc%fields(inc%delz)%field
   lm%pert%w    = inc%fields(inc%w)%field
endif

call d2a_ad(geom, ud, vd, ua, va)

lm%pert%u = ud(inc%isc:inc%iec,inc%jsc:inc%jec,:)
lm%pert%v = vd(inc%isc:inc%iec,inc%jsc:inc%jec,:)

lm%pert%ua(inc%isc:inc%iec,inc%jsc:inc%jec,:) = 0.0
lm%pert%va(inc%isc:inc%iec,inc%jsc:inc%jec,:) = 0.0

deallocate(ud,vd,ua,va)

end subroutine lm_to_inc_ad

! ------------------------------------------------------------------------------

subroutine inc_to_lm_ad(geom,inc,lm)

use wind_vt_mod, only: a2d_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_lm_type),   intent(in)    :: lm

integer :: k
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

allocate(ud(inc%isc:inc%iec  ,inc%jsc:inc%jec+1,1:inc%npz))
allocate(vd(inc%isc:inc%iec+1,inc%jsc:inc%jec  ,1:inc%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

inc%fields(inc%ua)%field   = 0.0
inc%fields(inc%va)%field   = 0.0
inc%fields(inc%t )%field   = 0.0
inc%fields(inc%ps)%field   = 0.0
inc%fields(inc%q )%field   = 0.0
inc%fields(inc%qi)%field   = 0.0
inc%fields(inc%ql)%field   = 0.0
inc%fields(inc%o3)%field   = 0.0

if (.not. inc%hydrostatic) then
   inc%fields(inc%delz)%field = 0.0
   inc%fields(inc%w)%field    = 0.0
endif

ud(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%u
vd(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%v
inc%fields(inc%t)%field    = lm%pert%t
inc%fields(inc%ps)%field = 0.0_kind_real
do k = 1,geom%npz
  inc%fields(inc%ps)%field(:,:,1) = inc%fields(inc%ps)%field(:,:,1) +  (geom%bk(k+1)-geom%bk(k))*lm%pert%delp(:,:,k)
enddo
inc%fields(inc%q)%field    = lm%pert%qv
inc%fields(inc%qi)%field   = lm%pert%qi
inc%fields(inc%ql)%field   = lm%pert%ql
inc%fields(inc%o3)%field   = lm%pert%o3
if (.not. inc%hydrostatic) then
  inc%fields(inc%delz)%field = lm%pert%delz
  inc%fields(inc%q)%field    = lm%pert%w
endif

!Convert A to D
call a2d_ad(geom, inc%fields(inc%ua)%field, inc%fields(inc%va)%field, ud, vd)

deallocate(ud,vd)

end subroutine inc_to_lm_ad

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
