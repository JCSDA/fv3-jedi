! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_tlm_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use duration_mod
use variables_mod

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment, create_inc=>create
use fv3jedi_traj_mod, only: fv3jedi_traj

use fv3jedi_lm_mod, only: fv3jedi_lm_type

use fv3jedi_linvarcha_a2m_mod

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
  type(fv3jedi_linvarcha_a2m) :: lvc
  type(fv3jedi_increment) :: xmod
end type fv3jedi_tlm

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine tlm_create(self, geom, c_conf, vars)

implicit none
type(fv3jedi_tlm),  intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom
type(c_ptr),        intent(in)    :: c_conf
type(oops_vars),    intent(in)    :: vars

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt
integer :: tmp

type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)

! Model time step
! ---------------
call f_conf%get_or_die("tstep",str)
ststep = str
deallocate(str)
dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)


! Model configuration and creation
! --------------------------------
call f_conf%get_or_die("lm_do_dyn",self%fv3jedi_lm%conf%do_dyn)
call f_conf%get_or_die("lm_do_trb",self%fv3jedi_lm%conf%do_phy_trb)
call f_conf%get_or_die("lm_do_mst",self%fv3jedi_lm%conf%do_phy_mst)

call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

call create_inc(self%xmod, geom, vars)

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
type(fv3jedi_increment), intent(inout) :: inc

call multiplyinverseadjoint(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyadjoint(self%lvc,geom,self%xmod,inc)

end subroutine tlm_initialize_ad

! ------------------------------------------------------------------------------

subroutine tlm_initialize_tl(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call multiply(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyinverse(self%lvc,geom,self%xmod,inc)

end subroutine tlm_initialize_tl

! ------------------------------------------------------------------------------

subroutine tlm_step_ad(geom, self, inc, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call multiplyinverseadjoint(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%step_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyadjoint(self%lvc,geom,self%xmod,inc)

end subroutine tlm_step_ad

! ------------------------------------------------------------------------------

subroutine tlm_step_tl(geom, self, inc, traj)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call multiply(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%step_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyinverse(self%lvc,geom,self%xmod,inc)

end subroutine tlm_step_tl

! ------------------------------------------------------------------------------

subroutine tlm_finalize_ad(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call multiplyinverseadjoint(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%final_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyadjoint(self%lvc,geom,self%xmod,inc)

end subroutine tlm_finalize_ad

! ------------------------------------------------------------------------------

subroutine tlm_finalize_tl(geom, self, inc)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),       intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call multiply(self%lvc,geom,inc,self%xmod)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%final_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

call multiplyinverse(self%lvc,geom,self%xmod,inc)

end subroutine tlm_finalize_tl

! ------------------------------------------------------------------------------

subroutine inc_to_lm(inc, lm)

implicit none
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

lm%pert%u    = inc%ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
lm%pert%v    = inc%vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
lm%pert%ua   = 0.0_kind_real
lm%pert%va   = 0.0_kind_real
lm%pert%t    = inc%t
lm%pert%delp = inc%delp
lm%pert%qv   = inc%q
lm%pert%qi   = inc%qi
lm%pert%ql   = inc%ql
lm%pert%o3   = inc%o3
if (.not. inc%hydrostatic) then
   lm%pert%delz = inc%delz
   lm%pert%w    = inc%w
endif

end subroutine inc_to_lm

! ------------------------------------------------------------------------------

subroutine lm_to_inc(lm, inc)

implicit none
type(fv3jedi_lm_type),   intent(in)    :: lm
type(fv3jedi_increment), intent(inout) :: inc

inc%ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%u
inc%vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%v
inc%t    = lm%pert%t
inc%delp = lm%pert%delp
inc%q    = lm%pert%qv
inc%qi   = lm%pert%qi
inc%ql   = lm%pert%ql
inc%o3   = lm%pert%o3
if (.not. inc%hydrostatic) then
  inc%delz = lm%pert%delz
  inc%w    = lm%pert%w
endif

end subroutine lm_to_inc

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
