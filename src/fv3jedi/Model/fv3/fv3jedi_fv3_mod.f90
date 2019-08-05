! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_io_gfs_mod, only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos
use fv3jedi_increment_mod, only: fv3jedi_increment

use fv3jedi_lm_mod, only: fv3jedi_lm_type

implicit none
private

public :: fv3_model
public :: fv3_create
public :: fv3_delete
public :: fv3_initialize
public :: fv3_step
public :: fv3_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: fv3_model
  type(fv3jedi_lm_type)                        :: fv3jedi_lm          !<Linearized model object
end type fv3_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3_create(self, geom, c_conf)

implicit none
type(c_ptr),         intent(in)  :: c_conf
type(fv3_model), intent(inout)   :: self
type(fv3jedi_geom),  intent(in)  :: geom

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt

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


! Safety checks
! -------------

!The full trajecotory of the tlm/adm is not output by this simplified model
!so if being used to generate the trajectry with physics the traj must be read
!from file or obtained by running GEOS or GFS.
if ((self%fv3jedi_lm%conf%do_phy_trb .ne. 0) .or. &
    (self%fv3jedi_lm%conf%do_phy_mst .ne. 0) ) then
   call abor1_ftn("fv3_model | FV3 : unless reading the trajecotory physics should be off")
endif

end subroutine fv3_create

! ------------------------------------------------------------------------------

subroutine fv3_delete(self)

implicit none
type(fv3_model), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine fv3_delete

! ------------------------------------------------------------------------------

subroutine fv3_initialize(self, state)

implicit none
type(fv3_model), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

call self%fv3jedi_lm%init_nl()

end subroutine fv3_initialize

! ------------------------------------------------------------------------------

subroutine fv3_step(self, state, geom, vdate)

implicit none
type(fv3_model),     intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_geom),  intent(inout) :: geom
type(datetime),      intent(inout) :: vdate !< Valid datetime after step

call state_to_lm(state,self%fv3jedi_lm)
call self%fv3jedi_lm%step_nl()
call lm_to_state(self%fv3jedi_lm,state)

end subroutine fv3_step

! ------------------------------------------------------------------------------

subroutine fv3_finalize(self, state)

implicit none
type(fv3_model), target :: self
type(fv3jedi_state)         :: state

call self%fv3jedi_lm%final_nl()

end subroutine fv3_finalize

! ------------------------------------------------------------------------------

subroutine state_to_lm( state, lm )

implicit none
type(fv3jedi_state),   intent(in)    :: state
type(fv3jedi_lm_type), intent(inout) :: lm

lm%traj%ua = 0.0_kind_real
lm%traj%va = 0.0_kind_real

lm%traj%u       = state%ud(state%isc:state%iec,state%jsc:state%jec,:)
lm%traj%v       = state%vd(state%isc:state%iec,state%jsc:state%jec,:)
if (associated(state%ua)) lm%traj%ua = state%ua
if (associated(state%va)) lm%traj%va = state%va
lm%traj%t       = state%t
lm%traj%delp    = state%delp
lm%traj%qv      = state%q
lm%traj%qi      = state%qi
lm%traj%ql      = state%ql
lm%traj%o3      = state%o3

if (.not. lm%conf%hydrostatic) then
lm%traj%w       = state%w
lm%traj%delz    = state%delz
endif

lm%traj%phis = state%phis(:,:,1)

end subroutine state_to_lm

! ------------------------------------------------------------------------------

subroutine lm_to_state( lm, state )

implicit none
type(fv3jedi_lm_type), intent(in)    :: lm
type(fv3jedi_state),   intent(inout) :: state

state%ud(state%isc:state%iec,state%jsc:state%jec,:)      = lm%traj%u
state%vd(state%isc:state%iec,state%jsc:state%jec,:)      = lm%traj%v
if (associated(state%ua)) state%ua = lm%traj%ua
if (associated(state%ua)) state%va = lm%traj%va
state%t       = lm%traj%t
state%delp    = lm%traj%delp
state%q       = lm%traj%qv
state%qi      = lm%traj%qi
state%ql      = lm%traj%ql
state%o3      = lm%traj%o3

if (.not. lm%conf%hydrostatic) then
state%w       = lm%traj%w
state%delz    = lm%traj%delz
endif

state%phis(:,:,1)    = lm%traj%phis

end subroutine lm_to_state

! ------------------------------------------------------------------------------

end module fv3jedi_fv3_mod
