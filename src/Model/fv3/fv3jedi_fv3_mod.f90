! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_fv3_mod

use iso_c_binding
use config_mod
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
  integer                                      :: readtraj            !<Read trajectory from file
  character(len=255)                           :: trajmodel           !<User specified model type for traj
  character(len=255)                           :: trajpath            !<User specified path to traj files
  character(len=255)                           :: trajfile            !<User specified path to traj files
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


! Model time step
! ---------------
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)


! Option to read traj from file instead of propagating model
! ----------------------------------------------------------
self%readtraj = config_get_int(c_conf,"readtraj")
if (self%readtraj == 1) then
  self%trajmodel = config_get_string(c_conf,len(self%trajmodel),"trajmodel")
  self%trajpath  = config_get_string(c_conf,len(self%trajpath ),"trajpath")
  self%trajfile  = config_get_string(c_conf,len(self%trajfile ),"trajfile")
endif


! Model configuration and creation
! --------------------------------
self%fv3jedi_lm%conf%do_dyn     = config_get_int(c_conf,"lm_do_dyn")
self%fv3jedi_lm%conf%do_phy_trb = config_get_int(c_conf,"lm_do_trb")
self%fv3jedi_lm%conf%do_phy_mst = config_get_int(c_conf,"lm_do_mst")

call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)


! Safety checks
! -------------

!The full trajecotory of the tlm/adm is not output by this simplified model
!so if being used to generate the trajectry with physics the traj must be read
!from file or obtained by running GEOS or GFS. 
if ((self%fv3jedi_lm%conf%do_phy_trb .ne. 0 .and. self%readtraj == 0) .or. &
    (self%fv3jedi_lm%conf%do_phy_mst .ne. 0 .and. self%readtraj == 0) ) then
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

if (self%readtraj == 0) then

  call state_to_lm(state,self%fv3jedi_lm)
  call self%fv3jedi_lm%step_nl()
  call lm_to_state(self%fv3jedi_lm,state)

else

  call psuedo_model( geom, self, state, vdate)

endif

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
 
lm%traj%u       = state%ud(state%isc:state%iec,state%jsc:state%jec,:)
lm%traj%v       = state%vd(state%isc:state%iec,state%jsc:state%jec,:)
lm%traj%ua      = state%ua
lm%traj%va      = state%va
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
state%ua      = lm%traj%ua
state%va      = lm%traj%va
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

subroutine psuedo_model( geom, self, state, vdate)

implicit none
type(fv3jedi_geom),  intent(inout) :: geom
type(fv3_model),     intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(inout) :: vdate !< Valid datetime after step

character(len=20)  :: vdatec
character(len=255) :: date, filename
character(len=4)   :: yyyy
character(len=2)   :: mm,dd,hh,mn,ss

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos

! Convert datetime to string
call datetime_to_string(vdate, vdatec)

! Write character form of date
write(date,*) vdatec(1:4),vdatec(6:7),vdatec(9:10),'_',vdatec(12:13),vdatec(15:16),'z.nc4'

! Date part of the filename
yyyy = vdatec(1 :4 )
mm   = vdatec(6 :7 )
dd   = vdatec(9 :10)
hh   = vdatec(12:13)
mn   = vdatec(15:16)
ss   = vdatec(18:19)

! File path/filename
if (trim(self%trajmodel) == "gfs") then

  gfs%datapath_ti = trim(self%trajpath)
  gfs%filename_core = yyyy//mm//dd//"."//hh//mn//ss//'.fv_core.res.nc'
  gfs%filename_trcr = yyyy//mm//dd//"."//hh//mn//ss//'.fv_tracer.res.nc'
  gfs%filename_sfcd = yyyy//mm//dd//"."//hh//mn//ss//'.sfc_data.nc'
  gfs%filename_sfcw = yyyy//mm//dd//"."//hh//mn//ss//'.srf_wnd.nc'
  gfs%filename_cplr = yyyy//mm//dd//"."//hh//mn//ss//'.coupler.res'
  gfs%datapath_sp = 'null'
  gfs%datapath_sp = 'null'

  call print_filename(self,gfs%filename_core)
  call gfs%read_meta  ( geom, vdate, state%calendar_type, state%date_init )
  call gfs%read_fields( geom, state%fields )

elseif (trim(self%trajmodel) == "geos") then

  filename = trim(self%trajpath)//trim(self%trajfile)//trim(yyyy)//trim(mm)//trim(dd)//"_"//trim(hh)//trim(mn)//'z.nc4'
  call print_filename(self,filename)
  call geos%create(geom, 'read', filename)
  call geos%read_fields(geom, state%fields)
  call geos%delete()

else

  call abor1_ftn("fv3jedi_fv3_mod: psuedo_model, model choice must be geos or gfs")

endif

end subroutine psuedo_model

! ------------------------------------------------------------------------------

subroutine print_filename(self,filename)

implicit none
type(fv3_model),  intent(in) :: self
character(len=*), intent(in) :: filename

! Print filename to the user
if (self%fv3jedi_lm%conf%rpe) print*, ' '
if (self%fv3jedi_lm%conf%rpe) print*, 'Psuedo model from file: ', trim(filename)
if (self%fv3jedi_lm%conf%rpe) print*, ' '

end subroutine print_filename

! ------------------------------------------------------------------------------

end module fv3jedi_fv3_mod
