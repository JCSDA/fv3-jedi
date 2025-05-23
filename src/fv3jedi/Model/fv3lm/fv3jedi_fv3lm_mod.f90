! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3lm_mod

! oops uses
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! Linear model
use fv3jedi_lm_mod,        only: fv3jedi_lm_type

! fv3-jedi uses
use fv3jedi_fmsnamelist_mod, only: fv3jedi_fmsnamelist
use fv3jedi_io_utils_mod,    only: vdate_to_datestring
use fv3jedi_kinds_mod,       only: kind_real
use fv3jedi_geom_mod,        only: fv3jedi_geom
use fv3jedi_state_mod,       only: fv3jedi_state

implicit none
private
public :: fv3lm_model

! --------------------------------------------------------------------------------------------------

type :: fv3lm_model
  type(fv3jedi_lm_type) :: fv3jedi_lm       !<Linearized model object
  logical :: model_has_been_initialized     !<Flag to indicate if the model has been initialized
  character(len=1024) :: datapath           !<Path to the restart file containing the D-Grid winds
  character(len=1024) :: filename_core      !<Filename for the D-Grid wind restart
  character(len=1024) :: datapath_out       !<Path to save the D-Grid winds for the next cycle
  character(len=1024) :: filename_core_out  !<Filename for the D-Grid wind restart for the next cyc
  logical :: a_to_d_on_init                 !<Flag to indicate if A -> D is used instead of files
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: initialize
    procedure, public :: step
    procedure, public :: finalize
end type fv3lm_model

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3lm_model),        intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt
character(len=:), allocatable :: str
type(fv3jedi_fmsnamelist) :: fmsnamelist

! Model time step
! ---------------
call conf%get_or_die("tstep",str)
ststep = str
deallocate(str)

dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)

! Config needs to provide the path to the FMS restart contining the D-Grid winds
! or specify that the internal D-Grid winds should be derived from A-Grid winds
! ------------------------------------------------------------------------------
if (conf%has("initialize model from A-Grid winds")) then
  ! For toy model cycling tests the A-Grid winds can be used to initialize the D-Grid winds
  ! this should never be used for a real forecasting system and is unphysical.
  call conf%get_or_die("initialize model from A-Grid winds", self%a_to_d_on_init)
endif

! In the normal case a_to_d_on_init is false the the user must provide paths for D-Grid wind
! input and output.
if (.not. self%a_to_d_on_init) then

  ! Path and file for input D-Grid wind files
  call conf%get_or_die("datapath", str)
  self%datapath = str
  deallocate(str)
  call conf%get_or_die("filename_core", str)
  self%filename_core = str
  deallocate(str)

  ! Path and file for saving D-Grid wind files for the next cycle
  self%datapath_out = self%datapath
  if (conf%has("datapath_out")) then
    call conf%get_or_die("datapath_out", str)
    self%datapath_out = str
    deallocate(str)
  end if
  call conf%get_or_die("filename_core_out", str)
  self%filename_core_out = str
  deallocate(str)

end if


! Model configuration and creation
! --------------------------------
call conf%get_or_die("lm_do_dyn",self%fv3jedi_lm%conf%do_dyn)
call conf%get_or_die("lm_do_trb",self%fv3jedi_lm%conf%do_phy_trb)
call conf%get_or_die("lm_do_mst",self%fv3jedi_lm%conf%do_phy_mst)

! Prepare namelist (nonlinear part)
! ---------------------------------
call fmsnamelist%replace_namelist(conf)

! Call model constructor
! ----------------------
call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

! Safety checks
! -------------

!The full trajecotory of the tlm/adm is not output by this simplified model
!so if being used to generate the trajectry with physics the traj must be read
!from file or obtained by running GEOS or GFS.
if ((self%fv3jedi_lm%conf%do_phy_trb .ne. 0) .or. &
    (self%fv3jedi_lm%conf%do_phy_mst .ne. 0) ) then
   call abor1_ftn("fv3lm_model | FV3LM : unless reading the trajectory physics should be off")
endif

! Set the flag to indicate that the model has not initialized
self%model_has_been_initialized = .false.

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3lm_model), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state)

class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

! Wind pointers
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()


! Make sure the tracers are allocated (false => trajectory only)
call self%fv3jedi_lm%allocate_tracers(state%ntracers, .false.)

! Copy the parts of the state that JEDI has into the state
call state_to_lm(state,self%fv3jedi_lm)

if (.not. self%model_has_been_initialized) then
  ! If this is the first time the model has been initialized then then initialize the internal grid
  ! staggered (D-Grid) winds of the model by reading the model restart
  if (.not. self%a_to_d_on_init) then
    call self%fv3jedi_lm%read_d_grid_winds(self%datapath, self%filename_core)
  else
    ! Option to intialize the model's internal D-Grid winds from the JEDI A-Grid winds.
    !!! THIS SHOULD NOT BE REPLICATED FOR A FULL FORECAST SYSTEM AS IT IS NOT A PHYSICALLY      !!!
    !!! CONSISTENT WAY TO INITIALIZE THE MODEL. IT IS ONLY USED FOR TOY MODEL TESTING PURPOSES. !!!
    call state%get_field('eastward_wind', ua)
    call state%get_field('northward_wind', va)
    call self%fv3jedi_lm%initialize_dwinds_from_awinds(ua, va)
  end if

  ! In case of outer loops, a copy of the winds (D and A) should be stored internally
  call self%fv3jedi_lm%store_winds()
else
  ! If the model has been initialized before then the internal D-Grid winds need to be reverted.
  ! If JEDI has changed the A-Grid winds then that contribution needs to be added to the reverted
  ! D-Grid winds
  call self%fv3jedi_lm%reinitialize_winds()

  ! Replace the stored wind values
  call self%fv3jedi_lm%store_winds()
endif

! Initialize the model
call self%fv3jedi_lm%init_nl()

! Set the flag to indicate that the model has been initialized
self%model_has_been_initialized = .true.

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, geom)

class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_geom),  intent(inout) :: geom

! Strings for datetime
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss
character(len=1024) :: filename_core_out

call self%fv3jedi_lm%step_nl()
call lm_to_state(self%fv3jedi_lm,state)

! Write out the D-Grid winds to a restart file
! --------------------------------------------
if (.not. self%a_to_d_on_init) then
  call vdate_to_datestring(state%time, yyyy=yyyy, mm=mm, dd=dd, hh=hh, min=min, ss=ss)
  filename_core_out = yyyy//mm//dd//'_'//hh//min//ss//'.'//trim(self%filename_core_out)
  call self%fv3jedi_lm%write_d_grid_winds(self%datapath_out, filename_core_out)
endif

end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine finalize(self, state)

class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state

call self%fv3jedi_lm%final_nl()

end subroutine finalize

! --------------------------------------------------------------------------------------------------

subroutine state_to_lm( state, lm )

! Arguments
type(fv3jedi_state),   intent(in)    :: state
type(fv3jedi_lm_type), intent(inout) :: lm

! Locals
integer :: ft, f, index
logical :: sphum_found = .false.
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: t  => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delp => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: q  => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: w  => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delz => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: phis => null()

! Assert that the state only has one flavor of pressure and temperature
! ---------------------------------------------------------------------
if (state%has_field('air_pressure_at_surface') .or. state%has_field('air_pressure') .or. &
    state%has_field('air_pressure_levels')) then
  call abor1_ftn("state_to_lm: When working in-core delp must be in the state and " // &
                 "other types of pressure must not be present. Otherwise it leads to " // &
                 "ambiguity in how increments are applied back to the model.")
end if

if (state%has_field('virtual_temperature') .or. state%has_field('pt')) then
  call abor1_ftn("state_to_lm: When working in-core temperature (t) must be in the state and " // &
                 "other types of temperature must not be present. Otherwise it leads to " // &
                 "ambiguity in how increments are applied back to the model.")
end if

! Required variables
! ------------------
call state%get_field('eastward_wind', ua)
call state%get_field('northward_wind', va)
call state%get_field('air_temperature', t)
call state%get_field('air_pressure_thickness', delp)
call state%get_field('geopotential_height_times_gravity_at_surface', phis)

lm%traj%ua   = ua
lm%traj%va   = va
lm%traj%t    = t
lm%traj%delp = delp
lm%traj%phis = phis(:,:,1)

! Tracer variables
! ----------------
ft = 1
do f = 1, state%nf
  if (state%fields(f)%tracer) then
    if (trim(state%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    lm%traj%tracers(:,:,:,index) = state%fields(f)%array
    lm%traj%tracer_names(index)  = trim(state%fields(f)%long_name)
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("state_to_lm: water_vapor_mixing_ratio_wrt_moist_air (sphum) not in tracer list")
end if

! Variables when non-hydrostatic
! ------------------------------
if (.not. lm%conf%hydrostatic) then
  call state%get_field('upward_air_velocity', w   )
  call state%get_field('layer_thickness', delz)
  lm%traj%w       = w
  lm%traj%delz    = delz
endif

end subroutine state_to_lm

! --------------------------------------------------------------------------------------------------

subroutine lm_to_state( lm, state )

! Arguments
type(fv3jedi_lm_type), intent(in)    :: lm
type(fv3jedi_state),   intent(inout) :: state

! Locals
integer :: ft, f, index
logical :: sphum_found = .false.
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: t => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delp => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: q => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: w => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delz => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: phis => null()

! Required variables
! ------------------
call state%get_field('eastward_wind', ua  )
call state%get_field('northward_wind', va  )
call state%get_field('air_temperature', t   )
call state%get_field('air_pressure_thickness', delp)
call state%get_field('geopotential_height_times_gravity_at_surface', phis)
ua          = lm%traj%ua
va          = lm%traj%va
t           = lm%traj%t
delp        = lm%traj%delp
phis(:,:,1) = lm%traj%phis

! Tracers
! -------
ft = 1
do f = 1, state%nf
  if (state%fields(f)%tracer) then
    if (trim(state%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    state%fields(f)%array = lm%traj%tracers(:,:,:,index)
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("lm_to_state: water_vapor_mixing_ratio_wrt_moist_air (sphum) not in tracer list")
end if

! Non-hydrostatic variables
! -------------------------
if (.not. lm%conf%hydrostatic) then
  call state%get_field('upward_air_velocity', w)
  call state%get_field('layer_thickness', delz)
  w       = lm%traj%w
  delz    = lm%traj%delz
endif

end subroutine lm_to_state

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fv3lm_mod
