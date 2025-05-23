! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_traj_mod

! fv3-jedi
use fv3jedi_kinds_mod,    only: kind_real
use fv3jedi_state_mod,    only: fv3jedi_state

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: fv3jedi_traj, set, wipe

! --------------------------------------------------------------------------------------------------

type :: fv3jedi_traj
  ! Very similar to the fv3-jedi-lm trajectory (model state), except without D-Grid winds.
  ! Note that D-Grid winds should not be interpolated directly and the traj is almost always
  ! interpolated so there is no benefit to having D-Grid winds here and doing so only serves to
  ! complicate the code.
  integer :: isc,iec,jsc,jec,npz,ntracers
  real(kind_real),     allocatable, dimension(:,:,:)   :: ua, va, t, delp
  real(kind_real),     allocatable, dimension(:,:,:,:) :: tracers
  character(len=2048), allocatable, dimension(:)       :: tracer_names
  real(kind_real),     allocatable, dimension(:,:,:)   :: w, delz
  real(kind_real),     allocatable, dimension(:,:,:)   :: cfcn
  real(kind_real),     allocatable, dimension(:,:,:)   :: qls, qcn
  real(kind_real),     allocatable, dimension(:,:)     :: phis
  real(kind_real),     allocatable, dimension(:,:)     :: frocean, frland
  real(kind_real),     allocatable, dimension(:,:)     :: varflt, ustar, bstar
  real(kind_real),     allocatable, dimension(:,:)     :: zpbl, cm, ct, cq
  real(kind_real),     allocatable, dimension(:,:)     :: kcbl, ts, khl, khu
end type fv3jedi_traj

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine set(traj, state)

type(fv3jedi_traj),  intent(inout) :: traj
type(fv3jedi_state), intent(in)    :: state

integer :: isc,iec,jsc,jec,npz,ft,f,index,number_tracers
logical :: sphum_found = .false.

real(kind=kind_real), pointer, dimension(:,:,:) :: phis => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: frocean => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: frland => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: varflt => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: ustar => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: bstar => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: zpbl => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: cm => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: ct => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: cq => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: kcbl => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: tsm => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: khl => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: khu => null()

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
npz = state%npz

traj%ntracers = state%ntracers

! Allocate traj
allocate(traj%ua     (isc:iec, jsc:jec, npz))
allocate(traj%va     (isc:iec, jsc:jec, npz))
allocate(traj%t      (isc:iec, jsc:jec, npz))
allocate(traj%delp   (isc:iec, jsc:jec, npz))
allocate(traj%w      (isc:iec, jsc:jec, npz))
allocate(traj%delz   (isc:iec, jsc:jec, npz))
allocate(traj%qls    (isc:iec, jsc:jec, npz))
allocate(traj%qcn    (isc:iec, jsc:jec, npz))
allocate(traj%cfcn   (isc:iec, jsc:jec, npz))
allocate(traj%phis   (isc:iec, jsc:jec))
allocate(traj%frocean(isc:iec, jsc:jec))
allocate(traj%frland (isc:iec, jsc:jec))
allocate(traj%varflt (isc:iec, jsc:jec))
allocate(traj%ustar  (isc:iec, jsc:jec))
allocate(traj%bstar  (isc:iec, jsc:jec))
allocate(traj%zpbl   (isc:iec, jsc:jec))
allocate(traj%cm     (isc:iec, jsc:jec))
allocate(traj%ct     (isc:iec, jsc:jec))
allocate(traj%cq     (isc:iec, jsc:jec))
allocate(traj%kcbl   (isc:iec, jsc:jec))
allocate(traj%ts     (isc:iec, jsc:jec))
allocate(traj%khl    (isc:iec, jsc:jec))
allocate(traj%khu    (isc:iec, jsc:jec))

! Initialize all to zero incase not in state
! ------------------------------------------
traj%ua      = 0.0_kind_real
traj%va      = 0.0_kind_real
traj%t       = 0.0_kind_real
traj%delp    = 0.0_kind_real
traj%w       = 0.0_kind_real
traj%delz    = 0.0_kind_real
traj%qls     = 0.0_kind_real
traj%qcn     = 0.0_kind_real
traj%cfcn    = 0.0_kind_real
traj%phis    = 0.0_kind_real
traj%frocean = 0.0_kind_real
traj%frland  = 0.0_kind_real
traj%varflt  = 0.0_kind_real
traj%ustar   = 0.0_kind_real
traj%bstar   = 0.0_kind_real
traj%zpbl    = 0.0_kind_real
traj%cm      = 0.0_kind_real
traj%ct      = 0.0_kind_real
traj%cq      = 0.0_kind_real
traj%kcbl    = 0.0_kind_real
traj%ts      = 0.0_kind_real
traj%khl     = 0.0_kind_real
traj%khu     = 0.0_kind_real


! Copy mandatory parts of the trajecotry
! --------------------------------------
call state%get_field('eastward_wind', traj%ua)
call state%get_field('northward_wind', traj%va)
call state%get_field('air_temperature', traj%t)
call state%get_field('air_pressure_thickness', traj%delp)


! Allocate all the tracers that will be part of what gets advected
! ----------------------------------------------------------------

! First allocate the number of tracers
number_tracers = state%ntracers

! Remove qls, qcn, cfcn
if (state%has_field('initial_mass_fraction_of_large_scale_cloud_condensate')) &
  number_tracers = number_tracers - 1
if (state%has_field('initial_mass_fraction_of_convective_cloud_condensate')) &
  number_tracers = number_tracers - 1
if (state%has_field('convective_cloud_area_fraction')) &
  number_tracers = number_tracers - 1

! Allocate the tracers

if (.not. allocated(traj%tracers) .or. size(traj%tracers, 4) .ne. number_tracers) then
  if (allocated(traj%tracers)) deallocate(traj%tracers)
  allocate(traj%tracers(state%isc:state%iec, state%jsc:state%jec, state%npz, number_tracers))
  if (allocated(traj%tracer_names)) deallocate(traj%tracer_names)
  allocate(traj%tracer_names(number_tracers))
end if

! Fill the tracers
ft = 1
do f = 1, state%nf
  if (state%fields(f)%tracer) then

    ! Skip parameterizaion specific cloud/moisture tracers
    if (trim(state%fields(f)%long_name) == &
       'initial_mass_fraction_of_large_scale_cloud_condensate') cycle
    if (trim(state%fields(f)%long_name) == &
       'initial_mass_fraction_of_convective_cloud_condensate') cycle
    if (trim(state%fields(f)%long_name) == &
       'convective_cloud_area_fraction') cycle

    ! Put specific humidity in the first spot
    if (trim(state%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    traj%tracers(:,:,:,index) = state%fields(f)%array
    traj%tracer_names(index)  = trim(state%fields(f)%long_name)
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("fv3jedi_traj_mod:set: (water_vapor_mixing_ratio_wrt_moist_air) sphum is not " // &
                 "listed in 'state variables'")
end if


! Copy all the optioanal parts of the trajectory
! ----------------------------------------------

! Copy optional parts of the trajecotry (Rank 3)
if (state%has_field('upward_air_velocity')) &
  call state%get_field('upward_air_velocity', traj%w)

if (state%has_field('layer_thickness')) &
  call state%get_field('layer_thickness', traj%delz)

if (state%has_field('initial_mass_fraction_of_large_scale_cloud_condensate')) &
  call state%get_field('initial_mass_fraction_of_large_scale_cloud_condensate', traj%qls)

if (state%has_field('initial_mass_fraction_of_convective_cloud_condensate')) &
  call state%get_field('initial_mass_fraction_of_convective_cloud_condensate', traj%qcn)

if (state%has_field('convective_cloud_area_fraction')) &
  call state%get_field('convective_cloud_area_fraction', traj%cfcn)

! Copy optional parts of the trajecotry (Rank 2)
if (state%has_field('geopotential_height_times_gravity_at_surface')) then
  call state%get_field('geopotential_height_times_gravity_at_surface', phis)
  traj%phis = phis(:,:,1)
endif
if (state%has_field('fraction_of_ocean')) then
  call state%get_field('fraction_of_ocean', frocean)
  traj%frocean = frocean(:,:,1)
endif
if (state%has_field('fraction_of_land')) then
  call state%get_field('fraction_of_land', frland)
  traj%frland = frland(:,:,1)
endif
if (state%has_field('isotropic_variance_of_filtered_topography')) then
  call state%get_field('isotropic_variance_of_filtered_topography', varflt)
  traj%varflt = varflt(:,:,1)
endif
if (state%has_field('surface_velocity_scale')) then
  call state%get_field('surface_velocity_scale', ustar)
  traj%ustar = ustar(:,:,1)
endif
if (state%has_field('surface_buoyancy_scale')) then
  call state%get_field('surface_buoyancy_scale', bstar)
  traj%bstar = bstar(:,:,1)
endif
if (state%has_field('planetary_boundary_layer_height')) then
  call state%get_field('planetary_boundary_layer_height', zpbl)
  traj%zpbl = zpbl(:,:,1)
endif
if (state%has_field('surface_exchange_coefficient_for_momentum')) then
  call state%get_field('surface_exchange_coefficient_for_momentum', cm)
  traj%cm = cm(:,:,1)
endif
if (state%has_field('surface_exchange_coefficient_for_heat')) then
  call state%get_field('surface_exchange_coefficient_for_heat', ct)
  traj%ct = ct(:,:,1)
endif
if (state%has_field('surface_exchange_coefficient_for_moisture')) then
  call state%get_field('surface_exchange_coefficient_for_moisture', cq)
  traj%cq = cq(:,:,1)
endif
if (state%has_field('KCBL_before_moist')) then
  call state%get_field('KCBL_before_moist', kcbl)
  traj%kcbl = kcbl(:,:,1)
endif
if (state%has_field('surface_temp_before_moist')) then
  call state%get_field('surface_temp_before_moist', tsm)
  traj%ts = tsm(:,:,1)
endif
if (state%has_field('lower_index_where_Kh_greater_than_2')) then
  call state%get_field('lower_index_where_Kh_greater_than_2', khl)
  traj%khl = khl(:,:,1)
endif
if (state%has_field('upper_index_where_Kh_greater_than_2')) then
  call state%get_field('upper_index_where_Kh_greater_than_2', khu)
  traj%khu = khu(:,:,1)
endif

end subroutine set

! --------------------------------------------------------------------------------------------------

subroutine wipe(traj)

! Arguments
type(fv3jedi_traj), intent(inout) :: traj

traj%isc = 0
traj%iec = 0
traj%jsc = 0
traj%jec = 0
traj%npz = 0
traj%ntracers = 0

if (allocated(traj%ua))           deallocate(traj%ua)
if (allocated(traj%va))           deallocate(traj%va)
if (allocated(traj%t))            deallocate(traj%t)
if (allocated(traj%delp))         deallocate(traj%delp)
if (allocated(traj%tracers))      deallocate(traj%tracers)
if (allocated(traj%tracer_names)) deallocate(traj%tracer_names)
if (allocated(traj%w))            deallocate(traj%w)
if (allocated(traj%delz))         deallocate(traj%delz)
if (allocated(traj%cfcn))         deallocate(traj%cfcn)
if (allocated(traj%qls))          deallocate(traj%qls)
if (allocated(traj%qcn))          deallocate(traj%qcn)
if (allocated(traj%phis))         deallocate(traj%phis)
if (allocated(traj%frocean))      deallocate(traj%frocean)
if (allocated(traj%frland))       deallocate(traj%frland)
if (allocated(traj%varflt))       deallocate(traj%varflt)
if (allocated(traj%ustar))        deallocate(traj%ustar)
if (allocated(traj%bstar))        deallocate(traj%bstar)
if (allocated(traj%zpbl))         deallocate(traj%zpbl)
if (allocated(traj%cm))           deallocate(traj%cm)
if (allocated(traj%ct))           deallocate(traj%ct)
if (allocated(traj%cq))           deallocate(traj%cq)
if (allocated(traj%kcbl))         deallocate(traj%kcbl)
if (allocated(traj%ts))           deallocate(traj%ts)
if (allocated(traj%khl))          deallocate(traj%khl)
if (allocated(traj%khu))          deallocate(traj%khu)

end subroutine wipe

! --------------------------------------------------------------------------------------------------

end module fv3jedi_traj_mod
