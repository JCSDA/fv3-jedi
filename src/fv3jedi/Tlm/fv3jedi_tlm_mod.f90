! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_tlm_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use duration_mod
use oops_variables_mod

! fv3-jedi-linearmodel
use fv3jedi_lm_mod, only: fv3jedi_lm_type

! fv3-jedi
use fv3jedi_fmsnamelist_mod, only: fv3jedi_fmsnamelist
use fv3jedi_geom_mod,        only: fv3jedi_geom
use fv3jedi_kinds_mod,       only: kind_real
use fv3jedi_increment_mod,   only: fv3jedi_increment
use fv3jedi_state_mod,       only: fv3jedi_state
use fv3jedi_traj_mod,        only: fv3jedi_traj
use wind_vt_mod,             only: a_to_d, d_to_a, a_to_d_ad, d_to_a_ad

implicit none
private
public :: fv3jedi_tlm

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold tlm definition
type:: fv3jedi_tlm
  type(fv3jedi_lm_type) :: fv3jedi_lm  !<Linearized model object
  contains
    procedure :: create
    procedure :: delete
    procedure :: initialize_tl
    procedure :: initialize_ad
    procedure :: step_tl
    procedure :: step_ad
    procedure :: finalize_tl
    procedure :: finalize_ad
end type fv3jedi_tlm

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_tlm),        intent(inout) :: self
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

! Model configuration and creation
! --------------------------------
call conf%get_or_die("lm_do_dyn",self%fv3jedi_lm%conf%do_dyn)
call conf%get_or_die("lm_do_trb",self%fv3jedi_lm%conf%do_phy_trb)
call conf%get_or_die("lm_do_mst",self%fv3jedi_lm%conf%do_phy_mst)

! Prepare namelist (nonlinear part)
! ---------------------------------
call fmsnamelist%replace_namelist(conf)

! Filename for namelist (tl/ad part)
! ----------------------------------
call conf%get_or_die("linear model namelist filename", str)
self%fv3jedi_lm%conf%inputpert_filename = str
deallocate(str)

! Call linear model constructor
! -----------------------------
call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_tlm), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize_ad(self, geom, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

! Make sure the tracers are allocated (true => both traj and pert tracers)
call self%fv3jedi_lm%allocate_tracers(inc%ntracers, .true.)
call traj_to_traj(geom, traj, self%fv3jedi_lm)

! The action of finalize_tl is {inc = LM; LM = 0}, the adjoint of which is {LM* = inc*; inc* = 0}.
! To keep the code simple, we first perform the extra step {LM* = 0}, so that we can then re-use
! the implementation of lm_to_inc_ad, which performs {LM* = LM* + inc*; inc* = 0}.
! Replacing lm_to_inc_ad with a method that just does the assignment may be a small optimization.

! Zero prior LM state
self%fv3jedi_lm%pert%u = 0.0_kind_real
self%fv3jedi_lm%pert%v = 0.0_kind_real
self%fv3jedi_lm%pert%t = 0.0_kind_real
self%fv3jedi_lm%pert%delp = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%tracers)) self%fv3jedi_lm%pert%tracers = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%w)) self%fv3jedi_lm%pert%w = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%delz)) self%fv3jedi_lm%pert%delz = 0.0_kind_real

call lm_to_inc_ad(geom, self%fv3jedi_lm,inc)
call self%fv3jedi_lm%init_ad()

end subroutine initialize_ad

! --------------------------------------------------------------------------------------------------

subroutine initialize_tl(self, geom, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

! Make sure the tracers are allocated (true => both traj and pert tracers)
call self%fv3jedi_lm%allocate_tracers(inc%ntracers, .true.)
call traj_to_traj(geom, traj, self%fv3jedi_lm)

call inc_to_lm(geom, inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_tl()
call lm_to_inc(geom, self%fv3jedi_lm,inc)

end subroutine initialize_tl

! --------------------------------------------------------------------------------------------------

subroutine step_ad(self, geom, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(geom, traj, self%fv3jedi_lm)

call lm_to_inc_ad(geom, self%fv3jedi_lm,inc)
call self%fv3jedi_lm%step_ad()

end subroutine step_ad

! --------------------------------------------------------------------------------------------------

subroutine step_tl(self, geom, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(geom, traj, self%fv3jedi_lm)

call self%fv3jedi_lm%step_tl()
call lm_to_inc(geom, self%fv3jedi_lm,inc)

end subroutine step_tl

! --------------------------------------------------------------------------------------------------

subroutine finalize_ad(self, geom, inc)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc

call lm_to_inc_ad(geom, self%fv3jedi_lm,inc)
call self%fv3jedi_lm%final_ad()
call inc_to_lm_ad(geom, inc,self%fv3jedi_lm)

end subroutine finalize_ad

! --------------------------------------------------------------------------------------------------

subroutine finalize_tl(self, geom, inc)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc

call self%fv3jedi_lm%final_tl()
call lm_to_inc(geom, self%fv3jedi_lm,inc)

! Destroy the LM state
self%fv3jedi_lm%pert%u = 0.0_kind_real
self%fv3jedi_lm%pert%v = 0.0_kind_real
self%fv3jedi_lm%pert%t = 0.0_kind_real
self%fv3jedi_lm%pert%delp = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%tracers)) self%fv3jedi_lm%pert%tracers = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%w)) self%fv3jedi_lm%pert%w = 0.0_kind_real
if (allocated(self%fv3jedi_lm%pert%delz)) self%fv3jedi_lm%pert%delz = 0.0_kind_real

end subroutine finalize_tl

! --------------------------------------------------------------------------------------------------

subroutine inc_to_lm(geom, inc, lm)

type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()


! Convert the increment A-Grid winds to model D-Grid winds
! --------------------------------------------------------

! Get pointers to A-Grid winds
call inc%get_field('eastward_wind', ua)
call inc%get_field('northward_wind', va)

! Allocate some temporary D-Grid winds with edges
allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))

! Convert
call a_to_d(geom, ua, va, ud, vd)

! Copy D-Grid to model internal
lm%pert%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = ud(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
lm%pert%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = vd(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)


! Temperature / pressure
! ----------------------
call inc%get_field('air_temperature', lm%pert%t)
call inc%get_field('air_pressure_thickness', lm%pert%delp)


! Tracers
! -------
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    lm%pert%tracers(:,:,:,index) = inc%fields(f)%array
    lm%pert%tracer_names(index) = inc%fields(f)%long_name
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("inc_to_lm: sphum is not listed in 'state variables'")
end if

! Optional fields
! ---------------
if (inc%has_field('upward_air_velocity')) call inc%get_field('upward_air_velocity', lm%pert%w)
if (inc%has_field('layer_thickness')) call inc%get_field('layer_thickness', lm%pert%delz)

end subroutine inc_to_lm

! --------------------------------------------------------------------------------------------------

subroutine lm_to_inc(geom, lm, inc)

type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_lm_type),   intent(in)    :: lm
type(fv3jedi_increment), intent(inout) :: inc

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()

! Convert internal D-Grid winds to A-Grid winds
! ---------------------------------------------

! Allocate D-Grid winds and copy to internal part
allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
ud = 0.0_kind_real
ud(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = lm%pert%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
vd = 0.0_kind_real
vd(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = lm%pert%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)

! Pointer to increment A-Grid winds
call inc%get_field('eastward_wind', ua)
call inc%get_field('northward_wind', va)

! Convert to A-Grid
call d_to_a(geom, ud, vd, ua, va)


! Temperature / pressure
! ----------------------
call inc%put_field('air_temperature', lm%pert%t)
call inc%put_field('air_pressure_thickness', lm%pert%delp)


! Tracers
! -------
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    endif
    inc%fields(f)%array = lm%pert%tracers(:,:,:,index)
    if (inc%fields(f)%long_name .ne. lm%pert%tracer_names(index)) then
      call abor1_ftn("lm_to_inc: Tracer names do not match")
    end if
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("lm_to_inc: sphum is not listed in 'state variables'")
end if


! Optional fields
! ---------------
if (inc%has_field('upward_air_velocity')) call inc%put_field('upward_air_velocity', lm%pert%w)
if (inc%has_field('layer_thickness')) call inc%put_field('layer_thickness', lm%pert%delz)


end subroutine lm_to_inc

! --------------------------------------------------------------------------------------------------

subroutine inc_to_lm_ad(geom, inc, lm)

implicit none
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: tmp => null()

! Optional fields
! ---------------
if (inc%has_field('upward_air_velocity')) then
  call inc%get_field('upward_air_velocity', tmp)
  tmp = tmp + lm%pert%w
end if
if (inc%has_field('layer_thickness')) then
  call inc%get_field('layer_thickness', tmp)
  tmp = tmp + lm%pert%delz
end if

! Tracers
! -------
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    endif
    inc%fields(f)%array = inc%fields(f)%array + lm%pert%tracers(:,:,:,index)
    if (inc%fields(f)%long_name .ne. lm%pert%tracer_names(index)) then
      call abor1_ftn("inc_to_lm_ad: Tracer names do not match")
    end if
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("inc_to_lm_ad: sphum is not listed in 'state variables'")
end if


! Temperature / pressure
! ----------------------
call inc%get_field('air_temperature', tmp)
tmp = tmp + lm%pert%t
call inc%get_field('air_pressure_thickness', tmp)
tmp = tmp + lm%pert%delp


! Winds (adjoint of A to D conversion)
! ------------------------------------

! Allocate some temporary D-Grid winds with edges
allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
! Copy D Grid from the model
ud = 0.0_kind_real
ud(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                         lm%pert%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
vd = 0.0_kind_real
vd(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                         lm%pert%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)

! Get pointers to A-Grid winds
call inc%get_field('eastward_wind', ua)
call inc%get_field('northward_wind', va)

! Convert (accumulation is internal)
call a_to_d_ad(geom, ua, va, ud, vd)

! Zero out the internal model increment
! -------------------------------------
lm%pert%u = 0.0_kind_real
lm%pert%v = 0.0_kind_real
lm%pert%t = 0.0_kind_real
lm%pert%delp = 0.0_kind_real
if (allocated(lm%pert%tracers)) lm%pert%tracers(:,:,:,:) = 0.0_kind_real
if (inc%has_field('upward_air_velocity')) lm%pert%w = 0.0_kind_real
if (inc%has_field('layer_thickness')) lm%pert%delz = 0.0_kind_real

end subroutine inc_to_lm_ad

! --------------------------------------------------------------------------------------------------

subroutine lm_to_inc_ad(geom, lm, inc)

implicit none
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_lm_type),   intent(inout) :: lm
type(fv3jedi_increment), intent(inout) :: inc

integer :: ft, f, index, k
logical :: sphum_found = .false.

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua   => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: va   => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: t    => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delp => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: w    => null()
real(kind=kind_real), pointer, dimension(:,:,:) :: delz => null()

! Optional fields
! ---------------
if (inc%has_field('upward_air_velocity')) then
  call inc%get_field('upward_air_velocity', w)
  lm%pert%w = lm%pert%w + w
end if
if (inc%has_field('layer_thickness')) then
  call inc%get_field('layer_thickness', delz)
  lm%pert%delz = lm%pert%delz + delz
end if


! Tracers
! -------
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%long_name) == 'water_vapor_mixing_ratio_wrt_moist_air') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    lm%pert%tracers(:,:,:,index) = lm%pert%tracers(:,:,:,index) + inc%fields(f)%array
    lm%pert%tracer_names(index) = inc%fields(f)%long_name
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("lm_to_inc_ad: sphum is not listed in 'state variables'")
end if


! Temperature / pressure
! ----------------------
call inc%get_field('air_temperature', t)
lm%pert%t = lm%pert%t + t
call inc%get_field('air_pressure_thickness', delp)
lm%pert%delp = lm%pert%delp + delp


! Adjoint of D-Grid winds to A-Grid winds
! ---------------------------------------

! Allocate D-Grid winds and copy to internal part
allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
ud = 0.0_kind_real  ! Set to zero because of internal accumulation
vd = 0.0_kind_real  ! Set to zero because of internal accumulation

! Pointer to increment A-Grid winds
call inc%get_field('eastward_wind', ua)
call inc%get_field('northward_wind', va)

! Convert to A-Grid
call d_to_a_ad(geom, ud, vd, ua, va)

! Copy temporary back into the model
lm%pert%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                     lm%pert%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) + &
                                            ud(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
lm%pert%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                     lm%pert%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) + &
                                            vd(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)


! Here we could call `inc%zero()` but that would zero additional increment fields that aren't
! participating in the LM forecast, and violate the adjoint. If in the future the increment fields
! are required to all be LM fields, then the code below could be replaced with just `inc%zero()`.
ua   = 0.0_kind_real
va   = 0.0_kind_real
t    = 0.0_kind_real
delp = 0.0_kind_real
if (associated(w)) w = 0.0_kind_real
if (associated(delz)) delz = 0.0_kind_real

! This assumes all tracers are evolved via LM
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    inc%fields(f)%array = 0.0_kind_real
  end if
end do

end subroutine lm_to_inc_ad

! --------------------------------------------------------------------------------------------------

subroutine traj_to_traj(geom, traj, lm)

type(fv3jedi_geom),    intent(in)    :: geom
type(fv3jedi_traj),    intent(in)    :: traj
type(fv3jedi_lm_type), intent(inout) :: lm

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd

! Convert the A-Grid trajectory to D-Grid trajectory needed by model
allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
call a_to_d(geom, traj%ua, traj%va, ud, vd)
lm%traj%u(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                                ud(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
lm%traj%v(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz) = &
                                                vd(geom%isc:geom%iec, geom%jsc:geom%jec, 1:geom%npz)
deallocate(ud, vd)

lm%traj%ua   = traj%ua
lm%traj%va   = traj%va
lm%traj%t    = traj%t
lm%traj%delp = traj%delp

! Copy the tracers
lm%traj%tracers = traj%tracers
lm%traj%tracer_names  = traj%tracer_names

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

! --------------------------------------------------------------------------------------------------

end module fv3jedi_tlm_mod
