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

subroutine initialize_ad(self, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

! Make sure the tracers are allocated (true => both traj and pert tracers)
call self%fv3jedi_lm%allocate_tracers(inc%isc, inc%iec, inc%jsc, inc%jec, &
                                      inc%npz, inc%ntracers, .true.)
call traj_to_traj(traj, self%fv3jedi_lm)

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

call lm_to_inc_ad(self%fv3jedi_lm,inc)
call self%fv3jedi_lm%init_ad()

end subroutine initialize_ad

! --------------------------------------------------------------------------------------------------

subroutine initialize_tl(self, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

! Make sure the tracers are allocated (true => both traj and pert tracers)
call self%fv3jedi_lm%allocate_tracers(inc%isc, inc%iec, inc%jsc, inc%jec, &
                                      inc%npz, inc%ntracers, .true.)
call traj_to_traj(traj, self%fv3jedi_lm)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine initialize_tl

! --------------------------------------------------------------------------------------------------

subroutine step_ad(self, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj, self%fv3jedi_lm)

call lm_to_inc_ad(self%fv3jedi_lm,inc)
call self%fv3jedi_lm%step_ad()

end subroutine step_ad

! --------------------------------------------------------------------------------------------------

subroutine step_tl(self, inc, traj)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj, self%fv3jedi_lm)

call self%fv3jedi_lm%step_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine step_tl

! --------------------------------------------------------------------------------------------------

subroutine finalize_ad(self, inc)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call lm_to_inc_ad(self%fv3jedi_lm,inc)
call self%fv3jedi_lm%final_ad()
call inc_to_lm_ad(inc,self%fv3jedi_lm)

end subroutine finalize_ad

! --------------------------------------------------------------------------------------------------

subroutine finalize_tl(self, inc)

class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call self%fv3jedi_lm%final_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

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

subroutine inc_to_lm(inc, lm)

type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd

! Bounds mismatch for D-Grid winds so first get pointer
call inc%get_field('ud', ud  )
call inc%get_field('vd', vd  )
lm%pert%u = ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
lm%pert%v = vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)

call inc%get_field('t'      , lm%pert%t   )
call inc%get_field('delp'   , lm%pert%delp)

! Tracers
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%short_name) == 'sphum') then
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
if (inc%has_field('w'      )) call inc%get_field('w'      , lm%pert%w   )
if (inc%has_field('delz'   )) call inc%get_field('delz'   , lm%pert%delz)

end subroutine inc_to_lm

! --------------------------------------------------------------------------------------------------

subroutine lm_to_inc(lm, inc)

type(fv3jedi_lm_type),   intent(in)    :: lm
type(fv3jedi_increment), intent(inout) :: inc

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd

! Bounds mismatch for D-Grid winds so first get pointer
call inc%get_field('ud'     , ud  )
call inc%get_field('vd'     , vd  )
ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%u
vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%v

call inc%put_field('t'      , lm%pert%t   )
call inc%put_field('delp'   , lm%pert%delp)

! Tracers
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%short_name) == 'sphum') then
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
if (inc%has_field('w'   )) call inc%put_field('w'   , lm%pert%w   )
if (inc%has_field('delz')) call inc%put_field('delz', lm%pert%delz)

end subroutine lm_to_inc

! --------------------------------------------------------------------------------------------------

subroutine inc_to_lm_ad(inc, lm)

implicit none
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), pointer, dimension(:,:,:) :: tmp

call inc%get_field('ud', tmp)
tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz) = tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz) + lm%pert%u
call inc%get_field('vd', tmp)
tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz) = tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz) + lm%pert%v

call inc%get_field('t', tmp)
tmp = tmp + lm%pert%t
call inc%get_field('delp', tmp)
tmp = tmp + lm%pert%delp

! Tracers
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%short_name) == 'sphum') then
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

! Optional fields
if (inc%has_field('w')) then
  call inc%get_field('w', tmp)
  tmp = tmp + lm%pert%w
end if
if (inc%has_field('delz')) then
  call inc%get_field('delz', tmp)
  tmp = tmp + lm%pert%delz
end if

lm%pert%u = 0.0_kind_real
lm%pert%v = 0.0_kind_real
lm%pert%t = 0.0_kind_real
lm%pert%delp = 0.0_kind_real
if (allocated(lm%pert%tracers)) lm%pert%tracers(:,:,:,:) = 0.0_kind_real
if (inc%has_field('w')) lm%pert%w = 0.0_kind_real
if (inc%has_field('delz')) lm%pert%delz = 0.0_kind_real

end subroutine inc_to_lm_ad

! --------------------------------------------------------------------------------------------------

subroutine lm_to_inc_ad(lm, inc)

implicit none
type(fv3jedi_lm_type),   intent(inout) :: lm
type(fv3jedi_increment), intent(inout) :: inc

integer :: ft, f, index
logical :: sphum_found = .false.

real(kind=kind_real), pointer, dimension(:,:,:) :: tmp

call inc%get_field('ud', tmp)
lm%pert%u = lm%pert%u + tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
call inc%get_field('vd', tmp)
lm%pert%v = lm%pert%v + tmp(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)

call inc%get_field('t', tmp)
lm%pert%t = lm%pert%t + tmp
call inc%get_field('delp', tmp)
lm%pert%delp = lm%pert%delp + tmp

! Tracers
ft = 1
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    if (trim(inc%fields(f)%short_name) == 'sphum') then
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

! Optional fields
if (inc%has_field('w')) then
  call inc%get_field('w', tmp)
  lm%pert%w = lm%pert%w + tmp
end if
if (inc%has_field('delz')) then
  call inc%get_field('delz', tmp)
  lm%pert%delz = lm%pert%delz + tmp
end if

! Here we could call `inc%zero()` but that would zero additional increment fields that aren't
! participating in the LM forecast, and violate the adjoint. If in the future the increment fields
! are required to all be LM fields, then the code below could be replaced with just `inc%zero()`.
call inc%get_field('ud', tmp)
tmp = 0.0_kind_real
call inc%get_field('vd', tmp)
tmp = 0.0_kind_real
call inc%get_field('t', tmp)
tmp = 0.0_kind_real
call inc%get_field('delp', tmp)
tmp = 0.0_kind_real

! This assumes all tracers are evolved via LM
do f = 1, inc%nf
  if (inc%fields(f)%tracer) then
    inc%fields(f)%array = 0.0_kind_real
  end if
end do

if (inc%has_field('w')) then
  call inc%get_field('w', tmp)
  tmp = 0.0_kind_real
end if
if (inc%has_field('delz')) then
  call inc%get_field('delz', tmp)
  tmp = 0.0_kind_real
end if

end subroutine lm_to_inc_ad

! --------------------------------------------------------------------------------------------------

subroutine traj_to_traj( traj, lm )

type(fv3jedi_traj),    intent(in)    :: traj
type(fv3jedi_lm_type), intent(inout) :: lm

lm%traj%u    = traj%u
lm%traj%v    = traj%v
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
