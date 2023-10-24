! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_traj_mod

! fv3-jedi-lm
use fv3jedi_lm_utils_mod, only: fv3jedi_traj => fv3jedi_lm_traj, deallocate_traj

! fv3-jedi
use fv3jedi_kinds_mod,    only: kind_real
use fv3jedi_state_mod,    only: fv3jedi_state

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: fv3jedi_traj
public :: set
public :: wipe

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine set(self, state)

implicit none
type(fv3jedi_traj),  intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

integer :: isc,iec,jsc,jec,npz,ft,f,index,number_tracers
logical :: sphum_found = .false.

! Pointers to rank 2 state
real(kind=kind_real), allocatable, dimension(:,:,:) :: u_tmp
real(kind=kind_real), allocatable, dimension(:,:,:) :: v_tmp

real(kind=kind_real), pointer, dimension(:,:,:) :: phis
real(kind=kind_real), pointer, dimension(:,:,:) :: frocean
real(kind=kind_real), pointer, dimension(:,:,:) :: frland
real(kind=kind_real), pointer, dimension(:,:,:) :: varflt
real(kind=kind_real), pointer, dimension(:,:,:) :: ustar
real(kind=kind_real), pointer, dimension(:,:,:) :: bstar
real(kind=kind_real), pointer, dimension(:,:,:) :: zpbl
real(kind=kind_real), pointer, dimension(:,:,:) :: cm
real(kind=kind_real), pointer, dimension(:,:,:) :: ct
real(kind=kind_real), pointer, dimension(:,:,:) :: cq
real(kind=kind_real), pointer, dimension(:,:,:) :: kcbl
real(kind=kind_real), pointer, dimension(:,:,:) :: tsm
real(kind=kind_real), pointer, dimension(:,:,:) :: khl
real(kind=kind_real), pointer, dimension(:,:,:) :: khu

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
npz = state%npz

self%ntracers = state%ntracers

! Allocate traj
allocate(self%u      (isc:iec, jsc:jec, npz))
allocate(self%v      (isc:iec, jsc:jec, npz))
allocate(self%ua     (isc:iec, jsc:jec, npz))
allocate(self%va     (isc:iec, jsc:jec, npz))
allocate(self%t      (isc:iec, jsc:jec, npz))
allocate(self%delp   (isc:iec, jsc:jec, npz))
allocate(self%w      (isc:iec, jsc:jec, npz))
allocate(self%delz   (isc:iec, jsc:jec, npz))
allocate(self%qls    (isc:iec, jsc:jec, npz))
allocate(self%qcn    (isc:iec, jsc:jec, npz))
allocate(self%cfcn   (isc:iec, jsc:jec, npz))
allocate(self%phis   (isc:iec, jsc:jec))
allocate(self%frocean(isc:iec, jsc:jec))
allocate(self%frland (isc:iec, jsc:jec))
allocate(self%varflt (isc:iec, jsc:jec))
allocate(self%ustar  (isc:iec, jsc:jec))
allocate(self%bstar  (isc:iec, jsc:jec))
allocate(self%zpbl   (isc:iec, jsc:jec))
allocate(self%cm     (isc:iec, jsc:jec))
allocate(self%ct     (isc:iec, jsc:jec))
allocate(self%cq     (isc:iec, jsc:jec))
allocate(self%kcbl   (isc:iec, jsc:jec))
allocate(self%ts     (isc:iec, jsc:jec))
allocate(self%khl    (isc:iec, jsc:jec))
allocate(self%khu    (isc:iec, jsc:jec))

! Initialize all to zero incase not in state
! ------------------------------------------
self%u       = 0.0_kind_real
self%v       = 0.0_kind_real
self%ua      = 0.0_kind_real
self%va      = 0.0_kind_real
self%t       = 0.0_kind_real
self%delp    = 0.0_kind_real
self%w       = 0.0_kind_real
self%delz    = 0.0_kind_real
self%qls     = 0.0_kind_real
self%qcn     = 0.0_kind_real
self%cfcn    = 0.0_kind_real
self%phis    = 0.0_kind_real
self%frocean = 0.0_kind_real
self%frland  = 0.0_kind_real
self%varflt  = 0.0_kind_real
self%ustar   = 0.0_kind_real
self%bstar   = 0.0_kind_real
self%zpbl    = 0.0_kind_real
self%cm      = 0.0_kind_real
self%ct      = 0.0_kind_real
self%cq      = 0.0_kind_real
self%kcbl    = 0.0_kind_real
self%ts      = 0.0_kind_real
self%khl     = 0.0_kind_real
self%khu     = 0.0_kind_real


! Copy mandatory parts of the trajecotry
! --------------------------------------
allocate(u_tmp(isc:iec  , jsc:jec+1, npz))
allocate(v_tmp(isc:iec+1, jsc:jec  , npz))

call state%get_field('ud'   , u_tmp     )
call state%get_field('vd'   , v_tmp     )
call state%get_field('t'    , self%t    )
call state%get_field('delp' , self%delp )

self%u = u_tmp(isc:iec, jsc:jec, :)
self%v = v_tmp(isc:iec, jsc:jec, :)

deallocate(u_tmp, v_tmp)


! Allocate all the tracers that will be part of what gets advected
! ----------------------------------------------------------------

! First allocate the number of tracers
number_tracers = state%ntracers

! Remove qls, qcn, cfcn
if (state%has_field('qls'))  number_tracers = number_tracers - 1
if (state%has_field('qcn'))  number_tracers = number_tracers - 1
if (state%has_field('cfcn')) number_tracers = number_tracers - 1

! Allocate the tracers

if (.not. allocated(self%tracers) .or. size(self%tracers, 4) .ne. number_tracers) then
  if (allocated(self%tracers)) deallocate(self%tracers)
  allocate(self%tracers(state%isc:state%iec, state%jsc:state%jec, state%npz, number_tracers))
  if (allocated(self%tracer_names)) deallocate(self%tracer_names)
  allocate(self%tracer_names(number_tracers))
end if

! Fill the tracers
ft = 1
do f = 1, state%nf
  if (state%fields(f)%tracer) then

    ! Skip if the short_name is qls, qcn or cfcn
    if (trim(state%fields(f)%short_name) == 'qls'  ) cycle
    if (trim(state%fields(f)%short_name) == 'qcn'  ) cycle
    if (trim(state%fields(f)%short_name) == 'cfcn' ) cycle

    ! Put specific humidity in the first spot
    if (trim(state%fields(f)%short_name) == 'sphum') then
      index = 1
      sphum_found = .true.
    else
      ft = ft + 1
      index = ft
    end if
    self%tracers(:,:,:,index) = state%fields(f)%array
    self%tracer_names(index)  = trim(state%fields(f)%long_name)
  end if
end do

if(.not.sphum_found) then
  call abor1_ftn("fv3jedi_traj_mod:set: sphum is not listed in 'state variables'")
end if


! Copy all the optioanal parts of the trajectory
! ----------------------------------------------

! Copy optional parts of the trajecotry (Rank 3)
if (state%has_field('ua'  )) call state%get_field('ua'  , self%ua  )
if (state%has_field('va'  )) call state%get_field('va'  , self%va  )
if (state%has_field('w'   )) call state%get_field('w'   , self%w   )
if (state%has_field('delz')) call state%get_field('delz', self%delz)
if (state%has_field('qls' )) call state%get_field('qls' , self%qls )
if (state%has_field('qcn' )) call state%get_field('qcn' , self%qcn )
if (state%has_field('cfcn')) call state%get_field('cfcn', self%cfcn)

! Copy optional parts of the trajecotry (Rank 2)
if (state%has_field('phis')) then
  call state%get_field('phis', phis)
  self%phis = phis(:,:,1)
endif
if (state%has_field('frocean')) then
  call state%get_field('frocean', frocean)
  self%frocean = frocean(:,:,1)
endif
if (state%has_field('frland')) then
  call state%get_field('frland', frland)
  self%frland = frland(:,:,1)
endif
if (state%has_field('varflt')) then
  call state%get_field('varflt', varflt)
  self%varflt = varflt(:,:,1)
endif
if (state%has_field('ustar')) then
  call state%get_field('ustar', ustar)
  self%ustar = ustar(:,:,1)
endif
if (state%has_field('bstar')) then
  call state%get_field('bstar', bstar)
  self%bstar = bstar(:,:,1)
endif
if (state%has_field('zpbl')) then
  call state%get_field('zpbl', zpbl)
  self%zpbl = zpbl(:,:,1)
endif
if (state%has_field('cm')) then
  call state%get_field('cm', cm)
  self%cm = cm(:,:,1)
endif
if (state%has_field('ct')) then
  call state%get_field('ct', ct)
  self%ct = ct(:,:,1)
endif
if (state%has_field('cq')) then
  call state%get_field('cq', cq)
  self%cq = cq(:,:,1)
endif
if (state%has_field('kcbl')) then
  call state%get_field('kcbl', kcbl)
  self%kcbl = kcbl(:,:,1)
endif
if (state%has_field('tsm')) then
  call state%get_field('tsm', tsm)
  self%ts = tsm(:,:,1)
endif
if (state%has_field('khl')) then
  call state%get_field('khl', khl)
  self%khl = khl(:,:,1)
endif
if (state%has_field('khu')) then
  call state%get_field('khu', khu)
  self%khu = khu(:,:,1)
endif

end subroutine set

! --------------------------------------------------------------------------------------------------

subroutine wipe(self)

implicit none
type(fv3jedi_traj), pointer :: self

call deallocate_traj(self)

end subroutine wipe

! --------------------------------------------------------------------------------------------------

end module fv3jedi_traj_mod
