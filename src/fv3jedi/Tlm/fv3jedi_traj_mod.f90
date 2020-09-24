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
use fv3jedi_field_mod,    only: has_field, copy_field_array, pointer_field_array

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

integer :: isc,iec,jsc,jec,npz

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

! Allocate traj
allocate(self%u      (isc:iec, jsc:jec, npz))
allocate(self%v      (isc:iec, jsc:jec, npz))
allocate(self%ua     (isc:iec, jsc:jec, npz))
allocate(self%va     (isc:iec, jsc:jec, npz))
allocate(self%t      (isc:iec, jsc:jec, npz))
allocate(self%delp   (isc:iec, jsc:jec, npz))
allocate(self%qv     (isc:iec, jsc:jec, npz))
allocate(self%qi     (isc:iec, jsc:jec, npz))
allocate(self%ql     (isc:iec, jsc:jec, npz))
allocate(self%o3     (isc:iec, jsc:jec, npz))
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

!Initialize all to zero incase not in state
self%u       = 0.0_kind_real
self%v       = 0.0_kind_real
self%ua      = 0.0_kind_real
self%va      = 0.0_kind_real
self%t       = 0.0_kind_real
self%delp    = 0.0_kind_real
self%qv      = 0.0_kind_real
self%qi      = 0.0_kind_real
self%ql      = 0.0_kind_real
self%o3      = 0.0_kind_real
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
allocate(u_tmp(isc:iec  , jsc:jec+1, npz))
allocate(v_tmp(isc:iec+1, jsc:jec  , npz))

call copy_field_array(state%fields, 'ud'  , u_tmp     )
call copy_field_array(state%fields, 'vd'  , v_tmp     )
call copy_field_array(state%fields, 't'   , self%t    )
call copy_field_array(state%fields, 'delp', self%delp )
call copy_field_array(state%fields, 'sphum'   , self%qv   )

self%u = u_tmp(isc:iec, jsc:jec, :)
self%v = v_tmp(isc:iec, jsc:jec, :)

deallocate(u_tmp, v_tmp)

! Copy optional parts of the trajecotry (Rank 3)
if (has_field(state%fields, 'ua'     )) call copy_field_array(state%fields, 'ua'     , self%ua  )
if (has_field(state%fields, 'va'     )) call copy_field_array(state%fields, 'va'     , self%va  )
if (has_field(state%fields, 'ice_wat')) call copy_field_array(state%fields, 'ice_wat', self%qi  )
if (has_field(state%fields, 'liq_wat')) call copy_field_array(state%fields, 'liq_wat', self%ql  )
if (has_field(state%fields, 'o3mr'   )) call copy_field_array(state%fields, 'o3mr'   , self%o3  )
if (has_field(state%fields, 'w'      )) call copy_field_array(state%fields, 'w'      , self%w   )
if (has_field(state%fields, 'delz'   )) call copy_field_array(state%fields, 'delz'   , self%delz)
if (has_field(state%fields, 'qls'    )) call copy_field_array(state%fields, 'qls'    , self%qls )
if (has_field(state%fields, 'qcn'    )) call copy_field_array(state%fields, 'qcn'    , self%qcn )
if (has_field(state%fields, 'cfcn'   )) call copy_field_array(state%fields, 'cfcn'   , self%cfcn)

! Copy optional parts of the trajecotry (Rank 2)
if (has_field(state%fields, 'phis')) then
  call pointer_field_array(state%fields, 'phis', phis)
  self%phis = phis(:,:,1)
endif
if (has_field(state%fields, 'frocean')) then
  call pointer_field_array(state%fields, 'frocean', frocean)
  self%frocean = frocean(:,:,1)
endif
if (has_field(state%fields, 'frland')) then
  call pointer_field_array(state%fields, 'frland', frland)
  self%frland = frland(:,:,1)
endif
if (has_field(state%fields, 'varflt')) then
  call pointer_field_array(state%fields, 'varflt', varflt)
  self%varflt = varflt(:,:,1)
endif
if (has_field(state%fields, 'ustar')) then
  call pointer_field_array(state%fields, 'ustar', ustar)
  self%ustar = ustar(:,:,1)
endif
if (has_field(state%fields, 'bstar')) then
  call pointer_field_array(state%fields, 'bstar', bstar)
  self%bstar = bstar(:,:,1)
endif
if (has_field(state%fields, 'zpbl')) then
  call pointer_field_array(state%fields, 'zpbl', zpbl)
  self%zpbl = zpbl(:,:,1)
endif
if (has_field(state%fields, 'cm')) then
  call pointer_field_array(state%fields, 'cm', cm)
  self%cm = cm(:,:,1)
endif
if (has_field(state%fields, 'ct')) then
  call pointer_field_array(state%fields, 'ct', ct)
  self%ct = ct(:,:,1)
endif
if (has_field(state%fields, 'cq')) then
  call pointer_field_array(state%fields, 'cq', cq)
  self%cq = cq(:,:,1)
endif
if (has_field(state%fields, 'kcbl')) then
  call pointer_field_array(state%fields, 'kcbl', kcbl)
  self%kcbl = kcbl(:,:,1)
endif
if (has_field(state%fields, 'tsm')) then
  call pointer_field_array(state%fields, 'tsm', tsm)
  self%ts = tsm(:,:,1)
endif
if (has_field(state%fields, 'khl')) then
  call pointer_field_array(state%fields, 'khl', khl)
  self%khl = khl(:,:,1)
endif
if (has_field(state%fields, 'khu')) then
  call pointer_field_array(state%fields, 'khu', khu)
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
