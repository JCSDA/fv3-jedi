! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_traj_mod

use fv3jedi_kinds_mod
use iso_c_binding
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_lm_utils_mod, only: fv3jedi_traj => fv3jedi_lm_traj, &
                                allocate_traj, deallocate_traj

implicit none
private

public :: fv3jedi_traj
public :: traj_prop
public :: traj_wipe
public :: traj_minmaxrms

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine traj_prop(state, self)

implicit none
type(fv3jedi_state) :: state
type(fv3jedi_traj)  :: self

integer :: isc,iec,jsc,jec,npz

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jsc
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
if (associated(state%ud)) then
  self%u = state%ud
else
  call abor1_ftn("fv3jedi_traj_mod.traj_prop: ud not found in state, minimally needed for TL/AD")
endif
if (associated(state%vd)) then
  self%v = state%vd
else
  call abor1_ftn("fv3jedi_traj_mod.traj_prop: vd not found in state, minimally needed for TL/AD")
endif
if (associated(state%t)) then
  self%t = state%t
else
  call abor1_ftn("fv3jedi_traj_mod.traj_prop: t not found in state, minimally needed for TL/AD")
endif
if (associated(state%delp)) then
  self%delp = state%delp
else
  call abor1_ftn("fv3jedi_traj_mod.traj_prop: delp not found in state, minimally needed for TL/AD")
endif
if (associated(state%q)) then
  self%qv = state%q
else
  call abor1_ftn("fv3jedi_traj_mod.traj_prop: q not found in state, minimally needed for TL/AD")
endif

! Copy optional parts of the trajecotry
if (associated(state%ua))      self%ua      = state%ua
if (associated(state%va))      self%va      = state%va
if (associated(state%qi))      self%qi      = state%qi
if (associated(state%ql))      self%ql      = state%ql
if (associated(state%o3))      self%o3      = state%o3
if (associated(state%w))       self%w       = state%w
if (associated(state%delz))    self%delz    = state%delz
if (associated(state%qls))     self%qls     = state%qls
if (associated(state%qcn))     self%qcn     = state%qcn
if (associated(state%cfcn))    self%cfcn    = state%cfcn
if (associated(state%phis))    self%phis    = state%phis(:,:,1)
if (associated(state%frocean)) self%frocean = state%frocean(:,:,1)
if (associated(state%frland))  self%frland  = state%frland(:,:,1)
if (associated(state%varflt))  self%varflt  = state%varflt(:,:,1)
if (associated(state%ustar))   self%ustar   = state%ustar(:,:,1)
if (associated(state%bstar))   self%bstar   = state%bstar(:,:,1)
if (associated(state%zpbl))    self%zpbl    = state%zpbl(:,:,1)
if (associated(state%cm))      self%cm      = state%cm(:,:,1)
if (associated(state%ct))      self%ct      = state%ct(:,:,1)
if (associated(state%cq))      self%cq      = state%cq(:,:,1)
if (associated(state%kcbl))    self%kcbl    = state%kcbl(:,:,1)
if (associated(state%tsm))     self%ts      = state%tsm(:,:,1)
if (associated(state%khl))     self%khl     = state%khl(:,:,1)
if (associated(state%khu))     self%khu     = state%khu(:,:,1)

end subroutine traj_prop

! ------------------------------------------------------------------------------

subroutine traj_wipe(self)

implicit none
type(fv3jedi_traj), pointer :: self

call deallocate_traj(self)

end subroutine traj_wipe

! ------------------------------------------------------------------------------

subroutine traj_minmaxrms(self, pminmax)

implicit none
type(fv3jedi_traj), intent(in)    :: self
real(c_double),     intent(inout) :: pminmax(3,11)

real(kind=kind_real) :: zz

pminmax = 0.0_kind_real

zz=real(size(self%u),kind_real)
pminmax(1,1)=minval(self%u(:,:,:))
pminmax(2,1)=maxval(self%u(:,:,:))
pminmax(3,1)=sqrt(sum(self%u(:,:,:)**2)/zz)

zz=real(size(self%v),kind_real)
pminmax(1,2)=minval(self%v(:,:,:))
pminmax(2,2)=maxval(self%v(:,:,:))
pminmax(3,2)=sqrt(sum(self%v(:,:,:)**2)/zz)

zz=real(size(self%t),kind_real)
pminmax(1,3)=minval(self%t(:,:,:))
pminmax(2,3)=maxval(self%t(:,:,:))
pminmax(3,3)=sqrt(sum(self%t(:,:,:)**2)/zz)

zz=real(size(self%delp),kind_real)
pminmax(1,4)=minval(self%delp(:,:,:))
pminmax(2,4)=maxval(self%delp(:,:,:))
pminmax(3,4)=sqrt(sum(self%delp(:,:,:)**2)/zz)

zz=real(size(self%qv),kind_real)
pminmax(1,5)=minval(self%qv(:,:,:))
pminmax(2,5)=maxval(self%qv(:,:,:))
pminmax(3,5)=sqrt(sum(self%qv(:,:,:)**2)/zz)

zz=real(size(self%qi),kind_real)
pminmax(1,6)=minval(self%qi(:,:,:))
pminmax(2,6)=maxval(self%qi(:,:,:))
pminmax(3,6)=sqrt(sum(self%qi(:,:,:)**2)/zz)

zz=real(size(self%ql),kind_real)
pminmax(1,7)=minval(self%ql(:,:,:))
pminmax(2,7)=maxval(self%ql(:,:,:))
pminmax(3,7)=sqrt(sum(self%ql(:,:,:)**2)/zz)

zz=real(size(self%o3),kind_real)
pminmax(1,8)=minval(self%o3(:,:,:))
pminmax(2,8)=maxval(self%o3(:,:,:))
pminmax(3,8)=sqrt(sum(self%o3(:,:,:)**2)/zz)

if (allocated(self%w)) then
   zz=real(size(self%w),kind_real)
   pminmax(1,9)=minval(self%w(:,:,:))
   pminmax(2,9)=maxval(self%w(:,:,:))
   pminmax(3,9)=sqrt(sum(self%w(:,:,:)**2)/zz)
endif

if (allocated(self%delz)) then
   zz=real(size(self%delz),kind_real)
   pminmax(1,10)=minval(self%delz(:,:,:))
   pminmax(2,10)=maxval(self%delz(:,:,:))
   pminmax(3,10)=sqrt(sum(self%delz(:,:,:)**2)/zz)
endif

zz=real(size(self%phis),kind_real)
pminmax(1,11)=minval(self%phis(:,:))
pminmax(2,11)=maxval(self%phis(:,:))
pminmax(3,11)=sqrt(sum(self%phis(:,:)**2)/zz)

end subroutine traj_minmaxrms

! ------------------------------------------------------------------------------

end module fv3jedi_traj_mod
