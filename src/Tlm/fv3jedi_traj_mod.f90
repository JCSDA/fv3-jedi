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

!Allocate trajectory field only if corresponding state field exists

if (associated(state%ud)) then
  allocate(self%u      (isc:iec, jsc:jec, npz))
  self%u       = state%ud
endif
if (associated(state%vd)) then
  allocate(self%v      (isc:iec, jsc:jec, npz))
  self%v       = state%vd
endif
if (associated(state%ua)) then
  allocate(self%ua     (isc:iec, jsc:jec, npz))
  self%ua      = state%ua
else
  self%ua = 0.0_kind_real
endif
if (associated(state%va)) then
  allocate(self%va     (isc:iec, jsc:jec, npz))
  self%va      = state%va
else
  self%va = 0.0_kind_real
endif
if (associated(state%t)) then
  allocate(self%t      (isc:iec, jsc:jec, npz))
  self%t       = state%t
endif
if (associated(state%delp)) then
  allocate(self%delp   (isc:iec, jsc:jec, npz))
  self%delp    = state%delp
endif
if (associated(state%q)) then
  allocate(self%qv     (isc:iec, jsc:jec, npz))
  self%qv      = state%q
endif
if (associated(state%qi)) then
  allocate(self%qi     (isc:iec, jsc:jec, npz))
  self%qi      = state%qi
endif
if (associated(state%ql)) then
  allocate(self%ql     (isc:iec, jsc:jec, npz))
  self%ql      = state%ql
endif
if (associated(state%o3)) then
  allocate(self%o3     (isc:iec, jsc:jec, npz))
  self%o3      = state%o3
endif
if (associated(state%w)) then
  allocate(self%w      (isc:iec, jsc:jec, npz))
  self%w       = state%w
endif
if (associated(state%delz)) then
  allocate(self%delz   (isc:iec, jsc:jec, npz))
  self%delz    = state%delz
endif
if (associated(state%qls)) then
  allocate(self%qls    (isc:iec, jsc:jec, npz))
  self%qls     = state%qls
endif
if (associated(state%qcn)) then
  allocate(self%qcn    (isc:iec, jsc:jec, npz))
  self%qcn     = state%qcn
endif
if (associated(state%cfcn)) then
  allocate(self%cfcn   (isc:iec, jsc:jec, npz)) 
  self%cfcn    = state%cfcn
endif
if (associated(state%phis)) then
  allocate(self%phis   (isc:iec, jsc:jec))
  self%phis    = state%phis(:,:,1)
endif
if (associated(state%frocean)) then
  allocate(self%frocean(isc:iec, jsc:jec))
  self%frocean = state%frocean(:,:,1)
endif
if (associated(state%frland)) then
  allocate(self%frland (isc:iec, jsc:jec))
  self%frland  = state%frland(:,:,1)
endif
if (associated(state%varflt)) then
  allocate(self%varflt (isc:iec, jsc:jec))
  self%varflt  = state%varflt(:,:,1)
endif
if (associated(state%ustar)) then
  allocate(self%ustar  (isc:iec, jsc:jec))
  self%ustar   = state%ustar(:,:,1)
endif
if (associated(state%bstar)) then
  allocate(self%bstar  (isc:iec, jsc:jec))
  self%bstar   = state%bstar(:,:,1)
endif
if (associated(state%zpbl)) then
  allocate(self%zpbl   (isc:iec, jsc:jec))
  self%zpbl    = state%zpbl(:,:,1)
endif
if (associated(state%cm)) then
  allocate(self%cm     (isc:iec, jsc:jec))
  self%cm      = state%cm(:,:,1)
endif
if (associated(state%ct)) then
  allocate(self%ct     (isc:iec, jsc:jec))
  self%ct      = state%ct(:,:,1)
endif
if (associated(state%cq)) then
  allocate(self%cq     (isc:iec, jsc:jec))
  self%cq      = state%cq(:,:,1)
endif
if (associated(state%kcbl)) then
  allocate(self%kcbl   (isc:iec, jsc:jec))
  self%kcbl    = state%kcbl(:,:,1)
endif
if (associated(state%ts)) then
  allocate(self%ts     (isc:iec, jsc:jec))
  self%ts      = state%ts(:,:,1)
endif
if (associated(state%khl)) then
  allocate(self%khl    (isc:iec, jsc:jec))
  self%khl     = state%khl(:,:,1)
endif
if (associated(state%khu)) then
  allocate(self%khu    (isc:iec, jsc:jec))
  self%khu     = state%khu(:,:,1)
endif

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
