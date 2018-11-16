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

call allocate_traj(self,state%isc,state%iec,state%jsc,state%jec,state%npz,&
                   state%hydrostatic,state%do_tlad_phymst)

self%u       = state%ud
self%v       = state%vd
self%ua      = state%ua
self%va      = state%va
self%t       = state%t
self%delp    = state%delp
self%qv      = state%q
self%ql      = state%qi
self%qi      = state%ql
self%o3      = state%o3

if (.not. state%hydrostatic) then
self%w       = state%w
self%delz    = state%delz
endif

if (state%do_tlad_phymst /= 0) then
self%qls     = state%qls
self%qcn     = state%qcn
self%cfcn    = state%cfcn
endif

self%phis    = state%phis
self%frocean = state%frocean
self%frland  = state%frland
self%varflt  = state%varflt
self%ustar   = state%ustar
self%bstar   = state%bstar
self%zpbl    = state%zpbl
self%cm      = state%cm
self%ct      = state%ct
self%cq      = state%cq
self%kcbl    = state%kcbl
self%ts      = state%ts
self%khl     = state%khl
self%khu     = state%khu

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
