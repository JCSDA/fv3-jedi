! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_traj_mod

use kinds
use iso_c_binding
use fv3jedi_model_mod, only: fv3jedi_model
use fv3jedi_state_mod, only: fv3jedi_state

implicit none
private

public :: fv3jedi_traj
public :: traj_prop
public :: traj_wipe
public :: traj_get
public :: traj_minmaxrms

! ------------------------------------------------------------------------------

!> Fortran derived type to hold the TL/AD traj
type :: fv3jedi_traj
  real(kind_real), allocatable, dimension(:,:,:) :: ud,vd,t,delp,w,delz
  real(kind_real), allocatable, dimension(:,:,:) :: q,qi,ql,o3
  real(kind_real), allocatable, dimension(:,:)   :: phis
end type fv3jedi_traj

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine traj_prop(model, state, self)

implicit none
type(fv3jedi_model) :: model
type(fv3jedi_state) :: state
type(fv3jedi_traj)  :: self

integer :: isc, iec, jsc, jec, npz

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
npz = state%npz

allocate(self%ud  (isc:iec,jsc:jec,npz))
allocate(self%vd  (isc:iec,jsc:jec,npz))
allocate(self%t   (isc:iec,jsc:jec,npz))
allocate(self%delp(isc:iec,jsc:jec,npz))
allocate(self%q   (isc:iec,jsc:jec,npz))
allocate(self%qi  (isc:iec,jsc:jec,npz))
allocate(self%ql  (isc:iec,jsc:jec,npz))
allocate(self%o3  (isc:iec,jsc:jec,npz))
if (.not. state%hydrostatic) then
   allocate(self%w   (isc:iec, jsc:jec, npz))
   allocate(self%delz(isc:iec, jsc:jec, npz))
endif
allocate(self%phis(isc:iec, jsc:jec))

self%ud  (isc:iec,jsc:jec,:) = state%ud  (isc:iec,jsc:jec,:)
self%vd  (isc:iec,jsc:jec,:) = state%vd  (isc:iec,jsc:jec,:)
self%t   (isc:iec,jsc:jec,:) = state%t   (isc:iec,jsc:jec,:)
self%delp(isc:iec,jsc:jec,:) = state%delp(isc:iec,jsc:jec,:)
self%q   (isc:iec,jsc:jec,:) = state%q   (isc:iec,jsc:jec,:)
self%qi  (isc:iec,jsc:jec,:) = state%qi  (isc:iec,jsc:jec,:)
self%ql  (isc:iec,jsc:jec,:) = state%ql  (isc:iec,jsc:jec,:)
self%o3  (isc:iec,jsc:jec,:) = state%o3  (isc:iec,jsc:jec,:)
if (.not.state%hydrostatic) then
   self%delz(isc:iec,jsc:jec,:) = state%delz(isc:iec,jsc:jec,:)
   self%w   (isc:iec,jsc:jec,:) = state%w   (isc:iec,jsc:jec,:)
endif
self%phis(isc:iec,jsc:jec) = state%phis(isc:iec,jsc:jec)

end subroutine traj_prop

! ------------------------------------------------------------------------------

subroutine traj_wipe(self)

implicit none
type(fv3jedi_traj), pointer :: self

if (allocated(self%ud  )) deallocate(self%ud  )
if (allocated(self%vd  )) deallocate(self%vd  )
if (allocated(self%t   )) deallocate(self%t   )
if (allocated(self%delp)) deallocate(self%delp)
if (allocated(self%q   )) deallocate(self%q   )
if (allocated(self%qi  )) deallocate(self%qi  )
if (allocated(self%ql  )) deallocate(self%ql  )
if (allocated(self%o3  )) deallocate(self%o3  )
if (allocated(self%w   )) deallocate(self%w   )
if (allocated(self%delz)) deallocate(self%delz)
if (allocated(self%phis)) deallocate(self%phis)

end subroutine traj_wipe

! ------------------------------------------------------------------------------

subroutine traj_get( self, &
                     isc,iec,jsc,jec, &
                     isd,ied,jsd,jed, &
                     npz, hydrostatic, &
                     ud,vd,t,delp,q,qi,ql,o3,w,delz,phis )

implicit none

type(fv3jedi_traj), intent(in) :: self
integer, intent(in) :: isd,ied,jsd,jed
integer, intent(in) :: isc,iec,jsc,jec
integer, intent(in) :: npz
logical, intent(in) :: hydrostatic
real(kind_real), intent(out) ::   ud(isd:ied  ,jsd:jed+1,npz)
real(kind_real), intent(out) ::   vd(isd:ied+1,jsd:jed  ,npz)
real(kind_real), intent(out) ::    t(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) :: delp(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) ::    q(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) ::   qi(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) ::   ql(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) ::   o3(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) ::    w(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) :: delz(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(out) :: phis(isd:ied  ,jsd:jed )

!Initialize to zero, including halo
ud   = 0.0
vd   = 0.0
t    = 0.0
delp = 0.0
q    = 0.0
qi   = 0.0
ql   = 0.0
o3   = 0.0
delz = 0.0
w    = 0.0
phis = 0.0

!Get compute domain from traj
ud  (isc:iec,jsc:jec,:) = self%ud  (isc:iec,jsc:jec,:)
vd  (isc:iec,jsc:jec,:) = self%vd  (isc:iec,jsc:jec,:)
t   (isc:iec,jsc:jec,:) = self%t   (isc:iec,jsc:jec,:)
delp(isc:iec,jsc:jec,:) = self%delp(isc:iec,jsc:jec,:)
q   (isc:iec,jsc:jec,:) = self%q   (isc:iec,jsc:jec,:)
qi  (isc:iec,jsc:jec,:) = self%qi  (isc:iec,jsc:jec,:)
ql  (isc:iec,jsc:jec,:) = self%ql  (isc:iec,jsc:jec,:)
o3  (isc:iec,jsc:jec,:) = self%o3  (isc:iec,jsc:jec,:)
if (.not. hydrostatic) then
   delz(isc:iec,jsc:jec,:) = self%delz(isc:iec,jsc:jec,:)
   w   (isc:iec,jsc:jec,:) = self%w   (isc:iec,jsc:jec,:)
endif
phis(isc:iec,jsc:jec) = self%phis(isc:iec,jsc:jec)

end subroutine traj_get

! ------------------------------------------------------------------------------

subroutine traj_minmaxrms(self, pminmax)

implicit none
type(fv3jedi_traj), intent(in)    :: self
real(c_double),     intent(inout) :: pminmax(3,11)

real(kind=kind_real) :: zz

pminmax = 0.0_kind_real

zz=real(size(self%ud),kind_real)
pminmax(1,1)=minval(self%ud(:,:,:))
pminmax(2,1)=maxval(self%ud(:,:,:))
pminmax(3,1)=sqrt(sum(self%ud(:,:,:)**2)/zz)

zz=real(size(self%vd),kind_real)
pminmax(1,2)=minval(self%vd(:,:,:))
pminmax(2,2)=maxval(self%vd(:,:,:))
pminmax(3,2)=sqrt(sum(self%vd(:,:,:)**2)/zz)

zz=real(size(self%t),kind_real)
pminmax(1,3)=minval(self%t(:,:,:))
pminmax(2,3)=maxval(self%t(:,:,:))
pminmax(3,3)=sqrt(sum(self%t(:,:,:)**2)/zz)

zz=real(size(self%delp),kind_real)
pminmax(1,4)=minval(self%delp(:,:,:))
pminmax(2,4)=maxval(self%delp(:,:,:))
pminmax(3,4)=sqrt(sum(self%delp(:,:,:)**2)/zz)

zz=real(size(self%q),kind_real)
pminmax(1,5)=minval(self%q(:,:,:))
pminmax(2,5)=maxval(self%q(:,:,:))
pminmax(3,5)=sqrt(sum(self%q(:,:,:)**2)/zz)

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
