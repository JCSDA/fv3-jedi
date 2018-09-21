! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_trajectories

use kinds

implicit none
private

public :: fv3jedi_trajectory, set_traj, get_traj, delete_traj
public :: fv3jedi_traj_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold the model trajectory
type :: fv3jedi_trajectory
  real(kind_real), allocatable, dimension(:,:,:) :: ud,vd,t,delp,w,delz
  real(kind_real), allocatable, dimension(:,:,:) :: q,qi,ql,o3
  real(kind_real), allocatable, dimension(:,:)   :: phis
end type fv3jedi_trajectory

#define LISTED_TYPE fv3jedi_trajectory

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_traj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine set_traj( self, &
                     isc,iec,jsc,jec, &
                     isd,ied,jsd,jed, &
                     npz, hydrostatic, &
                     ud,vd,t,delp,q,qi,ql,o3,w,delz,phis )

implicit none

type(fv3jedi_trajectory), intent(inout) :: self
integer, intent(in) :: isd,ied,jsd,jed
integer, intent(in) :: isc,iec,jsc,jec
integer, intent(in) :: npz
logical, intent(in) :: hydrostatic
real(kind_real), intent(in) ::   ud(isd:ied  ,jsd:jed+1,npz)
real(kind_real), intent(in) ::   vd(isd:ied+1,jsd:jed  ,npz)
real(kind_real), intent(in) ::    t(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) :: delp(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) ::    q(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) ::   qi(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) ::   ql(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) ::   o3(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) ::    w(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) :: delz(isd:ied  ,jsd:jed  ,npz)
real(kind_real), intent(in) :: phis(isd:ied  ,jsd:jed      )

allocate(self%ud  (isc:iec,jsc:jec,npz))
allocate(self%vd  (isc:iec,jsc:jec,npz))
allocate(self%t   (isc:iec,jsc:jec,npz))
allocate(self%delp(isc:iec,jsc:jec,npz))
allocate(self%q   (isc:iec,jsc:jec,npz))
allocate(self%qi  (isc:iec,jsc:jec,npz))
allocate(self%ql  (isc:iec,jsc:jec,npz))
allocate(self%o3  (isc:iec,jsc:jec,npz))
if (.not. hydrostatic) then
   allocate(self%w   (isc:iec, jsc:jec, npz))
   allocate(self%delz(isc:iec, jsc:jec, npz))
endif
allocate(self%phis(isc:iec, jsc:jec))


self%ud  (isc:iec,jsc:jec,:) = ud  (isc:iec,jsc:jec,:)
self%vd  (isc:iec,jsc:jec,:) = vd  (isc:iec,jsc:jec,:)
self%t   (isc:iec,jsc:jec,:) = t   (isc:iec,jsc:jec,:)
self%delp(isc:iec,jsc:jec,:) = delp(isc:iec,jsc:jec,:)
self%q   (isc:iec,jsc:jec,:) = q   (isc:iec,jsc:jec,:)
self%qi  (isc:iec,jsc:jec,:) = qi  (isc:iec,jsc:jec,:)
self%ql  (isc:iec,jsc:jec,:) = ql  (isc:iec,jsc:jec,:)
self%o3  (isc:iec,jsc:jec,:) = o3  (isc:iec,jsc:jec,:)
if (.not. hydrostatic) then
   self%delz(isc:iec,jsc:jec,:) = delz(isc:iec,jsc:jec,:)
   self%w   (isc:iec,jsc:jec,:) = w   (isc:iec,jsc:jec,:)
endif
self%phis(isc:iec,jsc:jec) = phis(isc:iec,jsc:jec)

end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine get_traj( self, &
                     isc,iec,jsc,jec, &
                     isd,ied,jsd,jed, &
                     npz, hydrostatic, &
                     ud,vd,t,delp,q,qi,ql,o3,w,delz,phis )

implicit none

type(fv3jedi_trajectory), intent(in) :: self
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

!Get compute domain from trajectory
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

end subroutine get_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(fv3jedi_trajectory), intent(inout) :: self

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

end subroutine delete_traj

! ------------------------------------------------------------------------------

subroutine c_minmax_traj(c_key_self, pminmax) bind(c,name='fv3jedi_traj_minmaxrms_f90')
use iso_c_binding
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: pminmax(3,7)
type(fv3jedi_trajectory), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_traj_registry%get(c_key_self,self)

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
pminmax(1,5)=minval(self%qi(:,:,:))
pminmax(2,5)=maxval(self%qi(:,:,:))
pminmax(3,5)=sqrt(sum(self%qi(:,:,:)**2)/zz)

zz=real(size(self%ql),kind_real)
pminmax(1,5)=minval(self%ql(:,:,:))
pminmax(2,5)=maxval(self%ql(:,:,:))
pminmax(3,5)=sqrt(sum(self%ql(:,:,:)**2)/zz)

zz=real(size(self%o3),kind_real)
pminmax(1,5)=minval(self%o3(:,:,:))
pminmax(2,5)=maxval(self%o3(:,:,:))
pminmax(3,5)=sqrt(sum(self%o3(:,:,:)**2)/zz)

if (allocated(self%w)) then
   zz=real(size(self%w),kind_real)
   pminmax(1,6)=minval(self%w(:,:,:))
   pminmax(2,6)=maxval(self%w(:,:,:))
   pminmax(3,6)=sqrt(sum(self%w(:,:,:)**2)/zz)
endif

if (allocated(self%delz)) then
   zz=real(size(self%delz),kind_real)
   pminmax(1,7)=minval(self%delz(:,:,:))
   pminmax(2,7)=maxval(self%delz(:,:,:))
   pminmax(3,7)=sqrt(sum(self%delz(:,:,:)**2)/zz)
endif

zz=real(size(self%phis),kind_real)
pminmax(1,5)=minval(self%phis(:,:))
pminmax(2,5)=maxval(self%phis(:,:))
pminmax(3,5)=sqrt(sum(self%phis(:,:)**2)/zz)

end subroutine c_minmax_traj

! ------------------------------------------------------------------------------

end module fv3jedi_trajectories
