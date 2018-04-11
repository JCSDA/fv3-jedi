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
  real(kind_real), allocatable, dimension(:,:,:)   :: u,v,pt,delp,w,delz
  real(kind_real), allocatable, dimension(:,:,:,:) :: q
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

subroutine set_traj(self,isd,ied,jsd,jed,npz,nq,hydrostatic,u,v,pt,delp,q,w,delz)

implicit none

type(fv3jedi_trajectory), intent(inout) :: self
integer, intent(in) :: isd,ied,jsd,jed,npz,nq
logical, intent(in) :: hydrostatic
real(kind_real), intent(in) ::    u(isd:ied,jsd:jed,npz)
real(kind_real), intent(in) ::    v(isd:ied,jsd:jed,npz)
real(kind_real), intent(in) ::   pt(isd:ied,jsd:jed,npz)
real(kind_real), intent(in) :: delp(isd:ied,jsd:jed,npz)
real(kind_real), intent(in) ::    q(isd:ied,jsd:jed,npz,nq)
real(kind_real), intent(in) ::    w(isd:ied,jsd:jed,npz)
real(kind_real), intent(in) :: delz(isd:ied,jsd:jed,npz)

allocate(self%u   (isd:ied,jsd:jed,npz))
allocate(self%v   (isd:ied,jsd:jed,npz))
allocate(self%pt  (isd:ied,jsd:jed,npz))
allocate(self%delp(isd:ied,jsd:jed,npz))
allocate(self%q   (isd:ied,jsd:jed,npz,nq))
!if (.not. hydrostatic) then
   allocate(self%w   (isd:ied,jsd:jed,npz))
   allocate(self%delz(isd:ied,jsd:jed,npz))
!endif

self%u    = u
self%v    = v
self%pt   = pt
self%delp = delp
self%q    = q
!if (.not. hydrostatic) then
   self%w    = w
   self%delz = delz
!endif

end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine get_traj(self,isd,ied,jsd,jed,npz,hydrostatic,nq,u,v,pt,delp,q,w,delz)

implicit none

type(fv3jedi_trajectory), intent(in) :: self
integer, intent(in) :: isd,ied,jsd,jed,npz,nq
logical, intent(in) :: hydrostatic
real(kind_real), intent(inout) ::    u(isd:ied,jsd:jed,npz)
real(kind_real), intent(inout) ::    v(isd:ied,jsd:jed,npz)
real(kind_real), intent(inout) ::   pt(isd:ied,jsd:jed,npz)
real(kind_real), intent(inout) :: delp(isd:ied,jsd:jed,npz)
real(kind_real), intent(inout) ::    q(isd:ied,jsd:jed,npz,nq)
real(kind_real), intent(inout) ::    w(isd:ied,jsd:jed,npz)
real(kind_real), intent(inout) :: delz(isd:ied,jsd:jed,npz)

u    = self%u
v    = self%v
pt   = self%pt
delp = self%delp
q    = self%q
!if (.not. hydrostatic) then
   w    = self%w
   delz = self%delz
!endif

end subroutine get_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(fv3jedi_trajectory), intent(inout) :: self

deallocate(self%u   )
deallocate(self%v   )
deallocate(self%pt  )
deallocate(self%delp)
deallocate(self%q   )
if (allocated(self%w   )) deallocate(self%w   )
if (allocated(self%delz)) deallocate(self%delz)

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

zz=real(size(self%u),kind_real)
pminmax(1,1)=minval(self%u(:,:,:))
pminmax(2,1)=maxval(self%u(:,:,:))
pminmax(3,1)=sqrt(sum(self%u(:,:,:)**2)/zz)

zz=real(size(self%v),kind_real)
pminmax(1,2)=minval(self%v(:,:,:))
pminmax(2,2)=maxval(self%v(:,:,:))
pminmax(3,2)=sqrt(sum(self%v(:,:,:)**2)/zz)

zz=real(size(self%pt),kind_real)
pminmax(1,3)=minval(self%pt(:,:,:))
pminmax(2,3)=maxval(self%pt(:,:,:))
pminmax(3,3)=sqrt(sum(self%pt(:,:,:)**2)/zz)

zz=real(size(self%delp),kind_real)
pminmax(1,4)=minval(self%delp(:,:,:))
pminmax(2,4)=maxval(self%delp(:,:,:))
pminmax(3,4)=sqrt(sum(self%delp(:,:,:)**2)/zz)

zz=real(size(self%q),kind_real)
pminmax(1,5)=minval(self%q(:,:,:,:))
pminmax(2,5)=maxval(self%q(:,:,:,:))
pminmax(3,5)=sqrt(sum(self%q(:,:,:,:)**2)/zz)

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

end subroutine c_minmax_traj

! ------------------------------------------------------------------------------

end module fv3jedi_trajectories
