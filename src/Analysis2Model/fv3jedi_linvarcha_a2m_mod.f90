! (C) Copyright 2018-2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_linvarcha_a2m_mod

use iso_c_binding
use config_mod
use fv3jedi_kinds_mod

use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment

use wind_vt_mod,     only: a2d, d2a, d2a_ad, a2d_ad
use mpp_domains_mod, only: mpp_get_boundary, DGRID_NE, mpp_get_boundary_ad

implicit none
private

public :: fv3jedi_linvarcha_a2m
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

type :: fv3jedi_linvarcha_a2m
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ua,va
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd
  real(kind=kind_real), allocatable, dimension(:,:)   :: ebuffery
  real(kind=kind_real), allocatable, dimension(:,:)   :: nbufferx
  real(kind=kind_real), allocatable, dimension(:,:)   :: wbuffery
  real(kind=kind_real), allocatable, dimension(:,:)   :: sbufferx
end type fv3jedi_linvarcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, c_conf)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_state),         intent(in)    :: bg
type(fv3jedi_state),         intent(in)    :: fg
type(c_ptr),                 intent(in)    :: c_conf

! Allocate A- and D-Grid wind holders
allocate(self%ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz))
allocate(self%va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz))
allocate(self%ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz))
allocate(self%vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz))

! Allocate D-Grid wind halo holders
allocate(self%wbuffery(geom%jsc:geom%jec,1:geom%npz))
allocate(self%sbufferx(geom%isc:geom%iec,1:geom%npz))
allocate(self%ebuffery(geom%jsc:geom%jec,1:geom%npz))
allocate(self%nbufferx(geom%isc:geom%iec,1:geom%npz))

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self

! Deallocate A- and D-grid wind holders
deallocate(self%ua,self%va)
deallocate(self%ud,self%vd)

! Deallocate D-Grid wind halo holders
deallocate(self%wbuffery,self%sbufferx,self%ebuffery,self%nbufferx)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

integer :: k

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real

call a2d(geom, xana%ua(xana%isc:xana%iec,xana%jsc:xana%jec,:), &
               xana%va(xana%isc:xana%iec,xana%jsc:xana%jec,:), &
               self%ud(xana%isc:xana%iec  ,xana%jsc:xana%jec+1,:), &
               self%vd(xana%isc:xana%iec+1,xana%jsc:xana%jec  ,:))

xmod%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:)
xmod%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:)

xmod%ua = 0.0_kind_real
xmod%va = 0.0_kind_real

xmod%t  = xana%t
do k = 1,geom%npz
  xmod%delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*xana%ps(:,:,1)
enddo

xmod%q  = xana%q
xmod%qi = xana%qi
xmod%ql = xana%ql
xmod%o3 = xana%o3

if (.not. xana%hydrostatic) then
  xmod%delz = xana%delz
  xmod%w    = xana%w
endif

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

integer :: k

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real

self%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xmod%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:)
self%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xmod%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:)

call a2d_ad(geom, xana%ua, xana%va, &
            self%ud(xana%isc:xana%iec  ,xana%jsc:xana%jec+1,:),&
            self%vd(xana%isc:xana%iec+1,xana%jsc:xana%jec  ,:))

xana%t    = xmod%t
xana%ps = 0.0_kind_real
do k = 1,geom%npz
  xana%ps(:,:,1) = xana%ps(:,:,1) +  (geom%bk(k+1)-geom%bk(k))*xmod%delp(:,:,k)
enddo

xana%q    = xmod%q
xana%qi   = xmod%qi
xana%ql   = xmod%ql
xana%o3   = xmod%o3

if (.not. xana%hydrostatic) then
  xana%delz = xmod%delz
  xana%q    = xmod%w
endif

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

integer :: i,j,k

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real
self%ua = 0.0_kind_real
self%va = 0.0_kind_real

self%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xmod%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:)
self%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xmod%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:)

call mpp_get_boundary( self%ud, self%vd, geom%domain, &
                       wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                       sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                       gridtype=DGRID_NE, complete=.true. )
do k=1,geom%npz
   do i=geom%isc,geom%iec
      self%ud(i,geom%jec+1,k) = self%nbufferx(i,k)
   enddo
enddo
do k=1,geom%npz
   do j=geom%jsc,geom%jec
      self%vd(geom%iec+1,j,k) = self%ebuffery(j,k)
   enddo
enddo

call d2a(geom, self%ud, self%vd, self%ua, self%va)

xana%ua = self%ua(xana%isc:xana%iec,xana%jsc:xana%jec,:)
xana%va = self%va(xana%isc:xana%iec,xana%jsc:xana%jec,:)

xana%t    = xmod%t
xana%ps(:,:,1)   = sum(xmod%delp,3)
xana%q    = xmod%q
xana%qi   = xmod%qi
xana%ql   = xmod%ql
xana%o3   = xmod%o3
if (.not. xana%hydrostatic) then
  xana%delz = xmod%delz
  xana%w    = xmod%w
endif

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

integer :: i,j,k

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real
self%ua = 0.0_kind_real
self%va = 0.0_kind_real

self%ua(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xana%ua
self%va(xana%isc:xana%iec,xana%jsc:xana%jec,:) = xana%va

call d2a_ad(geom, self%ud, self%vd, self%ua, self%va)

self%nbufferx = 0.0_kind_real
do k=1,geom%npz
   do i=geom%isc,geom%iec
      self%nbufferx(i,k) = self%ud(i,geom%jec+1,k)
   enddo
enddo
self%ebuffery = 0.0_kind_real
do k=1,geom%npz
   do j=geom%jsc,geom%jec
      self%ebuffery(j,k) = self%vd(geom%iec+1,j,k)
   enddo
enddo

call mpp_get_boundary_ad( self%ud, self%vd, geom%domain, &
                          wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                          sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                          gridtype=DGRID_NE, complete=.true. )

xmod%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:)
xmod%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:)

xmod%ua = 0.0_kind_real
xmod%va = 0.0_kind_real

xmod%t    = xana%t
do k = 1,geom%npz
  xmod%delp(:,:,k) = xana%ps(:,:,1)
enddo
xmod%q    = xana%q
xmod%qi   = xana%qi
xmod%ql   = xana%ql
xmod%o3   = xana%o3
if (.not. xana%hydrostatic) then
   xmod%delz = xana%delz
   xmod%w    = xana%w
endif

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_mod
