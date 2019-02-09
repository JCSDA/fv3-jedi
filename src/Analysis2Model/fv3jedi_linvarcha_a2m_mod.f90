! (C) Copyright 2018 UCAR
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

use wind_vt_mod, only: a2d, d2a, d2a_ad, a2d_ad

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
end type fv3jedi_linvarcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, bg, fg, geom, c_conf)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_state),         intent(in)    :: bg
type(fv3jedi_state),         intent(in)    :: fg
type(fv3jedi_geom),          intent(in)    :: geom
type(c_ptr),                 intent(in)    :: c_conf

! Allocate temporary A- and D-Grid wind holders
allocate(self%ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz))
allocate(self%va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz))
allocate(self%ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
allocate(self%vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self

! Deallocate temporary A- and D-grid wind holders
deallocate(self%ua,self%va)
deallocate(self%ud,self%vd)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

integer :: k

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real

call a2d(geom, xana%fields(xana%ua)%array(xana%isc:xana%iec,xana%jsc:xana%jec,:), &
               xana%fields(xana%va)%array(xana%isc:xana%iec,xana%jsc:xana%jec,:), &
               ud(xana%isc:xana%iec  ,xana%jsc:xana%jec+1,:), &
               vd(xana%isc:xana%iec+1,xana%jsc:xana%jec  ,:))

xmod%fields(xana%ud)%array(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%ud(xana%isc:xana%iec,xana%jsc:xana%jec,:)
xmod%fields(xana%vd)%array(xana%isc:xana%iec,xana%jsc:xana%jec,:) = self%vd(xana%isc:xana%iec,xana%jsc:xana%jec,:)

xmod%fields(xana%ua)%array   = xana%fields(xana%ua)%array
xmod%fields(xana%va)%array   = xana%fields(xana%va)%array

xmod%fields(xana%t)%array    = xana%fields(xana%t)%array

do k = 1,geom%npz
  xmod%fields(xana%delp)%array(:,:,k) = (geom%bk(k+1)-geom%bk(k))*xana%fields(xana%ps)%array(:,:,1)
enddo

xmod%fields(xana%qv)%array   = xana%fields(xana%q)%array
xmod%fields(xana%qi)%array   = xana%fields(xana%qi)%array
xmod%fields(xana%ql)%array   = xana%fields(xana%ql)%array
xmod%fields(xana%o3)%array   = xana%fields(xana%o3)%array

if (.not. xana%hydrostatic) then
   xmod%fields(xana%%delz)%array = xana%fields(xana%delz)%array
   xmod%fields(xana%%w)%array    = xana%fields(xana%w)%array
endif

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

self%ud = 0.0_kind_real
self%vd = 0.0_kind_real

xana%fields(xana%ua)%array   = 0.0
xana%fields(xana%va)%array   = 0.0
xana%fields(xana%t )%array   = 0.0
xana%fields(xana%ps)%array   = 0.0
xana%fields(xana%q )%array   = 0.0
xana%fields(xana%qi)%array   = 0.0
xana%fields(xana%ql)%array   = 0.0
xana%fields(xana%o3)%array   = 0.0

if (.not. xana%hydrostatic) then
   xana%fields(xana%delz)%array = 0.0
   xana%fields(xana%w)%array    = 0.0
endif

ud(xana%isc:xana%iec,xana%jsc:xana%jec,:) = lm%pert%u
vd(xana%isc:xana%iec,xana%jsc:xana%jec,:) = lm%pert%v
xana%fields(xana%t)%array    = lm%pert%t
xana%fields(xana%ps)%array = 0.0_kind_real
do k = 1,geom%npz
  xana%fields(xana%ps)%array(:,:,1) = xana%fields(xana%ps)%array(:,:,1) +  (geom%bk(k+1)-geom%bk(k))*lm%pert%delp(:,:,k)
enddo
xana%fields(xana%q)%array    = lm%pert%qv
xana%fields(xana%qi)%array   = lm%pert%qi
xana%fields(xana%ql)%array   = lm%pert%ql
xana%fields(xana%o3)%array   = lm%pert%o3
if (.not. xana%hydrostatic) then
  xana%fields(xana%delz)%array = lm%pert%delz
  xana%fields(xana%q)%array    = lm%pert%w
endif

!Convert A to D
call a2d_ad(geom, xana%fields(xana%ua)%array, xana%fields(xana%va)%array, self%ud, self%vd)

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

ud = 0.0_kind_real
vd = 0.0_kind_real
ua = 0.0_kind_real
va = 0.0_kind_real

ud(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%u
vd(inc%isc:inc%iec,inc%jsc:inc%jec,:) = lm%pert%v

call d2a(geom, ud, vd, ua, va)

inc%fields(inc%ua)%array   = ua(inc%isc:inc%iec,inc%jsc:inc%jec,:)
inc%fields(inc%va)%array   = va(inc%isc:inc%iec,inc%jsc:inc%jec,:)
inc%fields(inc%t)%array    = lm%pert%t
inc%fields(inc%ps)%array(:,:,1)   = sum(lm%pert%delp,3)
inc%fields(inc%q)%array    = lm%pert%qv
inc%fields(inc%qi)%array   = lm%pert%qi
inc%fields(inc%ql)%array   = lm%pert%ql
inc%fields(inc%o3)%array   = lm%pert%o3
if (.not. inc%hydrostatic) then
  inc%fields(inc%delz)%array = lm%pert%delz
  inc%fields(inc%w)%array    = lm%pert%w
endif

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

ud = 0.0_kind_real
vd = 0.0_kind_real

ua = 0.0_kind_real
va = 0.0_kind_real

ua(inc%isc:inc%iec,inc%jsc:inc%jec,:) = inc%fields(inc%ua)%array
va(inc%isc:inc%iec,inc%jsc:inc%jec,:) = inc%fields(inc%va)%array

lm%pert%t    = inc%fields(inc%t)%array
do k = 1,geom%npz
  lm%pert%delp(:,:,k) = inc%fields(inc%ps)%array(:,:,1)
enddo
lm%pert%qv   = inc%fields(inc%q)%array
lm%pert%qi   = inc%fields(inc%qi)%array
lm%pert%ql   = inc%fields(inc%ql)%array
lm%pert%o3   = inc%fields(inc%o3)%array
if (.not. inc%hydrostatic) then
   lm%pert%delz = inc%fields(inc%delz)%array
   lm%pert%w    = inc%fields(inc%w)%array
endif

call d2a_ad(geom, ud, vd, ua, va)

lm%pert%u = ud(inc%isc:inc%iec,inc%jsc:inc%jec,:)
lm%pert%v = vd(inc%isc:inc%iec,inc%jsc:inc%jec,:)

lm%pert%ua(inc%isc:inc%iec,inc%jsc:inc%jec,:) = 0.0
lm%pert%va(inc%isc:inc%iec,inc%jsc:inc%jec,:) = 0.0

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_mod
