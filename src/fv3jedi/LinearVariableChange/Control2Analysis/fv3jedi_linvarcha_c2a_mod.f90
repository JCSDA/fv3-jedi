! (C) Copyright 2018-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_linvarcha_c2a_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! fv3jedi
use fv3jedi_constants_mod, only: constoz
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_field_mod,     only: copy_subset, field_clen
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state

use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use wind_vt_mod

implicit none
private

public :: fv3jedi_linvarcha_c2a
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_linvarcha_c2a
 real(kind=kind_real), allocatable :: ttraj(:,:,:)
 real(kind=kind_real), allocatable :: tvtraj(:,:,:)
 real(kind=kind_real), allocatable :: qtraj(:,:,:)
 real(kind=kind_real), allocatable :: qsattraj(:,:,:)
end type fv3jedi_linvarcha_c2a

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, conf)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(fv3jedi_state), target, intent(in)    :: bg
type(fv3jedi_state), target, intent(in)    :: fg
type(fckit_configuration),   intent(in)    :: conf

real(kind=kind_real), pointer :: t   (:,:,:)=>NULL()
real(kind=kind_real), pointer :: tv  (:,:,:)=>NULL()
real(kind=kind_real), pointer :: q   (:,:,:)=>NULL()
real(kind=kind_real), allocatable :: delp(:,:,:)
real(kind=kind_real), pointer :: ps  (:,:,:)=>NULL()

!> Pointers to the background state
if ( bg%has_field('t') ) then
  call bg%get_field('t', t)
endif
if ( bg%has_field('sphum') ) then
  call bg%get_field('sphum', q)
endif

!> Pressure
if (bg%has_field('delp')) then
  call bg%get_field('delp', delp)
elseif (bg%has_field('ps')) then
  call bg%get_field('ps', ps)
  allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call ps_to_delp(geom, ps, delp)
else
  call abor1_ftn("fv3jedi_linvarcha_c2a_mod.create : delp or ps should be present")
endif

!> Virtual temperature trajectory
allocate(self%tvtraj(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
if ( bg%has_field('tv')) then
  call bg%get_field('tv', tv)
  self%tvtraj = tv
else
  if (associated(t).and.associated(q)) then
    call T_to_Tv(geom,t,q,self%tvtraj)
  endif
endif

!> Temperature trajectory
if(associated(t)) then
  allocate(self%ttraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  self%ttraj = t
endif

!> Specific humidity trajecotory
if(associated(t)) then
  allocate(self%qtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  self%qtraj = q
endif

if(associated(t).and.associated(q).and.allocated(delp)) then
  !> Compute saturation specific humidity for q to RH transform
  allocate(self%qsattraj(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

  !> Compute saturation specific humidity
  call get_qsat(geom,delp,t,q,self%qsattraj)
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self

if (allocated(self%tvtraj)) deallocate(self%tvtraj)
if (allocated(self%ttraj)) deallocate(self%ttraj)
if (allocated(self%qtraj)) deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, geom, dxc, dxa)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: dxc
type(fv3jedi_increment),     intent(inout) :: dxa

integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Winds
logical :: have_uava, have_udvd
real(kind=kind_real), allocatable, dimension(:,:,:) :: psi
real(kind=kind_real), allocatable, dimension(:,:,:) :: chi
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud
real(kind=kind_real), allocatable, dimension(:,:,:) :: vd
real(kind=kind_real), allocatable, dimension(:,:,:) :: ua
real(kind=kind_real), allocatable, dimension(:,:,:) :: va

! Specific humidity
logical :: have_q
real(kind=kind_real), pointer,     dimension(:,:,:) :: q

! Temperature
logical :: have_t
real(kind=kind_real), pointer,     dimension(:,:,:) :: tv
real(kind=kind_real), allocatable, dimension(:,:,:) :: t

! Surface Pressure
logical :: have_ps
real(kind=kind_real), allocatable, dimension(:,:,:) :: ps

! Ozone
logical :: have_o3ppmv
logical :: have_o3mr
real(kind=kind_real), pointer,     dimension(:,:,:) :: o3ctl
real(kind=kind_real), allocatable, dimension(:,:,:) :: o3ana

! Identity part of the change of fields
! -------------------------------------
call copy_subset(dxc%fields, dxa%fields, fields_to_do)

! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return

! Winds
! -----
have_udvd = .false.
if (dxc%has_field('psi') .and. dxc%has_field('chi')) then
  call dxc%get_field('psi', psi)
  call dxc%get_field('chi', chi)
  allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,geom%npz))
  allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz))
  call psichi_to_udvd(geom, psi, chi, ud, vd)
  have_udvd = .true.
endif

! A-Grid winds
! ------------
have_uava = .false.
if (have_udvd) then
  allocate(ua(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
  allocate(va(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
  call d_to_a(geom, ud, vd, ua, va)
  have_uava=.true.
endif

! Temperature
! -----------
have_t = .false.
have_q = dxc%has_field('sphum')
if (dxc%has_field('t')) then
  allocate(t(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
  call dxc%get_field('t', t)
  have_t = .true.
elseif (dxc%has_field('tv') .and. have_q) then
  allocate(t(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
  call dxc%get_field('sphum', q)
  call dxc%get_field('tv'  , tv)
  call Tv_to_T_tl(geom, self%tvtraj, tv, self%qtraj, q, t)
  have_t = .true.
endif

! Surface pressure
! ----------------
have_ps = .false.
if (dxc%has_field('ps')) then
  call dxc%get_field('ps', ps)
! allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
! call ps_to_delp_tl(geom, ps, delp)
  have_ps = .true.
endif

! Ozone
! -----
have_o3mr = .false.
have_o3ppmv = .false.
if (dxc%has_field('o3ppmv').and.dxa%has_field('o3mr')) then
   call dxc%get_field('o3ppmv', o3ctl)
   allocate(o3ana(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
   o3ana = o3ctl / constoz
   have_o3mr=.true.
endif
if (dxc%has_field('o3mr').and.dxa%has_field('o3ppmv')) then
   call dxc%get_field('o3mr', o3ctl)
   allocate(o3ana(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
   o3ana = o3ctl * constoz
   have_o3ppmv=.true.
endif

! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call dxa%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ps")

    if (.not. have_ps) call field_fail('tl_'//fields_to_do(f))
    field_ptr = ps

  case ("ud")

    if (.not. have_uava) call field_fail('tl_'//fields_to_do(f))
    field_ptr = ud

  case ("vd")

    if (.not. have_uava) call field_fail('tl_'//fields_to_do(f))
    field_ptr = vd

  case ("ua")

    if (.not. have_uava) call field_fail('tl_'//fields_to_do(f))
    field_ptr = ua

  case ("va")

    if (.not. have_uava) call field_fail('tl_'//fields_to_do(f))
    field_ptr = va

  case ("t")

    if (.not. have_t) call field_fail('tl_'//fields_to_do(f))
    field_ptr = t

  case ("o3mr")

    if (.not. have_o3mr) call field_fail('tl_'//fields_to_do(f))
    field_ptr = o3ana

  case ("o3ppmv")

    if (.not. have_o3ppmv) call field_fail('tl_'//fields_to_do(f))
    field_ptr = o3ana

  case default

    call abor1_ftn("fv3jedi_linvarcha_c2a_mod.multiply unknown field: "//trim(fields_to_do(f)) &
                   //". Not in input field and no transform case specified.")

  end select

enddo

if(allocated(psi)) deallocate(psi)
if(allocated(chi)) deallocate(chi)
if(allocated(ud )) deallocate(ud )
if(allocated(vd )) deallocate(vd )
if(allocated(ua )) deallocate(ua )
if(allocated(va )) deallocate(va )
if(allocated(t  )) deallocate(t  )
if(allocated(ps )) deallocate(ps )
if(allocated(o3ana)) deallocate(o3ana)

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,dxa,dxc)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(inout) :: dxa
type(fv3jedi_increment),     intent(inout) :: dxc

integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Winds
logical :: have_psichi, have_udvd, have_uava
real(kind=kind_real), pointer,     dimension(:,:,:) :: ua
real(kind=kind_real), pointer,     dimension(:,:,:) :: va
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud
real(kind=kind_real), allocatable, dimension(:,:,:) :: vd
real(kind=kind_real), allocatable, dimension(:,:,:) :: psi
real(kind=kind_real), allocatable, dimension(:,:,:) :: chi

! Specfic humidity
real(kind=kind_real), pointer,     dimension(:,:,:) :: q

! Virtual temperature
logical :: have_tv
real(kind=kind_real), pointer,     dimension(:,:,:) :: t
real(kind=kind_real), allocatable, dimension(:,:,:) :: tv

! Surface pressure
logical :: have_ps
real(kind=kind_real), pointer,     dimension(:,:,:) :: ps

! Ozone
logical :: have_o3mr,have_o3ppmv
real(kind=kind_real), pointer,     dimension(:,:,:) :: o3ana
real(kind=kind_real), allocatable, dimension(:,:,:) :: o3ctl


! Zero output
! -----------
call dxc%zero()

! Identity part of the change of fields
! -------------------------------------
call copy_subset(dxa%fields, dxc%fields, fields_to_do)

! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return

! Virtual temperature
! -------------------
have_tv = .false.
if (dxa%has_field('t') .and. dxa%has_field('sphum')) then
  call dxa%get_field('t', t)
  call dxa%get_field('sphum', q)
  allocate(tv(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
  tv = 0.0_kind_real
  call Tv_to_T_ad(geom, self%tvtraj, tv, self%qtraj, q, t)
  have_tv = .true.
endif

! A-Grid winds
! ------------
have_udvd = .false.
have_uava = .false.
if (dxa%has_field('ud') .and. dxa%has_field('vd')) then
  call dxa%get_field('ud', ud)
  call dxa%get_field('vd', vd)
  have_udvd = .true.
elseif (dxa%has_field('ua') .and. dxa%has_field('va')) then
  call dxa%get_field('ua', ua)
  call dxa%get_field('va', va)
  allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
  allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
  ud = 0.0_kind_real
  vd = 0.0_kind_real
  call d_to_a_ad(geom, ud, vd, ua, va)
  have_udvd = .true.
  have_uava = .true.
endif

! Winds
! -----
have_psichi = .false.
if (have_udvd) then
  allocate(psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  psi = 0.0_kind_real
  chi = 0.0_kind_real
  call psichi_to_udvd_adm(geom, psi, chi, ud, vd)
  have_psichi = .true.
endif

! Surface pressure
! ----------------
have_ps = .false.
if (dxa%has_field('ps')) then
  call dxa%get_field('ps', ps)
endif

! Ozone
! -----
have_o3mr = .false.
have_o3ppmv = .false.
if (dxc%has_field('o3ppmv').and.dxa%has_field('o3mr')) then
   call dxa%get_field('o3mr', o3ana)
   allocate(o3ctl(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
   o3ctl = o3ana * constoz
   have_o3ppmv = .true.
endif
if (dxc%has_field('o3mr').and.dxa%has_field('o3ppmv')) then
   call dxa%get_field('o3ppmv', o3ana)
   allocate(o3ctl(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
   o3ctl = o3ana / constoz
   have_o3mr = .true.
endif

! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call dxc%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ps")

    if (.not. have_ps) call field_fail('ad_'//fields_to_do(f))
    field_ptr = ps

  case ("ua")

    if (.not. have_uava) call field_fail('ad_'//fields_to_do(f))
    field_ptr = ua

  case ("va")

    if (.not. have_uava) call field_fail('ad_'//fields_to_do(f))
    field_ptr = va

  case ("psi")

    if (.not. have_psichi) call field_fail('ad_'//fields_to_do(f))
    field_ptr = psi

  case ("chi")

    if (.not. have_psichi) call field_fail('ad_'//fields_to_do(f))
    field_ptr = chi

  case ("tv")

    if (.not. have_tv) call field_fail('ad_'//fields_to_do(f))
    field_ptr = tv

  case ("o3mr")

    if (.not. have_o3mr) call field_fail('ad_'//fields_to_do(f))
    field_ptr = o3ctl

  case ("o3ppmv")

    if (.not. have_o3ppmv) call field_fail('ad_'//fields_to_do(f))
    field_ptr = o3ctl

  case default

    call abor1_ftn("fv3jedi_linvarcha_c2a_mod.multiplyadjoint unknown field: "//trim(fields_to_do(f)) &
                   //". Not in input field and no transform case specified.")

  end select

enddo

if(allocated(psi)) deallocate(psi)
if(allocated(chi)) deallocate(chi)
if(allocated(ud )) deallocate(ud )
if(allocated(vd )) deallocate(vd )
if(allocated(tv )) deallocate(tv )
if(allocated(o3ctl)) deallocate(o3ctl)


end subroutine multiplyadjoint

! --------------------------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,dxa,dxc)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: dxa
type(fv3jedi_increment),     intent(inout) :: dxc

integer :: f

! Forced identity
! ---------------
do f = 1, size(dxc%fields)
  dxc%fields(f)%array = dxa%fields(f)%array
enddo

end subroutine multiplyinverse

! --------------------------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,dxc,dxa)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: dxc
type(fv3jedi_increment),     intent(inout) :: dxa

integer :: f

! Forced identity
! ---------------
do f = 1, size(dxc%fields)
  dxa%fields(f)%array = dxc%fields(f)%array
enddo

end subroutine multiplyinverseadjoint

! --------------------------------------------------------------------------------------------------

end module fv3jedi_linvarcha_c2a_mod
