! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_linvarcha_c2a_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fv3jedi_kinds_mod

use fv3jedi_field_mod, only: copy_subset, has_field, pointer_field_array

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

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, conf)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(fv3jedi_state), target, intent(in)    :: bg
type(fv3jedi_state), target, intent(in)    :: fg
type(fckit_configuration),   intent(in)    :: conf

real(kind=kind_real), pointer :: t   (:,:,:)
real(kind=kind_real), pointer :: q   (:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)

!> Pointers to the background state
call pointer_field_array(bg%fields, 't'   , t)
call pointer_field_array(bg%fields, 'q'   , q)
call pointer_field_array(bg%fields, 'delp', delp)

!> Virtual temperature trajectory
allocate(self%tvtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
call T_to_Tv(geom,t,q,self%tvtraj)

!> Temperature trajectory
allocate(self%ttraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%ttraj = t

!> Specific humidity trajecotory
allocate(self%qtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%qtraj = q

!> Compute saturation specific humidity for q to RH transform
allocate(self%qsattraj(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

!> Compute saturation specific humidity
call get_qsat(geom,delp,t,q,self%qsattraj)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self

if (allocated(self%tvtraj)) deallocate(self%tvtraj)
if (allocated(self%ttraj)) deallocate(self%ttraj)
if (allocated(self%qtraj)) deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xctl,xana)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xctl
type(fv3jedi_increment),     intent(inout) :: xana

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_q

real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_rh
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_ps

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xctl%fields,xana%fields)

!Tangent linear of control to analysis variables
call pointer_field_array(xctl%fields, 'psi', xctl_psi)
call pointer_field_array(xctl%fields, 'chi', xctl_chi)
call pointer_field_array(xctl%fields, 'tv' , xctl_tv)
call pointer_field_array(xctl%fields, 'rh' , xctl_rh)
call pointer_field_array(xana%fields, 'ua' , xana_ua)
call pointer_field_array(xana%fields, 'va' , xana_va)
call pointer_field_array(xana%fields, 't'  , xana_t)
call pointer_field_array(xana%fields, 'q'  , xana_q)

call control_to_analysis_tlm(geom, xctl_psi, xctl_chi, xctl_tv, xctl_rh, &
                                   xana_ua,  xana_va,  xana_t,  xana_q,  &
                                   self%tvtraj,self%qtraj,self%qsattraj )

! Copy calendar infomation
xana%calendar_type = xctl%calendar_type
xana%date_init = xctl%date_init

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xana,xctl)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(inout) :: xana
type(fv3jedi_increment),     intent(inout) :: xctl

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_q

real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_rh
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_ps

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xana%fields,xctl%fields)

!Adjoint of control to analysis variables
call pointer_field_array(xctl%fields, 'psi', xctl_psi)
call pointer_field_array(xctl%fields, 'chi', xctl_chi)
call pointer_field_array(xctl%fields, 'tv' , xctl_tv)
call pointer_field_array(xctl%fields, 'rh' , xctl_rh)
call pointer_field_array(xana%fields, 'ua' , xana_ua)
call pointer_field_array(xana%fields, 'va' , xana_va)
call pointer_field_array(xana%fields, 't'  , xana_t)
call pointer_field_array(xana%fields, 'q'  , xana_q)
call control_to_analysis_adm(geom, xctl_psi, xctl_chi, xctl_tv, xctl_rh, &
                                   xana_ua,  xana_va,  xana_t,  xana_q,  &
                                   self%tvtraj,self%qtraj,self%qsattraj )

! Copy calendar infomation
xctl%calendar_type = xana%calendar_type
xctl%date_init = xana%date_init

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xana,xctl)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xctl

real(kind=kind_real), allocatable, dimension(:,:,:) :: vort, divg, ua, va

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xana%fields,xctl%fields)

! Copy calendar infomation
xctl%calendar_type = xana%calendar_type
xctl%date_init = xana%date_init

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xctl,xana)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xctl
type(fv3jedi_increment),     intent(inout) :: xana

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xctl%fields,xana%fields)

! Copy calendar infomation
xana%calendar_type = xctl%calendar_type
xana%date_init = xctl%date_init

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

subroutine control_to_analysis_tlm(geom,psi, chi, tv, rh, &
                                        ua , va , t , q, &
                                   tvt, qt, qsat)

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Input: control variables
 real(kind=kind_real), intent(in)    ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real), intent(in)    ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(in)    ::   tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(in)    ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Output: analysis variables
 real(kind=kind_real), intent(inout) ::   ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !A-grid winds (ua)
 real(kind=kind_real), intent(inout) ::   va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !A-grid winds (va)
 real(kind=kind_real), intent(inout) ::    t(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Trajectory for virtual temperature to temperature
 real(kind=kind_real), intent(in)    ::  tvt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in)    ::   qt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Sat spec hum

 real(kind=kind_real), allocatable, dimension(:,:,:) :: psi_dom, chi_dom

 ua = 0.0_kind_real
 va = 0.0_kind_real
 t  = 0.0_kind_real
 q  = 0.0_kind_real

 !psi and chi to A-grid u and v
 !-----------------------------
 allocate(psi_dom(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 allocate(chi_dom(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 psi_dom = 0.0_kind_real
 chi_dom = 0.0_kind_real

 psi_dom(geom%isc:geom%iec,geom%jsc:geom%jec,:) = psi
 chi_dom(geom%isc:geom%iec,geom%jsc:geom%jec,:) = chi

 call psichi_to_uava(geom,psi_dom,chi_dom,ua,va)

 deallocate(psi_dom, chi_dom)

 !Relative humidity to specific humidity
 !--------------------------------------
 call rh_to_q_tl(geom,qsat,rh,q)

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_tl(geom,Tvt,Tv,qt,q,T)

endsubroutine control_to_analysis_tlm

! ------------------------------------------------------------------------------

!> Control variables to state variables - Adjoint

subroutine control_to_analysis_adm(geom,psi, chi, tv, rh, &
                                        ua , va , t , q, &
                                   tvt, qt, qsat)

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Output: control variables
 real(kind=kind_real), intent(inout) ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout) ::   tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Input: analysis variables
 real(kind=kind_real), intent(inout) ::   ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::   va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::    t(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Trajectory for virtual temperature to temperaturc
 real(kind=kind_real), intent(in)    ::  tvt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in)    ::   qt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Sat spec hum

 real(kind=kind_real), allocatable, dimension(:,:,:) :: psi_dom, chi_dom

 psi = 0.0_kind_real
 chi = 0.0_kind_real
 tv  = 0.0_kind_real
 rh  = 0.0_kind_real

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_ad(geom,Tvt,Tv,qt,q,T)

 !Relative humidity to specific humidity
 !--------------------------------------
 call rh_to_q_ad(geom,qsat,rh,q)

 !psi and chi to D-grid u and v
 !-----------------------------
 allocate(psi_dom(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 allocate(chi_dom(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 psi_dom = 0.0_kind_real
 chi_dom = 0.0_kind_real

 call psichi_to_uava_adm(geom,psi_dom,chi_dom,ua,va)

 psi = psi_dom(geom%isc:geom%iec,geom%jsc:geom%jec,:)
 chi = chi_dom(geom%isc:geom%iec,geom%jsc:geom%jec,:)

 deallocate(psi_dom, chi_dom)

endsubroutine control_to_analysis_adm

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_c2a_mod
