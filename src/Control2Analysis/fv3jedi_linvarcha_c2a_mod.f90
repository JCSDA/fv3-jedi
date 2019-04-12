! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_linvarcha_c2a_mod

use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom
use iso_c_binding
use config_mod
use fv3jedi_kinds_mod

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
 integer :: degsubs   = 100
 real(8) :: tmintbl   = 150.0_8, tmaxtbl = 333.0_8
 integer :: tablesize
 real(8), allocatable :: estblx(:)
 real(kind=kind_real), allocatable :: ttraj(:,:,:)
 real(kind=kind_real), allocatable :: tvtraj(:,:,:)
 real(kind=kind_real), allocatable :: qtraj(:,:,:)
 real(kind=kind_real), allocatable :: qsattraj(:,:,:)
end type fv3jedi_linvarcha_c2a

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, c_conf)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(fv3jedi_state), target, intent(in)    :: bg
type(fv3jedi_state), target, intent(in)    :: fg
type(c_ptr),                 intent(in)    :: c_conf

real(kind=kind_real), allocatable :: pe(:,:,:)
real(kind=kind_real), allocatable :: pm(:,:,:)
real(kind=kind_real), allocatable :: dqsatdt(:,:,:)

!Create lookup table for computing saturation specific humidity
self%tablesize = nint(self%tmaxtbl-self%tmintbl)*self%degsubs + 1
allocate(self%estblx(self%tablesize))
call esinit(self%tablesize,self%degsubs,self%tmintbl,self%tmaxtbl,self%estblx)

!> Virtual temperature trajectory
allocate(self%tvtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
call T_to_Tv(geom,bg%t,bg%q,self%tvtraj)

!> Temperature trajectory
allocate(self%ttraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%ttraj = bg%t

!> Specific humidity trajecotory
allocate(self%qtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%qtraj = bg%q

!> Compute saturation specific humidity for q to RH transform
allocate(self%qsattraj(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

allocate(pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
allocate(pm(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz  ))
allocate(dqsatdt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz  ))

call delp_to_pe_p_logp(geom,bg%delp,pe,pm)
call dqsat_calc( geom,bg%t,pm,self%degsubs,self%tmintbl,self%tmaxtbl,&
                 self%tablesize,self%estblx,dqsatdt,self%qsattraj)

deallocate(dqsatdt)
deallocate(pm)
deallocate(pe)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_c2a), intent(inout) :: self

if (allocated(self%estblx)) deallocate(self%estblx)
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

!Ps (identity)
xana%ps = xctl%ps

!Tracers (identity)
xana%qi = xctl%qi
xana%ql = xctl%ql
xana%o3 = xctl%o3

!Tangent linear of control to analysis variables
call control_to_analysis_tlm(geom, xctl%psi, xctl%chi, xctl%tv, xctl%rh, &
                                   xana%ua,  xana%va,  xana%t,  xana%q,  &
                                   self%tvtraj,self%qtraj,self%qsattraj )

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xana,xctl)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(inout) :: xana
type(fv3jedi_increment),     intent(inout) :: xctl

!Ps (identity)
xctl%ps = xana%ps

!Tracers (identity)
xctl%qi = xana%qi
xctl%ql = xana%ql
xctl%o3 = xana%o3

!Adjoint of control to analysis variables
call control_to_analysis_adm(geom, xctl%psi, xctl%chi, xctl%tv, xctl%rh, &
                                   xana%ua,  xana%va,  xana%t,  xana%q,  &
                                   self%tvtraj,self%qtraj,self%qsattraj )

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xana,xctl)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xctl

real(kind=kind_real), allocatable, dimension(:,:,:) :: vort, divg, ua, va

!Tangent linear inverse (analysis to control)

xctl%psi = xana%ua
xctl%chi = xana%va
xctl%tv  = xana%t 
xctl%ps  = xana%ps
xctl%rh  = xana%q 
xctl%qi  = xana%qi
xctl%ql  = xana%ql
xctl%o3  = xana%o3

!allocate (vort(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (divg(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (  ua(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (  va(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!
!!> Convert u,v to vorticity and divergence
!call uv_to_vortdivg(geom,xana%ua,xana%va,ua,va,vort,divg)
!
!!> Poisson solver for vorticity and divergence to psi and chi
!call vortdivg_to_psichi(geom,vort,divg,xctl%psi,xctl%chi)
!
!!> T to Tv
!call t_to_tv_tl(geom,self%ttraj,xana%t,self%qtraj,xana%q)
!xctl%tv = xana%t
!
!!> q to RH
!call q_to_rh_tl(geom,self%qsattraj,xana%q,xctl%rh)
!
!!Deallocate
!deallocate(vort,divg,ua,va)

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xctl,xana)

implicit none
type(fv3jedi_linvarcha_c2a), intent(in)    :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xctl
type(fv3jedi_increment),     intent(inout) :: xana

xana%ua = xctl%psi
xana%va = xctl%chi
xana%t  = xctl%tv 
xana%ps = xctl%ps 
xana%q  = xctl%rh 
xana%qi = xctl%qi 
xana%ql = xctl%ql 
xana%o3 = xctl%o3 

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
