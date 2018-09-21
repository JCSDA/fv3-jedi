! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_varcha_c2m_mod

use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom
use iso_c_binding
use config_mod
use kinds

use moisture_vt_mod, only: esinit, dqsat

implicit none

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_varcha_c2m
 integer :: degsubs   = 100
 real(8) :: tmintbl   = 150.0_8, tmaxtbl = 333.0_8
 integer :: tablesize
 real(8), allocatable :: estblx(:)
 real(kind=kind_real), allocatable :: ttraj(:,:,:)
 real(kind=kind_real), allocatable :: tvtraj(:,:,:)
 real(kind=kind_real), allocatable :: qtraj(:,:,:)
 real(kind=kind_real), allocatable :: qsattraj(:,:,:)
end type fv3jedi_varcha_c2m

#define LISTED_TYPE fv3jedi_varcha_c2m

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_varcha_c2m_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_setup(self, bg, fg, geom, c_conf)

use tmprture_vt_mod, only: T_to_Tv
use pressure_vt_mod, only: delp_to_pe_p_logp

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self    !< Change variable structure
type(fv3jedi_state), target, intent(in) :: bg
type(fv3jedi_state), target, intent(in) :: fg
type(fv3jedi_geom), target,  intent(in) :: geom
type(c_ptr),                 intent(in) :: c_conf  !< Configuration

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
call dqsat( geom,bg%t,pm,self%degsubs,self%tmintbl,self%tmaxtbl,&
            self%tablesize,self%estblx,dqsatdt,self%qsattraj)

deallocate(dqsatdt)
deallocate(pm)
deallocate(pe)

end subroutine fv3jedi_varcha_c2m_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_delete(self)

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self

if (allocated(self%estblx)) deallocate(self%estblx)
if (allocated(self%tvtraj)) deallocate(self%tvtraj)
if (allocated(self%ttraj)) deallocate(self%ttraj)
if (allocated(self%qtraj)) deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

end subroutine fv3jedi_varcha_c2m_delete

! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_multiply(self,geom,xctl,xmod)

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self
type(fv3jedi_geom), target,  intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: xctl
type(fv3jedi_increment), intent(inout) :: xmod

!Tangent linear of analysis (control) to model variables
call control_to_model_tlm(geom,xctl%psi,xctl%chi,xctl%tv,xctl%ps  ,xctl%qc, &
                                    xmod%ua  ,xmod%va  ,xmod%t,xmod%delp,xmod%q, &
                                    self%tvtraj,self%qtraj,self%qsattraj)

end subroutine fv3jedi_varcha_c2m_multiply

! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_multiplyadjoint(self,geom,xmod,xctl)

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self
type(fv3jedi_geom), target,  intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: xmod
type(fv3jedi_increment), intent(inout) :: xctl

!Adjoint of analysis (control) to model variables
call control_to_model_adm(geom,xctl%psi,xctl%chi,xctl%tv,xctl%ps  ,xctl%qc, &
                                    xmod%ua  ,xmod%va  ,xmod%t,xmod%delp,xmod%q, &
                                    self%tvtraj,self%qtraj,self%qsattraj)

end subroutine fv3jedi_varcha_c2m_multiplyadjoint

! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_multiplyinverse(self,geom,xmod,xctr)

!use wind_vt_mod, only: uv_to_vortdivg, vortdivg_to_psichi
use tmprture_vt_mod, only: t_to_tv_tl
use moisture_vt_mod, only: q_to_rh_tl

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self
type(fv3jedi_geom), target,  intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: xmod
type(fv3jedi_increment), intent(inout) :: xctr

real(kind=kind_real), allocatable, dimension(:,:,:) :: vort, divg, ua, va

!Tangent linear inverse (model to control)

 xctr%psi = 0.0_kind_real
 xctr%chi = 0.0_kind_real
 xctr%tv  = 0.0_kind_real
 xctr%ps  = 0.0_kind_real
 xctr%qc  = 0.0_kind_real
 xctr%qic = 0.0_kind_real
 xctr%qlc = 0.0_kind_real
 xctr%o3c = 0.0_kind_real

 xctr%psi =  xmod%ua
 xctr%chi =  xmod%va
 xctr%tv  =  xmod%t
 xctr%ps  =  xmod%delp(:,:,geom%npz)
 xctr%qc  =  xmod%q
 xctr%qic =  xmod%qi
 xctr%qlc =  xmod%ql
 xctr%o3c =  xmod%o3

!allocate (vort(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (divg(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (  ua(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!allocate (  va(geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
!
!!> Convert u,v to vorticity and divergence
!call uv_to_vortdivg(geom,xmod%ua,xmod%va,ua,va,vort,divg)
!
!!> Poisson solver for vorticity and divergence to psi and chi
!call vortdivg_to_psichi(geom,vort,divg,xctr%psi,xctr%chi)
!
!!> T to Tv
!call t_to_tv_tl(geom,self%ttraj,xmod%t,self%qtraj,xmod%q)
!xctr%tv = xmod%t
!
!!> q to RH
!call q_to_rh_tl(geom,self%qsattraj,xmod%q,xctr%qc)
!
!!Deallocate
!deallocate(vort,divg,ua,va)

end subroutine fv3jedi_varcha_c2m_multiplyinverse

! ------------------------------------------------------------------------------

subroutine fv3jedi_varcha_c2m_multiplyinverseadjoint(self,geom,xctr,xmod)

implicit none
type(fv3jedi_varcha_c2m), intent(inout) :: self
type(fv3jedi_geom), target,  intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: xctr
type(fv3jedi_increment), intent(inout) :: xmod

!> Not implemented

end subroutine fv3jedi_varcha_c2m_multiplyinverseadjoint

! ------------------------------------------------------------------------------

subroutine control_to_model_tlm(geom,psi,chi,tv,ps,qc,ua,va,t,delp,qs,tvt,qt,qsat)

! use wind_vt_mod, only: psichi_to_udvd
 use tmprture_vt_mod, only: tv_to_t_tl
 use pressure_vt_mod, only: ps_to_delp_tl
 use moisture_vt_mod, only: rh_to_q_tl 

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Input: control vector
 real(kind=kind_real), intent(inout) ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout) ::   tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(inout) ::   ps(geom%isc:geom%iec,geom%jsc:geom%jec           ) !Surface pressure
 real(kind=kind_real), intent(inout) ::   qc(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Output: state/model vector
 real(kind=kind_real), intent(inout) ::   ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::   va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::    t(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Pressure thickness
 real(kind=kind_real), intent(inout) ::   qs(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity
 
 !Trajectory for virtual temperature to temperature
 real(kind=kind_real), intent(in   ) ::  tvt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in   ) ::   qt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in   ) :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Sat spec hum

 ua = 0.0_kind_real
 va = 0.0_kind_real
 t = 0.0_kind_real
 delp = 0.0_kind_real
 qs = 0.0_kind_real

 !psi and chi to D-grid u and v 
 !-----------------------------
! call psichi_to_udvd(geom,psi,chi,ua,va)
 
 !ps to delp
 !----------
 call ps_to_delp_tl(geom,ps,delp)

 !Relative humidity to specific humidity
 !--------------------------------------
 call rh_to_q_tl(geom,qsat,qc,qs)

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_tl(geom,Tvt,Tv,qt,qs,T)

endsubroutine control_to_model_tlm

! ------------------------------------------------------------------------------

!> Control variables to state variables - Adjoint

subroutine control_to_model_adm(geom,psi,chi,tv,ps,qc,ua,va,t,delp,qs,tvt,qt,qsat)

! use wind_vt_mod, only: psichi_to_udvd_adm
 use tmprture_vt_mod, only: tv_to_t_ad
 use pressure_vt_mod, only: ps_to_delp_ad
 use moisture_vt_mod, only: rh_to_q_ad 

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Input: control vector
 real(kind=kind_real), intent(inout) ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout) ::   tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(inout) ::   ps(geom%isc:geom%iec,geom%jsc:geom%jec           ) !Surface pressure
 real(kind=kind_real), intent(inout) ::   qc(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity

 !Output: state/model vector
 real(kind=kind_real), intent(inout) ::   ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::   va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::    t(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Pressure thickness
 real(kind=kind_real), intent(inout) ::   qs(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity
 
 !Trajectory for virtual temperature to temperaturc
 real(kind=kind_real), intent(in   ) ::  tvt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in   ) ::   qt(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in   ) :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Sat spec hum

 psi = 0.0_kind_real
 chi = 0.0_kind_real
 tv  = 0.0_kind_real
 ps  = 0.0_kind_real
 qc  = 0.0_kind_real

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_ad(geom,Tvt,Tv,qt,qs,T)

 !Relative humidity to specific humidity
 !--------------------------------------
 call rh_to_q_ad(geom,qsat,qc,qs)

 !Ps to delp
 !----------
 call ps_to_delp_ad(geom,ps,delp)

 !psi and chi to D-grid u and v 
 !-----------------------------
! call psichi_to_udvd_adm(geom,psi,chi,ua,va)

endsubroutine control_to_model_adm

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_c2m_mod
