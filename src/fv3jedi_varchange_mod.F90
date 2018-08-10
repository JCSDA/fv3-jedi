! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_varchange_mod

use fv3jedi_fields_mod, only: fv3jedi_field
use fv3jedi_geom_mod,   only: fv3jedi_geom
use iso_c_binding
use config_mod
use kinds

use moisture_vt_mod, only: esinit, dqsat
use pressure_vt_mod, only: delp_to_pe_p_logp

implicit none

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_varchange
 integer :: degsubs   = 100
 real(8) :: tmintbl   = 150.0, tmaxtbl = 333.0
 integer :: tablesize
 real(8), allocatable :: estblx(:)
 real(kind=kind_real), allocatable :: ttraj(:,:,:)
 real(kind=kind_real), allocatable :: qtraj(:,:,:)
 real(kind=kind_real), allocatable :: qsattraj(:,:,:)
end type fv3jedi_varchange

#define LISTED_TYPE fv3jedi_varchange

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_varchange_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_setup(self, c_conf)

implicit none
type(fv3jedi_varchange), intent(inout) :: self    !< Change variable structure
type(c_ptr),             intent(in)    :: c_conf  !< Configuration

!Create lookup table for computing saturation specific humidity
self%tablesize = nint(self%tmaxtbl-self%tmintbl)*self%degsubs + 1
allocate(self%estblx(self%tablesize))
call esinit(self%tablesize,self%degsubs,self%tmintbl,self%tmaxtbl,self%estblx)

end subroutine fv3jedi_varchange_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_delete(self)

implicit none
type(fv3jedi_varchange), intent(inout) :: self

if (allocated(self%estblx)) deallocate(self%estblx)
if (allocated(self%ttraj)) deallocate(self%ttraj)
if (allocated(self%qtraj)) deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

end subroutine fv3jedi_varchange_delete

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_linearize(self,geom,traj)

implicit none
type(fv3jedi_varchange), intent(inout) :: self
type(fv3jedi_geom),  intent(inout) :: geom
type(fv3jedi_field), intent(inout) :: traj

real(kind=kind_real), allocatable :: pe(:,:,:)
real(kind=kind_real), allocatable :: pm(:,:,:)
real(kind=kind_real), allocatable :: dqsatdt(:,:,:)

!> Compute saturation specific humidity for q to RH transform

allocate(self%ttraj   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz))
allocate(self%qtraj   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz))
allocate(self%qsattraj(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz))

self%ttraj = traj%Atm%pt
self%qtraj = traj%Atm%q(:,:,:,1)

allocate(pe(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz+1))
allocate(pm(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz  ))

call delp_to_pe_p_logp(geom,traj%Atm%delp,pe,pm)

allocate(dqsatdt(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz  ))

call dqsat( geom,traj%Atm%pt,pm,self%degsubs,self%tmintbl,self%tmaxtbl,&
            self%tablesize,self%estblx,dqsatdt,self%qsattraj)

deallocate(dqsatdt)
deallocate(pm)
deallocate(pe)

end subroutine fv3jedi_varchange_linearize

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_transform(self,xctl,xmod)

implicit none
type(fv3jedi_varchange), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xctl
type(fv3jedi_field), intent(inout) :: xmod

!Tangent linear of analysis (control) to model variables
xctl%Atm%qct = 0.0_kind_real

call control_to_state_tlm(xctl%geom,xctl%Atm%psi,xctl%Atm%chi,xctl%Atm%tv,xctl%Atm%ps  ,xctl%Atm%qct(:,:,:,1), &
                                    xmod%Atm%u  ,xmod%Atm%v  ,xmod%Atm%pt,xmod%Atm%delp,xmod%Atm%q  (:,:,:,1), &
                                    self%ttraj,self%qtraj,self%qsattraj)

end subroutine fv3jedi_varchange_transform

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_transformadjoint(self,xmod,xctl)

implicit none
type(fv3jedi_varchange), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xmod
type(fv3jedi_field), intent(inout) :: xctl

!Adjoint of analysis (control) to model variables
call control_to_state_adm(xctl%geom,xctl%Atm%psi,xctl%Atm%chi,xctl%Atm%tv,xctl%Atm%ps  ,xctl%Atm%qct(:,:,:,1), &
                                    xmod%Atm%u  ,xmod%Atm%v  ,xmod%Atm%pt,xmod%Atm%delp,xmod%Atm%q  (:,:,:,1), &
                                    self%ttraj,self%qtraj,self%qsattraj)

end subroutine fv3jedi_varchange_transformadjoint

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_transforminverse(self,xinc,xctr)

implicit none
type(fv3jedi_varchange), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

!> Not implemented

end subroutine fv3jedi_varchange_transforminverse

! ------------------------------------------------------------------------------

subroutine fv3jedi_varchange_transforminverseadjoint(self,xinc,xctr)

implicit none
type(fv3jedi_varchange), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

!> Not implemented

end subroutine fv3jedi_varchange_transforminverseadjoint

! ------------------------------------------------------------------------------

subroutine control_to_state_tlm(geom,psi,chi,tv,ps,qc,u,v,t,delp,qs,tvt,qt,qsat)

! use wind_vt_mod, only: psichi_to_udvd
 use wind_vt_mod, only: psichi_to_udvd
 use tmprture_vt_mod, only: tv_to_t_tl
 use pressure_vt_mod, only: ps_to_delp_tl
 use moisture_vt_mod, only: q_to_rh_tl 

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Input: control vector
 real(kind=kind_real), intent(inout) ::  psi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) ::  chi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout) ::   tv(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(inout) ::   ps(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed             ) !Surface pressure
 real(kind=kind_real), intent(inout) ::   qc(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity

 !Output: state/model vector
 real(kind=kind_real), intent(inout) ::    u(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::    v(geom%bd%isd:geom%bd%ied+1,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::    t(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) :: delp(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Pressure thickness
 real(kind=kind_real), intent(inout) ::   qs(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity
 
 !Trajectory for virtual temperature to temperature
 real(kind=kind_real), intent(in   ) ::  tvt(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in   ) ::   qt(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in   ) :: qsat(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Sat spec hum

 u    = 0.0_kind_real
 v    = 0.0_kind_real
 t    = 0.0_kind_real
 delp = 0.0_kind_real
 qs   = 0.0_kind_real

 !psi and chi to D-grid u and v 
 !-----------------------------
 call psichi_to_udvd(geom,psi,chi,u,v)

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_tl(geom,Tvt,Tv,qt,qc,T)

 !Ps to delp
 !----------
 call ps_to_delp_tl(geom,ps,delp)

 !Specific to relative humidity
 !-----------------------------
 call q_to_rh_tl(geom,qsat,qc,qs)

endsubroutine control_to_state_tlm

! ------------------------------------------------------------------------------

!> Control variables to state variables - Adjoint

subroutine control_to_state_adm(geom,psi,chi,tv,ps,qc,u,v,t,delp,qs,tvt,qt,qsat)

! use wind_vt_mod, only: psichi_to_udvd_adm
 use wind_vt_mod, only: psichi_to_udvd_adm
 use tmprture_vt_mod, only: tv_to_t_ad
 use pressure_vt_mod, only: ps_to_delp_ad
 use moisture_vt_mod, only: q_to_rh_ad 

 implicit none
 type(fv3jedi_geom), intent(inout) :: geom

 !Input: control vector
 real(kind=kind_real), intent(inout) ::  psi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) ::  chi(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout) ::   tv(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Virtual temp
 real(kind=kind_real), intent(inout) ::   ps(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed             ) !Surface pressure
 real(kind=kind_real), intent(inout) ::   qc(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity

 !Output: state/model vector
 real(kind=kind_real), intent(inout) ::    u(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::    v(geom%bd%isd:geom%bd%ied+1,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::    t(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Dry temperature
 real(kind=kind_real), intent(inout) :: delp(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Pressure thickness
 real(kind=kind_real), intent(inout) ::   qs(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity
 
 !Trajectory for virtual temperature to temperature
 real(kind=kind_real), intent(in   ) ::  tvt(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !VTemperature traj
 real(kind=kind_real), intent(in   ) ::   qt(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Specific humidity traj
 real(kind=kind_real), intent(in   ) :: qsat(geom%bd%isd:geom%bd%ied  ,geom%bd%jsd:geom%bd%jed  ,1:geom%npz) !Sat spec hum

 psi = 0.0_kind_real
 chi = 0.0_kind_real
 tv  = 0.0_kind_real
 ps  = 0.0_kind_real
 qc  = 0.0_kind_real

 !Specific to relative humidity
 !-----------------------------
 call q_to_rh_ad(geom,qsat,qc,qs)

 !Ps to delp
 !----------
 call ps_to_delp_ad(geom,ps,delp)

 !Virtual temperature to temperature
 !----------------------------------
 call Tv_to_T_ad(geom,Tvt,Tv,qt,qc,T)

 !psi and chi to D-grid u and v 
 !-----------------------------
 call psichi_to_udvd_adm(geom,psi,chi,u,v)

endsubroutine control_to_state_adm

! ------------------------------------------------------------------------------

end module fv3jedi_varchange_mod
