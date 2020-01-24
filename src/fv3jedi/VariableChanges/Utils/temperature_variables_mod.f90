! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module temperature_vt_mod

use fv3jedi_constants_mod, only: kappa, epsilon
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

implicit none
private

public t_to_tv, t_to_tv_tl, t_to_tv_ad
public tv_to_t, tv_to_t_tl, tv_to_t_ad
public pt_to_t, pt_to_t_tl, pt_to_t_ad
public t_to_pt

contains

!----------------------------------------------------------------------------
! Temperature to Virtual Temperature ----------------------------------------
!----------------------------------------------------------------------------

subroutine t_to_tv(geom,t,q,tv)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: t (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Virtual temperature (K)

 tv = t*(1.0_kind_real + epsilon*q)

end subroutine t_to_tv

!----------------------------------------------------------------------------

subroutine t_to_tv_tl(geom,t,t_tl,q,q_tl,tv_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom
 real(kind=kind_real), intent(in ) :: t   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: t_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: q   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: q_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(out) :: tv_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 tv_tl = t_tl*(1.0_kind_real + epsilon*q) + t*epsilon*q_tl

end subroutine t_to_tv_tl

!----------------------------------------------------------------------------

subroutine t_to_tv_ad(geom,t,t_ad,q,q_ad,tv_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: t    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: t_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: tv_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 t_ad = t_ad + tv_ad * (1.0_kind_real + epsilon*q)
 q_ad = q_ad + tv_ad *                  epsilon*t
 tv_ad = 0.0_kind_real

end subroutine t_to_tv_ad

!----------------------------------------------------------------------------
! Virtual Temperature to Temperature ----------------------------------------
!----------------------------------------------------------------------------

subroutine tv_to_t(geom,tv,q,t)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: t (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Virtual temperature (K)

 t = tv/(1.0_kind_real + epsilon*q)

end subroutine tv_to_t

!----------------------------------------------------------------------------

subroutine tv_to_t_tl(geom,tv,tv_tl,q,q_tl,t_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: tv   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: tv_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q_tl (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: t_tl (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 t_tl = (tv_tl*(1.0_kind_real+epsilon*q)-tv*epsilon*q_tl)/(1.0_kind_real+epsilon*q)**2

end subroutine tv_to_t_tl

!----------------------------------------------------------------------------

subroutine tv_to_t_ad(geom,tv,tv_ad,q,q_ad,t_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: tv   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: tv_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: t_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 real(kind=kind_real) :: temp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 temp = t_ad/(epsilon*q+1.0_kind_real)

 tv_ad = tv_ad + temp
 q_ad  = q_ad  - tv*epsilon*temp/(epsilon*q+1.0_kind_real)
 t_ad = 0.0_kind_real

end subroutine tv_to_t_ad

!----------------------------------------------------------------------------
! Potential Temperature to Temperature --------------------------------------
!----------------------------------------------------------------------------

subroutine pt_to_t(geom,pkz,pt,t)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: pkz(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Pressure to the kappa
 real(kind=kind_real), intent(in ) :: pt (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Potential temperature (K)
 real(kind=kind_real), intent(out) :: t  (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)

 t = pt*pkz

end subroutine pt_to_t

!----------------------------------------------------------------------------

subroutine pt_to_t_tl(geom,pkz,pkz_tl,pt,pt_tl,t_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: pkz   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: pkz_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: pt    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: pt_tl (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: t_tl  (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 t_tl = pt_tl*pkz + pt*pkz_tl

end subroutine pt_to_t_tl

!----------------------------------------------------------------------------

subroutine pt_to_t_ad(geom,pkz,pkz_ad,pt,pt_ad,t_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: pkz   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: pkz_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: pt    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: pt_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: t_ad  (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 pt_ad = pt_ad + pkz*t_ad
 pkz_ad = pkz_ad + pt*t_ad
 t_ad = 0.0

end subroutine pt_to_t_ad

!----------------------------------------------------------------------------

subroutine t_to_pt(geom,pkz,t,pt)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: pkz(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Pressure to the kappa
 real(kind=kind_real), intent(in ) :: t  (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Potential temperature (K)
 real(kind=kind_real), intent(out) :: pt (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)

 pt = t/pkz

end subroutine t_to_pt

!----------------------------------------------------------------------------

end module temperature_vt_mod
