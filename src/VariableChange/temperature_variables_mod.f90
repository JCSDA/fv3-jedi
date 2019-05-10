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

public T_to_Tv, T_to_Tv_tl, T_to_Tv_ad
public Tv_to_T, Tv_to_T_tl, Tv_to_T_ad

contains

!----------------------------------------------------------------------------
! Temperature to Virtual Temperature ----------------------------------------
!----------------------------------------------------------------------------

subroutine T_to_Tv(geom,T,q,Tv)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: T (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: Tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Virtual temperature (K)

 Tv = T*(1.0_kind_real + epsilon*q)

end subroutine T_to_Tv

!----------------------------------------------------------------------------

subroutine T_to_Tv_tl(geom,T,T_tl,q,q_tl,Tv_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom
 real(kind=kind_real), intent(in ) :: T   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: T_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: q   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in ) :: q_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(out) :: Tv_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 Tv_tl = T_tl*(1.0_kind_real + epsilon*q) + T*epsilon*q_tl

end subroutine T_to_Tv_tl

!----------------------------------------------------------------------------

subroutine T_to_Tv_ad(geom,T,T_ad,q,q_ad,Tv_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: T    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: Tv_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 T_ad = T_ad + Tv_ad * (1.0_kind_real + epsilon*q)
 q_ad = q_ad + Tv_ad *                  epsilon*T
 tv_ad = 0.0_kind_real

end subroutine T_to_Tv_ad

!----------------------------------------------------------------------------
! Virtual Temperature to Temperature ----------------------------------------
!----------------------------------------------------------------------------

subroutine Tv_to_T(geom,Tv,q,T)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: Tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: T (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)  !Virtual temperature (K)

 T = Tv/(1.0_kind_real + epsilon*q)

end subroutine Tv_to_T

!----------------------------------------------------------------------------

subroutine Tv_to_T_tl(geom,Tv,Tv_tl,q,q_tl,T_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: Tv   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: Tv_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q_tl (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_tl (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 t_tl = (tv_tl*(1.0_kind_real+epsilon*q)-tv*epsilon*q_tl)/(1.0_kind_real+epsilon*q)**2

end subroutine Tv_to_T_tl

!----------------------------------------------------------------------------

subroutine Tv_to_T_ad(geom,Tv,Tv_ad,q,q_ad,T_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: Tv   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: Tv_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_ad (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 real(kind=kind_real) :: temp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 temp = t_ad/(epsilon*q+1.0_kind_real)

 tv_ad = tv_ad + temp
 q_ad  = q_ad  - tv*epsilon*temp/(epsilon*q+1.0_kind_real)
 t_ad = 0.0_kind_real

end subroutine Tv_to_T_ad

!----------------------------------------------------------------------------

end module temperature_vt_mod
