! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on temperature variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module tmprture_vt_mod

use fv3jedi_constants, only: kappa, epsilon
use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real

implicit none
public

contains

!----------------------------------------------------------------------------
! Temperature to Virtual Temperature ----------------------------------------
!----------------------------------------------------------------------------

subroutine T_to_Tv(geom,T,q,Tv)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: T (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: Tv(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Virtual temperature (K)

 integer :: isc,iec,jsc,jec

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 Tv(isc:iec,jsc:jec,:) = T(isc:iec,jsc:jec,:)*(1.0 + epsilon*q(isc:iec,jsc:jec,:))

end subroutine T_to_Tv

subroutine T_to_Tv_tl(geom,T,T_tl,q,q_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: T   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_tl(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q_tl(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)

 real(kind=kind_real) :: Tv_tl(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)

 integer :: isc,iec,jsc,jec

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 Tv_tl(isc:iec,jsc:jec,:) = T_tl(isc:iec,jsc:jec,:)*(1.0 + epsilon*q(isc:iec,jsc:jec,:)) + T(isc:iec,jsc:jec,:)*epsilon*q_tl(isc:iec,jsc:jec,:)

 !Replace temperature with virtual temperature
 T_tl(isc:iec,jsc:jec,:) = Tv_tl(isc:iec,jsc:jec,:)

end subroutine T_to_Tv_tl

subroutine T_to_Tv_ad(geom,T,T_ad,q,q_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in )   :: geom 
 real(kind=kind_real), intent(in )   :: T   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in )   :: q   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_ad(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)

 real(kind=kind_real) :: TV_ad(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)

 integer :: isc,iec,jsc,jec

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 Tv_ad(isc:iec,jsc:jec,:) = T_ad(isc:iec,jsc:jec,:)
 T_ad(isc:iec,jsc:jec,:)  =                    (1.0 + epsilon*q(isc:iec,jsc:jec,:))*Tv_ad(isc:iec,jsc:jec,:)
 q_ad(isc:iec,jsc:jec,:)  = q_ad(isc:iec,jsc:jec,:) + epsilon*T(isc:iec,jsc:jec,:) *Tv_ad(isc:iec,jsc:jec,:)
 Tv_ad(isc:iec,jsc:jec,:) = 0.0

end subroutine T_to_Tv_ad

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine Tv_to_T(geom,Tv,q,T)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: Tv(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Temperature (K)
 real(kind=kind_real), intent(in ) :: q (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Specific humidity (kg/kg)
 real(kind=kind_real), intent(out) :: T (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)  !Virtual temperature (K)

 integer :: isc,iec,jsc,jec, i, j, k

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 do k = 1,geom%npz
   do j = jsc,jec
     do i = isc,iec

       T(i,j,k) = Tv(i,j,k)/(1.0_kind_real + epsilon*q(i,j,k))

     enddo
   enddo
 enddo

end subroutine Tv_to_T

subroutine Tv_to_T_tl(geom,Tv,Tv_tl,q,q_tl,T_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: Tv   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: Tv_tl(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q_tl (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_tl (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 
 integer :: isc,iec,jsc,jec,npz,i,j,k

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec
 npz = geom%npz

 t_tl = 0.0
 do k=1,npz
   do j=jsc,jec
     do i=isc,iec
       t_tl(i,j,k) = (tv_tl(i,j,k)*(1.0_kind_real+epsilon*q(i,j,k))-tv(i,j,k)*epsilon*q_tl(i,j,k))/(1.0+epsilon*q(i,j,k))**2
     end do
   end do
 end do

end subroutine Tv_to_T_tl

subroutine Tv_to_T_ad(geom,Tv,Tv_ad,q,q_ad,T_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom
 real(kind=kind_real), intent(in   ) :: Tv   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: Tv_ad(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: q    (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: q_ad (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 real(kind=kind_real), intent(inout) :: T_ad (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
 
 integer :: isc,iec,jsc,jec,npz,i,j,k
 real(kind=kind_real) :: temp

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec
 npz = geom%npz

 do k=npz,1,-1
   do j=jec,jsc,-1
     do i=iec,isc,-1
       temp = epsilon*q(i,j,k) + 1.0_kind_real
       tv_ad(i,j,k) = tv_ad(i,j,k) + t_ad(i,j,k)/temp
       q_ad(i,j,k) = q_ad(i,j,k) - tv(i,j,k)*epsilon*t_ad(i,j,k)/temp**2
       t_ad(i,j,k) = 0.0
     end do
   end do
 end do

end subroutine Tv_to_T_ad

end module tmprture_vt_mod
