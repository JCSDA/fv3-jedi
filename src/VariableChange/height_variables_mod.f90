! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module height_vt_mod

use fv3jedi_constants_mod, only: grav, rvap, tice,rdry, zvir
use fv3jedi_geom_mod,  only: fv3jedi_geom
use fv3jedi_kinds_mod,             only: kind_real

implicit none
private
! Constants from GSI
! Constants for compressibility factor (Davis et al 1992)
real(kind_real),parameter ::  cpf_a0 =  1.58123e-6_kind_real ! K/Pa
real(kind_real),parameter ::  cpf_a1 = -2.9331e-8_kind_real  ! 1/Pa
real(kind_real),parameter ::  cpf_a2 =  1.1043e-10_kind_real ! 1/K 1/Pa
real(kind_real),parameter ::  cpf_b0 =  5.707e-6_kind_real   ! K/Pa
real(kind_real),parameter ::  cpf_b1 = -2.051e-8_kind_real   ! 1/Pa
real(kind_real),parameter ::  cpf_c0 =  1.9898e-4_kind_real  ! K/Pa
real(kind_real),parameter ::  cpf_c1 = -2.376e-6_kind_real   ! 1/Pa
real(kind_real),parameter ::  cpf_d  =  1.83e-11_kind_real   ! K2/Pa2
real(kind_real),parameter ::  cpf_e  = -0.765e-8_kind_real   ! K2/Pa2

real(kind_real),parameter ::  psv_a =  1.2378847e-5_kind_real       !  (1/K2)
real(kind_real),parameter ::  psv_b = -1.9121316e-2_kind_real       !  (1/K)
real(kind_real),parameter ::  psv_c = 33.93711047_kind_real         !
real(kind_real),parameter ::  psv_d = -6.3431645e+3_kind_real       !  (K)
! Constants for enhancement factor to calculating the mole fraction of water vapor
real(kind_real),parameter ::  ef_alpha = 1.00062_kind_real           !
real(kind_real),parameter ::  ef_beta  = 3.14e-8_kind_real           !  (1/Pa)
real(kind_real),parameter ::  ef_gamma = 5.6e-7_kind_real

public geop_height
public geop_height_levels

contains

subroutine geop_height(geom,prs,prsi,T,q,phis,use_compress,gph)

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
real(kind_real), intent(in ) :: prs(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)    !mid layerpressure
real(kind_real), intent(in ) :: prsi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !interface pressure
real(kind_real), intent(in ) :: phis(geom%isc:geom%iec,geom%jsc:geom%jec)              !Surface geopotential (grav*Z_sfc)
real(kind_real), intent(in ) :: T(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
real(kind_real), intent(in ) :: q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)     ! specific humidity
real(kind_real), intent(out) :: gph(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !geopotential height (meters)

!locals
real(kind_real)       :: Tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
real(kind_real)       :: qmr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) ! mixing ratio|kg/kg
logical               :: use_compress
integer               :: isc,iec,jsc,jec,npz,i,j,k
real(kind=kind_real)  :: Tkk,Tvk,Tc, qmk,Pak,dpk,dz
real(kind=kind_real)  :: prs_sv, prs_v
real(kind=kind_real)  :: ehn_fct,x_v,cmpr

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


!get qmr--mixing ratio and virtual temeprature
qmr(isc:iec,jsc:jec,:) = q(isc:iec,jsc:jec,:)/(1.0 - q(isc:iec,jsc:jec,:))
Tv(isc:iec,jsc:jec,:)  = T(isc:iec,jsc:jec,:)*(1.0 + zvir*qmr(isc:iec,jsc:jec,:))

if (use_compress) then

!  Compute compressibility factor (Picard et al 2008) and geopotential heights at midpoint
  do j = jsc,jec
  do i = isc,iec

     do k = geom%npz, 1, -1
        if ( k == geom%npz) then
           Tkk  = T(i,j,k)
           Tvk  = Tv(i,j,k)
           Pak  = exp(0.5_kind_real*(log(prsi(i,j,k+1)*0.01)+log(prs(i,j,k)*0.01)))
           dpk  = prsi(i,j,k+1)/prs(i,j,k)
        else
           Tkk  = 0.5_kind_real * ( T(i,j,k+1) +  T(i,j,k) )
           Tvk  = 0.5_kind_real * (Tv(i,j,k+1) + Tv(i,j,k) )
           Pak  = exp(0.5_kind_real*(log(prs(i,j,k+1)*0.01)+log(prs(i,j,k)*0.01)))
           dpk  = prs(i,j,k+1)/prs(i,j,k)
        end if

        Tc   = Tkk - tice
        qmk  = qmr(i,j,k)
        prs_sv  = exp(psv_a*Tkk**2 + psv_b*Tkk + psv_c + psv_d/Tkk ) ! Pvap sat, eq A1.1 (Pa)
        ehn_fct = ef_alpha + ef_beta*Pak + ef_gamma*Tc**2 ! enhancement factor (eq. A1.2)
        prs_v   = qmk* Pak/(1.0+qmk*rdry/rvap)            ! vapor pressure (Pa)
        x_v     = prs_v/prs_sv * ehn_fct * prs_sv/Pak     ! molar fraction of water vapor (eq. A1.3)
!       Compressibility factor (eq A1.4 from Picard et al 2008)
        cmpr    = 1.0_kind_real - (Pak/Tkk) * (cpf_a0 + cpf_a1*Tc + cpf_a2*Tc**2 &
                    + (cpf_b0 + cpf_b1*Tc)*x_v + (cpf_c0 + cpf_c1*Tc)*x_v**2 ) &
                    + (Pak**2/Tkk**2) * (cpf_d + cpf_e*x_v**2)
        dz      = rdry/grav * Tvk * cmpr * log(dpk)
        if ( k == geom%npz) then
          gph(i,j,k) = phis(i,j)/grav + dz
        else
          gph(i,j,k) = gph(i,j,k+1) + dz
        end if
      end do ! end k loop

    end do ! end i loop
    end do ! end j loop

else  ! not use compressivity

  do j = jsc,jec
  do i = isc,iec

     k  = geom%npz
     dz         = rdry/grav * Tv(i,j,k) * log(prsi(i,j,k+1)/prs(i,j,k))
     gph(i,j,k) = phis(i,j)/grav + dz

     do k = geom%npz-1, 1, -1
        dz         = rdry/grav * 0.5_kind_real * (Tv(i,j,k+1)+Tv(i,j,k)) * log(prs(i,j,k+1)/prs(i,j,k))
        gph(i,j,k) = gph(i,j,k+1) + dz
     end do

   end do
   end do

end if

end subroutine geop_height

subroutine geop_height_levels(geom,prs,prsi,T,q,phis,use_compress,gphi)

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model

real(kind_real), intent(in ) :: prs(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)    !mid layerpressure
real(kind_real), intent(in ) :: prsi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !interface pressure
real(kind_real), intent(in ) :: phis(geom%isc:geom%iec,geom%jsc:geom%jec)              !Surface geopotential (grav*Z_sfc)
real(kind_real), intent(in ) :: T(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
real(kind_real), intent(in ) :: q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)      !specific humidity
real(kind_real), intent(out) :: gphi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !geopotential height at interface levels (m)

!locals
real(kind_real)       :: Tv(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
real(kind_real)       :: qmr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) ! mixing ratio|kg/kg
logical               :: use_compress

integer               :: isc,iec,jsc,jec,npz,i,j,k
real(kind=kind_real)  :: Tkk,Tvk,Tc, qmk,Pak,dpk,dz
real(kind=kind_real)  :: prs_sv, prs_v
real(kind=kind_real)  :: ehn_fct,x_v,cmpr

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz
!get qmr--mixing ratio and virtual temeprature
qmr(isc:iec,jsc:jec,:) = q(isc:iec,jsc:jec,:)/(1.0 - q(isc:iec,jsc:jec,:))
Tv(isc:iec,jsc:jec,:)  = T(isc:iec,jsc:jec,:)*(1.0 + zvir*qmr(isc:iec,jsc:jec,:))

if (use_compress) then

!  Compute compressibility factor (Picard et al 2008) and geopotential heights at midpoint
  do j = jsc,jec
  do i = isc,iec

     gphi(i,j,geom%npz+1) = phis(i,j)/grav !phis is gh? or gopoential
     do k = geom%npz, 1, -1
        if ( k == 1) then
           Pak  = exp(0.5_kind_real*(log(prsi(i,j,k+1)*0.01)+log(prs(i,j,k)*0.01)))
           dpk  = prsi(i,j,k+1)/prs(i,j,k)
        else
           Pak  = exp(0.5_kind_real*(log(prsi(i,j,k+1)*0.01)+log(prsi(i,j,k)*0.01)))
           dpk  = prsi(i,j,k+1)/prsi(i,j,k)
        end if
        Tkk  = T(i,j,k)
        Tvk  = Tv(i,j,k)
        Tc   = Tkk - tice
        qmk  = qmr(i,j,k)
        prs_sv  = exp(psv_a*Tkk**2 + psv_b*Tkk + psv_c + psv_d/Tkk ) ! Pvap sat, eq A1.1 (Pa)
        ehn_fct = ef_alpha + ef_beta*Pak + ef_gamma*Tc**2 ! enhancement factor (eq. A1.2)
        prs_v   = qmk* Pak/(1.0+qmk*rdry/rvap)            ! vapor pressure (Pa)
        x_v     = prs_v/prs_sv * ehn_fct * prs_sv/Pak     ! molar fraction of water vapor (eq. A1.3)
!       Compressibility factor (eq A1.4 from Picard et al 2008)

        cmpr    = 1.0_kind_real - (Pak/Tkk) * (cpf_a0 + cpf_a1*Tc + cpf_a2*Tc**2 &
                    + (cpf_b0 + cpf_b1*Tc)*x_v + (cpf_c0 + cpf_c1*Tc)*x_v**2 ) &
                    + (Pak**2/Tkk**2) * (cpf_d + cpf_e*x_v**2)
        dz      = rdry/grav * Tvk * cmpr * log(dpk)
        gphi(i,j,k) = gphi(i,j,k+1) + dz
      end do ! end k loop

  end do ! end i loop
  end do ! end j loop

else  ! not use compressivity
  do j = jsc,jec
  do i = isc,iec

     gphi(i,j, geom%npz+1) = phis(i,j)/grav

     do k = geom%npz, 1, -1
        if(k ==1) then
          dz         = rdry/grav * Tv(i,j,k) * log(prsi(i,j,k+1)/prs(i,j,k))
        else
          dz         = rdry/grav * Tv(i,j,k) * log(prsi(i,j,k+1)/prsi(i,j,k))
        end if
        gphi(i,j,k) = gphi(i,j,k+1) + dz
     end do

   end do
   end do

end if

end subroutine geop_height_levels
end module height_vt_mod
