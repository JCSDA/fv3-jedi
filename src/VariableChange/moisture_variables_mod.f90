! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on moisture variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module moisture_vt_mod

use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_constants_mod, only: rdry,grav,tice,zvir

implicit none
private

public crtm_ade_efr
public crtm_mixratio
public crtm_mixratio_tl
public crtm_mixratio_ad
public rh_to_q
public rh_to_q_tl
public rh_to_q_ad
public q_to_rh
public q_to_rh_tl
public q_to_rh_ad
public ESINIT
public dqsat

contains

!>----------------------------------------------------------------------------
!> Compute cloud area density and effective radius for the crtm --------------
!>----------------------------------------------------------------------------

subroutine crtm_ade_efr( geom,p,T,delp,sea_frac,q,ql,qi,ql_ade,qi_ade,ql_efr,qi_efr)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom
real(kind=kind_real), intent(in)  :: p(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Pressure | Pa
real(kind=kind_real), intent(in)  :: t(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Temperature | K
real(kind=kind_real), intent(in)  :: delp(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Layer thickness | Pa
real(kind=kind_real), intent(in)  :: sea_frac(geom%isc:geom%iec,geom%jsc:geom%jec)          !Sea fraction | 1
real(kind=kind_real), intent(in)  :: q(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Specific humidity | kg/kg
real(kind=kind_real), intent(in)  :: ql(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of cloud liquid water | kg/kg
real(kind=kind_real), intent(in)  :: qi(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of cloud ice water | kg/kg

real(kind=kind_real), intent(out) :: ql_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for cloud liquid water | kg/m^2
real(kind=kind_real), intent(out) :: qi_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for cloud ice water | kg/m^2
real(kind=kind_real), intent(out) :: ql_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for cloud liquid water | microns
real(kind=kind_real), intent(out) :: qi_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for cloud ice | microns

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
logical, allocatable :: seamask(:,:)
real(kind=kind_real) :: tem1, tem2, tem3, kgkg_to_kgm2


! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


! Set outputs to zero
! -------------------
ql_ade = 0.0_kind_real
qi_ade = 0.0_kind_real
ql_efr = 0.0_kind_real
qi_efr = 0.0_kind_real


! Sea mask
! --------
allocate(seamask(isc:iec,jsc:jec))
seamask = .false.
do j = jsc,jec
  do i = isc,iec
     seamask(i,j) = min(max(0.0_kind_real,sea_frac(i,j)),1.0_kind_real)  >= 0.99_kind_real
  enddo
enddo


! Convert clouds from kg/kg into kg/m^2
! -------------------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
       if (seamask(i,j)) then

         kgkg_to_kgm2 = delp(i,j,k) / grav
         ql_ade(i,j,k) = max(ql(i,j,k),0.0_kind_real) * kgkg_to_kgm2
         qi_ade(i,j,k) = max(qi(i,j,k),0.0_kind_real) * kgkg_to_kgm2

         if (t(i,j,k) - tice > -20.0_kind_real) then
            ql_ade(i,j,k) = max(1.001_kind_real*1.0E-6_kind_real,ql_ade(i,j,k))
         endif

         if (t(i,j,k) < tice) then
            qi_ade(i,j,k) = max(1.001_kind_real*1.0E-6_kind_real,qi_ade(i,j,k))
         endif

       endif
    enddo
  enddo
enddo

! Cloud liquid water effective radius
! -----------------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
       if (seamask(i,j)) then
         tem1 = max(0.0_kind_real,(tice-t(i,j,k))*0.05_kind_real)
         ql_efr(i,j,k) = 5.0_kind_real + 5.0_kind_real * min(1.0_kind_real, tem1)
       endif
    enddo
  enddo
enddo


! Cloud ice water effective radius
! ---------------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec

      if (seamask(i,j)) then

        tem2 = t(i,j,k) - tice
        tem1 = grav/rdry
        tem3 = tem1 * qi_ade(i,j,k) * (p(i,j,k)/delp(i,j,k))/t(i,j,k) * (1.0_kind_real + zvir * q(i,j,k))

        if (tem2 < -50.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.917_kind_real)*tem3**0.109_kind_real
        elseif (tem2 < -40.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.337_kind_real)*tem3**0.08_kind_real
        elseif (tem2 < -30.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.208_kind_real)*tem3**0.055_kind_real
        else
           qi_efr(i,j,k) =  (1250._kind_real/9.387_kind_real)*tem3**0.031_kind_real
        endif

      endif

    enddo
  enddo
enddo 


ql_efr = max(0.0_kind_real,ql_efr)
qi_efr = max(0.0_kind_real,qi_efr)

deallocate(seamask)

end subroutine crtm_ade_efr

!----------------------------------------------------------------------------
! Compute mixing ratio from specific humidity -------------------------------
!----------------------------------------------------------------------------

subroutine crtm_mixratio(geom,q,qmr)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in ) :: geom
real(kind=kind_real), intent(in ) :: q  (geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(out) :: qmr(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Mixing ratio | 1

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
real(kind=kind_real) :: c3


! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


! Convert specific humidity
! -------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
       c3 = 1.0_kind_real / (1.0_kind_real - q(i,j,k))
       qmr(i,j,k) = 1000.0_kind_real * q(i,j,k) * c3
    enddo
  enddo
enddo 


end subroutine crtm_mixratio

subroutine crtm_mixratio_tl(geom,q,q_tl,qmr_tl)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in ) :: geom
real(kind=kind_real), intent(in ) :: q     (geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(in ) :: q_tl  (geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(out) :: qmr_tl(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Mixing ratio | 1

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
real(kind=kind_real) :: c3, c3_tl


! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


! Convert specific humidity
! -------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
       c3_tl = -((-q_tl(i,j,k))/(1.0_kind_real-q(i,j,k))**2)
       c3 = 1.0_kind_real / (1.0_kind_real - q(i,j,k))
       qmr_tl(i,j,k) = 1000.0_kind_real*(q_tl(i,j,k)*c3+q(i,j,k)*c3_tl)
    enddo
  enddo
enddo 

end subroutine crtm_mixratio_tl

subroutine crtm_mixratio_ad(geom,q,q_ad,qmr_ad)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in )   :: geom
real(kind=kind_real), intent(in )   :: q     (geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(inout) :: q_ad  (geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(inout) :: qmr_ad(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Mixing ratio | 1

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
real(kind=kind_real) :: c3, c3_ad


! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz


! Convert specific humidity
! -------------------------
q_ad = 0.0_kind_real
do k=npz,1,-1
  do j=jec,jsc,-1
    do i=iec,isc,-1
     c3_ad = 1000.0_kind_real*q(i,j,k)*qmr_ad(i,j,k)
     c3 = 1.0_kind_real / (1.0_kind_real - q(i,j,k))
     q_ad(i,j,k) = q_ad(i,j,k) + c3_ad/(1.0_kind_real-q(i, j, k))**2 + 1000.0_kind_real*c3*qmr_ad(i,j,k)
     qmr_ad(i, j, k) = 0.0_kind_real
    end do
  end do
end do

end subroutine crtm_mixratio_ad

!----------------------------------------------------------------------------

subroutine rh_to_q(geom,qsat,rh,q)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       q(i,j,k) = rh(i,j,k) * qsat(i,j,k)
     end do
   end do
 end do

end subroutine rh_to_q

!----------------------------------------------------------------------------

subroutine rh_to_q_tl(geom,qsat,rh,q)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       q(i,j,k) = rh(i,j,k) * qsat(i,j,k)
     end do
   end do
 end do

end subroutine rh_to_q_tl

!----------------------------------------------------------------------------

subroutine rh_to_q_ad(geom,qsat,rh,q)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=geom%npz,1,-1
   do j=geom%jec,geom%jsc,-1
     do i=geom%iec,geom%isc,-1
       rh(i,j,k) = rh(i,j,k) + q(i,j,k) * qsat(i,j,k)
     end do
   end do
 end do

end subroutine rh_to_q_ad

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine q_to_rh(geom,qsat,q,rh)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in)    ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       rh(i,j,k) = q(i,j,k) / qsat(i,j,k)
     end do
   end do
 end do

end subroutine q_to_rh

!----------------------------------------------------------------------------

subroutine q_to_rh_tl(geom,qsat,q,rh)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in)    ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       rh(i,j,k) = q(i,j,k) / qsat(i,j,k)
     end do
   end do
 end do

end subroutine q_to_rh_tl

!----------------------------------------------------------------------------

subroutine q_to_rh_ad(geom,qsat,q,rh)

 implicit none
 type(fv3jedi_geom),   intent(in)    :: geom
 real(kind=kind_real), intent(in)    :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::    q(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) ::   rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 
 integer :: i,j,k
 
 do k=geom%npz,1,-1
   do j=geom%jec,geom%jsc,-1
     do i=geom%iec,geom%isc,-1
       q(i,j,k) = rh(i,j,k) / qsat(i,j,k)
     end do
   end do
 end do

end subroutine q_to_rh_ad

!----------------------------------------------------------------------------

subroutine ESINIT(TABLESIZE,DEGSUBS,TMINTBL,TMAXTBL,ESTBLX)

 IMPLICIT NONE

 integer, intent(in) :: TABLESIZE,DEGSUBS
 real(8), intent(in) :: TMINTBL,TMAXTBL

 !OUTPUT
 real(8), dimension(TABLESIZE) :: ESTBLX

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMIX = -20.0

 real(8), dimension(TABLESIZE) :: ESTBLE, ESTBLW

 integer :: I
 real(8)    :: T, DELTA_T

 DELTA_T = 1.0/DEGSUBS

 do I=1,TABLESIZE

    T = (I-1)*DELTA_T + TMINTBL

    if(T>ZEROC) then
       call QSATLQU0(ESTBLE(I),T,TMAXTBL)
    else
       call QSATICE0(ESTBLE(I),T)
    end if

    call QSATLQU0(ESTBLW(I),T,TMAXTBL)

    T = T-ZEROC
    if(T>=TMIX .and. T<0.0) then
       ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
    else
       ESTBLX(I) = ESTBLE(I)
    end if

 end do

 end subroutine ESINIT


subroutine QSATLQU0(QS,TL,TMAXTBL)
!SUPERSATURATED AS LIQUID

 IMPLICIT NONE

 !INPUTS
 real(8) :: TL, TMAXTBL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC   = 273.16
 real(8), parameter :: TMINLQU = ZEROC - 40.0

 real(8),  parameter :: B6 = 6.136820929E-11*100.0
 real(8),  parameter :: B5 = 2.034080948E-8 *100.0
 real(8),  parameter :: B4 = 3.031240396E-6 *100.0
 real(8),  parameter :: B3 = 2.650648471E-4 *100.0
 real(8),  parameter :: B2 = 1.428945805E-2 *100.0
 real(8),  parameter :: B1 = 4.436518521E-1 *100.0
 real(8),  parameter :: B0 = 6.107799961E+0 *100.0

 real(8) :: TX, EX, TI, TT

 TX = TL

 if    (TX<TMINLQU) then
    TI = TMINLQU
 elseif(TX>TMAXTBL) then
    TI = TMAXTBL
 else
    TI = TX
 end if

 TT = TI-ZEROC  !Starr polynomial fit
 EX = (TT*(TT*(TT*(TT*(TT*(TT*B6+B5)+B4)+B3)+B2)+B1)+B0)

 TL = TX
 QS = EX

 return

end subroutine QSATLQU0


subroutine QSATICE0(QS,TL)
!SUPERSATURATED AS ICE

 IMPLICIT NONE   
      
 !INPUTS
 real(8) :: TL

 !OUTPUTS
 real(8) :: QS

 !LOCALS
 real(8), parameter :: ZEROC = 273.16, TMINSTR = -95.0
 real(8), parameter :: TMINICE = ZEROC + TMINSTR

 real(8), parameter :: TSTARR1 = -75.0, TSTARR2 = -65.0, TSTARR3 = -50.0,  TSTARR4 = -40.0

 real(8),  parameter :: BI6= 1.838826904E-10*100.0
 real(8),  parameter :: BI5= 4.838803174E-8 *100.0
 real(8),  parameter :: BI4= 5.824720280E-6 *100.0
 real(8),  parameter :: BI3= 4.176223716E-4 *100.0
 real(8),  parameter :: BI2= 1.886013408E-2 *100.0
 real(8),  parameter :: BI1= 5.034698970E-1 *100.0
 real(8),  parameter :: BI0= 6.109177956E+0 *100.0
 real(8),  parameter :: S16= 0.516000335E-11*100.0
 real(8),  parameter :: S15= 0.276961083E-8 *100.0
 real(8),  parameter :: S14= 0.623439266E-6 *100.0
 real(8),  parameter :: S13= 0.754129933E-4 *100.0
 real(8),  parameter :: S12= 0.517609116E-2 *100.0
 real(8),  parameter :: S11= 0.191372282E+0 *100.0
 real(8),  parameter :: S10= 0.298152339E+1 *100.0
 real(8),  parameter :: S26= 0.314296723E-10*100.0
 real(8),  parameter :: S25= 0.132243858E-7 *100.0
 real(8),  parameter :: S24= 0.236279781E-5 *100.0
 real(8),  parameter :: S23= 0.230325039E-3 *100.0
 real(8),  parameter :: S22= 0.129690326E-1 *100.0
 real(8),  parameter :: S21= 0.401390832E+0 *100.0
 real(8),  parameter :: S20= 0.535098336E+1 *100.0

 real(8) :: TX, TI, TT, W, EX

 TX = TL

 if (TX<TMINICE) then
    TI = TMINICE
 elseif(TX>ZEROC  ) then
    TI = ZEROC
 else
    TI = TX
 end if

 TT = TI - ZEROC
 if (TT < TSTARR1 ) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10)
 elseif(TT >= TSTARR1 .and. TT < TSTARR2) then
     W = (TSTARR2 - TT)/(TSTARR2-TSTARR1)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S16+S15)+S14)+S13)+S12)+S11)+S10) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR2 .and. TT < TSTARR3) then
     EX = (TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20)
 elseif(TT >= TSTARR3 .and. TT < TSTARR4) then
     W = (TSTARR4 - TT)/(TSTARR4-TSTARR3)
     EX =       W *(TT*(TT*(TT*(TT*(TT*(TT*S26+S25)+S24)+S23)+S22)+S21)+S20) &
              + (1.-W)*(TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 else
     EX = (TT*(TT*(TT*(TT*(TT*(TT*BI6+BI5)+BI4)+BI3)+BI2)+BI1)+BI0)
 endif

 QS = EX

 return
 
end subroutine QSATICE0

!----------------------------------------------------------------------------

subroutine dqsat(geom,temp,pmid,degsubs,tmintbl,tmaxtbl,tablesize,estblx,dqsi,qssi)

!computes saturation vapour pressure qssi and gradient w.r.t temperature dqsi.
!inputs are temperature and plo (pressure at t-levels)
!vales are computed from look-up talbe (piecewise linear)

 use fv3jedi_constants_mod, only: h2omw, airmw

 implicit none

 !inputs
 type(fv3jedi_geom),   intent(in) :: geom
 integer,              intent(in) :: degsubs
 real(8),              intent(in) :: tmintbl, tmaxtbl
 integer,              intent(in) :: tablesize
 real(kind=kind_real), intent(in) :: temp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in) :: pmid(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(8),              intent(in) :: estblx(tablesize)

 !outputs
 real(kind=kind_real), intent(inout) :: dqsi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(inout) :: qssi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 !locals
 real(8), parameter :: max_mixing_ratio = 1.0_8
 real(8) :: esfac

 integer :: i, j, k
 real(8) :: tt, ti, dqq, qq, dd
 integer :: it
 real(8) :: temp8, pmid8, qssi8, dqsi8

 dqsi = 0.0_kind_real
 qssi = 0.0_kind_real

 esfac = real(h2omw,8)/real(airmw,8)

 do k=1,geom%npz
    do j=geom%jsc,geom%jec
      do i=geom%isc,geom%iec

          temp8 = real(temp(i,j,k),8)
          pmid8 = real(pmid(i,j,k),8)

          if (temp8<=tmintbl) then
             ti = tmintbl
          elseif(temp8>=tmaxtbl-.001_8) then
             ti = tmaxtbl-.001_8
          else
             ti = temp8
          end if

          tt = (ti - tmintbl)*real(degsubs,8)+1.0_8
          it = int(tt)

          dqq =  estblx(it+1) - estblx(it)
          qq  =  (tt-real(it,8))*dqq + estblx(it)

          if (pmid8 <= qq) then
             qssi8 = max_mixing_ratio
             dqsi8 = 0.0
          else
             dd = 1.0/(pmid8 - (1.0-esfac)*qq)
             qssi8 = esfac*qq*dd
             dqsi8 = (esfac*degsubs)*dqq*pmid8*(dd*dd)
          end if

          dqsi(i,j,k) = real(dqsi8,kind_real)
          qssi(i,j,k) = real(qssi8,kind_real)

       end do
    end do
 end do

end subroutine dqsat

!----------------------------------------------------------------------------

end module moisture_vt_mod
