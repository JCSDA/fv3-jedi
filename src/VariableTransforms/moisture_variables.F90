! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on moisture variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module moisture_vt_mod

use kinds, only: kind_real
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_constants, only: rdry,grav,tice,zvir

implicit none
public

contains

!----------------------------------------------------------------------------
! Compute cloud area density and effective radius for the crtm --------------
!----------------------------------------------------------------------------

subroutine crtm_ade_efr( geom,p,T,delp,sea_frac,q,ql,qi,ql_ade,qi_ade,ql_efr,qi_efr)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom
real(kind=kind_real), intent(in)  :: p(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)     !Pressure | Pa
real(kind=kind_real), intent(in)  :: t(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)     !Temperature | K
real(kind=kind_real), intent(in)  :: delp(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)  !Layer thickness | Pa
real(kind=kind_real), intent(in)  :: sea_frac(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)          !Sea fraction | 1
real(kind=kind_real), intent(in)  :: q(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)     !Specific humidity | kg/kg
real(kind=kind_real), intent(in)  :: ql(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)    !Mixing ratio of cloud liquid water | kg/kg
real(kind=kind_real), intent(in)  :: qi(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)    !Mixing ratio of cloud ice water | kg/kg

real(kind=kind_real), intent(out) :: ql_ade(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz) !area density for cloud liquid water | kg/m^2
real(kind=kind_real), intent(out) :: qi_ade(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz) !area density for cloud ice water | kg/m^2
real(kind=kind_real), intent(out) :: ql_efr(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz) !efr for cloud liquid water | microns
real(kind=kind_real), intent(out) :: qi_efr(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz) !efr for cloud ice | microns

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
logical, allocatable :: seamask(:,:)
real(kind=kind_real) :: tem1, tem2, tem3, kgkg_to_kgm2


! Grid convenience
! ----------------
isc = geom%bd%isc
iec = geom%bd%iec
jsc = geom%bd%jsc
jec = geom%bd%jec
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
         ql_ade(i,j,k) = ql_ade(i,j,k) * kgkg_to_kgm2
         qi_ade(i,j,k) = qi_ade(i,j,k) * kgkg_to_kgm2
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

ql_efr = max(0.0_kind_real,ql_efr)
   

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
real(kind=kind_real), intent(in ) :: q  (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)  !Specific humidity | kg/kg
real(kind=kind_real), intent(out) :: qmr(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed, 1:geom%npz)  !Mixing ratio | 1

!Locals
integer :: isc,iec,jsc,jec,npz
integer :: i,j,k
real(kind=kind_real) :: c3


! Grid convenience
! ----------------
isc = geom%bd%isc
iec = geom%bd%iec
jsc = geom%bd%jsc
jec = geom%bd%jec
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

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module moisture_vt_mod
