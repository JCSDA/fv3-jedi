! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on wind variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module wind_vt_mod

use fv3jedi_constants, only: pi, rad2deg
use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real

implicit none
public

contains

!----------------------------------------------------------------------------
! Lowest model level winds to 10m speed and direction -----------------------
!----------------------------------------------------------------------------

subroutine sfc_10m_winds(geom,usrf,vsrf,f10r,spd10m,dir10m)

 implicit none

 !Arguments
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: usrf(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Lowest model level u m/s
 real(kind=kind_real), intent(in ) :: vsrf(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Lowest model level v m/s
 real(kind=kind_real), intent(in ) :: f10r(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Ratio of lowest level to 10m
 real(kind=kind_real), intent(out) :: spd10m(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !10m wind speed u m/s
 real(kind=kind_real), intent(out) :: dir10m(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !10m model wind direction

 !Locals
 integer :: isc,iec,jsc,jec,i,j
 integer :: iquad
 real(kind=kind_real) :: windangle, windratio
 real(kind=kind_real), parameter :: windscale = 999999.0_kind_real
 real(kind=kind_real), parameter :: windlimit = 0.0001_kind_real
 real(kind=kind_real), parameter :: quadcof(4,2) = reshape((/ 0.0_kind_real,  1.0_kind_real, 1.0_kind_real,  2.0_kind_real,  &
                                                              1.0_kind_real, -1.0_kind_real, 1.0_kind_real, -1.0_kind_real /), (/4, 2/))

 !In GSI these calculations are done after interpolation to obs location

 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 !10m wind speed
 spd10m(isc:iec,jsc:jec) = f10r(isc:iec,jsc:jec)*sqrt( usrf(isc:iec,jsc:jec)*usrf(isc:iec,jsc:jec) + &
                                                      vsrf(isc:iec,jsc:jec)*vsrf(isc:iec,jsc:jec) ) 

 !10m wind direction
 do j = jsc,jec
   do i = isc,iec
     
     if (usrf(i,j)*f10r(i,j) >= 0.0_kind_real .and. vsrf(i,j)*f10r(i,j) >= 0.0_kind_real) iquad = 1
     if (usrf(i,j)*f10r(i,j) >= 0.0_kind_real .and. vsrf(i,j)*f10r(i,j) <  0.0_kind_real) iquad = 2
     if (usrf(i,j)*f10r(i,j) <  0.0_kind_real .and. vsrf(i,j)*f10r(i,j) >= 0.0_kind_real) iquad = 3
     if (usrf(i,j)*f10r(i,j) <  0.0_kind_real .and. vsrf(i,j)*f10r(i,j) <  0.0_kind_real) iquad = 4

     if (abs(vsrf(i,j)*f10r(i,j)) >= windlimit) then
        windratio = (usrf(i,j)*f10r(i,j)) / (vsrf(i,j)*f10r(i,j))
     else
        windratio = 0.0_kind_real 
        if (abs(usrf(i,j)*f10r(i,j)) > windlimit) then 
          windratio = windscale * usrf(i,j)*f10r(i,j) 
        endif 
     endif

     windangle = atan(abs(windratio))

     dir10m(i,j) = rad2deg*(quadcof(iquad,1) * pi + windangle * quadcof(iquad, 2))

   enddo
 enddo

end subroutine sfc_10m_winds

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module wind_vt_mod
