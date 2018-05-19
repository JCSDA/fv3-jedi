! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on wind variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real

implicit none
public

contains

!----------------------------------------------------------------------------
! Surface quantities in the form needed by the crtm ------------------------- 
!----------------------------------------------------------------------------

subroutine crtm_surface( geom, slmsk, sheleg, tsea, vtype, & 
                         stype, vfrac, stc, smc, &
                         land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                         snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                         snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth )

 implicit none

 !Arguments
 type(fv3jedi_geom)  , intent(in )                 :: geom !Geometry for the model
 integer             , intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: slmsk
 real(kind=kind_real), intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: sheleg
 real(kind=kind_real), intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: tsea
 integer             , intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vtype
 integer             , intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: stype
 real(kind=kind_real), intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vfrac
 real(kind=kind_real), intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4) :: stc
 real(kind=kind_real), intent(in ), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4) :: smc
 integer             , intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vegetation_type
 integer             , intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: land_type
 integer             , intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: soil_type
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: water_coverage
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: land_coverage
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: ice_coverage
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: snow_coverage
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: lai
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: water_temperature
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: land_temperature
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: ice_temperature
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: snow_temperature
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: soil_moisture_content
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vegetation_fraction
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: soil_temperature
 real(kind=kind_real), intent(out), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: snow_depth

 !Locals
 integer :: wlis_msk(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Water, land, ice + snow
 integer :: lai_type(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)

 integer :: isc,iec,jsc,jec,i,j

 ! Indices for the CRTM NPOESS EmisCoeff file
 integer, parameter :: invalid_land = 0
 integer, parameter :: compacted_soil = 1
 integer, parameter :: tilled_soil = 2
 integer, parameter :: irrigated_low_vegetation = 5
 integer, parameter :: meadow_grass = 6
 integer, parameter :: scrub = 7
 integer, parameter :: broadleaf_forest = 8
 integer, parameter :: pine_forest = 9
 integer, parameter :: tundra = 10
 integer, parameter :: grass_soil = 11
 integer, parameter :: broadleaf_pine_forest = 12
 integer, parameter :: grass_scrub = 13
 integer, parameter :: urban_concrete = 15
 integer, parameter :: broadleaf_brush = 17
 integer, parameter :: wet_soil = 18
 integer, parameter :: scrub_soil = 19
 integer, parameter :: vegetation_n_types = 13
 integer, parameter :: soil_n_types = 9
 integer, parameter, dimension(0:13) :: jedi_to_crtm=(/compacted_soil, broadleaf_forest, broadleaf_forest, &
                                                       broadleaf_pine_forest, pine_forest, &
                                                       pine_forest, broadleaf_brush, scrub, scrub, scrub_soil, tundra, &
                                                       compacted_soil, tilled_soil, compacted_soil/)
 integer :: itype, istype


!*!*!*!*! TODO TODO TODO
! THE BELOW NEEDS TO BE REPLACED WITH THE PROPER WAY, BY TAKING FRACTIONS OF THE SURROUNDING GRID BOXES
!*!*!*!*! TODO TODO TODO


 ! Convenience
 ! -----------
 isc = geom%bd%isc
 iec = geom%bd%iec
 jsc = geom%bd%jsc
 jec = geom%bd%jec

 ! Zero outputs
 ! ------------
 vegetation_type = 0
 land_type = 0
 soil_type = 0
 water_coverage = 0.0_kind_real
 land_coverage = 0.0_kind_real
 ice_coverage = 0.0_kind_real
 snow_coverage = 0.0_kind_real
 lai = 0.0_kind_real
 water_temperature = 0.0_kind_real
 land_temperature = 0.0_kind_real
 ice_temperature = 0.0_kind_real
 snow_temperature = 0.0_kind_real
 soil_moisture_content = 0.0_kind_real
 vegetation_fraction = 0.0_kind_real
 soil_temperature = 0.0_kind_real
 snow_depth = 0.0_kind_real


 ! Vegation, land and soil types
 ! -----------------------------
 do j = jsc,jec
   do i = isc,iec

     itype  = min(max(0,vtype(i,j)),vegetation_n_types)
     istype = min(max(1,stype(i,j)),soil_n_types)

     vegetation_type(i,j) = jedi_to_crtm(itype)
     land_type(i,j) = max(1,itype)
     soil_type(i,j) = istype
     lai_type(i,j) = itype

   enddo
 enddo


 ! Surface types and temperatures
 ! ------------------------------

 !Surface mask - 0: water | 1: land | 2: ice
 wlis_msk(isc:iec,jsc:jec) = slmsk(isc:iec,jsc:jec)

 do j = jsc,jec
   do i = isc,iec

      !Add if snow | 3: snow
      if(wlis_msk(i,j) >=1 .and. sheleg(i,j) > 1.0_kind_real/10.0_kind_real) wlis_msk(i,j) = 3

      !Set percentages (wrong way to do this but done to get working!)
      if (wlis_msk(i,j) == 0) then
        water_coverage(i,j) = 1.0_kind_real
      elseif  (wlis_msk(i,j) == 1) then
        land_coverage(i,j) = 1.0_kind_real
      elseif (wlis_msk(i,j) == 2) then
        ice_coverage(i,j) = 1.0_kind_real
      elseif (wlis_msk(i,j) == 3) then
        snow_coverage(i,j) = 1.0_kind_real
      endif

      !Set percentages (wrong way to do this but done to get working!)
      if (wlis_msk(i,j) == 0) then
        water_temperature(i,j) = tsea(i,j)
      elseif  (wlis_msk(i,j) == 1) then
        land_temperature(i,j) = tsea(i,j)
      elseif (wlis_msk(i,j) == 2) then
        ice_temperature(i,j) = tsea(i,j)
      elseif (wlis_msk(i,j) == 3) then
        snow_temperature(i,j) = tsea(i,j)
      endif

   enddo
 enddo


 !TODO BY NOW WE SHOULD WE WORKING WITH OBSLOCATION DATA


 ! Bounds on coverage
 ! ------------------
 water_coverage = min(max(0.0_kind_real,water_coverage),1.0_kind_real)
 land_coverage  = min(max(0.0_kind_real,land_coverage ),1.0_kind_real)
 ice_coverage   = min(max(0.0_kind_real,ice_coverage  ),1.0_kind_real)
 snow_coverage  = min(max(0.0_kind_real,snow_coverage ),1.0_kind_real)
     
 ! Get vegetation lai from summer and winter values
 ! ------------------------------------------------
 Lai  = 0.0_kind_real
 do j = jsc,jec
   do i = isc,iec

     if (land_coverage(i,j) > 0.0_kind_real) then

       if (lai_type(i,j) > 0) then
         !call get_lai(data_s,nchanl,nreal,itime,ilate,lai_type,lai(i,j))
       endif     
     
       !For Glacial land ice soil type and vegetation type
       if (Soil_Type(i,j) == 9 .OR. Vegetation_Type(i,j) == 13) then
         ice_coverage(i,j) = min(ice_coverage(i,j) + land_coverage(i,j), 1.0_kind_real)
         land_coverage = 0.0_kind_real
       endif

     endif

   enddo
 enddo

 ! Set some bounds on temperatures
 ! -------------------------------
 water_temperature = max(water_temperature,270._kind_real)
 ice_temperature   = min(ice_temperature  ,280._kind_real)
 snow_temperature  = min(snow_temperature ,280._kind_real)


 soil_moisture_content = smc(:,:,1)
 vegetation_fraction = vfrac
 soil_temperature = stc(:,:,1)
 snow_depth = sheleg

end subroutine crtm_surface

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module surface_vt_mod
