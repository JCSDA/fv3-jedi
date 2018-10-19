module fv3jedi_getvalues_mod

use fckit_mpi_module, only: fckit_mpi_comm
use type_bump, only: bump_type
use ioda_locs_mod, only: ioda_locs
use ufo_vars_mod, only: ufo_vars
use ufo_geovals_mod, only: ufo_geovals

use fv3jedi_constants_mod, only: rad2deg, constoz, grav
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_state_utils_mod, only: fv3jedi_state
use fv3jedi_increment_utils_mod, only: fv3jedi_increment

use surface_vt_mod
use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use height_vt_mod

implicit none

private
public :: getvalues, getvalues_tl, getvalues_ad

!This module contains the routines for interfacing with UFO
!This is where the GeoVaLs get set for the state and increment.

contains

! ------------------------------------------------------------------------------

subroutine getvalues(geom, state, locs, vars, gom, traj)

implicit none
type(fv3jedi_geom),                         intent(inout) :: geom 
type(fv3jedi_state),                        intent(inout) :: state 
type(ioda_locs),                            intent(in)    :: locs 
type(ufo_vars),                             intent(in)    :: vars
type(ufo_geovals),                          intent(inout) :: gom
type(fv3jedi_getvaltraj), optional, target, intent(inout) :: traj

character(len=*), parameter :: myname = 'interp'

type(bump_type), target  :: bump
type(bump_type), pointer :: pbump
logical, target :: bump_alloc
logical, pointer :: pbumpa

integer :: ii, jj, ji, jvar, jlev, ngrid, nobs
real(kind=kind_real), allocatable :: mod_state(:,:)
real(kind=kind_real), allocatable :: obs_state(:,:)
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
integer :: nvl
logical :: do_interp

integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,npz,i,j

integer :: nt, trcount
character(len=20) :: trname

!Local pressure variables
real(kind=kind_real), allocatable :: prsi(:,:,:) !Pressure Pa, interfaces
real(kind=kind_real), allocatable :: prs (:,:,:) !Pressure Pa, midpoint
real(kind=kind_real), allocatable :: logp(:,:,:) !Log(pressue), (Pa) midpoint

!Local CRTM moisture variables
real(kind=kind_real), allocatable :: ql_ade(:,:,:) !Cloud liq water kgm^2
real(kind=kind_real), allocatable :: qi_ade(:,:,:) !Cloud ice water kgm^2
real(kind=kind_real), allocatable :: ql_efr(:,:,:) !Cloud effective radius microns
real(kind=kind_real), allocatable :: qi_efr(:,:,:) !Cloud effective radium microns
real(kind=kind_real), allocatable :: qmr(:,:,:)    !Moisture mixing ratio
real(kind=kind_real), allocatable :: water_coverage_m(:,:) !Water coverage, model grid

!Local CRTM surface variables
integer             , allocatable :: vegetation_type(:)          !Index of vege type              | surface(1)%Vegetation_Type
integer             , allocatable :: land_type(:)                !Index of land type              | surface(1)%Land_Type
integer             , allocatable :: soil_type(:)                !Index of soil type              | surface(1)%Soil_Type
real(kind=kind_real), allocatable :: wind_speed(:)               !10 meter wind speed m/s         | surface(1)%wind_speed
real(kind=kind_real), allocatable :: wind_direction(:)           !10 meter wind direction degrees | surface(1)%wind_direction
real(kind=kind_real), allocatable :: water_coverage(:)           !Fraction of water coverage      | surface(1)%water_coverage
real(kind=kind_real), allocatable :: land_coverage(:)            !Fraction of land coverage       | surface(1)%land_coverage
real(kind=kind_real), allocatable :: ice_coverage(:)             !Fraction of ice coverage        | surface(1)%ice_coverage
real(kind=kind_real), allocatable :: snow_coverage(:)            !Fraction of snow coverage       | surface(1)%snow_coverage
real(kind=kind_real), allocatable :: lai(:)                      !Leaf area index                 ! surface(1)%lai
real(kind=kind_real), allocatable :: water_temperature(:)        !Water temp (K)                  | surface(1)%water_temperature    
real(kind=kind_real), allocatable :: land_temperature(:)         !Land temp (K)                   | surface(1)%land_temperature     
real(kind=kind_real), allocatable :: ice_temperature(:)          !Ice temp (K)                    | surface(1)%ice_temperature      
real(kind=kind_real), allocatable :: snow_temperature(:)         !Snow temp (K)                   | surface(1)%snow_temperature     
real(kind=kind_real), allocatable :: soil_moisture_content(:)    !Soil moisture content           | surface(1)%soil_moisture_content
real(kind=kind_real), allocatable :: vegetation_fraction(:)      !Vegetation fraction             | surface(1)%vegetation_fraction  
real(kind=kind_real), allocatable :: soil_temperature(:)         !Soil temperature                | surface(1)%soil_temperature     
real(kind=kind_real), allocatable :: snow_depth(:)               !Snow depth                      | surface(1)%snow_depth           
logical,  parameter                ::use_compress = .true.  !!could be a fv3 namelist option?


! Grid convenience
! ----------------
isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
isd = state%isd
ied = state%ied
jsd = state%jsd
jed = state%jed
npz = state%npz

ngrid = (iec-isc+1)*(jec-jsc+1)
nobs = locs%nlocs 

! Initialize the interpolation
! ----------------------------
if (present(traj)) then

  pbump => traj%bump

  if (.not. traj%lalloc) then
  
     traj%ngrid = ngrid
     traj%nobs = nobs
   
     if (.not.allocated(traj%t)) allocate(traj%t(isc:iec,jsc:jec,1:npz))
     if (.not.allocated(traj%q)) allocate(traj%q(isc:iec,jsc:jec,1:npz))
  
     traj%t = state%t
     traj%q = state%q
 
     pbumpa => traj%lalloc

  endif

else

  pbump => bump
  bump_alloc = .false.
  pbumpa => bump_alloc

endif

if (.not. pbumpa) then
   call initialize_bump(geom, locs, pbump)
   pbumpa = .true.
endif


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_state(ngrid,1))
allocate(obs_state(nobs,1))

! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))

! Get pressures at edge, center & log center
! ------------------------------------------
allocate(prsi(isc:iec,jsc:jec,npz+1))
allocate(prs (isc:iec,jsc:jec,npz  ))
allocate(logp(isc:iec,jsc:jec,npz  ))

call delp_to_pe_p_logp(geom,state%delp,prsi,prs,logp)

! Get CRTM surface variables
! ----------------------
allocate(wind_speed(nobs))
allocate(wind_direction(nobs))
allocate(land_type(nobs))
allocate(vegetation_type(nobs))
allocate(soil_type(nobs))
allocate(water_coverage(nobs))
allocate(land_coverage(nobs))
allocate(ice_coverage(nobs))
allocate(snow_coverage(nobs))
allocate(lai(nobs))
allocate(water_temperature(nobs))
allocate(land_temperature(nobs))
allocate(ice_temperature(nobs))
allocate(snow_temperature(nobs))
allocate(soil_moisture_content(nobs))
allocate(vegetation_fraction(nobs))
allocate(soil_temperature(nobs))
allocate(snow_depth(nobs))

wind_speed = 0.0_kind_real
wind_direction = 0.0_kind_real
land_type = 0
vegetation_type = 0
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

if (state%havecrtmfields) then
  !TODO only if a radiance
  call crtm_surface( geom, nobs, ngrid, locs%lat(:), locs%lon(:), &
                     state%slmsk, state%sheleg, state%tsea, state%vtype, &
                     state%stype, state%vfrac, state%stc, state%smc, state%snwdph, &
                     state%u_srf,state%v_srf,state%f10m, &
                     land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                     snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                     snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                     wind_speed, wind_direction )
endif


! Get CRTM moisture variables
! ---------------------------
allocate(ql_ade(isc:iec,jsc:jec,npz))
allocate(qi_ade(isc:iec,jsc:jec,npz))
allocate(ql_efr(isc:iec,jsc:jec,npz))
allocate(qi_efr(isc:iec,jsc:jec,npz))
allocate(qmr(isc:iec,jsc:jec,npz))
allocate(water_coverage_m(isc:iec,jsc:jec))

ql_ade = 0.0_kind_real
qi_ade = 0.0_kind_real
ql_efr = 0.0_kind_real
qi_efr = 0.0_kind_real

if (state%havecrtmfields) then

  !TODO Is it water_coverage or sea_coverage fed in here?
  water_coverage_m = 0.0_kind_real
  do j = jsc,jec
    do i = isc,iec
      if (state%slmsk(i,j) == 0) water_coverage_m(i,j) = 1.0_kind_real
    enddo
  enddo
  
  call crtm_ade_efr( geom,prsi,state%t,state%delp,water_coverage_m,state%q, &
                     state%ql,state%qi,ql_ade,qi_ade,ql_efr,qi_efr )
  
  call crtm_mixratio(geom,state%q,qmr)

endif


! Variable transforms and interpolate to obs locations
! ----------------------------------------------------

do jvar = 1, vars%nv

  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  ! Convert to observation variables/units
  ! --------------------------------------
  select case (trim(vars%fldnames(jvar)))

  case ("upper_air_u_component")

    nvl = npz
    do_interp = .true.
    geovalm = state%ua
    geoval => geovalm

  case ("upper_air_v_component")

    nvl = npz
    do_interp = .true.
    geovalm = state%va
    geoval => geovalm

  case ("temperature")

    nvl = npz
    do_interp = .true.
    geovalm = state%t
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geovalm = state%q
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    call T_to_Tv(geom,state%t,state%q,geovalm)
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

    nvl = npz
    do_interp = .true.
    geovalm = log(0.001_kind_real) + logp !to kPa
    geoval => geovalm

  case ("humidity_mixing_ratio")
    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm

  case ("air_pressure")

    nvl = npz
    do_interp = .true.
    geovalm = prs / 100.0_kind_real !to hPa
    geoval => geovalm

  case ("air_pressure_levels")

    nvl = npz + 1
    do_interp = .true.
    geovale = prsi / 100.0_kind_real !to hPa
    geoval => geovale

  case ("geopotential_height")

    call geop_height(geom,prs,prsi,state%t,state%q,state%phis,use_compress,geovalm)
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("geopotential_height_levels")

    call geop_height_levels(geom,prs,prsi,state%t,state%q,state%phis,use_compress,geovale)
    nvl = npz + 1
    do_interp = .true.
    geoval => geovale

  case ("sfc_geopotential_height")

    nvl = 1
    do_interp = .true.
    geovalm(:,:,1) = state%phis / grav
    geoval => geovalm

  case ("mass_concentration_of_ozone_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = state%o3 * constoz
   geoval => geovalm

  case ("mass_concentration_of_carbon_dioxide_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = 407.0_kind_real !Just a constant for now
   geoval => geovalm

  case ("atmosphere_mass_content_of_cloud_liquid_water")

   nvl = npz
   do_interp = .true.
   geovalm = ql_ade
   geoval => geovalm

  case ("atmosphere_mass_content_of_cloud_ice")

   nvl = npz
   do_interp = .true.
   geovalm = qi_ade
   geoval => geovalm

  case ("effective_radius_of_cloud_liquid_water_particle")

   nvl = npz
   do_interp = .true.
   geovalm = ql_efr
   geoval => geovalm

  case ("effective_radius_of_cloud_ice_particle")

   nvl = npz
   do_interp = .true.
   geovalm = qi_efr
   geoval => geovalm

  case ("Water_Fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_coverage

  case ("Land_Fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_coverage

  case ("Ice_Fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_coverage

  case ("Snow_Fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_coverage

  case ("Water_Temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_temperature

  case ("Land_Temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_temperature

  case ("Ice_Temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_temperature

  case ("Snow_Temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_temperature

  case ("Snow_Depth")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_depth

  case ("Vegetation_Fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = vegetation_fraction

  case ("Sfc_Wind_Speed")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_speed

  case ("Sfc_Wind_Direction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_direction

  case ("Lai")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = lai

  case ("Soil_Moisture")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_moisture_content

  case ("Soil_Temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_temperature

  case ("Land_Type_Index")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(land_type,kind_real)

  case ("Vegetation_Type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(vegetation_type,kind_real)

  case ("Soil_Type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(soil_type,kind_real)

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select


  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nobs,nvl)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call interp_checks(myname, locs, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = 1, nvl
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          mod_state(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call pbump%apply_obsop(mod_state,obs_state)
      gom%geovals(jvar)%vals(jlev,:) = obs_state(:,1)
    enddo
  else
    gom%geovals(jvar)%vals(nvl,:) = obs_state(:,1)
  endif

  nullify(geoval)

enddo

if (.not. present(traj)) then
  call pbump%dealloc
endif

nullify(pbump)

deallocate(mod_state)
deallocate(obs_state)
deallocate(geovale)
deallocate(geovalm)
deallocate(prsi)
deallocate(prs )
deallocate(logp)
deallocate(wind_speed)
deallocate(wind_direction)
deallocate(land_type)
deallocate(vegetation_type)
deallocate(soil_type)
deallocate(water_coverage)
deallocate(land_coverage)
deallocate(ice_coverage)
deallocate(snow_coverage)
deallocate(lai)
deallocate(water_temperature)
deallocate(land_temperature)
deallocate(ice_temperature)
deallocate(snow_temperature)
deallocate(soil_moisture_content)
deallocate(vegetation_fraction)
deallocate(soil_temperature)
deallocate(snow_depth)
deallocate(ql_ade)
deallocate(qi_ade)
deallocate(ql_efr)
deallocate(qi_efr)
deallocate(qmr)
deallocate(water_coverage_m)

!write(*,*)'interp geovals t min, max= ',minval(gom%geovals(1)%vals(:,:)),maxval(gom%geovals(1)%vals(:,:))
!write(*,*)'interp geovals p min, max= ',minval(gom%geovals(2)%vals(:,:)),maxval(gom%geovals(2)%vals(:,:))

end subroutine getvalues

! ------------------------------------------------------------------------------

subroutine getvalues_tl(geom, inc, locs, vars, gom, traj)

implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
type(fv3jedi_increment),      intent(inout) :: inc 
type(ioda_locs),          intent(in)    :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_tl'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


! Grid convenience
! ----------------
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
isd = inc%isd
ied = inc%ied
jsd = inc%jsd
jed = inc%jed
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))


! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv
 
  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  select case (trim(vars%fldnames(jvar)))
   
  case ("upper_air_u_component")
  
    nvl = npz
    do_interp = .true.
    geovalm = inc%ua
    geoval => geovalm

  case ("upper_air_v_component")
  
    nvl = npz
    do_interp = .true.
    geovalm = inc%va
    geoval => geovalm

  case ("temperature")
  
    nvl = npz
    do_interp = .true.
    geovalm = inc%t
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geovalm = inc%q
    geoval => geovalm

  case ("virtual_temperature")

    nvl = inc%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, inc%t, traj%q, inc%q )
    geovalm = inc%t
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = inc%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, inc%q, geovalm)
    geoval => geovalm  

  case ("air_pressure")

  case ("air_pressure_levels")
 
  case ("geopotential_height")

  case ("geopotential_height_levels")

  case ("sfc_geopotential_height")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,traj%nobs,nvl)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call interp_checks(myname, locs, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = 1, nvl
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          mod_increment(ii, 1) = geoval(ji, jj, jlev)
        enddo
      enddo
      call traj%bump%apply_obsop(mod_increment,obs_increment)
      gom%geovals(jvar)%vals(jlev,:) = obs_increment(:,1)
    enddo
  else
    gom%geovals(jvar)%vals(nvl,:) = obs_increment(:,1)
  endif

  nullify(geoval)

enddo

deallocate(geovalm,geovale)

deallocate(mod_increment)
deallocate(obs_increment)

end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(geom, inc, locs, vars, gom, traj)

implicit none
type(fv3jedi_geom),       intent(inout) :: geom 
type(fv3jedi_increment),      intent(inout) :: inc 
type(ioda_locs),           intent(in)   :: locs 
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvaltraj), intent(in)    :: traj

character(len=*), parameter :: myname = 'interp_ad'

integer :: ii, jj, ji, jvar, jlev
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


! Grid convenience
! ----------------
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
isd = inc%isd
ied = inc%ied
jsd = inc%jsd
jed = inc%jed
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(traj%nobs,1))


! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))

geovale = 0.0_kind_real
geovalm = 0.0_kind_real

! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nv

  ! PART 1, do_interp flag
  ! ----------------------
  do_interp = .false.

  select case (trim(vars%fldnames(jvar)))
   
  case ("upper_air_u_component")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("upper_air_v_component")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("temperature")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("geopotential_height")

  case ("geopotential_height_levels")

  case ("sfc_geopotential_height")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  !Part 2, apply adjoint of interp
  !-------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = nvl, 1, -1
      obs_increment(:,1) = gom%geovals(jvar)%vals(jlev,:)
      gom%geovals(jvar)%vals(jlev,:) = 0.0_kind_real
      call traj%bump%apply_obsop_ad(obs_increment,mod_increment)
      ii = 0
      do jj = jsc, jec
        do ji = isc, iec
          ii = ii + 1
          geoval(ji, jj, jlev) = mod_increment(ii, 1)
        enddo
      enddo
    enddo
  else
    obs_increment(:,1) = gom%geovals(jvar)%vals(nvl,:)
  endif

  !Part 3, back to state variables
  !-------------------------------
 
  select case (trim(vars%fldnames(jvar)))
 
  case ("upper_air_u_component")

    inc%ua = geovalm

  case ("upper_air_v_component")

    inc%va = geovalm

  case ("temperature")

    inc%t = geovalm

  case ("specific_humidity")

    inc%q = geovalm

  case ("virtual_temperature")
    
    inc%t = geovalm
    call T_to_Tv_ad(geom, traj%t, inc%t, traj%q, inc%q )

  case ("atmosphere_ln_pressure_coordinate")

  case ("humidity_mixing_ratio")
  
    call crtm_mixratio_ad(geom, traj%q, inc%q, geovalm)

  case ("air_pressure")

  case ("air_pressure_levels")
 
  case ("geopotential_height")

  case ("geopotential_height_levels")

  case ("sfc_geopotential_height")

  case ("mass_concentration_of_ozone_in_air")

  case ("mass_concentration_of_carbon_dioxide_in_air")

  case ("atmosphere_mass_content_of_cloud_liquid_water")

  case ("atmosphere_mass_content_of_cloud_ice")

  case ("effective_radius_of_cloud_liquid_water_particle")

  case ("effective_radius_of_cloud_ice_particle")

  case ("Water_Fraction")

  case ("Land_Fraction")
 
  case ("Ice_Fraction")
 
  case ("Snow_Fraction")
 
  case ("Water_Temperature")
 
  case ("Land_Temperature")
 
  case ("Ice_Temperature")

  case ("Snow_Temperature")

  case ("Snow_Depth")

  case ("Vegetation_Fraction")

  case ("Sfc_Wind_Speed")

  case ("Sfc_Wind_Direction")

  case ("Lai")

  case ("Soil_Moisture")

  case ("Soil_Temperature")

  case ("Land_Type_Index")

  case ("Vegetation_Type")

  case ("Soil_Type")

  case default

    call abor1_ftn(trim(myname)//"unknown variable")

  end select

  geovale = 0.0_kind_real
  geovalm = 0.0_kind_real

enddo

deallocate(mod_increment)
deallocate(obs_increment)

end subroutine getvalues_ad

! ------------------------------------------------------------------------------

subroutine allocate_geovals_vals(gom,jvar,nobs,gvlev)

implicit none
integer, intent(in) :: jvar, nobs, gvlev
type(ufo_geovals), intent(inout) :: gom

! Allocate geovals for this jvar
if (allocated(gom%geovals(jvar)%vals)) deallocate(gom%geovals(jvar)%vals)

allocate(gom%geovals(jvar)%vals(gvlev,nobs))

gom%geovals(jvar)%nval = gvlev
gom%geovals(jvar)%nobs = nobs
gom%geovals(jvar)%vals = 0.0_kind_real

gom%linit  = .true.

end subroutine allocate_geovals_vals

! ------------------------------------------------------------------------------

subroutine initialize_bump(geom, locs, bump)

implicit none

!Arguments
type(fv3jedi_geom), intent(in)    :: geom
type(ioda_locs),    intent(in)    :: locs
type(bump_type),    intent(inout) :: bump

!Locals
integer :: mod_num
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

integer, save :: bumpcount = 0
character(len=5) :: cbumpcount
character(len=16) :: bump_nam_prefix

type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

! Each bump%nam%prefix must be distinct
! -------------------------------------
bumpcount = bumpcount + 1
write(cbumpcount,"(I0.5)") bumpcount
bump_nam_prefix = 'fv3jedi_bump_data_'//cbumpcount


!Get the Solution dimensions
!---------------------------
mod_num = (geom%iec - geom%isc + 1) * (geom%jec - geom%jsc + 1)


!Calculate interpolation weight using BUMP
!-----------------------------------------
allocate( mod_lat(mod_num), mod_lon(mod_num) )
mod_lat = reshape( rad2deg*geom%grid_lat(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] )  
mod_lon = reshape( rad2deg*geom%grid_lon(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] ) - 180.0_kind_real

!Important namelist options
call bump%nam%init
bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)

!Less important namelist options (should not be changed)
bump%nam%prefix = bump_nam_prefix   ! Prefix for files output
bump%nam%default_seed = .true.
bump%nam%new_obsop = .true.

!Initialize geometry
allocate(area(mod_num))
allocate(vunit(mod_num,1))
allocate(lmask(mod_num,1))
area = 1.0           ! Dummy area
vunit = 1.0          ! Dummy vertical unit
lmask = .true.       ! Mask

!Initialize BUMP
call bump%setup_online( f_comm%communicator(),mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                                nobs=locs%nlocs,lonobs=locs%lon(:)-180.0_kind_real,latobs=locs%lat(:) )

!Release memory
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate( mod_lat, mod_lon )

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, locs, vars, gom, jvar)
implicit none
character(len=*), intent(in) :: cop
type(ioda_locs), intent(in)     :: locs
type(ufo_vars), intent(in)      :: vars
type(ufo_geovals), intent(in)   :: gom
integer, intent(in)             :: jvar

character(len=255) :: cinfo

cinfo="fv3jedi_state:checks "//trim(cop)//" : "

!Check things are the sizes we expect
!------------------------------------
if (gom%nobs /= locs%nlocs ) then
   call abor1_ftn(trim(cinfo)//"geovals wrong size")
endif
if( gom%nvar .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//"nvar wrong size")
endif
if( .not. allocated(gom%geovals) )then
   call abor1_ftn(trim(cinfo)//"geovals not allocated")
endif
if( size(gom%geovals) .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//"geovals wrong size")
endif
if (.not.gom%linit) then
   call abor1_ftn(trim(cinfo)//"geovals not initialized")
endif
if (allocated(gom%geovals(jvar)%vals)) then  
   if( gom%geovals(jvar)%nobs .ne. locs%nlocs )then
      call abor1_ftn(trim(cinfo)//"nobs wrong size")
   endif
   if( size(gom%geovals(jvar)%vals, 2) .ne. locs%nlocs )then
      call abor1_ftn(trim(cinfo)//"vals wrong size 2")
   endif       
else
  call abor1_ftn(trim(cinfo)//"vals not allocated")
endif 

end subroutine interp_checks

! ------------------------------------------------------------------------------

end module fv3jedi_getvalues_mod
