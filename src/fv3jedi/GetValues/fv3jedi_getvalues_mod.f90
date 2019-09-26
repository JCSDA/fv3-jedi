module fv3jedi_getvalues_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use type_bump, only: bump_type
use ufo_geovals_mod, only: ufo_geovals, ufo_geovals_write_netcdf
use ufo_locs_mod, only: ufo_locs
use variables_mod, only: oops_vars

use fv3jedi_constants_mod, only: rad2deg, constoz, grav
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_getvalues_traj_mod, only: fv3jedi_getvalues_traj
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

type(fv3jedi_geom),                             intent(in)    :: geom
type(fv3jedi_state),                            intent(in)    :: state
type(ufo_locs),                                 intent(in)    :: locs
type(oops_vars),                                intent(in)    :: vars
type(ufo_geovals),                              intent(inout) :: gom
type(fv3jedi_getvalues_traj), optional, target, intent(inout) :: traj

character(len=*), parameter :: myname = 'getvalues'

type(fckit_mpi_comm) :: f_comm
type(bump_type), target  :: bump
type(bump_type), pointer :: pbump => null()
logical, target  :: bump_alloc = .false.
logical, pointer :: pbump_alloc => null()
integer, target, save :: bumpid = 1000
integer, pointer      :: pbumpid => null()

integer :: ii, jj, ji, jvar, jlev, jloc, ngrid, nlocs, nlocsg
real(kind=kind_real), allocatable :: mod_state(:,:)
real(kind=kind_real), allocatable :: obs_state(:,:)
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:), geovals(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
integer :: nvl
logical :: do_interp

integer :: isc,iec,jsc,jec,npz,i,j

!Local pressure variables
real(kind=kind_real), allocatable :: delp(:,:,:) !Pressure thickness Pa
real(kind=kind_real), allocatable :: prsi(:,:,:) !Pressure Pa, interfaces
real(kind=kind_real), allocatable :: prs (:,:,:) !Pressure Pa, midpoint
real(kind=kind_real), allocatable :: logp(:,:,:) !Log(pressue), (Pa) midpoint

real(kind=kind_real), allocatable :: t   (:,:,:) !Temperature
real(kind=kind_real), allocatable :: qsat(:,:,:) !Saturation specific humidity
real(kind=kind_real), allocatable :: rh  (:,:,:) !Relative humidity

!Local CRTM moisture variables
real(kind=kind_real), allocatable :: ql_ade(:,:,:) !Cloud liq water kgm^2
real(kind=kind_real), allocatable :: qi_ade(:,:,:) !Cloud ice water kgm^2
real(kind=kind_real), allocatable :: ql_efr(:,:,:) !Cloud effective radius microns
real(kind=kind_real), allocatable :: qi_efr(:,:,:) !Cloud effective radium microns
real(kind=kind_real), allocatable :: qmr(:,:,:)    !Moisture mixing ratio
real(kind=kind_real), allocatable :: water_coverage_m(:,:) !Water coverage, model grid

!Flags on variables
logical :: have_t, have_pressures, have_crtm_srf, have_crtm_cld, have_rh, have_qmr

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
logical,  parameter               :: use_compress = .true.       !Could be a fv3 namelist option?

!For writing GeoVaLs foe debugging
!type(fckit_mpi_comm) :: comm
!character(len=800)   :: gomfilename
!character(len=10)    :: cproc

! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

ngrid = (iec-isc+1)*(jec-jsc+1)
nlocs = locs%nlocs


!If no observations can early exit
!---------------------------------
f_comm = fckit_mpi_comm()
call f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
if (nlocsg == 0) then
  if (present(traj)) then
     traj%lalloc = .true.
     traj%noobs = .true.
  endif
  return
endif


! Variable transforms
! -------------------

! Temperature
have_t = .false.

if (associated(state%t)) then
  allocate(t(isc:iec,jsc:jec,npz))
  t = state%t
  have_t = .true.
elseif (associated(state%pt)) then
  if (.not. associated(state%pkz)) then
    call abor1_ftn("fv3jedi_getvalues_mod.getvalues: A state with potential temperature needs pressure to the kappa")
  endif
  allocate(t(isc:iec,jsc:jec,npz))
  call pt_to_t(geom,state%pkz,state%pt,t)
  have_t = .true.
endif

! Initialize the interpolation trajectory
! ---------------------------------------
if (present(traj)) then

  pbump => traj%bump

  if (.not. traj%lalloc) then

     traj%ngrid = ngrid

     if (have_t) then
       if (.not.allocated(traj%t)) allocate(traj%t(isc:iec,jsc:jec,1:npz))
       traj%t = t
     endif

     if (associated(state%q)) then
       if (.not.allocated(traj%q)) allocate(traj%q(isc:iec,jsc:jec,1:npz))
       traj%q = state%q
     endif

     if (associated(state%o3)) then
       if (.not.allocated(traj%o3)) allocate(traj%o3(isc:iec,jsc:jec,1:npz))
       traj%o3 = state%o3
     endif

     pbump_alloc => traj%lalloc
     pbumpid => traj%bumpid

  endif

else

  pbump => bump
  bump_alloc = .false.
  pbump_alloc => bump_alloc
  bumpid = bumpid + 1
  pbumpid => bumpid

endif

if (.not. pbump_alloc) then
   call initialize_bump(geom, locs, pbump, pbumpid)
   pbump_alloc = .true.
endif

! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_state(ngrid,1))
allocate(obs_state(nlocs,1))


! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))
allocate(geovals(isc:iec,jsc:jec,1))

! Get pressures at edge, center & log center
! ------------------------------------------
have_pressures = .false.

if (associated(state%delp)) then
  allocate(delp(isc:iec,jsc:jec,npz))
  delp = state%delp
  have_pressures = .true.
elseif (associated(state%ps)) then
  allocate(delp(isc:iec,jsc:jec,npz  ))
  do jlev = 1,geom%npz
    delp(:,:,jlev) = (geom%ak(jlev+1)-geom%ak(jlev))+(geom%bk(jlev+1)-geom%bk(jlev))*state%ps(:,:,1)
  enddo
  have_pressures = .true.
elseif (associated(state%pe)) then
  allocate(delp(isc:iec,jsc:jec,npz  ))
  do jlev = 1,geom%npz
    delp(:,:,jlev) = state%pe(:,:,jlev+1) - state%pe(:,:,jlev)
  enddo
  have_pressures = .true.
endif

if (have_pressures) then
  allocate(prsi(isc:iec,jsc:jec,npz+1))
  allocate(prs (isc:iec,jsc:jec,npz  ))
  allocate(logp(isc:iec,jsc:jec,npz  ))
  call delp_to_pe_p_logp(geom,delp,prsi,prs,logp)
endif

! Compute relative humidity
! -------------------------
have_rh = .false.
if (have_t .and. have_pressures .and. associated(state%q)) then

  allocate(qsat(isc:iec,jsc:jec,npz))
  allocate(rh(isc:iec,jsc:jec,npz))

  call get_qsat(geom,delp,t,state%q,qsat)
  call q_to_rh(geom,qsat,state%q,rh)

  deallocate(qsat)
  have_rh = .true.
endif

! Get CRTM surface variables
! ----------------------
have_crtm_srf =.false.
if (associated(state%slmsk ) .and. associated(state%slmsk ) .and. &
    associated(state%sheleg) .and. associated(state%tsea  ) .and. &
    associated(state%vtype ) .and. associated(state%stype ) .and. &
    associated(state%vfrac ) .and. associated(state%stc   ) .and. &
    associated(state%smc   ) .and. associated(state%snwdph) .and. &
    associated(state%u_srf ) .and. associated(state%v_srf ) .and. &
    associated(state%f10m  )) then

  allocate(wind_speed(nlocs))
  allocate(wind_direction(nlocs))
  allocate(land_type(nlocs))
  allocate(vegetation_type(nlocs))
  allocate(soil_type(nlocs))
  allocate(water_coverage(nlocs))
  allocate(land_coverage(nlocs))
  allocate(ice_coverage(nlocs))
  allocate(snow_coverage(nlocs))
  allocate(lai(nlocs))
  allocate(water_temperature(nlocs))
  allocate(land_temperature(nlocs))
  allocate(ice_temperature(nlocs))
  allocate(snow_temperature(nlocs))
  allocate(soil_moisture_content(nlocs))
  allocate(vegetation_fraction(nlocs))
  allocate(soil_temperature(nlocs))
  allocate(snow_depth(nlocs))

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

  call crtm_surface( geom, nlocs, ngrid, locs%lat(:), locs%lon(:), &
                     state%slmsk,  state%sheleg, &
                     state%tsea,   state%vtype, &
                     state%stype,  state%vfrac, &
                     state%stc,    state%smc, &
                     state%snwdph, state%u_srf, &
                     state%v_srf,  state%f10m, &
                     land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                     snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                     snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                     wind_speed, wind_direction )

  have_crtm_srf = .true.

endif


! Get CRTM moisture variables
! ---------------------------

have_crtm_cld = .false.
if (associated(state%slmsk) .and. have_t .and. have_pressures &
  .and. associated(state%q) .and. associated(state%qi) .and. associated(state%ql)) then

  allocate(ql_ade(isc:iec,jsc:jec,npz))
  allocate(qi_ade(isc:iec,jsc:jec,npz))
  allocate(ql_efr(isc:iec,jsc:jec,npz))
  allocate(qi_efr(isc:iec,jsc:jec,npz))
  allocate(water_coverage_m(isc:iec,jsc:jec))

  ql_ade = 0.0_kind_real
  qi_ade = 0.0_kind_real
  ql_efr = 0.0_kind_real
  qi_efr = 0.0_kind_real

  !TODO Is it water_coverage or sea_coverage fed in here?
  water_coverage_m = 0.0_kind_real
  do j = jsc,jec
    do i = isc,iec
      if (state%slmsk(i,j,1) == 0) water_coverage_m(i,j) = 1.0_kind_real
    enddo
  enddo

  call crtm_ade_efr( geom,prsi,t,delp, &
                     water_coverage_m,state%q, &
                     state%ql,state%qi, &
                     ql_ade,qi_ade,ql_efr,qi_efr )

  have_crtm_cld = .true.

endif

! CRTM mixing ratio
! -----------------
have_qmr = .false.
if (associated(state%q)) then
  allocate(qmr(isc:iec,jsc:jec,npz))
  call crtm_mixratio(geom,state%q,qmr)
  have_qmr = .true.
endif

! Variable transforms and interpolate to obs locations
! ----------------------------------------------------

do jvar = 1, vars%nv

  geovalm = 0.0_kind_real
  geovals = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  ! Convert to observation variables/units
  ! --------------------------------------
  select case (trim(vars%fldnames(jvar)))

  case ("eastward_wind")

    if (.not. associated(state%ua)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%ua")

    nvl = npz
    do_interp = .true.
    geovalm = state%ua
    geoval => geovalm

  case ("northward_wind")

    if (.not. associated(state%va)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%va")

    nvl = npz
    do_interp = .true.
    geovalm = state%va
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not. have_t) &
      call variable_fail(trim(vars%fldnames(jvar)),"t")

    nvl = npz
    do_interp = .true.
    geovalm = t
    geoval => geovalm

  case ("specific_humidity")

    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%q")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%q,0.0_kind_real)
    geoval => geovalm

  case ("virtual_temperature")

    if (.not. have_t) &
      call variable_fail(trim(vars%fldnames(jvar)),"t")
    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%q")

    nvl = npz
    do_interp = .true.
    call T_to_Tv(geom,t,state%q,geovalm)
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) &
      call variable_fail(trim(vars%fldnames(jvar)),"qmr")

    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm

  case ("relative_humidity")

    if (.not. have_rh) &
      call variable_fail(trim(vars%fldnames(jvar)),"rh")

    nvl = npz
    do_interp = .true.
    geovalm = max(rh,0.0_kind_real)
    geoval => geovalm

  case ("air_pressure")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%fldnames(jvar)),"prs")

    nvl = npz
    do_interp = .true.
    geovalm = prs
    geoval => geovalm

  case ("air_pressure_levels")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%fldnames(jvar)),"prsi")

    nvl = npz + 1
    do_interp = .true.
    geovale = prsi
    geoval => geovale

  case ("geopotential_height")

    if (.not. have_t) &
      call variable_fail(trim(vars%fldnames(jvar)),"t")
    if (.not. have_pressures) &
      call variable_fail(trim(vars%fldnames(jvar)),"prs,prsi")
    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%q")
    if (.not. associated(state%phis)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%phis")

    call geop_height(geom,prs,prsi,t,state%q,&
                     state%phis(:,:,1),use_compress,geovalm)
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("geopotential_height_levels")

    if (.not. have_t) &
      call variable_fail(trim(vars%fldnames(jvar)),"t")
    if (.not. have_pressures) &
      call variable_fail(trim(vars%fldnames(jvar)),"prs,prsi")
    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%q")
    if (.not. associated(state%phis)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%phis")

    call geop_height_levels(geom,prs,prsi,t,state%q,&
                            state%phis(:,:,1),use_compress,geovale)
    nvl = npz + 1
    do_interp = .true.
    geoval => geovale

  case ("surface_geopotential_height")

    if (.not. associated(state%phis)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%phis")

    nvl = 1
    do_interp = .true.
    geovalm(:,:,1) = state%phis(:,:,1) / grav
    geoval => geovalm

  case ("mole_fraction_of_ozone_in_air")

    if (.not. associated(state%o3)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%o3")

   nvl = npz
   do_interp = .true.
   geovalm = max(0.0_kind_real,state%o3) * constoz
   geoval => geovalm

  case ("mole_fraction_of_carbon_dioxide_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = 407.0_kind_real !Just a constant for now
   geoval => geovalm

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%fldnames(jvar)),"ql_ade")

   nvl = npz
   do_interp = .true.
   geovalm = ql_ade
   geoval => geovalm

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%fldnames(jvar)),"qi_ade")

   nvl = npz
   do_interp = .true.
   geovalm = qi_ade
   geoval => geovalm

  case ("effective_radius_of_cloud_liquid_water_particle")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%fldnames(jvar)),"ql_efr")

   nvl = npz
   do_interp = .true.
   geovalm = ql_efr
   geoval => geovalm

  case ("effective_radius_of_cloud_ice_particle")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%fldnames(jvar)),"qi_efr")

   nvl = npz
   do_interp = .true.
   geovalm = qi_efr
   geoval => geovalm

  case ("water_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"water_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_coverage

  case ("land_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"land_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_coverage

  case ("ice_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"ice_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_coverage

  case ("surface_snow_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"snow_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_coverage

  case ("surface_temperature_where_sea")

   if (.not. have_crtm_srf) &
     call variable_fail(trim(vars%fldnames(jvar)),"water_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_temperature


  case ("sea_surface_temperature")

    if (.not. associated(state%tsea)) &
      call variable_fail(trim(vars%fldnames(jvar)),"state%tsea")

   nvl = 1
   do_interp = .true.
   geovals = state%tsea
   geoval => geovals

  case ("surface_temperature_where_land")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"land_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_temperature

  case ("surface_temperature_where_ice")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"ice_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_temperature

  case ("surface_temperature_where_snow")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"snow_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_temperature

  case ("surface_snow_thickness")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"snow_depth")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_depth

  case ("vegetation_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"vegetation_fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = vegetation_fraction

  case ("surface_wind_speed")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"wind_speed")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_speed

  case ("surface_wind_from_direction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"wind_direction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_direction

  case ("leaf_area_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"lai")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = lai

  case ("volume_fraction_of_condensed_water_in_soil")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"soil_moisture_content")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_moisture_content

  case ("soil_temperature")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"soil_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_temperature

  case ("land_type_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"land_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(land_type,kind_real)

  case ("vegetation_type_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"vegetation_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(vegetation_type,kind_real)

  case ("soil_type")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%fldnames(jvar)),"soil_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(soil_type,kind_real)

  case ("sulf","so4")

   if (.not. associated(state%so4)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%so4")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%so4,0.0_kind_real)
   geoval => geovalm

  case ("bc1","bcphobic")

   if (.not. associated(state%bcphobic)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%bcphobic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%bcphobic,0.0_kind_real)
   geoval => geovalm

  case ("bc2","bcphilic")

   if (.not. associated(state%bcphilic)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%bcphilic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%bcphilic,0.0_kind_real)
   geoval => geovalm

  case ("oc1","ocphobic")

   if (.not. associated(state%ocphobic)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ocphobic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ocphobic,0.0_kind_real)
   geoval => geovalm

  case ("oc2","ocphilic")

   if (.not. associated(state%ocphilic)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ocphilic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ocphilic,0.0_kind_real)
   geoval => geovalm

  case ("dust1","du001")

   if (.not. associated(state%du001)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%du001")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du001,0.0_kind_real)
   geoval => geovalm

  case ("dust2","du002")

   if (.not. associated(state%du002)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%du002")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du002,0.0_kind_real)
   geoval => geovalm

  case ("dust3","du003")

   if (.not. associated(state%du003)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%du003")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du003,0.0_kind_real)
   geoval => geovalm

  case ("dust4","du004")

   if (.not. associated(state%du004)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%du004")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du004,0.0_kind_real)
   geoval => geovalm

  case ("dust5","du005")

   if (.not. associated(state%du005)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%du005")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du005,0.0_kind_real)
   geoval => geovalm

  case ("seas1","ss001")

   if (.not. associated(state%ss001)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ss001")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss001,0.0_kind_real)
   geoval => geovalm

  case ("seas2","ss002")

   if (.not. associated(state%ss002)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ss002")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss002,0.0_kind_real)
   geoval => geovalm

  case ("seas3","ss003")

   if (.not. associated(state%ss003)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ss003")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss003,0.0_kind_real)
   geoval => geovalm

  case ("seas4","ss004")

   if (.not. associated(state%ss004)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ss004")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss004,0.0_kind_real)
   geoval => geovalm

  case ("seas5","ss005")

   if (.not. associated(state%ss005)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%ss005")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%ss005,0.0_kind_real)
    geoval => geovalm

  case ("no3an1")

   if (.not. associated(state%no3an1)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%no3an1")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an1,0.0_kind_real)
    geoval => geovalm

  case ("no3an2")

   if (.not. associated(state%no3an2)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%no3an2")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an2,0.0_kind_real)
    geoval => geovalm

  case ("no3an3")

   if (.not. associated(state%no3an3)) &
     call variable_fail(trim(vars%fldnames(jvar)),"state%no3an3")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an3,0.0_kind_real)
    geoval => geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select

  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nv)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, vars, gom, jvar)


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
      do jloc = 1,locs%nlocs
        gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_state(jloc,1)
      enddo
    enddo
  else
    do jloc = 1,locs%nlocs
      gom%geovals(jvar)%vals(nvl,locs%indx(jloc)) = obs_state(jloc,1)
    enddo
  endif

  nullify(geoval)

enddo


!For debugging we can write the geovals
!comm = fckit_mpi_comm()
!write(cproc,fmt='(i4.4)') comm%rank()
!gomfilename = 'Data/fv3jedi_geovals_'//trim(cproc)//'.nc4'
!call ufo_geovals_write_netcdf(gom, trim(gomfilename), locs%lat(:), locs%lon(:))


if (.not. present(traj)) then
  call pbump%dealloc
endif

nullify(pbump)
nullify(pbump_alloc)
nullify(pbumpid)

! Deallocate local memory
! -----------------------
if (allocated(mod_state            )) deallocate(mod_state            )
if (allocated(obs_state            )) deallocate(obs_state            )
if (allocated(geovale              )) deallocate(geovale              )
if (allocated(geovalm              )) deallocate(geovalm              )
if (allocated(geovals              )) deallocate(geovals              )
if (allocated(delp                 )) deallocate(delp                 )
if (allocated(prsi                 )) deallocate(prsi                 )
if (allocated(prs                  )) deallocate(prs                  )
if (allocated(logp                 )) deallocate(logp                 )
if (allocated(t                    )) deallocate(t                    )
if (allocated(rh                   )) deallocate(rh                   )
if (allocated(wind_speed           )) deallocate(wind_speed           )
if (allocated(wind_direction       )) deallocate(wind_direction       )
if (allocated(land_type            )) deallocate(land_type            )
if (allocated(vegetation_type      )) deallocate(vegetation_type      )
if (allocated(soil_type            )) deallocate(soil_type            )
if (allocated(water_coverage       )) deallocate(water_coverage       )
if (allocated(land_coverage        )) deallocate(land_coverage        )
if (allocated(ice_coverage         )) deallocate(ice_coverage         )
if (allocated(snow_coverage        )) deallocate(snow_coverage        )
if (allocated(lai                  )) deallocate(lai                  )
if (allocated(water_temperature    )) deallocate(water_temperature    )
if (allocated(land_temperature     )) deallocate(land_temperature     )
if (allocated(ice_temperature      )) deallocate(ice_temperature      )
if (allocated(snow_temperature     )) deallocate(snow_temperature     )
if (allocated(soil_moisture_content)) deallocate(soil_moisture_content)
if (allocated(vegetation_fraction  )) deallocate(vegetation_fraction  )
if (allocated(soil_temperature     )) deallocate(soil_temperature     )
if (allocated(snow_depth           )) deallocate(snow_depth           )
if (allocated(ql_ade               )) deallocate(ql_ade               )
if (allocated(qi_ade               )) deallocate(qi_ade               )
if (allocated(ql_efr               )) deallocate(ql_efr               )
if (allocated(qi_efr               )) deallocate(qi_efr               )
if (allocated(qmr                  )) deallocate(qmr                  )
if (allocated(water_coverage_m     )) deallocate(water_coverage_m     )

end subroutine getvalues

! ------------------------------------------------------------------------------

subroutine getvalues_tl(geom, inc, locs, vars, gom, traj)

implicit none
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_increment),  intent(in)    :: inc
type(ufo_locs),           intent(in)    :: locs
type(oops_vars),          intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_tl'

integer :: ii, jj, ji, jvar, jlev, jloc, i, j, k
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


!If no observations can early exit
!---------------------------------
if (traj%noobs)  return


! Grid convenience
! ----------------
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(locs%nlocs,1))


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

  case ("eastward_wind")

    nvl = npz
    do_interp = .true.
    geovalm = inc%ua
    geoval => geovalm

  case ("northward_wind")

    nvl = npz
    do_interp = .true.
    geovalm = inc%va
    geoval => geovalm

  case ("air_temperature","temperature")

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

    if (.not.allocated(traj%t)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%t")
    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%q")

    nvl = inc%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, inc%t, traj%q, inc%q, geovalm )
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%q")

    nvl = inc%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, inc%q, geovalm)
    geoval => geovalm

  case ("mole_fraction_of_ozone_in_air")

    if (.not.allocated(traj%o3)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%o3")

    nvl = npz
    do_interp = .true.
	  do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
	        if (traj%o3(i,j,k) .ge. 0.0_kind_real) then
	          geovalm = inc%o3
	        else
	          geovalm = 0.0_kind_real
	        endif
        enddo
      enddo
    enddo

    geoval => geovalm

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4")

   nvl = npz
   do_interp = .true.
   geovalm = inc%so4
   geoval => geovalm

  case ("bc1","bcphobic")

   nvl = npz
   do_interp = .true.
   geovalm = inc%bcphobic
   geoval => geovalm

  case ("bc2","bcphilic")

   nvl = npz
   do_interp = .true.
   geovalm = inc%bcphilic
   geoval => geovalm

  case ("oc1","ocphobic")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ocphobic
   geoval => geovalm

  case ("oc2","ocphilic")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ocphilic
   geoval => geovalm

  case ("dust1","du001")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du001
   geoval => geovalm

  case ("dust2","du002")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du002
   geoval => geovalm

  case ("dust3","du003")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du003
   geoval => geovalm

  case ("dust4","du004")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du004
   geoval => geovalm

  case ("dust5","du005")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du005
   geoval => geovalm

  case ("seas1","ss001")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss001
   geoval => geovalm

  case ("seas2","ss002")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss002
   geoval => geovalm

  case ("seas3","ss003")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss003
   geoval => geovalm

  case ("seas4","ss004")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss004
   geoval => geovalm

  case ("seas5","ss005")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss005
   geoval => geovalm

  case ("no3an1")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an1
    geoval => geovalm

  case ("no3an2")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an2
    geoval => geovalm

  case ("no3an3")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an3
    geoval => geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select


  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nv)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, vars, gom, jvar)


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
      do jloc = 1,locs%nlocs
        gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_increment(jloc,1)
      enddo
    enddo
  else
    do jloc = 1,locs%nlocs
      gom%geovals(jvar)%vals(nvl,locs%indx(jloc)) = obs_increment(jloc,1)
    enddo
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
type(fv3jedi_increment),  intent(inout) :: inc
type(ufo_locs),           intent(in)    :: locs
type(oops_vars),          intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_ad'

integer :: ii, jj, ji, jvar, jlev, jloc
real(kind=kind_real), allocatable :: mod_increment(:,:)
real(kind=kind_real), allocatable :: obs_increment(:,:)

integer :: nvl, i, j, k
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, npz


! Check traj is implemented
! -------------------------
if (.not.traj%lalloc) &
call abor1_ftn(trim(myname)//" trajectory for this obs op not found")


!If no observations can early exit
!---------------------------------
if (traj%noobs)  return


! Grid convenience
! ----------------
isc = inc%isc
iec = inc%iec
jsc = inc%jsc
jec = inc%jec
npz = inc%npz


! Create Buffer for interpolated values
! --------------------------------------
allocate(mod_increment(traj%ngrid,1))
allocate(obs_increment(locs%nlocs,1))


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

  case ("eastward_wind")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("northward_wind")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_temperature","temperature")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("specific_humidity")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("mole_fraction_of_ozone_in_air")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("humidity_mixing_ratio")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("bc1","bcphobic")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("bc2","bcphilic")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("oc1","ocphobic")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("oc2","ocphilic")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust1","du001")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust2","du002")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust3","du003")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust4","du004")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust5","du005")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas1","ss001")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas2","ss002")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas3","ss003")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas4","ss004")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas5","ss005")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("no3an1","no3an2","no3an3")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select

  !Part 2, apply adjoint of interp
  !-------------------------------
  if (do_interp) then
    !Perform level-by-level interpolation using BUMP
    do jlev = nvl, 1, -1
      do jloc = 1,locs%nlocs
        obs_increment(jloc,1) = gom%geovals(jvar)%vals(jlev,locs%indx(jloc))
      enddo
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
    do jloc = 1,locs%nlocs
      obs_increment(jloc,1) = gom%geovals(jvar)%vals(nvl,locs%indx(jloc))
    enddo
  endif

  !Part 3, back to increment variables
  !-----------------------------------

  select case (trim(vars%fldnames(jvar)))

  case ("eastward_wind")

    inc%ua = inc%ua + geovalm

  case ("northward_wind")

    inc%va = inc%va + geovalm

  case ("air_temperature","temperature")

    inc%t = inc%t + geovalm

  case ("mole_fraction_of_ozone_in_air")

    if (.not.allocated(traj%o3)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%o3")

  	do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
	        if (traj%o3(i,j,k) .ge. 0.0_kind_real) then
	          inc%o3 = inc%o3 + geovalm
	        else
	          inc%o3 = 0.0_kind_real
	        endif
        enddo
      enddo
    enddo

  case ("specific_humidity")

    inc%q = inc%q + geovalm

  case ("virtual_temperature")

    if (.not.allocated(traj%t)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%t")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%q")

    call T_to_Tv_ad(geom, traj%t, inc%t, traj%q, inc%q, geovalm )

  case ("humidity_mixing_ratio")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%fldnames(jvar)),"traj%q")

    call crtm_mixratio_ad(geom, traj%q, inc%q, geovalm)

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4")

   inc%so4 = inc%so4 + geovalm

  case ("bc1","bcphobic")

   inc%bcphobic = inc%bcphobic + geovalm

  case ("bc2","bcphilic")

   inc%bcphilic = inc%bcphilic + geovalm

  case ("oc1","ocphobic")

   inc%ocphobic = inc%ocphobic + geovalm

  case ("oc2","ocphilic")

   inc%ocphilic = inc%ocphilic + geovalm

  case ("dust1","du001")

   inc%du001 = inc%du001 + geovalm

  case ("dust2","du002")

   inc%du002 = inc%du002 + geovalm

  case ("dust3","du003")

   inc%du003 = inc%du003 + geovalm

  case ("dust4","du004")

   inc%du004 = inc%du004 + geovalm

  case ("dust5","du005")

   inc%du005 = inc%du005 + geovalm

  case ("seas1","ss001")

   inc%ss001 = inc%ss001 + geovalm

  case ("seas2","ss002")

   inc%ss002 = inc%ss002 + geovalm

  case ("seas3","ss003")

   inc%ss003 = inc%ss003 + geovalm

  case ("seas4","ss004")

   inc%ss004 = inc%ss004 + geovalm

  case ("seas5","ss005")

    inc%ss005 = inc%ss005 + geovalm

  case ("no3an1")

    inc%no3an1 = inc%no3an1 + geovalm

  case ("no3an2")

    inc%no3an2 = inc%no3an2 + geovalm

  case ("no3an3")

    inc%no3an3 = inc%no3an3 + geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select

  geovale = 0.0_kind_real
  geovalm = 0.0_kind_real

enddo

deallocate(mod_increment)
deallocate(obs_increment)

end subroutine getvalues_ad

! ------------------------------------------------------------------------------

subroutine allocate_geovals_vals(gom,jvar,gvlev,lastvar)

implicit none
type(ufo_geovals), intent(inout) :: gom     !GeoVaL
integer,           intent(in)    :: jvar    !Current variable
integer,           intent(in)    :: gvlev   !Number of model levels
logical,           intent(in)    :: lastvar !Logical true if on last var to go into gom

if (.not.allocated(gom%geovals(jvar)%vals)) then

  ! Set number of levels, nobs already set from total locs
  gom%geovals(jvar)%nval = gvlev

  ! Allocate %vals
  allocate(gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nlocs))
  gom%geovals(jvar)%vals = 0.0_kind_real

  ! Set flag for internal data arrays having been set
  if (lastvar) gom%linit  = .true.

endif

end subroutine allocate_geovals_vals

! ------------------------------------------------------------------------------

subroutine initialize_bump(geom, locs, bump, bumpid)

implicit none

!Arguments
type(fv3jedi_geom), intent(in)    :: geom
type(ufo_locs),     intent(in)    :: locs
type(bump_type),    intent(inout) :: bump
integer,            intent(in)    :: bumpid

!Locals
integer :: mod_num
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:)
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

character(len=5)    :: cbumpcount
character(len=1024) :: bump_nam_prefix

type(fckit_mpi_comm) :: f_comm


! Communicator from OOPS
! ----------------------
f_comm = fckit_mpi_comm()


! Each bump%nam%prefix must be distinct
! -------------------------------------
write(cbumpcount,"(I0.5)") bumpid
bump_nam_prefix = 'fv3jedi_bumpobsop_data_'//cbumpcount

!Get the Solution dimensions
!---------------------------
mod_num = (geom%iec - geom%isc + 1) * (geom%jec - geom%jsc + 1)


!Calculate interpolation weight using BUMP
!-----------------------------------------
allocate(mod_lat(mod_num))
allocate(mod_lon(mod_num))
mod_lat = reshape( rad2deg*geom%grid_lat(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] )
mod_lon = reshape( rad2deg*geom%grid_lon(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] ) - 180.0_kind_real


! Namelist options
! ----------------

!Important namelist options
call bump%nam%init

!Less important namelist options (should not be changed)
bump%nam%prefix = trim(bump_nam_prefix)   ! Prefix for files output
bump%nam%default_seed = .true.
bump%nam%new_obsop = .true.

bump%nam%write_obsop = .false.
bump%nam%verbosity = "none"

! Initialize geometry
! -------------------
allocate(area(mod_num))
allocate(vunit(mod_num,1))
allocate(lmask(mod_num,1))
area = 1.0           ! Dummy area
vunit = 1.0          ! Dummy vertical unit
lmask = .true.       ! Mask

! Initialize BUMP
! ---------------
call bump%setup_online( mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                        nobs=locs%nlocs,lonobs=locs%lon(:)-180.0_kind_real,latobs=locs%lat(:))

!Run BUMP drivers
call bump%run_drivers

!Partial deallocate option
call bump%partial_dealloc

! Release memory
! --------------
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate(mod_lat)
deallocate(mod_lon)

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine getvalues_checks(cop, vars, gom, jvar)
implicit none
character(len=*),   intent(in) :: cop
type(oops_vars),    intent(in) :: vars
type(ufo_geovals),  intent(in) :: gom
integer,            intent(in) :: jvar

character(len=255) :: cinfo

cinfo="fv3jedi_"//trim(cop)//" checks:"

!Check things are the sizes we expect
!------------------------------------
if( gom%nvar .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//" nvar wrong size")
endif
if( .not. allocated(gom%geovals) )then
   call abor1_ftn(trim(cinfo)//" geovals not allocated")
endif
if( size(gom%geovals) .ne. vars%nv )then
   call abor1_ftn(trim(cinfo)//" size geovals does not match number of vars from UFo")
endif
if (.not.gom%linit) then
!   call abor1_ftn(trim(cinfo)//" geovals initialization flag not set")
endif
if (.not. allocated(gom%geovals(jvar)%vals)) then
   call abor1_ftn(trim(cinfo)//"vals not allocated")
endif

end subroutine getvalues_checks

! ------------------------------------------------------------------------------

subroutine variable_fail(ufo_var,fv3_var)

implicit none
character(len=*), intent(in) :: ufo_var
character(len=*), intent(in) :: fv3_var

call abor1_ftn("GetValues.variable_fail: ufo variable "//trim(ufo_var)//&
               " needs fv3-jedi variable "//trim(fv3_var)//", but it is not available.")

end subroutine variable_fail

! ------------------------------------------------------------------------------

end module fv3jedi_getvalues_mod
