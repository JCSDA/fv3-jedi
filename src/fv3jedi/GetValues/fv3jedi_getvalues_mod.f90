module fv3jedi_getvalues_mod

use atlas_module
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use type_bump, only: bump_type
use ufo_geovals_mod, only: ufo_geovals, ufo_geovals_write_netcdf
use ufo_locs_mod, only: ufo_locs
use oops_variables_mod, only: oops_variables
use unstructured_interpolation_mod, only: unstrc_interp

use fv3jedi_bump_mod,      only: bump_init, bump_apply, bump_apply_ad
use fv3jedi_constants_mod, only: rad2deg, constoz, grav
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_getvalues_traj_mod, only: fv3jedi_getvalues_traj
use fv3jedi_kinds_mod, only: kind_int, kind_real
use fv3jedi_state_utils_mod, only: fv3jedi_state
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_field_mod, only: get_field_array

use surface_vt_mod
use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use height_vt_mod
use wind_vt_mod, only: d2a, d2a_ad

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
type(oops_variables),                           intent(in)    :: vars
type(ufo_geovals),                              intent(inout) :: gom
type(fv3jedi_getvalues_traj), optional, target, intent(inout) :: traj

character(len=*), parameter :: myname = 'getvalues'

type(fckit_mpi_comm) :: f_comm

! Interpolation
logical, target  :: interp_alloc = .false.
logical, pointer :: pinterp_alloc => null()
type(bump_type), target  :: bump
type(bump_type), pointer :: pbump => null()
integer, target, save :: bumpid = 1000
integer, pointer      :: pbumpid => null()
type(unstrc_interp), target  :: unsinterp
type(unstrc_interp), pointer :: punsinterp => null()
real(kind=kind_real), allocatable :: lats_in(:), lons_in(:)

integer :: ii, jj, ji, jvar, jlev, jloc, ngrid, nlocs, nlocsg
real(kind=kind_real), allocatable :: mod_state(:), obs_state(:,:)
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

! Height variables
logical :: have_heights
real(kind=kind_real), allocatable :: height (:,:,:) !Height m
real(kind=kind_real), allocatable :: heighti(:,:,:) !Height m, interfaces

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
real(kind=kind_real), allocatable :: salinity(:)                 !Sea Surface Salinity            | surface(1)%salinity
logical,  parameter               :: use_compress = .true.       !Could be a fv3 namelist option?

! Pointers to needed fields
real(kind=kind_real), pointer :: slmsk   (:,:,:)
real(kind=kind_real), pointer :: frocean (:,:,:)
real(kind=kind_real), pointer :: frlake  (:,:,:)
real(kind=kind_real), pointer :: frseaice(:,:,:)
real(kind=kind_real), pointer :: tsea    (:,:,:)
real(kind=kind_real), pointer :: f10m    (:,:,:)
real(kind=kind_real), pointer :: ua      (:,:,:)
real(kind=kind_real), pointer :: va      (:,:,:)
real(kind=kind_real), pointer :: u_srf   (:,:,:)
real(kind=kind_real), pointer :: v_srf   (:,:,:)

real(kind=kind_real), allocatable :: state_f10m(:,:,:), state_slmsk(:,:,:)
real(kind=kind_real) :: wspd

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
f_comm = geom%f_comm

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
  punsinterp => traj%unsinterp

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

     pinterp_alloc => traj%lalloc
     pbumpid => traj%bumpid

  endif

else

  pbump => bump
  bumpid = bumpid + 1
  pbumpid => bumpid

  punsinterp => unsinterp

  interp_alloc = .false.
  pinterp_alloc => interp_alloc

endif

if (.not. pinterp_alloc) then
  if (trim(geom%interp_method) == 'bump') then
    call bump_init(geom, locs%nlocs, locs%lat, locs%lon, pbump, pbumpid)
  elseif (trim(geom%interp_method) == 'barycent') then
    allocate(lats_in(ngrid))
    allocate(lons_in(ngrid))
    jj = 0
    do j = jsc,jec
      do i = isc,iec
         jj = jj + 1
         lats_in(jj) = rad2deg*geom%grid_lat(i,j)
         lons_in(jj) = rad2deg*geom%grid_lon(i,j)
      enddo
    enddo
    call punsinterp%create( f_comm, 4, trim(geom%interp_method), &
                           ngrid, lats_in, lons_in, &
                           locs%nlocs, locs%lat, locs%lon )
  endif
  pinterp_alloc = .true.
endif

! Create Buffer for interpolated values
! --------------------------------------
if (trim(geom%interp_method) == 'barycent') allocate(mod_state(ngrid))
allocate(obs_state(nlocs,npz+1))


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

if (associated(state%rh)) then
  allocate(rh(isc:iec,jsc:jec,npz))
  rh = state%rh
  have_rh = .true.
elseif (have_t .and. have_pressures .and. associated(state%q)) then

  allocate(qsat(isc:iec,jsc:jec,npz))
  allocate(rh(isc:iec,jsc:jec,npz))

  call get_qsat(geom,delp,t,state%q,qsat)
  call q_to_rh(geom,qsat,state%q,rh)

  deallocate(qsat)
  have_rh = .true.
endif

! Compute heights
! ---------------
have_heights = .false.

if (have_pressures .and. have_t .and. associated(state%q) .and. associated(state%phis)) then

  have_heights = .true.

  allocate(height(isc:iec,jsc:jec,npz))
  allocate(heighti(isc:iec,jsc:jec,npz+1))

  call geop_height(geom,prs,prsi,t,state%q,&
                   state%phis(:,:,1),use_compress,height)

  call geop_height_levels(geom,prs,prsi,t,state%q,&
                          state%phis(:,:,1),use_compress,heighti)

endif

! Get CRTM surface variables
! ----------------------

! Model may not have slmsk
if (get_field_array(state%fields,'slmsk',slmsk)) then

  allocate(state_slmsk(isc:iec,jsc:jec,1))
  state_slmsk(isc:iec,jsc:jec,1) = slmsk(isc:iec,jsc:jec,1)

elseif ( get_field_array(state%fields,'frocean' ,frocean ) .and. &
         get_field_array(state%fields,'frlake'  ,frlake  ) .and. &
         get_field_array(state%fields,'frseaice',frseaice) .and. &
         get_field_array(state%fields,'tsea'    ,tsea    ) ) then

  allocate(state_slmsk(isc:iec,jsc:jec,1))

  state_slmsk = 1.0_kind_real !Land
  do j = jsc,jec
    do i = isc,iec
      if ( frocean(i,j,1) + frlake(i,j,1) >= 0.6_kind_real) then
        state_slmsk(i,j,1) = 0.0_kind_real ! Water
      endif
      if ( state_slmsk(i,j,1) == 0.0_kind_real .and. frseaice(i,j,1) > 0.5_kind_real) then
        state_slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
      if ( state_slmsk(i,j,1) == 0.0_kind_real .and. tsea(i,j,1) < 271.4_kind_real ) then
        state_slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
    enddo
  enddo

endif

! Model may not have f10m
if (get_field_array(state%fields,'f10m',f10m)) then

  allocate(state_f10m(isc:iec,jsc:jec,1))
  state_f10m(isc:iec,jsc:jec,1) = f10m(isc:iec,jsc:jec,1)

elseif ( get_field_array(state%fields,'ua'    , ua    ) .and. &
         get_field_array(state%fields,'va'    , va    ) .and. &
         get_field_array(state%fields,'u_srf' , u_srf ) .and. &
         get_field_array(state%fields,'v_srf' , v_srf ) ) then

  allocate(state_f10m(isc:iec,jsc:jec,1))

  state_f10m(isc:iec,jsc:jec,1) = sqrt(u_srf(isc:iec,jsc:jec,1)**2 + v_srf(isc:iec,jsc:jec,1)**2)

  do j = jsc,jec
    do i = isc,iec
      wspd = sqrt(ua(i,j,geom%npz)**2 +  va(i,j,geom%npz)**2)
      if (state_f10m(i,j,1) > 0.0_kind_real) then
        state_f10m(i,j,1) = state_f10m(i,j,1)/wspd
      else
        state_f10m(i,j,1) = 1.0_kind_real
      endif
    enddo
  enddo

endif

have_crtm_srf =.false.
if ( associated(state%sheleg) .and. associated(state%tsea  ) .and. &
     associated(state%vtype ) .and. associated(state%stype ) .and. &
     associated(state%vfrac ) .and. associated(state%stc   ) .and. &
     associated(state%smc   ) .and. associated(state%u_srf ) .and. &
     associated(state%v_srf ) .and. allocated(state_slmsk  ) .and. &
     allocated(state_f10m) ) then

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

  if ( associated(state%sss)) then
    allocate(salinity(nlocs))
    salinity = 0.0_kind_real
    call crtm_surface( geom, nlocs, ngrid, locs%lat(:), locs%lon(:), &
                       state_slmsk,  state%sheleg, &
                       state%tsea,   state%vtype, &
                       state%stype,  state%vfrac, &
                       state%stc(:,:,1), state%smc(:,:,1), &
                       state%u_srf, state%v_srf, state_f10m, &
                       land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                       snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                       snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                       wind_speed, wind_direction, state%sss, salinity)
  else
    call crtm_surface( geom, nlocs, ngrid, locs%lat(:), locs%lon(:), &
                       state_slmsk,  state%sheleg, &
                       state%tsea,   state%vtype, &
                       state%stype,  state%vfrac, &
                       state%stc(:,:,1), state%smc(:,:,1), &
                       state%u_srf, state%v_srf, state_f10m, &
                       land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                       snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                       snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                       wind_speed, wind_direction)
  endif

  have_crtm_srf = .true.

endif

! Get CRTM moisture variables
! ---------------------------

have_crtm_cld = .false.
if (allocated(state_slmsk) .and. have_t .and. have_pressures &
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
      if (state_slmsk(i,j,1) == 0) water_coverage_m(i,j) = 1.0_kind_real
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

if (allocated(state_slmsk)) deallocate(state_slmsk)
if (allocated(state_f10m)) deallocate(state_f10m)

! Variable transforms and interpolate to obs locations
! ----------------------------------------------------

do jvar = 1, vars%nvars()

  geovalm = 0.0_kind_real
  geovals = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  ! Convert to observation variables/units
  ! --------------------------------------
  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    if (.not. associated(state%ua)) &
      call variable_fail(trim(vars%variable(jvar)),"state%ua")

    nvl = npz
    do_interp = .true.
    geovalm = state%ua
    geoval => geovalm

  case ("northward_wind")

    if (.not. associated(state%va)) &
      call variable_fail(trim(vars%variable(jvar)),"state%va")

    nvl = npz
    do_interp = .true.
    geovalm = state%va
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not. have_t) &
      call variable_fail(trim(vars%variable(jvar)),"t")

    nvl = npz
    do_interp = .true.
    geovalm = t
    geoval => geovalm

  case ("specific_humidity")

    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%variable(jvar)),"state%q")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%q,0.0_kind_real)
    geoval => geovalm

  case ("virtual_temperature")

    if (.not. have_t) &
      call variable_fail(trim(vars%variable(jvar)),"t")
    if (.not. associated(state%q)) &
      call variable_fail(trim(vars%variable(jvar)),"state%q")

    nvl = npz
    do_interp = .true.
    call T_to_Tv(geom,t,state%q,geovalm)
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) &
      call variable_fail(trim(vars%variable(jvar)),"qmr")

    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm

  case ("relative_humidity")

    if (.not. have_rh) &
      call variable_fail(trim(vars%variable(jvar)),"rh")

    nvl = npz
    do_interp = .true.
    geovalm = max(rh,0.0_kind_real)
    geoval => geovalm

  case ("surface_pressure")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%variable(jvar)),"ps")

      nvl = 1
      do_interp = .true.

      if (associated(state%ps)) then
        geovals(:,:,1) = state%ps(:,:,1)
      elseif (associated(state%delp)) then
        geovals(:,:,1) = sum(state%delp,3)
      else
        call abor1_ftn(trim(myname)//" no way to get surface pressure ")
      endif

      geoval => geovals

  case ("air_pressure")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%variable(jvar)),"prs")

    nvl = npz
    do_interp = .true.
    geovalm = prs
    geoval => geovalm

  case ("air_pressure_levels")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%variable(jvar)),"prsi")

    nvl = npz + 1
    do_interp = .true.
    geovale = prsi
    geoval => geovale

  case ("air_pressure_thickness")

    if (.not. have_pressures) &
      call variable_fail(trim(vars%variable(jvar)),"delp")

    nvl = npz
    do_interp = .true.
    geovalm = delp
    geoval => geovalm


  case ("geopotential_height","height")

    if (.not. have_heights) &
      call variable_fail(trim(vars%variable(jvar)),"have_heights")

    nvl = npz
    do_interp = .true.
    geovalm = height
    geoval => geovalm

  case ("geopotential_height_levels")

    if (.not. have_heights) &
      call variable_fail(trim(vars%variable(jvar)),"have_heights")

    nvl = npz
    do_interp = .true.
    geovale = heighti
    geoval => geovale

  case ("surface_geopotential_height","surface_altitude")

    if (.not. associated(state%phis)) &
      call variable_fail(trim(vars%variable(jvar)),"state%phis")

    nvl = 1
    do_interp = .true.
    geovalm(:,:,1) = state%phis(:,:,1) / grav
    geoval => geovalm

  case ("mole_fraction_of_ozone_in_air")

    if (.not. associated(state%o3)) &
      call variable_fail(trim(vars%variable(jvar)),"state%o3")

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
      call variable_fail(trim(vars%variable(jvar)),"ql_ade")

   nvl = npz
   do_interp = .true.
   geovalm = ql_ade
   geoval => geovalm

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%variable(jvar)),"qi_ade")

   nvl = npz
   do_interp = .true.
   geovalm = qi_ade
   geoval => geovalm

  case ("effective_radius_of_cloud_liquid_water_particle")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%variable(jvar)),"ql_efr")

   nvl = npz
   do_interp = .true.
   geovalm = ql_efr
   geoval => geovalm

  case ("effective_radius_of_cloud_ice_particle")

    if (.not. have_crtm_cld) &
      call variable_fail(trim(vars%variable(jvar)),"qi_efr")

   nvl = npz
   do_interp = .true.
   geovalm = qi_efr
   geoval => geovalm

  case ("water_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"water_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_coverage

  case ("land_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"land_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_coverage

  case ("ice_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"ice_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_coverage

  case ("surface_snow_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"snow_coverage")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_coverage

  case ("surface_temperature_where_sea")

   if (.not. have_crtm_srf) &
     call variable_fail(trim(vars%variable(jvar)),"water_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = water_temperature


  case ("sea_surface_temperature")

    if (.not. associated(state%tsea)) &
      call variable_fail(trim(vars%variable(jvar)),"state%tsea")

   nvl = 1
   do_interp = .true.
   geovals = state%tsea
   geoval => geovals

   case ("sea_surface_salinity")

   if (.not. associated(state%sss)) &
      call variable_fail(trim(vars%variable(jvar)),"salinity")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = salinity



  case ("surface_temperature_where_land")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"land_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = land_temperature

  case ("surface_temperature_where_ice")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"ice_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = ice_temperature

  case ("surface_temperature_where_snow")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"snow_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_temperature

  case ("surface_snow_thickness")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"snow_depth")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = snow_depth

  case ("vegetation_area_fraction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"vegetation_fraction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = vegetation_fraction

  case ("surface_wind_speed")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"wind_speed")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_speed

  case ("surface_wind_from_direction")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"wind_direction")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = wind_direction

  case ("leaf_area_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"lai")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = lai

  case ("volume_fraction_of_condensed_water_in_soil")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"soil_moisture_content")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_moisture_content

  case ("soil_temperature")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"soil_temperature")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = soil_temperature

  case ("land_type_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"land_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(land_type,kind_real)

  case ("vegetation_type_index")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"vegetation_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(vegetation_type,kind_real)

  case ("soil_type")

    if (.not. have_crtm_srf) &
      call variable_fail(trim(vars%variable(jvar)),"soil_type")

   nvl = 1
   do_interp = .false.
   obs_state(:,1) = real(soil_type,kind_real)

  case ("sulf","so4","mass_fraction_of_sulfate_in_air")

   if (.not. associated(state%so4)) &
     call variable_fail(trim(vars%variable(jvar)),"state%so4")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%so4,0.0_kind_real)
   geoval => geovalm

  case ("bc1","bcphobic","mass_fraction_of_hydrophobic_black_carbon_in_air")

   if (.not. associated(state%bcphobic)) &
     call variable_fail(trim(vars%variable(jvar)),"state%bcphobic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%bcphobic,0.0_kind_real)
   geoval => geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

   if (.not. associated(state%bcphilic)) &
     call variable_fail(trim(vars%variable(jvar)),"state%bcphilic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%bcphilic,0.0_kind_real)
   geoval => geovalm

  case ("oc1","ocphobic","mass_fraction_of_hydrophobic_organic_carbon_in_air")

   if (.not. associated(state%ocphobic)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ocphobic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ocphobic,0.0_kind_real)
   geoval => geovalm

  case ("oc2","ocphilic","mass_fraction_of_hydrophilic_organic_carbon_in_air")

   if (.not. associated(state%ocphilic)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ocphilic")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ocphilic,0.0_kind_real)
   geoval => geovalm

  case ("dust1","du001","mass_fraction_of_dust001_in_air")

   if (.not. associated(state%du001)) &
     call variable_fail(trim(vars%variable(jvar)),"state%du001")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du001,0.0_kind_real)
   geoval => geovalm

  case ("dust2","du002","mass_fraction_of_dust002_in_air")

   if (.not. associated(state%du002)) &
     call variable_fail(trim(vars%variable(jvar)),"state%du002")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du002,0.0_kind_real)
   geoval => geovalm

  case ("dust3","du003","mass_fraction_of_dust003_in_air")

   if (.not. associated(state%du003)) &
     call variable_fail(trim(vars%variable(jvar)),"state%du003")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du003,0.0_kind_real)
   geoval => geovalm

  case ("dust4","du004","mass_fraction_of_dust004_in_air")

   if (.not. associated(state%du004)) &
     call variable_fail(trim(vars%variable(jvar)),"state%du004")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du004,0.0_kind_real)
   geoval => geovalm

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

   if (.not. associated(state%du005)) &
     call variable_fail(trim(vars%variable(jvar)),"state%du005")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%du005,0.0_kind_real)
   geoval => geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

   if (.not. associated(state%ss001)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ss001")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss001,0.0_kind_real)
   geoval => geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

   if (.not. associated(state%ss002)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ss002")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss002,0.0_kind_real)
   geoval => geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

   if (.not. associated(state%ss003)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ss003")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss003,0.0_kind_real)
   geoval => geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

   if (.not. associated(state%ss004)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ss004")

   nvl = npz
   do_interp = .true.
   geovalm = max(state%ss004,0.0_kind_real)
   geoval => geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

   if (.not. associated(state%ss005)) &
     call variable_fail(trim(vars%variable(jvar)),"state%ss005")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%ss005,0.0_kind_real)
    geoval => geovalm

  case ("no3an1","mass_fraction_of_nitrate001_in_air")

   if (.not. associated(state%no3an1)) &
     call variable_fail(trim(vars%variable(jvar)),"state%no3an1")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an1,0.0_kind_real)
    geoval => geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

   if (.not. associated(state%no3an2)) &
     call variable_fail(trim(vars%variable(jvar)),"state%no3an2")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an2,0.0_kind_real)
    geoval => geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

   if (.not. associated(state%no3an3)) &
     call variable_fail(trim(vars%variable(jvar)),"state%no3an3")

    nvl = npz
    do_interp = .true.
    geovalm = max(state%no3an3,0.0_kind_real)
    geoval => geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%variable(jvar)))

  end select

  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nvars())


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    if (trim(geom%interp_method) == 'bump') then
      call bump_apply(nvl, geom, geoval(:,:,1:nvl), locs%nlocs, obs_state(:,1:nvl), pbump)
      do jlev = 1, nvl
        do jloc = 1,locs%nlocs
          gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_state(jloc,jlev)
        enddo
      enddo
    elseif (trim(geom%interp_method) == 'barycent') then
      do jlev = 1, nvl
        ii = 0
        do jj = jsc, jec
          do ji = isc, iec
            ii = ii + 1
            mod_state(ii) = geoval(ji, jj, jlev)
          enddo
        enddo
        call punsinterp%apply(mod_state,obs_state(:,1))
        do jloc = 1,locs%nlocs
          gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_state(jloc,1)
        enddo
      enddo
    endif
  else
    do jloc = 1,locs%nlocs
      gom%geovals(jvar)%vals(nvl,locs%indx(jloc)) = obs_state(jloc,1)
    enddo
  endif

  nullify(geoval)

enddo


! For debugging we can write the geovals
! comm = geom%f_comm
! write(cproc,fmt='(i4.4)') comm%rank()
! gomfilename = 'Data/fv3jedi_geovals_'//trim(cproc)//'.nc4'
! call ufo_geovals_write_netcdf(gom, trim(gomfilename), locs%lat(:), locs%lon(:))


if (.not. present(traj)) then
  if (trim(geom%interp_method) == 'bump') then
    call pbump%dealloc
  elseif (trim(geom%interp_method) == 'barycent') then
    call punsinterp%delete
  endif
endif

nullify(pbump)
nullify(pinterp_alloc)
nullify(pbumpid)
nullify(punsinterp)

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
type(oops_variables),     intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_tl'

integer :: ii, jj, ji, jvar, jlev, jloc, i, j, k
real(kind=kind_real), allocatable :: mod_increment(:), obs_increment(:,:)

integer :: nvl
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:), geovals(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, npz, ngrid

logical :: have_winds
real(kind=kind_real), pointer :: u(:,:,:), v(:,:,:)
real(kind=kind_real), allocatable :: inc_ua(:,:,:), inc_va(:,:,:)

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
if (trim(geom%interp_method) == 'barycent') then
  ngrid = (iec-isc+1)*(jec-jsc+1)
  allocate(mod_increment(ngrid))
endif
allocate(obs_increment(locs%nlocs,npz+1))


! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))
allocate(geovals(isc:iec,jsc:jec,1))


! Wind transforms
! ---------------
have_winds = .false.
if (get_field_array(inc%fields,'ua',u) .and. get_field_array(inc%fields,'va',v)) then

  have_winds = .true.
  allocate(inc_ua(isc:iec,jsc:jec,npz))
  allocate(inc_va(isc:iec,jsc:jec,npz))
  inc_ua(isc:iec,jsc:jec,:) = u(isc:iec,jsc:jec,:)
  inc_va(isc:iec,jsc:jec,:) = v(isc:iec,jsc:jec,:)

elseif (get_field_array(inc%fields,'ud',u) .and. get_field_array(inc%fields,'vd',v)) then

  have_winds = .true.
  allocate(inc_ua(isc:iec,jsc:jec,npz))
  allocate(inc_va(isc:iec,jsc:jec,npz))
  call d2a(geom, u, v, inc_ua, inc_va)

endif


! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nvars()

  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    if (.not. have_winds) &
      call variable_fail(trim(vars%variable(jvar)),"eastward_wind (tl)")

    nvl = npz
    do_interp = .true.
    geovalm = inc_ua
    geoval => geovalm

  case ("northward_wind")

    if (.not. have_winds) &
      call variable_fail(trim(vars%variable(jvar)),"northward_wind (tl)")

    nvl = npz
    do_interp = .true.
    geovalm = inc_va
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not. associated(inc%t)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%t (tl)")

    nvl = npz
    do_interp = .true.
    geovalm = inc%t
    geoval => geovalm

  case ("specific_humidity")

    if (.not. associated(inc%q)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%q (tl)")

    nvl = npz
    do_interp = .true.
    geovalm = inc%q
    geoval => geovalm

  case ("virtual_temperature")

    if (.not.allocated(traj%t)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%t (tl)")
    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%q (tl)")
    if (.not. associated(inc%t)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%t (tl)")
    if (.not. associated(inc%q)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%q (tl)")

    nvl = inc%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, inc%t, traj%q, inc%q, geovalm )
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%q (tl)")
    if (.not.associated(inc%q)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%q (tl)")

    nvl = inc%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, inc%q, geovalm)
    geoval => geovalm

  case ("surface_pressure")

    nvl = 1
    do_interp = .true.

    if (associated(inc%ps)) then
      geovals = inc%ps
    elseif (associated(inc%delp)) then
      geovals(:,:,1) = sum(inc%delp,3)
    else
      call abor1_ftn(trim(myname)//" no way to get surface pressure ")
    endif

    geoval => geovals

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("mole_fraction_of_ozone_in_air")

    if (.not.allocated(traj%o3)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%o3")
    if (.not.associated(inc%o3)) &
      call variable_fail(trim(vars%variable(jvar)),"inc%o3 (tl)")

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

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%so4
   geoval => geovalm

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%bcphobic
   geoval => geovalm

  case ("bc2","bcphilic", "mass_fraction_of_hydrophilic_black_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%bcphilic
   geoval => geovalm

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ocphobic
   geoval => geovalm

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ocphilic
   geoval => geovalm

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du001
   geoval => geovalm

  case ("dust2","du002", "mass_fraction_of_dust002_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du002
   geoval => geovalm

  case ("dust3","du003", "mass_fraction_of_dust003_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du003
   geoval => geovalm

  case ("dust4","du004", "mass_fraction_of_dust004_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du004
   geoval => geovalm

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%du005
   geoval => geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss001
   geoval => geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss002
   geoval => geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss003
   geoval => geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss004
   geoval => geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

   nvl = npz
   do_interp = .true.
   geovalm = inc%ss005
   geoval => geovalm

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an1
    geoval => geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an2
    geoval => geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    nvl = npz
    do_interp = .true.
    geovalm = inc%no3an3
    geoval => geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%variable(jvar)))

  end select


  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nvars())


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, vars, gom, jvar)


  ! Find observation location equivlent
  ! -----------------------------------
  if (do_interp) then
    if (trim(geom%interp_method) == 'bump') then
      call bump_apply(nvl, geom, geoval(:,:,1:nvl), locs%nlocs, obs_increment(:,1:nvl), traj%bump)
      do jlev = 1, nvl
        do jloc = 1,locs%nlocs
          gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_increment(jloc,jlev)
        enddo
      enddo
    elseif (trim(geom%interp_method) == 'barycent') then
      do jlev = 1, nvl
        ii = 0
        do jj = jsc, jec
          do ji = isc, iec
            ii = ii + 1
            mod_increment(ii) = geoval(ji, jj, jlev)
          enddo
        enddo
        call traj%unsinterp%apply(mod_increment,obs_increment(:,1))
        do jloc = 1,locs%nlocs
          gom%geovals(jvar)%vals(jlev,locs%indx(jloc)) = obs_increment(jloc,1)
        enddo
      enddo
    endif
  else
    do jloc = 1,locs%nlocs
      gom%geovals(jvar)%vals(nvl,locs%indx(jloc)) = obs_increment(jloc,1)
    enddo
  endif

  nullify(geoval)

enddo

if (trim(geom%interp_method) == 'barycent') deallocate(mod_increment)
deallocate(obs_increment)
deallocate(geovalm,geovale,geovals)
if (allocated(inc_ua) .and. allocated(inc_va)) deallocate(inc_ua,inc_va)

end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(geom, inc, locs, vars, gom, traj)

implicit none
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_increment),  intent(inout) :: inc
type(ufo_locs),           intent(in)    :: locs
type(oops_variables),     intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_ad'

integer :: ii, jj, ji, jvar, jlev, jloc
real(kind=kind_real), allocatable :: mod_increment(:), obs_increment(:,:)

integer :: nvl, i, j, k
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:), geovals(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
logical :: do_interp

integer :: isc, iec, jsc, jec, npz, ngrid

logical :: have_winds
real(kind=kind_real), pointer :: u(:,:,:), v(:,:,:)
real(kind=kind_real), allocatable :: inc_ua(:,:,:), inc_va(:,:,:)
real(kind=kind_real), allocatable :: inc_u(:,:,:), inc_v(:,:,:)

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
if (trim(geom%interp_method) == 'barycent') then
  ngrid = (iec-isc+1)*(jec-jsc+1)
  allocate(mod_increment(ngrid))
endif
allocate(obs_increment(locs%nlocs,npz+1))

! Local GeoVals
! -------------
allocate(geovale(isc:iec,jsc:jec,npz+1))
allocate(geovalm(isc:iec,jsc:jec,npz))
allocate(geovals(isc:iec,jsc:jec,1))

geovale = 0.0_kind_real
geovalm = 0.0_kind_real
geovals = 0.0_kind_real

! Flag on whether wind observations are used
! ------------------------------------------
have_winds = .false.

! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nvars()

  ! PART 1, do_interp flag
  ! ----------------------
  do_interp = .false.

  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    have_winds = .true.
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("northward_wind")

    have_winds = .true.
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

  case ("surface_pressure")

   nvl = 1
   do_interp = .true.
   geoval => geovals

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

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust2","du002","mass_fraction_of_dust002_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust3","du003","mass_fraction_of_dust003_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust4","du004","mass_fraction_of_dust004_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

   nvl = npz
   do_interp = .true.
   geoval => geovalm

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    nvl = npz
    do_interp = .true.
    geoval => geovalm


  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    nvl = npz
    do_interp = .true.
    geoval => geovalm


  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%variable(jvar)))

  end select

  !Part 2, apply adjoint of interp
  !-------------------------------
  if (do_interp) then
    if (trim(geom%interp_method) == 'bump') then
      do jlev = 1, nvl
        do jloc = 1,locs%nlocs
          obs_increment(jloc,jlev) = gom%geovals(jvar)%vals(jlev,locs%indx(jloc))
        enddo
      enddo
      call bump_apply_ad(nvl, geom, geoval(:,:,1:nvl), locs%nlocs, obs_increment(:,1:nvl), traj%bump)
    elseif (trim(geom%interp_method) == 'barycent') then
      do jlev = nvl, 1, -1
        do jloc = 1,locs%nlocs
          obs_increment(jloc,1) = gom%geovals(jvar)%vals(jlev,locs%indx(jloc))
        enddo
        call traj%unsinterp%apply_ad(mod_increment,obs_increment(:,1))
        ii = 0
        do jj = jsc, jec
          do ji = isc, iec
            ii = ii + 1
            geoval(ji, jj, jlev) = mod_increment(ii)
          enddo
        enddo
      enddo
    endif
  else
    do jloc = 1,locs%nlocs
      obs_increment(jloc,1) = gom%geovals(jvar)%vals(nvl,locs%indx(jloc))
    enddo
  endif

  !Part 3, back to increment variables
  !-----------------------------------

  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    allocate(inc_ua(isc:iec,jsc:jec,npz))
    inc_ua = geovalm

  case ("northward_wind")

    allocate(inc_va(isc:iec,jsc:jec,npz))
    inc_va = geovalm

  case ("air_temperature","temperature")

    inc%t = inc%t + geovalm

  case ("mole_fraction_of_ozone_in_air")

    if (.not.allocated(traj%o3)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%o3")

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
      call variable_fail(trim(vars%variable(jvar)),"traj%t")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%q")

    call T_to_Tv_ad(geom, traj%t, inc%t, traj%q, inc%q, geovalm )

  case ("humidity_mixing_ratio")

    if (.not.allocated(traj%q)) &
      call variable_fail(trim(vars%variable(jvar)),"traj%q")

    call crtm_mixratio_ad(geom, traj%q, inc%q, geovalm)

  case ("surface_pressure")

    if (associated(inc%ps)) then
      inc%ps = inc%ps + geovals
    elseif (associated(inc%delp)) then
      do k = 1,geom%npz
        inc%delp(:,:,k) = inc%delp(:,:,k) + geovals(:,:,1)
      enddo
    else
      call abor1_ftn(trim(myname)//" no way to set surface pressure ")
    endif

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

   inc%so4 = inc%so4 + geovalm

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

   inc%bcphobic = inc%bcphobic + geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

   inc%bcphilic = inc%bcphilic + geovalm

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

   inc%ocphobic = inc%ocphobic + geovalm

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

   inc%ocphilic = inc%ocphilic + geovalm

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

   inc%du001 = inc%du001 + geovalm

  case ("dust2","du002", "mass_fraction_of_dust002_in_air")

   inc%du002 = inc%du002 + geovalm

  case ("dust3","du003", "mass_fraction_of_dust003_in_air")

   inc%du003 = inc%du003 + geovalm

  case ("dust4","du004", "mass_fraction_of_dust004_in_air")

   inc%du004 = inc%du004 + geovalm

  case ("dust5","du005", "mass_fraction_of_dust005_in_air")

   inc%du005 = inc%du005 + geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

   inc%ss001 = inc%ss001 + geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

   inc%ss002 = inc%ss002 + geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

   inc%ss003 = inc%ss003 + geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

   inc%ss004 = inc%ss004 + geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

    inc%ss005 = inc%ss005 + geovalm

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    inc%no3an1 = inc%no3an1 + geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    inc%no3an2 = inc%no3an2 + geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    inc%no3an3 = inc%no3an3 + geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%variable(jvar)))

  end select

  geovale = 0.0_kind_real
  geovalm = 0.0_kind_real

enddo

! Conversion of wind fields
! -------------------------
if (have_winds) then

  if (get_field_array(inc%fields,'ua',u) .and. get_field_array(inc%fields,'va',v)) then

    u = u + inc_ua
    v = v + inc_va

  elseif (get_field_array(inc%fields,'ud',u) .and. get_field_array(inc%fields,'vd',v)) then

    allocate(inc_u(isc:iec  ,jsc:jec+1,npz))
    allocate(inc_v(isc:iec+1,jsc:jec  ,npz))
    call d2a_ad(geom, inc_u, inc_v, inc_ua, inc_va)
    u = u + inc_u
    v = v + inc_u

  else

    call abor1_ftn("getvalues_ad: problem getting winds from increment")

  endif

endif

if (trim(geom%interp_method) == 'barycent') deallocate(mod_increment)
deallocate(obs_increment)
deallocate(geovalm,geovale,geovals)
if (allocated(inc_ua) .and. allocated(inc_va)) deallocate(inc_ua,inc_va)
if (allocated(inc_u) .and. allocated(inc_v)) deallocate(inc_u,inc_v)

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

subroutine getvalues_checks(cop, vars, gom, jvar)
implicit none
character(len=*),     intent(in) :: cop
type(oops_variables), intent(in) :: vars
type(ufo_geovals),    intent(in) :: gom
integer,              intent(in) :: jvar

character(len=255) :: cinfo

cinfo="fv3jedi_"//trim(cop)//" checks:"

!Check things are the sizes we expect
!------------------------------------
if( gom%nvar .ne. vars%nvars() )then
   call abor1_ftn(trim(cinfo)//" nvar wrong size")
endif
if( .not. allocated(gom%geovals) )then
   call abor1_ftn(trim(cinfo)//" geovals not allocated")
endif
if( size(gom%geovals) .ne. vars%nvars() )then
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
