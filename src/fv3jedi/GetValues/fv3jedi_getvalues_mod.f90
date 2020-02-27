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
use fv3jedi_field_mod, only: has_field, copy_field_array, allocate_copy_field_array, &
                             pointer_field_array, pointer_field

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

type(fv3jedi_geom),                             intent(inout) :: geom
type(fv3jedi_state),                            intent(in)    :: state
type(ufo_locs),                                 intent(in)    :: locs
type(oops_variables),                           intent(in)    :: vars
type(ufo_geovals),                              intent(inout) :: gom
type(fv3jedi_getvalues_traj), optional, target, intent(inout) :: traj

character(len=*), parameter :: myname = 'getvalues'

! Interpolation objects
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

integer :: isc,iec,jsc,jec,npz,i,j,k


! Transform flags
! ---------------
logical :: have_winds, have_t, have_pressures, have_crtm_srf, have_crtm_cld, have_q, have_rh, have_qmr
logical :: have_tv, have_salinity, have_geoph, have_o3, have_ps
logical,  parameter :: use_compress = .true.       !Could be a fv3 namelist option?

! States
! ------
real(kind=kind_real), allocatable :: ua    (:,:,:)         ! A-grid winds, u component
real(kind=kind_real), allocatable :: va    (:,:,:)         ! A-grid winds, v component
real(kind=kind_real), allocatable :: t     (:,:,:)         !Temperature
real(kind=kind_real), allocatable :: tv    (:,:,:)         !Virtual temperature
real(kind=kind_real), allocatable :: qsat  (:,:,:)         !Saturation specific humidity
real(kind=kind_real), allocatable :: rh    (:,:,:)         !Relative humidity
real(kind=kind_real), allocatable :: ps    (:,:,:)         !Pressure thickness Pa
real(kind=kind_real), allocatable :: delp  (:,:,:)         !Pressure thickness Pa
real(kind=kind_real), allocatable :: o3    (:,:,:)         !Ozone
real(kind=kind_real), allocatable :: prsi  (:,:,:)         !Pressure Pa, interfaces
real(kind=kind_real), allocatable :: prs   (:,:,:)         !Pressure Pa, midpoint
real(kind=kind_real), allocatable :: geophi(:,:,:)         !Pressure Pa, interfaces
real(kind=kind_real), allocatable :: geoph (:,:,:)         !Pressure Pa, midpoint
real(kind=kind_real), allocatable :: logp  (:,:,:)         !Log(pressue), (Pa) midpoint
real(kind=kind_real), allocatable :: ql_ade(:,:,:)       !Cloud liq water kgm^2
real(kind=kind_real), allocatable :: qi_ade(:,:,:)       !Cloud ice water kgm^2
real(kind=kind_real), allocatable :: ql_efr(:,:,:)       !Cloud effective radius microns
real(kind=kind_real), allocatable :: qi_efr(:,:,:)       !Cloud effective radium microns
real(kind=kind_real), allocatable :: qmr   (:,:,:)          !Moisture mixing ratio
real(kind=kind_real), allocatable :: water_coverage_m(:,:) !Water coverage, model grid

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


! State pointers
! --------------
real(kind=kind_real), pointer :: ud      (:,:,:)
real(kind=kind_real), pointer :: vd      (:,:,:)
real(kind=kind_real), pointer :: pt      (:,:,:)
real(kind=kind_real), pointer :: pkz     (:,:,:)
real(kind=kind_real), pointer :: phis    (:,:,:)
real(kind=kind_real), pointer :: q       (:,:,:)
real(kind=kind_real), pointer :: qi      (:,:,:)
real(kind=kind_real), pointer :: ql      (:,:,:)
real(kind=kind_real), pointer :: slmsk   (:,:,:)
real(kind=kind_real), pointer :: frocean (:,:,:)
real(kind=kind_real), pointer :: frlake  (:,:,:)
real(kind=kind_real), pointer :: frseaice(:,:,:)
real(kind=kind_real), pointer :: tsea    (:,:,:)
real(kind=kind_real), pointer :: f10m    (:,:,:)
real(kind=kind_real), pointer :: uap     (:,:,:)
real(kind=kind_real), pointer :: vap     (:,:,:)
real(kind=kind_real), pointer :: u_srf   (:,:,:)
real(kind=kind_real), pointer :: v_srf   (:,:,:)
real(kind=kind_real), pointer :: sheleg  (:,:,:)
real(kind=kind_real), pointer :: vtype   (:,:,:)
real(kind=kind_real), pointer :: stype   (:,:,:)
real(kind=kind_real), pointer :: vfrac   (:,:,:)
real(kind=kind_real), pointer :: stc     (:,:,:)
real(kind=kind_real), pointer :: smc     (:,:,:)
real(kind=kind_real), pointer :: sss     (:,:,:)


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


! ------------------------------- !
! Compute all variable transforms !
! ------------------------------- !

! Temperature
! -----------
have_t = .false.

if (has_field(state%fields,'t')) then
  call allocate_copy_field_array(state%fields,'t',t)
  have_t = .true.
elseif (has_field(state%fields,'pt')) then
  if (.not. has_field(state%fields, 'pkz')) &
    call abor1_ftn("fv3jedi_getvalues_mod.getvalues: A state with pt needs pkz")
  allocate(t(isc:iec,jsc:jec,npz))
  call pointer_field_array(state%fields, 'pt',  pt)
  call pointer_field_array(state%fields, 'pkz', pkz)
  call pt_to_t(geom,pkz,pt,t)
  have_t = .true.
endif

! Specific humidity
! -----------------
have_q = .false.
if (has_field(state%fields,'q')) then
  call pointer_field_array(state%fields, 'q',  q)
  have_q = .true.
endif

! Virtual temperature
! -------------------
have_tv = .false.
if (have_t .and. have_q) then
  allocate(tv(isc:iec,jsc:jec,npz))
  call t_to_tv(geom, t, q, tv)
  have_tv = .true.
endif

! Ozone
! -----
have_o3 = .false.
if (has_field(state%fields,'o3')) then
  call allocate_copy_field_array(state%fields, 'o3', o3)
  have_o3 = .true.
  do k = 1, npz
    do j = jsc, jec
      do i = isc, iec
        if (o3(i,j,k) >= 0.0_kind_real) then
          o3(i,j,k) = o3(i,j,k) * constoz
        else
          o3(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif

! Wind transforms
! ---------------
have_winds = .false.
if (state%have_agrid) then

  call allocate_copy_field_array(state%fields,'ua',ua)
  call allocate_copy_field_array(state%fields,'va',va)

  have_winds = .true.

elseif (state%have_dgrid) then

  call pointer_field_array(state%fields, 'ud', ud)
  call pointer_field_array(state%fields, 'vd', vd)

  allocate(ua(isc:iec,jsc:jec,npz))
  allocate(va(isc:iec,jsc:jec,npz))
  call d2a(geom, ud, vd, ua, va)

  have_winds = .true.

endif

! Get pressures at edge, center & log center
! ------------------------------------------
have_pressures = .false.

if (has_field(state%fields,'delp')) then
  call allocate_copy_field_array(state%fields, 'delp', delp)
  allocate(ps(isc:iec,jsc:jec,1))
  ps(:,:,1) = sum(delp,3)
  have_pressures = .true.
elseif (has_field(state%fields,'ps')) then
  call allocate_copy_field_array(state%fields, 'ps', ps)
  allocate(delp(isc:iec,jsc:jec,npz))
  do jlev = 1,geom%npz
    delp(:,:,jlev) = (geom%ak(jlev+1)-geom%ak(jlev))+(geom%bk(jlev+1)-geom%bk(jlev))*ps(:,:,1)
  enddo
  have_pressures = .true.
elseif (has_field(state%fields,'pe')) then
  call allocate_copy_field_array(state%fields, 'pe', prsi)
  allocate(ps(isc:iec,jsc:jec,1))
  ps(:,:,1) = prsi(:,:,npz+1)
  allocate(delp(isc:iec,jsc:jec,npz))
  do jlev = 1,geom%npz
    delp(:,:,jlev) = prsi(:,:,jlev+1) - prsi(:,:,jlev)
  enddo
  have_pressures = .true.
endif

if (have_pressures) then
  if (.not.allocated(prsi)) allocate(prsi(isc:iec,jsc:jec,npz+1))
  if (.not.allocated(prs )) allocate(prs (isc:iec,jsc:jec,npz  ))
  if (.not.allocated(logp)) allocate(logp(isc:iec,jsc:jec,npz  ))
  call delp_to_pe_p_logp(geom, delp, prsi, prs, logp)
endif

! Geopotential height
! -------------------
have_geoph = .false.
if (have_t .and. have_pressures .and. have_q .and. has_field(state%fields,'phis')) then
  call pointer_field_array(state%fields, 'phis',  phis)
  if (.not.allocated(geophi)) allocate(geophi(isc:iec,jsc:jec,npz+1))
  if (.not.allocated(geoph )) allocate(geoph (isc:iec,jsc:jec,npz  ))
  call geop_height(geom, prs, prsi, t, q, phis(:,:,1), use_compress, geoph)
  call geop_height_levels(geom, prs, prsi, t, q, phis(:,:,1), use_compress, geophi)
  have_geoph = .true.
endif

! Relative humidity
! -----------------
have_rh = .false.
if (has_field(state%fields,'rh')) then
  call allocate_copy_field_array(state%fields, 'rh', rh)
  have_rh = .true.
elseif (have_t .and. have_pressures .and. have_q) then

  allocate(qsat(isc:iec,jsc:jec,npz))
  allocate(rh  (isc:iec,jsc:jec,npz))

  call get_qsat(geom,delp,t,q,qsat)
  call q_to_rh(geom,qsat,q,rh)

  deallocate(qsat)
  have_rh = .true.
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

     if (have_q) then
       if (.not.allocated(traj%q)) allocate(traj%q(isc:iec,jsc:jec,1:npz))
       traj%q = q
     endif

     if (has_field(state%fields,'o3')) then
       call allocate_copy_field_array(state%fields, 'o3', traj%o3)
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


! Get CRTM surface variables
! ----------------------

! Model may not have slmsk
if (has_field(state%fields,'slmsk')) then

  call allocate_copy_field_array(state%fields, 'slmsk', state_slmsk)

elseif ( has_field(state%fields,'frocean' ) .and. &
         has_field(state%fields,'frlake'  ) .and. &
         has_field(state%fields,'frseaice') .and. &
         has_field(state%fields,'tsea'    ) ) then

  call pointer_field_array(state%fields, 'frocean' , frocean )
  call pointer_field_array(state%fields, 'frlake'  , frlake  )
  call pointer_field_array(state%fields, 'frseaice', frseaice)
  call pointer_field_array(state%fields, 'tsea'    , tsea    )

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
if (has_field(state%fields,'f10m')) then

  call allocate_copy_field_array(state%fields, 'f10m', state_f10m)

elseif ( has_field(state%fields, 'ua'   ) .and. &
         has_field(state%fields, 'va'   ) .and. &
         has_field(state%fields, 'u_srf') .and. &
         has_field(state%fields, 'v_srf') ) then

  call pointer_field_array(state%fields, 'ua'    , uap  )
  call pointer_field_array(state%fields, 'va'    , vap  )
  call pointer_field_array(state%fields, 'u_srf' , u_srf)
  call pointer_field_array(state%fields, 'v_srf' , v_srf)

  allocate(state_f10m(isc:iec,jsc:jec,1))

  state_f10m(isc:iec,jsc:jec,1) = sqrt(u_srf(isc:iec,jsc:jec,1)**2 + v_srf(isc:iec,jsc:jec,1)**2)

  do j = jsc,jec
    do i = isc,iec
      wspd = sqrt(uap(i,j,geom%npz)**2 +  vap(i,j,geom%npz)**2)
      if (state_f10m(i,j,1) > 0.0_kind_real) then
        state_f10m(i,j,1) = state_f10m(i,j,1)/wspd
      else
        state_f10m(i,j,1) = 1.0_kind_real
      endif
    enddo
  enddo

endif

have_crtm_srf =.false.
if ( has_field(state%fields, 'sheleg') .and. has_field(state%fields, 'tsea'  ) .and. &
     has_field(state%fields, 'vtype' ) .and. has_field(state%fields, 'stype' ) .and. &
     has_field(state%fields, 'vfrac' ) .and. has_field(state%fields, 'stc'   ) .and. &
     has_field(state%fields, 'smc'   ) .and. has_field(state%fields, 'u_srf' ) .and. &
     has_field(state%fields, 'v_srf' ) .and. allocated(state_slmsk  ) .and. &
     allocated(state_f10m) ) then

  call pointer_field_array(state%fields, 'sheleg', sheleg )
  call pointer_field_array(state%fields, 'tsea'  , tsea   )
  call pointer_field_array(state%fields, 'vtype' , vtype  )
  call pointer_field_array(state%fields, 'stype' , stype  )
  call pointer_field_array(state%fields, 'vfrac' , vfrac  )
  call pointer_field_array(state%fields, 'stc'   , stc    )
  call pointer_field_array(state%fields, 'smc'   , smc    )
  call pointer_field_array(state%fields, 'u_srf' , u_srf  )
  call pointer_field_array(state%fields, 'v_srf' , v_srf  )

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

  if (has_field(state%fields, 'sss' )) then

    call pointer_field_array(state%fields, 'sss', sss )

    allocate(salinity(nlocs))
    salinity = 0.0_kind_real

    call crtm_surface( geom, nlocs, ngrid, locs%lat(:), locs%lon(:), &
                       state_slmsk,  sheleg, &
                       tsea,   vtype, &
                       stype,  vfrac, &
                       stc(:,:,1), smc(:,:,1), &
                       u_srf, v_srf, state_f10m, &
                       land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                       snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                       snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                       wind_speed, wind_direction, sss, salinity)

    have_salinity = .true.

  else

    call crtm_surface( geom, nlocs, ngrid, locs%lat(:), locs%lon(:), &
                       state_slmsk,  sheleg, &
                       tsea,   vtype, &
                       stype,  vfrac, &
                       stc(:,:,1), smc(:,:,1), &
                       u_srf, v_srf, state_f10m, &
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
  .and. have_q .and. has_field(state%fields, 'qi') .and. has_field(state%fields, 'ql')) then

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

  call pointer_field_array(state%fields, 'ql', ql)
  call pointer_field_array(state%fields, 'qi', qi)

  call crtm_ade_efr( geom, prsi, t, delp, water_coverage_m, q, ql, qi, &
                     ql_ade, qi_ade, ql_efr, qi_efr )

  have_crtm_cld = .true.

endif

! CRTM mixing ratio
! -----------------
have_qmr = .false.
if (have_q) then
  allocate(qmr(isc:iec,jsc:jec,npz))
  call crtm_mixratio(geom,q,qmr)
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

    if (.not.have_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = ua
    geoval => geovalm

  case ("northward_wind")

    if (.not.have_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = va
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not.has_field(state%fields, 't')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = t
    geoval => geovalm

  case ("specific_humidity")

    if (.not.has_field(state%fields, 'q')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = max(q, 0.0_kind_real)
    geoval => geovalm

  case ("virtual_temperature")

    if (.not. have_tv) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = tv
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm

  case ("relative_humidity")

    if (.not. have_rh) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = max(rh, 0.0_kind_real)
    geoval => geovalm

  case ("surface_pressure")

    if (.not. have_pressures) call variable_fail(vars%variable(jvar))

      nvl = 1
      do_interp = .true.
      geovals = ps
      geoval => geovals

  case ("air_pressure")

    if (.not. have_pressures) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = prs
    geoval => geovalm

  case ("air_pressure_levels")

    if (.not. have_pressures) call variable_fail(vars%variable(jvar))

    nvl = npz + 1
    do_interp = .true.
    geovale = prsi
    geoval => geovale

  case ("air_pressure_thickness")

    if (.not. have_pressures) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = delp
    geoval => geovalm

  case ("geopotential_height","height")

    if (.not. have_geoph) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = geoph
    geoval => geovalm

  case ("geopotential_height_levels")

    if (.not. have_geoph) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovale = geophi
    geoval => geovale

  case ("surface_geopotential_height","surface_altitude")

    if (.not.has_field(state%fields, 'phis')) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .true.
    call pointer_field_array(state%fields, 'phis', geoval)
    geoval = geoval / grav

  case ("mole_fraction_of_ozone_in_air")

    if (.not.have_o3) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = o3
    geoval => geovalm

  case ("mole_fraction_of_carbon_dioxide_in_air")

    nvl = npz
    do_interp = .true.
    geovalm = 407.0_kind_real !Just a constant for now
    geoval => geovalm

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

    if (.not. have_crtm_cld) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = ql_ade
    geoval => geovalm

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

    if (.not. have_crtm_cld) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = qi_ade
    geoval => geovalm

  case ("effective_radius_of_cloud_liquid_water_particle")

    if (.not. have_crtm_cld) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = ql_efr
    geoval => geovalm

  case ("effective_radius_of_cloud_ice_particle")

    if (.not. have_crtm_cld) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = qi_efr
    geoval => geovalm

  case ("water_area_fraction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = water_coverage

  case ("land_area_fraction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = land_coverage

  case ("ice_area_fraction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = ice_coverage

  case ("surface_snow_area_fraction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = snow_coverage

  case ("surface_temperature_where_sea")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = water_temperature

  case ("sea_surface_temperature")

    if (.not.has_field(state%fields, 'tsea')) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .true.
    call pointer_field_array(state%fields, 'tsea', geoval)

  case ("sea_surface_salinity")

    if (.not. have_salinity) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = salinity

  case ("surface_temperature_where_land")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = land_temperature

  case ("surface_temperature_where_ice")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = ice_temperature

  case ("surface_temperature_where_snow")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = snow_temperature

  case ("surface_snow_thickness")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = snow_depth

  case ("vegetation_area_fraction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = vegetation_fraction

  case ("surface_wind_speed")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = wind_speed

  case ("surface_wind_from_direction")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = wind_direction

  case ("leaf_area_index")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = lai

  case ("volume_fraction_of_condensed_water_in_soil")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = soil_moisture_content

  case ("soil_temperature")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = soil_temperature

  case ("land_type_index")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = real(land_type,kind_real)

  case ("vegetation_type_index")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = real(vegetation_type,kind_real)

  case ("soil_type")

    if (.not. have_crtm_srf) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .false.
    obs_state(:,1) = real(soil_type,kind_real)

  case ("sulf","so4","mass_fraction_of_sulfate_in_air")

   if (.not.has_field(state%fields,'so4')) call variable_fail(vars%variable(jvar))

   nvl = npz
   do_interp = .true.
   call copy_field_array(state%fields, 'so4', geovalm)
   geovalm = max(geovalm, 0.0_kind_real)
   geoval => geovalm

  case ("bc1","bcphobic","mass_fraction_of_hydrophobic_black_carbon_in_air")

    if (.not.has_field(state%fields,'bcphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'bcphobic', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

    if (.not.has_field(state%fields,'bcphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'bcphilic', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("oc1","ocphobic","mass_fraction_of_hydrophobic_organic_carbon_in_air")

    if (.not.has_field(state%fields,'ocphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ocphobic', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("oc2","ocphilic","mass_fraction_of_hydrophilic_organic_carbon_in_air")

    if (.not.has_field(state%fields,'ocphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ocphilic', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("dust1","du001","mass_fraction_of_dust001_in_air")

    if (.not.has_field(state%fields,'du001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'du001', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("dust2","du002","mass_fraction_of_dust002_in_air")

    if (.not.has_field(state%fields,'du002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'du002', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("dust3","du003","mass_fraction_of_dust003_in_air")

    if (.not.has_field(state%fields,'du003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'du003', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("dust4","du004","mass_fraction_of_dust004_in_air")

    if (.not.has_field(state%fields,'du004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'du004', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

    if (.not.has_field(state%fields,'du005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'du005', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

    if (.not.has_field(state%fields,'ss001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ss001', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

    if (.not.has_field(state%fields,'ss002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ss002', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

    if (.not.has_field(state%fields,'ss003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ss003', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

    if (.not.has_field(state%fields,'ss004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ss004', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

    if (.not.has_field(state%fields,'ss005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'ss005', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("no3an1","mass_fraction_of_nitrate001_in_air")

    if (.not.has_field(state%fields,'no3an1')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'no3an1', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    if (.not.has_field(state%fields,'no3an2')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'no3an2', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
    geoval => geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    if (.not.has_field(state%fields,'no3an3')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call copy_field_array(state%fields, 'no3an3', geovalm)
    geovalm = max(geovalm, 0.0_kind_real)
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

logical :: have_winds, have_tv, have_qmr, have_o3, have_ps

real(kind=kind_real), pointer :: ud  (:,:,:)
real(kind=kind_real), pointer :: vd  (:,:,:)
real(kind=kind_real), pointer :: t   (:,:,:)
real(kind=kind_real), pointer :: q   (:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)

real(kind=kind_real), allocatable :: ua (:,:,:)
real(kind=kind_real), allocatable :: va (:,:,:)
real(kind=kind_real), allocatable :: tv (:,:,:)
real(kind=kind_real), allocatable :: qmr(:,:,:)
real(kind=kind_real), allocatable :: o3 (:,:,:)
real(kind=kind_real), allocatable :: ps (:,:,:)

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
if (inc%have_agrid) then

  call allocate_copy_field_array(inc%fields,'ua',ua)
  call allocate_copy_field_array(inc%fields,'va',va)

  have_winds = .true.

elseif (inc%have_dgrid) then

  call pointer_field_array(inc%fields, 'ud', ud)
  call pointer_field_array(inc%fields, 'vd', vd)

  allocate(ua(isc:iec,jsc:jec,npz))
  allocate(va(isc:iec,jsc:jec,npz))
  call d2a(geom, ud, vd, ua, va)

  have_winds = .true.

endif

! Virtual temperature
! -------------------
have_tv = .false.
if (allocated(traj%t) .and. allocated(traj%t) .and. &
    has_field(inc%fields,'t') .and. has_field(inc%fields,'q')) then
  call pointer_field_array(inc%fields, 't', t)
  call pointer_field_array(inc%fields, 'q', q)
  allocate(tv(isc:iec,jsc:jec,npz))
  call T_to_Tv_tl(geom, traj%t, t, traj%q, q, tv )
  have_tv = .true.
endif

! Humidity mixing ratio
! ---------------------
have_qmr = .false.
if (allocated(traj%q) .and. has_field(inc%fields,'q')) then
  call pointer_field_array(inc%fields, 'q', q)
  allocate(qmr(isc:iec,jsc:jec,npz))
  call crtm_mixratio_tl(geom, traj%q, q, qmr)
  have_qmr = .true.
endif

! Ozone
! -----
have_o3 = .false.
if (allocated(traj%o3) .and. has_field(inc%fields, 'o3')) then
  call allocate_copy_field_array(inc%fields,'o3',o3)
  have_o3 = .true.
  do k = 1, npz
    do j = jsc, jec
      do i = isc, iec
        if (traj%o3(i,j,k) >= 0.0_kind_real) then
          o3(i,j,k) = o3(i,j,k)  * constoz
        else
          o3(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif

! Surface pressure
! ----------------
have_ps = .false.
if (has_field(inc%fields,'ps')) then
  call allocate_copy_field_array(inc%fields,'ps',ps)
  have_ps = .true.
elseif (has_field(inc%fields,'delp')) then
  call pointer_field_array(inc%fields,'delp',delp)
  allocate(ps(isc:iec,jsc:jec,1))
  ps(:,:,1) = sum(delp,3)
  have_ps = .true.
endif


! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nvars()

  geovalm = 0.0_kind_real
  geovale = 0.0_kind_real

  do_interp = .false.

  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    if (.not. have_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = ua
    geoval => geovalm

  case ("northward_wind")

    if (.not. have_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = va
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not. has_field(inc%fields,'t')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 't', geoval)

  case ("specific_humidity")

    if (.not. has_field(inc%fields,'q')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'q', geoval)

  case ("virtual_temperature")

    if (.not. have_tv) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = tv
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = qmr
    geoval => geovalm

  case ("surface_pressure")

    if (.not. have_ps) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .true.
    geovals = ps
    geoval => geovals

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("mole_fraction_of_ozone_in_air")

    if (.not. have_o3) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geovalm = o3
    geoval => geovalm

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

    if (.not.has_field(inc%fields,'so4')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'so4', geoval)

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

    if (.not.has_field(inc%fields,'bcphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'bcphobic', geoval)

  case ("bc2","bcphilic", "mass_fraction_of_hydrophilic_black_carbon_in_air")

    if (.not.has_field(inc%fields,'bcphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'bcphilic', geoval)

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

    if (.not.has_field(inc%fields,'ocphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ocphobic', geoval)

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

    if (.not.has_field(inc%fields,'ocphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ocphilic', geoval)

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

    if (.not.has_field(inc%fields,'du001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'du001', geoval)

  case ("dust2","du002", "mass_fraction_of_dust002_in_air")

    if (.not.has_field(inc%fields,'du002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'du002', geoval)

  case ("dust3","du003", "mass_fraction_of_dust003_in_air")

    if (.not.has_field(inc%fields,'du003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'du003', geoval)

  case ("dust4","du004", "mass_fraction_of_dust004_in_air")

    if (.not.has_field(inc%fields,'du004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'du004', geoval)

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

    if (.not.has_field(inc%fields,'du005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'du005', geoval)

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

    if (.not.has_field(inc%fields,'ss001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ss001', geoval)

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

    if (.not.has_field(inc%fields,'ss002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ss002', geoval)

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

    if (.not.has_field(inc%fields,'ss003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ss003', geoval)

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

    if (.not.has_field(inc%fields,'ss004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ss004', geoval)

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

    if (.not.has_field(inc%fields,'ss005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'ss005', geoval)

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    if (.not.has_field(inc%fields,'no3an1')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'no3an1', geoval)

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    if (.not.has_field(inc%fields,'no3an2')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'no3an2', geoval)

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    if (.not.has_field(inc%fields,'no3an3')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    call pointer_field_array(inc%fields, 'no3an3', geoval)

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

if (allocated(mod_increment)) deallocate(mod_increment)
deallocate(obs_increment)
deallocate(geovalm,geovale,geovals)

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

logical :: have_inc_winds
logical :: have_inc_tv
logical :: have_inc_qmr
logical :: have_inc_o3
logical :: have_inc_ps
logical :: have_gom_winds
logical :: have_gom_tv
logical :: have_gom_qmr
logical :: have_gom_o3
logical :: have_gom_ps

real(kind=kind_real), pointer :: uap (:,:,:)
real(kind=kind_real), pointer :: vap (:,:,:)
real(kind=kind_real), pointer :: ud  (:,:,:)
real(kind=kind_real), pointer :: vd  (:,:,:)
real(kind=kind_real), pointer :: t   (:,:,:)
real(kind=kind_real), pointer :: q   (:,:,:)
real(kind=kind_real), pointer :: o3p (:,:,:)
real(kind=kind_real), pointer :: psp (:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)

real(kind=kind_real), pointer :: so4     (:,:,:)
real(kind=kind_real), pointer :: bcphobic(:,:,:)
real(kind=kind_real), pointer :: bcphilic(:,:,:)
real(kind=kind_real), pointer :: ocphobic(:,:,:)
real(kind=kind_real), pointer :: ocphilic(:,:,:)
real(kind=kind_real), pointer :: du001   (:,:,:)
real(kind=kind_real), pointer :: du002   (:,:,:)
real(kind=kind_real), pointer :: du003   (:,:,:)
real(kind=kind_real), pointer :: du004   (:,:,:)
real(kind=kind_real), pointer :: du005   (:,:,:)
real(kind=kind_real), pointer :: ss001   (:,:,:)
real(kind=kind_real), pointer :: ss002   (:,:,:)
real(kind=kind_real), pointer :: ss003   (:,:,:)
real(kind=kind_real), pointer :: ss004   (:,:,:)
real(kind=kind_real), pointer :: ss005   (:,:,:)
real(kind=kind_real), pointer :: no3an1  (:,:,:)
real(kind=kind_real), pointer :: no3an2  (:,:,:)
real(kind=kind_real), pointer :: no3an3  (:,:,:)

real(kind=kind_real), allocatable :: ua (:,:,:)
real(kind=kind_real), allocatable :: va (:,:,:)
real(kind=kind_real), allocatable :: tv (:,:,:)
real(kind=kind_real), allocatable :: o3 (:,:,:)
real(kind=kind_real), allocatable :: qmr(:,:,:)
real(kind=kind_real), allocatable :: ps (:,:,:)


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


! Transform flags to false
! ------------------------
have_inc_winds = .false.
have_inc_tv    = .false.
have_inc_qmr   = .false.
have_inc_o3    = .false.
have_gom_winds = .false.
have_gom_tv    = .false.
have_gom_qmr   = .false.
have_gom_o3    = .false.
have_gom_ps    = .false.


! Wind transforms
! ---------------
have_inc_winds = .false.
if (inc%have_agrid .or. inc%have_dgrid) then
  have_inc_winds = .true.
endif

! Virtual temperature
! -------------------
have_inc_tv = .false.
if (allocated(traj%t) .and. allocated(traj%q) .and. &
    has_field(inc%fields,'t') .and. has_field(inc%fields,'q')) then
  have_inc_tv = .true.
endif

! Humidity mixing ratio
! ---------------------
have_inc_qmr = .false.
if (allocated(traj%q) .and. has_field(inc%fields,'q')) then
  have_inc_qmr = .true.
endif

! Ozone
! -----
have_inc_o3 = .false.
if (allocated(traj%o3) .and. has_field(inc%fields,'o3')) then
  have_inc_o3 = .true.
endif

! Surface pressure
! ----------------
have_inc_ps = .false.
if (has_field(inc%fields,'ps')) then
  have_inc_ps = .true.
elseif (has_field(inc%fields,'delp')) then
  have_inc_ps = .true.
endif


! Interpolate increment to obs locations using pre-calculated weights
! ----------------------------------------------------------------
do jvar = 1, vars%nvars()

  ! PART 1, do_interp flag
  ! ----------------------
  do_interp = .false.

  select case (trim(vars%variable(jvar)))

  case ("eastward_wind")

    if (.not. have_inc_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("northward_wind")

    if (.not. have_inc_winds) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_temperature","temperature")

    if (.not. has_field(inc%fields,'t')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("specific_humidity")

    if (.not. has_field(inc%fields,'q')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("surface_pressure")

    if (.not. have_inc_ps) call variable_fail(vars%variable(jvar))

    nvl = 1
    do_interp = .true.
    geoval => geovals

  case ("mole_fraction_of_ozone_in_air")

    if (.not. have_inc_o3) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("virtual_temperature")

    if (.not. have_inc_tv) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("humidity_mixing_ratio")

    if (.not. have_inc_qmr) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

    if (.not.has_field(inc%fields,'so4')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

    if (.not.has_field(inc%fields,'bcphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

    if (.not.has_field(inc%fields,'bcphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

    if (.not.has_field(inc%fields,'ocphobic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

    if (.not.has_field(inc%fields,'ocphilic')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

    if (.not.has_field(inc%fields,'du001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("dust2","du002","mass_fraction_of_dust002_in_air")

    if (.not.has_field(inc%fields,'du002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("dust3","du003","mass_fraction_of_dust003_in_air")

    if (.not.has_field(inc%fields,'du003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("dust4","du004","mass_fraction_of_dust004_in_air")

    if (.not.has_field(inc%fields,'du004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("dust5","du005","mass_fraction_of_dust005_in_air")

    if (.not.has_field(inc%fields,'du005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

    if (.not.has_field(inc%fields,'ss001')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

    if (.not.has_field(inc%fields,'ss002')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

    if (.not.has_field(inc%fields,'ss003')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

    if (.not.has_field(inc%fields,'ss004')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

    if (.not.has_field(inc%fields,'ss005')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    if (.not.has_field(inc%fields,'no3an1')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    if (.not.has_field(inc%fields,'no3an2')) call variable_fail(vars%variable(jvar))

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    if (.not.has_field(inc%fields,'no3an3')) call variable_fail(vars%variable(jvar))

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

    have_gom_winds = .true.
    allocate(ua(isc:iec,jsc:jec,npz))
    ua = geovalm

  case ("northward_wind")

    have_gom_winds = .true.
    allocate(va(isc:iec,jsc:jec,npz))
    va = geovalm

  case ("air_temperature","temperature")

    call pointer_field_array(inc%fields, 't', t)
    t = t + geovalm

  case ("mole_fraction_of_ozone_in_air")

    have_gom_o3 = .true.
    allocate(o3(isc:iec,jsc:jec,npz))
    o3 = geovalm

  case ("specific_humidity")

    call pointer_field_array(inc%fields, 'q', q)
    q = q + geovalm

  case ("virtual_temperature")

    have_gom_tv = .true.
    allocate(tv(isc:iec,jsc:jec,npz))
    tv = geovalm

  case ("humidity_mixing_ratio")

    have_gom_qmr = .true.
    allocate(qmr(isc:iec,jsc:jec,npz))
    qmr = geovalm

  case ("surface_pressure")

    have_gom_ps = .true.
    allocate(ps(isc:iec,jsc:jec,npz))
    ps = geovals

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

  case ("air_pressure")

  case ("air_pressure_levels")

  case ("sulf","so4", "mass_fraction_of_sulfate_in_air")

    call pointer_field_array(inc%fields, 'so4', so4)
    so4 = so4 + geovalm

  case ("bc1","bcphobic", "mass_fraction_of_hydrophobic_black_carbon_in_air")

    call pointer_field_array(inc%fields, 'bcphobic', bcphobic)
    bcphobic = bcphobic + geovalm

  case ("bc2","bcphilic","mass_fraction_of_hydrophilic_black_carbon_in_air")

    call pointer_field_array(inc%fields, 'bcphilic', bcphilic)
    bcphilic = bcphilic + geovalm

  case ("oc1","ocphobic", "mass_fraction_of_hydrophobic_organic_carbon_in_air")

    call pointer_field_array(inc%fields, 'ocphobic', ocphobic)
    ocphobic = ocphobic + geovalm

  case ("oc2","ocphilic", "mass_fraction_of_hydrophilic_organic_carbon_in_air")

    call pointer_field_array(inc%fields, 'ocphilic', ocphilic)
    ocphilic = ocphilic + geovalm

  case ("dust1","du001", "mass_fraction_of_dust001_in_air")

    call pointer_field_array(inc%fields, 'du001', du001)
    du001 = du001 + geovalm

  case ("dust2","du002", "mass_fraction_of_dust002_in_air")

    call pointer_field_array(inc%fields, 'du002', du002)
    du002 = du002 + geovalm

  case ("dust3","du003", "mass_fraction_of_dust003_in_air")

    call pointer_field_array(inc%fields, 'du003', du003)
    du003 = du003 + geovalm

  case ("dust4","du004", "mass_fraction_of_dust004_in_air")

    call pointer_field_array(inc%fields, 'du004', du004)
    du004 = du004 + geovalm

  case ("dust5","du005", "mass_fraction_of_dust005_in_air")

    call pointer_field_array(inc%fields, 'du005', du005)
    du005 = du005 + geovalm

  case ("seas1","ss001","mass_fraction_of_sea_salt001_in_air")

    call pointer_field_array(inc%fields, 'ss001', ss001)
    ss001 = ss001 + geovalm

  case ("seas2","ss002","mass_fraction_of_sea_salt002_in_air")

    call pointer_field_array(inc%fields, 'ss002', ss002)
    ss002 = ss002 + geovalm

  case ("seas3","ss003","mass_fraction_of_sea_salt003_in_air")

    call pointer_field_array(inc%fields, 'ss003', ss003)
    ss003 = ss003 + geovalm

  case ("seas4","ss004","mass_fraction_of_sea_salt004_in_air")

    call pointer_field_array(inc%fields, 'ss004', ss004)
    ss004 = ss004 + geovalm

  case ("seas5","ss005","mass_fraction_of_sea_salt005_in_air")

    call pointer_field_array(inc%fields, 'ss005', ss005)
    ss005 = ss005 + geovalm

  case ("no3an1", "mass_fraction_of_nitrate001_in_air")

    call pointer_field_array(inc%fields, 'no3an1', no3an1)
    no3an1 = no3an1 + geovalm

  case ("no3an2", "mass_fraction_of_nitrate002_in_air")

    call pointer_field_array(inc%fields, 'no3an2', no3an2)
    no3an2 = no3an2 + geovalm

  case ("no3an3", "mass_fraction_of_nitrate003_in_air")

    call pointer_field_array(inc%fields, 'no3an3', no3an3)
    no3an3 = no3an3 + geovalm

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%variable(jvar)))

  end select

  geovale = 0.0_kind_real
  geovalm = 0.0_kind_real

enddo

! Conversion of wind fields
! -------------------------
if (have_gom_winds) then

  if (inc%have_agrid) then

    call pointer_field_array(inc%fields, 'ua', uap)
    call pointer_field_array(inc%fields, 'va', vap)

    uap = uap + ua
    vap = vap + va

  elseif (inc%have_dgrid) then

    call pointer_field_array(inc%fields, 'ud', ud)
    call pointer_field_array(inc%fields, 'vd', vd)

    call d2a_ad(geom, ud, vd, ua, va)

  endif
endif

! Virtual temperature
! -------------------
if (have_gom_tv) then
  call pointer_field_array(inc%fields, 't', t)
  call pointer_field_array(inc%fields, 'q', q)
  call T_to_Tv_ad(geom, traj%t, t, traj%q, q, tv )
endif

! Moisture mixing ratio
! ---------------------
if (have_gom_qmr) then
  call pointer_field_array(inc%fields, 'q', q)
  call crtm_mixratio_ad(geom, traj%q, q, qmr)
endif

! Ozone
! -----
if (have_gom_o3) then
  do k = 1, npz
    do j = jsc, jec
      do i = isc, iec
        if (traj%o3(i,j,k) >= 0.0_kind_real) then
          o3(i,j,k) = o3(i,j,k) * constoz
        else
          o3(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
  call pointer_field_array(inc%fields, 'o3', o3p)
  o3p = o3p + o3
endif

! Surface pressure
! ----------------
if (have_gom_ps) then
  if (has_field(inc%fields,'ps')) then
    call pointer_field_array(inc%fields, 'ps', psp)
    psp = psp + ps
  elseif (has_field(inc%fields,'delp')) then
    call pointer_field_array(inc%fields, 'delp', delp)
    do k = 1,npz
      delp(:,:,k) = delp(:,:,k) + ps(:,:,1)
    enddo
  endif
endif

if (allocated(mod_increment)) deallocate(mod_increment)
deallocate(obs_increment)
deallocate(geovalm,geovale,geovals)

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

! --------------------------------------------------------------------------------------------------

subroutine variable_fail(ufo_var)

implicit none
character(len=*), intent(in) :: ufo_var

call abor1_ftn("GetValues.variable_fail: ufo variable "//trim(ufo_var)//" cannot be obtained from available model variables.")

end subroutine variable_fail

! --------------------------------------------------------------------------------------------------

end module fv3jedi_getvalues_mod
