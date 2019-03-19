module fv3jedi_getvalues_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use type_bump, only: bump_type
use ufo_geovals_mod, only: ufo_geovals
use ufo_locs_mod, only: ufo_locs
use ufo_vars_mod, only: ufo_vars

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
type(ufo_vars),                                 intent(in)    :: vars
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
real(kind=kind_real), target, allocatable :: geovale(:,:,:), geovalm(:,:,:)
real(kind=kind_real), pointer :: geoval(:,:,:)
integer :: nvl
logical :: do_interp

integer :: isc,iec,jsc,jec,npz,i,j

!Local pressure variables
real(kind=kind_real), allocatable :: delp(:,:,:) !Pressure thickness Pa
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
logical,  parameter               :: use_compress = .true.       !Could be a fv3 namelist option?



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

! Initialize the interpolation trajectory
! ---------------------------------------
if (present(traj)) then

  pbump => traj%bump

  if (.not. traj%lalloc) then
  
     traj%ngrid = ngrid
   
     if (.not.allocated(traj%t)) allocate(traj%t(isc:iec,jsc:jec,1:npz))
     if (.not.allocated(traj%q)) allocate(traj%q(isc:iec,jsc:jec,1:npz))
  
     traj%t = state%t
     traj%q = state%q
 
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

! Get pressures at edge, center & log center
! ------------------------------------------
allocate(delp(isc:iec,jsc:jec,npz  ))
allocate(prsi(isc:iec,jsc:jec,npz+1))
allocate(prs (isc:iec,jsc:jec,npz  ))
allocate(logp(isc:iec,jsc:jec,npz  ))

if (associated(state%delp)) then
  delp = state%delp
elseif (associated(state%ps)) then
  do jlev = 1,geom%npz
    delp(:,:,jlev) = (geom%ak(jlev+1)-geom%ak(jlev))+(geom%bk(jlev+1)-geom%bk(jlev))*state%ps(:,:,1)
  enddo
else
  call abor1_ftn("fv3jedi_getvalues_mod.getvalues: No way to compute delp from the state")
endif

call delp_to_pe_p_logp(geom,delp,prsi,prs,logp)

! Get CRTM surface variables
! ----------------------
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

if (associated(state%slmsk)) then
  !TODO only if a radiance
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

if (associated(state%slmsk)) then

  !TODO Is it water_coverage or sea_coverage fed in here?
  water_coverage_m = 0.0_kind_real
  do j = jsc,jec
    do i = isc,iec
      if (state%slmsk(i,j,1) == 0) water_coverage_m(i,j) = 1.0_kind_real
    enddo
  enddo
  
  call crtm_ade_efr( geom,prsi,state%t,delp, &
                     water_coverage_m,state%q, &
                     state%ql,state%qi, &
                     ql_ade,qi_ade,ql_efr,qi_efr )
  
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

  case ("eastward_wind")

    nvl = npz
    do_interp = .true.
    geovalm = state%ua
    geoval => geovalm

  case ("northward_wind")

    nvl = npz
    do_interp = .true.
    geovalm = state%va
    geoval => geovalm

  case ("air_temperature","temperature")

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

    call geop_height(geom,prs,prsi,state%t,state%q,&
                     state%phis(:,:,1),use_compress,geovalm)
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("geopotential_height_levels")

    call geop_height_levels(geom,prs,prsi,state%t,state%q,&
                            state%phis(:,:,1),use_compress,geovale)
    nvl = npz + 1
    do_interp = .true.
    geoval => geovale

  case ("sfc_geopotential_height")

    nvl = 1
    do_interp = .true.
    geovalm(:,:,1) = state%phis(:,:,1) / grav
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

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select

  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nv)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, locs, vars, gom, jvar)


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
if (allocated(delp                 )) deallocate(delp                 )
if (allocated(prsi                 )) deallocate(prsi                 )
if (allocated(prs                  )) deallocate(prs                  )
if (allocated(logp                 )) deallocate(logp                 )
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
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_tl'

integer :: ii, jj, ji, jvar, jlev, jloc
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

    nvl = inc%npz
    do_interp = .true.
    call T_to_Tv_tl(geom, traj%t, inc%t, traj%q, inc%q, geovalm )
    geoval => geovalm

  case ("humidity_mixing_ratio")
  
    nvl = inc%npz
    do_interp = .true.
    call crtm_mixratio_tl(geom, traj%q, inc%q, geovalm)
    geoval => geovalm  

  case ("air_pressure")

  case default

    call abor1_ftn(trim(myname)//"unknown variable: "//trim(vars%fldnames(jvar)))

  end select


  ! Allocate geovals%val for this jvars
  ! -----------------------------------
  call allocate_geovals_vals(gom,jvar,nvl,jvar==vars%nv)


  !Run some basic checks on the interpolation
  !------------------------------------------
  call getvalues_checks(myname, locs, vars, gom, jvar)


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
type(ufo_vars),           intent(in)    :: vars
type(ufo_geovals),        intent(inout) :: gom
type(fv3jedi_getvalues_traj), intent(inout)    :: traj

character(len=*), parameter :: myname = 'getvalues_ad'

integer :: ii, jj, ji, jvar, jlev, jloc
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

  case ("virtual_temperature")

    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("humidity_mixing_ratio")
  
    nvl = npz
    do_interp = .true.
    geoval => geovalm

  case ("air_pressure")
 
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

    inc%ua = geovalm

  case ("northward_wind")

    inc%va = geovalm

  case ("air_temperature","temperature")

    inc%t = geovalm

  case ("specific_humidity")

    inc%q = geovalm

  case ("virtual_temperature")
    
    call T_to_Tv_ad(geom, traj%t, inc%t, traj%q, inc%q, geovalm )

  case ("humidity_mixing_ratio")
  
    call crtm_mixratio_ad(geom, traj%q, inc%q, geovalm)

  case ("air_pressure")

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
  allocate(gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nobs))
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
bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)

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

subroutine getvalues_checks(cop, locs, vars, gom, jvar)
implicit none
character(len=*),  intent(in) :: cop
type(ufo_locs),    intent(in) :: locs
type(ufo_vars),    intent(in) :: vars
type(ufo_geovals), intent(in) :: gom
integer, intent(in)           :: jvar

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

end module fv3jedi_getvalues_mod
