program cold_to_warm

! fckit uses
use fckit_module,               only: fckit_main, fckit_configuration, fckit_pathname, &
                                      fckit_YAMLConfiguration
use fckit_mpi_module,           only: fckit_mpi_comm, fckit_mpi_comm, fckit_mpi_sum

! fms uses
use constants_mod,              only: grav, rdgas, rvgas
use fms_io_mod,                 only: nullify_domain
use fms_mod,                    only: fms_init
use mpp_mod,                    only: mpp_exit, mpp_pe, mpp_npes, mpp_error, FATAL, NOTE, &
                                      mpp_root_pe
use mpp_domains_mod,            only: domain2D, mpp_deallocate_domain, mpp_define_layout, &
                                      mpp_define_mosaic, mpp_define_io_domain, mpp_domains_exit, &
                                      mpp_domains_set_stack_size, mpp_update_domains
use field_manager_mod,          only: fm_string_len, field_manager_init, MODEL_ATMOS
use tracer_manager_mod,         only: get_number_tracers, get_tracer_names, get_tracer_index, &
                                      NO_TRACER, set_tracer_profile
use fms_io_mod,                 only: restart_file_type, register_restart_field
use fms_io_mod,                 only: free_restart_type, restore_state, save_restart
use fms_io_mod,                 only: set_domain, nullify_domain
use mpp_domains_mod,            only: east, north, center

! fv3 uses
use fv3jedi_fv3_arrays_mod,     only: fv_atmos_type, R_GRID
use fv3jedi_fv3_grid_utils_mod, only: mid_pt_sphere, get_unit_vect2, get_latlon_vector, &
                                      inner_prod, g_sum
use fv3jedi_fv3_control_mod,    only: fv_control_init

! fv3jedi uses ***this code should not make use of state or increment type***
use fv3jedi_geom_mod,           only: initialize_fms => initialize, fv3jedi_geom
use fv3jedi_fmsnamelist_mod,    only: fv3jedi_fmsnamelist
use fv3jedi_kinds_mod,          only: kind_int, kind_real


! Nothing implicit
implicit none

type :: state_cold_type
  integer :: isc, iec, jsc, jec, npz
  real(kind=kind_real), dimension(:,:,:),   allocatable :: u_w_cold, v_w_cold, u_s_cold, v_s_cold
  real(kind=kind_real), dimension(:,:,:),   allocatable :: ud_cold
  real(kind=kind_real), dimension(:,:,:),   allocatable :: vd_cold
  real(kind=kind_real), dimension(:,:,:),   allocatable :: t_cold
  real(kind=kind_real), dimension(:,:,:,:), allocatable :: q_cold
  real(kind=kind_real), dimension(:,:,:),   allocatable :: zh_cold
  real(kind=kind_real), dimension(:,:,:),   allocatable :: w_cold
  real(kind=kind_real), dimension(:,:),     allocatable :: ps_cold
  real(kind=kind_real), dimension(:,:),     allocatable :: orog_filt
end type state_cold_type

! Local variables
type(fv3jedi_geom) :: geom
integer :: arg_count, gtile, ptile=1, nlevs, ntracers, ntprog, ierr, isc, iec, jsc, jec, npx, npy, npz
integer :: this_grid
character(len=256) :: yaml_file
type(fckit_configuration) :: config, config_fms, config_geom, config_input, config_remap, &
                             config_output, config_test
type(fckit_mpi_comm)      :: comm
type(fv_atmos_type), allocatable :: Atm(:)
logical, allocatable :: grids_on_this_pe(:)
type(fv3jedi_fmsnamelist) :: fmsnamelist

type(state_cold_type)       :: state_cold

real(kind=kind_real), dimension(:), allocatable :: phis_flat

! Parameters
logical :: data_source_fv3gfs
real(kind=kind_real), parameter:: zvir = rvgas/rdgas - 1.
real(kind=kind_real), parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.

! Initialize fckit and communicator
! ---------------------------------
call fckit_main%initialise()
comm = fckit_mpi_comm("world")

! Check if a YAML file was provided as an argument
! ------------------------------------------------
arg_count = command_argument_count()
if (arg_count < 1) then
    print *, 'Usage: yaml_config_example <yaml_file>'
    stop 1
end if

! Initialize config
! -----------------
call get_command_argument(1, yaml_file)
config = fckit_YAMLConfiguration( fckit_pathname(yaml_file))

! Split into sub configs
! ----------------------
call config%get_or_die("fms",      config_fms)
call config%get_or_die("geometry", config_geom)
call config%get_or_die("input",    config_input)
call config%get_or_die("remap",    config_remap)
call config%get_or_die("output",   config_output)

! Initialize fms
! --------------
call initialize_fms(config_fms, comm)

! Intialize fv3 geometry
! ---------------------
call fmsnamelist%replace_namelist(config_geom)
call geom%create(config_geom, comm, npx, npy, npz)

! Create fv3 object
! -----------------
call fmsnamelist%replace_namelist(config_remap)
call fv_control_init(Atm, 300.0_kind_real, this_grid, grids_on_this_pe, ptile, &
                     skip_nml_read_in=.true.)
call fmsnamelist%revert_namelist

Atm(1)%flagstruct%nggps_ic = .true.
Atm(1)%ak = real(geom%ak,kind_real)
Atm(1)%bk = real(geom%bk,kind_real)
Atm%ptop = real(geom%ak(1),kind_real)

Atm(1)%u    = 0.0
Atm(1)%v    = 0.0
Atm(1)%pt   = 0.0
Atm(1)%delp = 0.0
Atm(1)%phis = 0.0
Atm(1)%q    = 0.0

! Aliases
isc = Atm(1)%bd%isc
iec = Atm(1)%bd%iec
jsc = Atm(1)%bd%jsc
jec = Atm(1)%bd%jec
npz = Atm(1)%npz

! Delete geom
! -----------
call geom%delete()

! Allocate the cold state
! -----------------------
allocate(state_cold% u_w_cold(isc:iec+1, jsc:jec  , 1:npz+1))
allocate(state_cold% v_w_cold(isc:iec+1, jsc:jec  , 1:npz+1))
allocate(state_cold% u_s_cold(isc:iec  , jsc:jec+1, 1:npz+1))
allocate(state_cold% v_s_cold(isc:iec  , jsc:jec+1, 1:npz+1))
allocate(state_cold%  ud_cold(isc:iec  , jsc:jec+1, 1:npz+1))
allocate(state_cold%  vd_cold(isc:iec+1, jsc:jec  , 1:npz+1))
allocate(state_cold%   t_cold(isc:iec  , jsc:jec  , 1:npz+1))
allocate(state_cold%  zh_cold(isc:iec  , jsc:jec  , 1:npz+2))
allocate(state_cold%   w_cold(isc:iec  , jsc:jec  , 1:npz+1))
allocate(state_cold%  ps_cold(isc:iec  , jsc:jec           ))
allocate(state_cold%orog_filt(isc:iec  , jsc:jec           ))

call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)
allocate(state_cold%q_cold(isc:iec,jsc:jec,1:npz+1,ntracers))

! Read the restarts
! -----------------

call read_restarts(Atm, config_input, state_cold)

! Convert the winds
! -----------------
call convert_winds(Atm, state_cold)

! Perform remapping
! -----------------

call remap(config_remap, state_cold, Atm)

! Write the restarts
! ------------------
call write_restarts(config_output, Atm)

! Check output numbers
! --------------------
if (config%has("testing")) then

  ! Get testing config
  call config%get_or_die("testing", config_test)

  ! Check output
  call testing(comm, geom, config_test, Atm)

endif


! --------------------------------------------------------------------------------------------------


contains


! --------------------------------------------------------------------------------------------------


subroutine read_restarts(Atm, conf, state_cold)

! Arguments
type(fv_atmos_type),       intent(in)    :: Atm(:)
type(fckit_configuration), intent(in)    :: conf
type(state_cold_type),     intent(inout) :: state_cold

! Locals
character(len=1024) :: datapath, fname_cold, fname_orog
character(len=:), allocatable :: str
integer :: spechum, liq_wat, ice_wat, rainwat, snowwat, graupel, ntclamt, o3mr
integer :: idrst
type(restart_file_type) :: rstc, rsto

! Read file path and names from config
! ------------------------------------
call conf%get_or_die("datapath", str)
datapath = str
deallocate(str)
call conf%get_or_die("filename_cold", str)
fname_cold = str
deallocate(str)
call conf%get_or_die("filename_orog", str)
fname_orog = str
deallocate(str)

! Winds
idrst = register_restart_field(rstc, trim(fname_cold), 'u_w_cold', state_cold%u_w_cold, &
                               domain=Atm(1)%domain, position=east  )
idrst = register_restart_field(rstc, trim(fname_cold), 'v_w_cold', state_cold%v_w_cold, &
                               domain=Atm(1)%domain, position=east  )
idrst = register_restart_field(rstc, trim(fname_cold), 'u_s_cold', state_cold%u_s_cold, &
                               domain=Atm(1)%domain, position=north )
idrst = register_restart_field(rstc, trim(fname_cold), 'v_s_cold', state_cold%v_s_cold, &
                               domain=Atm(1)%domain, position=north )
idrst = register_restart_field(rstc, trim(fname_cold), 'w_cold'  , state_cold%w_cold,   &
                               domain=Atm(1)%domain, position=center)

! Thermo
idrst = register_restart_field(rstc, trim(fname_cold), 't_cold'  , state_cold%t_cold,   &
                               domain=Atm(1)%domain, position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'ps_cold' , state_cold%ps_cold,  &
                               domain=Atm(1)%domain, position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'zh_cold' , state_cold%zh_cold,  &
                               domain=Atm(1)%domain, position=center)

! Tracers
spechum = get_tracer_index(MODEL_ATMOS, 'sphum'  )
liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
rainwat = get_tracer_index(MODEL_ATMOS, 'rainwat')
snowwat = get_tracer_index(MODEL_ATMOS, 'snowwat')
graupel = get_tracer_index(MODEL_ATMOS, 'graupel')

idrst = register_restart_field(rstc, trim(fname_cold), 'sphum_cold'  , &
                               state_cold%q_cold(:,:,:,spechum), domain=Atm(1)%domain, &
                               position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'liq_wat_cold', &
                               state_cold%q_cold(:,:,:,liq_wat), domain=Atm(1)%domain, &
                               position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'ice_wat_cold', &
                               state_cold%q_cold(:,:,:,ice_wat), domain=Atm(1)%domain, &
                               position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'rainwat_cold', &
                               state_cold%q_cold(:,:,:,rainwat), domain=Atm(1)%domain, &
                               position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'snowwat_cold', &
                               state_cold%q_cold(:,:,:,snowwat), domain=Atm(1)%domain, &
                               position=center)
idrst = register_restart_field(rstc, trim(fname_cold), 'graupel_cold', &
                               state_cold%q_cold(:,:,:,graupel), domain=Atm(1)%domain, &
                               position=center)

o3mr = get_tracer_index(MODEL_ATMOS, 'o3mr')
if (o3mr > 0) &
idrst = register_restart_field(rstc, trim(fname_cold), 'o3mr_cold', state_cold%q_cold(:,:,:,o3mr), &
                               domain=Atm(1)%domain, position=center)

ntclamt = get_tracer_index(MODEL_ATMOS, 'cld_amt')
if (ntclamt > 0) &
state_cold%q_cold(:,:,:,ntclamt) = 0.0
!idrst = register_restart_field(rstc, trim(fname_cold), 'cldamt_cold', &
!                               state_cold%q_cold(:,:,:,ntclamt), domain=Atm(1)%domain, &
!                               position=center)

! Read
call restore_state(rstc, directory=trim(adjustl(datapath)))
call free_restart_type(rstc)

! Orography
idrst = register_restart_field(rsto, trim(fname_orog), 'orog_filt', state_cold%orog_filt, &
                               domain=Atm(1)%domain, position=center)

call restore_state(rsto, directory=trim(adjustl(datapath)))
call free_restart_type(rsto)


end subroutine read_restarts


! --------------------------------------------------------------------------------------------------


subroutine write_restarts(conf, Atm)

! Arguments
type(fckit_configuration), intent(in) :: conf
type(fv_atmos_type),       intent(in) :: Atm(:)

! Locals
character(len=:), allocatable :: str
character(len=20) :: isodatetime
character(len=1024) :: prefix, dpath, fcore, ftrac
integer :: idrst, date(6)
type(restart_file_type) :: rcore, rtrac
integer :: nt, ntracers, ntprog
character(len=64) :: tracer_name
integer :: isc, iec, jsc, jec, npz

! Path to write files to
call conf%get_or_die("datapath", str)
dpath = str
deallocate(str)

! Check for prefix in the filename
prefix = ""
if (conf%has("prefix")) then
  call conf%get_or_die("prefix", str)
  prefix = str
  deallocate(str)
endif

! Indices
isc = Atm(1)%bd%isc
iec = Atm(1)%bd%iec
jsc = Atm(1)%bd%jsc
jec = Atm(1)%bd%jec
npz = Atm(1)%npz

! Filenames
fcore = trim(adjustl(prefix))//"fv_core.res.nc"
ftrac = trim(adjustl(prefix))//"fv_tracer.res.nc"

idrst = register_restart_field( rcore, trim(fcore), 'u', Atm(1)%u(isc:iec, jsc:jec+1, 1:npz), &
                                domain=Atm(1)%domain, &
                                position=north, longname = 'u_component_of_native_D_grid_wind', &
                                units = 'ms-1' )
idrst = register_restart_field( rcore, trim(fcore), 'v', Atm(1)%v(isc:iec+1, jsc:jec, 1:npz), &
                                domain=Atm(1)%domain, &
                                position=east, longname = 'v_component_of_native_D_grid_wind', &
                                units = 'ms-1' )
idrst = register_restart_field( rcore, trim(fcore), 'T', Atm(1)%pt(isc:iec,jsc:jec,1:npz), &
                                domain=Atm(1)%domain, &
                                position=center, longname = 'air_temperature', &
                                units = 'K' )
idrst = register_restart_field( rcore, trim(fcore), 'w', Atm(1)%w(isc:iec,jsc:jec,1:npz), &
                                domain=Atm(1)%domain, &
                                position=center, longname = 'upward_air_velocity', &
                                units = 'ms-1' )
idrst = register_restart_field( rcore, trim(fcore), 'DELP', Atm(1)%delp(isc:iec,jsc:jec,1:npz), &
                                domain=Atm(1)%domain, &
                                position=center, longname = 'air_pressure_thickness', &
                                units = 'Pa' )
idrst = register_restart_field( rcore, trim(fcore), 'delz', Atm(1)%delz(isc:iec,jsc:jec,1:npz),&
                                domain=Atm(1)%domain, &
                                position=center, longname = 'layer_thickness', &
                                units = 'm' )
idrst = register_restart_field( rcore, trim(fcore), 'phis', Atm(1)%phis, domain=Atm(1)%domain, &
                                position=center, longname = 'sfc_geopotential_height_times_grav', &
                                units = 'm' )

! Register tracers
call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)
do nt = 1, ntracers
  call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
  idrst = register_restart_field( rtrac, trim(ftrac), tracer_name, &
                                  Atm(1)%q(isc:iec,jsc:jec,1:npz,nt), &
                                  domain=Atm(1)%domain, position=center, longname = tracer_name, &
                                  units = 'kgkg-1' )
enddo

! Save the restarts
call save_restart(rcore, directory=trim(adjustl(dpath)))
call save_restart(rtrac, directory=trim(adjustl(dpath)))

! Free up
call free_restart_type(rcore)
call free_restart_type(rtrac)

! Get datetime from config (yyyy-mm-ddThh:mm:ss)
call conf%get_or_die("datetime", str)
isodatetime = str
deallocate(str)

! Convert first for digits to integer year
read(isodatetime(1:4),*) date(1)
read(isodatetime(6:7),*) date(2)
read(isodatetime(9:10),*) date(3)
read(isodatetime(12:13),*) date(4)
read(isodatetime(15:16),*) date(5)
read(isodatetime(18:19),*) date(6)

!Write date/time info in coupler.res
!-----------------------------------
if (mpp_pe() == mpp_root_pe()) then
   open(101, file = trim(adjustl(dpath))//'/'//trim(prefix)//'coupler.res', form='formatted')
   write( 101, '(i6,8x,a)' ) 2, &
        '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
   write( 101, '(6i6,8x,a)') date, 'Model start time:   year, month, day, hour, minute, second'
   write( 101, '(6i6,8x,a)') date, 'Current model time: year, month, day, hour, minute, second'
   close(101)
endif

end subroutine write_restarts


! --------------------------------------------------------------------------------------------------


subroutine convert_winds(Atm, state_cold)

! Arguments
type(fv_atmos_type),   intent(in)    :: Atm(:)
type(state_cold_type), intent(inout) :: state_cold

! Locals
integer :: i, j, k, levp
real(kind=R_GRID), dimension(2) :: p1, p2, p3
real(kind=R_GRID), dimension(3) :: e1, e2, ex, ey

! Potential for different levels coming in
levp = size(state_cold%u_w_cold,3)

do k = 1, levp
do j = Atm(1)%bd%jsc, Atm(1)%bd%jec+1
    do i = Atm(1)%bd%isc, Atm(1)%bd%iec
    p1(:) = Atm(1)%gridstruct%grid(i,  j,1:2)
    p2(:) = Atm(1)%gridstruct%grid(i+1,j,1:2)
    call  mid_pt_sphere(p1, p2, p3)
    call get_unit_vect2(p1, p2, e1)
    call get_latlon_vector(p3, ex, ey)
    state_cold%ud_cold(i,j,k) = state_cold%u_s_cold(i,j,k)*inner_prod(e1, ex) + &
                                state_cold%v_s_cold(i,j,k)*inner_prod(e1, ey)
    enddo
enddo
do j = Atm(1)%bd%jsc, Atm(1)%bd%jec
    do i = Atm(1)%bd%isc, Atm(1)%bd%iec+1
    p1(:) = Atm(1)%gridstruct%grid(i,j  ,1:2)
    p2(:) = Atm(1)%gridstruct%grid(i,j+1,1:2)
    call  mid_pt_sphere(p1, p2, p3)
    call get_unit_vect2(p1, p2, e2)
    call get_latlon_vector(p3, ex, ey)
    state_cold%vd_cold(i,j,k) = state_cold%u_w_cold(i,j,k)*inner_prod(e2, ex) + &
                                state_cold%v_w_cold(i,j,k)*inner_prod(e2, ey)
    enddo
enddo
enddo

end subroutine convert_winds


! --------------------------------------------------------------------------------------------------


subroutine remap(conf, state_cold, Atm)

! Arguments
type(fckit_configuration), intent(in)    :: conf
type(state_cold_type),     intent(in)    :: state_cold
type(fv_atmos_type),       intent(inout) :: Atm(:)

! Config
logical :: from_cold_start, checker_tr
integer :: nt_checker
character(len=:), allocatable :: str
integer:: i, j, k, nt, ntracers, ntprog, itoa, levp, isc, iec, jsc, jec, npz, nts
integer:: liq_wat, ice_wat, rainwat, snowwat, graupel, ntclamt
character(len=64) :: tracer_name
real(kind=kind_real), allocatable :: ak(:), bk(:)
real(kind=kind_real) :: wt, qt, m_fac


! Parse the config
! ----------------
if( .not. conf%get('input is cold starts', from_cold_start) ) from_cold_start = .true.
if( .not. conf%get('check tracers',    checker_tr) ) checker_tr = .false.
if( .not. conf%get('check tracers nt', nt_checker) ) nt_checker = 0


! Remap to new Lagrangian vertical coordinate
! -------------------------------------------
  ! Shortcuts
  isc = Atm(1)%bd%is
  iec = Atm(1)%bd%ie
  jsc = Atm(1)%bd%js
  jec = Atm(1)%bd%je
  npz = Atm(1)%npz

  ! Orography
  Atm(1)%phis(isc:iec,jsc:jec) = state_cold%orog_filt(isc:iec,jsc:jec) * grav  ! Convert to phis

  ! akbk from cold starts have different levels (for now extra levels are zero)
  levp = size(state_cold%w_cold,3)
  itoa = levp - npz + 1

  allocate(ak(levp+1))
  allocate(bk(levp+1))
  ak = 0.0_kind_real
  bk = 0.0_kind_real
  ak(itoa:levp+1) = Atm(1)%ak(1:npz+1)
  bk(itoa:levp+1) = Atm(1)%bk(1:npz+1)
  ak(1) = max(1.e-9_kind_real, ak(1))

  ! Tracers
  ! -------
  ! Number of tracers in Atm
  call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)

  ! initialize all tracers to default values prior to being input
  do nt = 1, ntprog
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    ! set all tracers to an initial profile value
    call set_tracer_profile (MODEL_ATMOS, nt, Atm(1)%q(:,:,:,nt)  )
  enddo
  do nt = ntprog+1, ntracers
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    ! set all tracers to an initial profile value
    call set_tracer_profile (MODEL_ATMOS, nt, Atm(1)%qdiag(:,:,:,nt)  )
  enddo

  ! Call remapping non-wind variables
  ! ---------------------------------
  if (allocated(state_cold%t_cold)) then
    call remap_scalar(Atm(1), levp, npz, ntracers, ak, bk, state_cold%ps_cold, &
                      state_cold%q_cold, state_cold%zh_cold, state_cold%w_cold, state_cold%t_cold)
  else
    call remap_scalar(Atm(1), levp, npz, ntracers, ak, bk, state_cold%ps_cold, &
                      state_cold%q_cold, state_cold%zh_cold, state_cold%w_cold)
  endif

  ! Call remapping wind variables
  ! -----------------------------
  call remap_dwinds(levp, npz, ak, bk, state_cold%ps_cold, state_cold%ud_cold, &
                    state_cold%vd_cold, Atm(1))

  ! Tracer weighting
  ! ----------------
  liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
  ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
  rainwat = get_tracer_index(MODEL_ATMOS, 'rainwat')
  snowwat = get_tracer_index(MODEL_ATMOS, 'snowwat')
  graupel = get_tracer_index(MODEL_ATMOS, 'graupel')
  ntclamt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

  if (from_cold_start) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          wt = Atm(1)%delp(i,j,k)
          if ( Atm(1)%flagstruct%nwat == 6 ) then
            qt = wt*(1.0_kind_real + Atm(1)%q(i,j,k,liq_wat) + Atm(1)%q(i,j,k,ice_wat) + &
                                    Atm(1)%q(i,j,k,rainwat) + Atm(1)%q(i,j,k,snowwat) + &
                                    Atm(1)%q(i,j,k,graupel))
          else
            qt = wt*(1.0_kind_real + sum(Atm(1)%q(i,j,k,2:Atm(1)%flagstruct%nwat)))
          endif
          Atm(1)%delp(i,j,k) = qt
          if (ntclamt > 0) Atm(1)%q(i,j,k,ntclamt) = 0.0
        enddo
      enddo
    enddo
  else
    ! TODO Question is do we need this if just remapping after adding the increment, say?
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          wt = Atm(1)%delp(i,j,k)
          if ( Atm(1)%flagstruct%nwat == 6 ) then
            qt = wt*(1.0_kind_real + Atm(1)%q(i,j,k,liq_wat) + Atm(1)%q(i,j,k,ice_wat) + &
                                    Atm(1)%q(i,j,k,rainwat) + Atm(1)%q(i,j,k,snowwat) + &
                                    Atm(1)%q(i,j,k,graupel))
          else
             qt = wt*(1.0_kind_real + sum(Atm(1)%q(i,j,k,2:Atm(1)%flagstruct%nwat)))
          endif
          m_fac = wt / qt
          do nt=1,ntracers
            Atm(1)%q(i,j,k,nt) = m_fac * Atm(1)%q(i,j,k,nt)
          enddo
          Atm(1)%delp(i,j,k) = qt
          if (ntclamt > 0) Atm(1)%q(i,j,k,ntclamt) = 0.0
        enddo
      enddo
    enddo
  endif

if (checker_tr) then
  nts = ntracers - nt_checker+1
  call checker_tracers(isc, iec, jsc, jec, Atm(1)%bd%isd, Atm(1)%bd%ied, Atm(1)%bd%jsd, &
                       Atm(1)%bd%jed, nt_checker, npz, Atm(1)%q(:,:,:,nts:ntracers), &
                       Atm(1)%gridstruct%agrid_64(isc:iec,jsc:jec,1),     &
                       Atm(1)%gridstruct%agrid_64(isc:iec,jsc:jec,2), &
                       9.0_kind_real, 9.0_kind_real)
endif


end subroutine remap


! --------------------------------------------------------------------------------------------------


subroutine testing(comm, geom, conf, Atm)

type(fckit_mpi_comm),      intent(in) :: comm
type(fv3jedi_geom),        intent(in) :: geom
type(fckit_configuration), intent(in) :: conf
type(fv_atmos_type),       intent(in) :: Atm(:)

real(kind=kind_real) :: rms_u_c, rms_v_c, rms_t_c, rms_d_c, rms_p_c, tol
real(kind=kind_real) :: rms_u_r, rms_v_r, rms_t_r, rms_d_r, rms_p_r
integer :: isc, iec, jsc, jec, npz
real(kind=kind_real) :: tmp(5), gs3, gs3g

! Get min/max/rms from config
call conf%get_or_die("rms_u", rms_u_c)
call conf%get_or_die("rms_v", rms_v_c)
call conf%get_or_die("rms_t", rms_t_c)
call conf%get_or_die("rms_delp", rms_d_c)
call conf%get_or_die("rms_phis", rms_p_c)
call conf%get_or_die("tol", tol)

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

! Compute global sum over the field
gs3 = real((geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)*geom%npz, kind_real)
call comm%allreduce(gs3,gs3g,fckit_mpi_sum())

! Min/Max/SumSquares (u)
tmp(1) = sum(Atm(1)%u(isc:iec,jsc:jec,1:npz)**2)
tmp(2) = sum(Atm(1)%v(isc:iec,jsc:jec,1:npz)**2)
tmp(3) = sum(Atm(1)%pt(isc:iec,jsc:jec,1:npz)**2)
tmp(4) = sum(Atm(1)%delp(isc:iec,jsc:jec,1:npz)**2)
tmp(5) = sum(Atm(1)%phis(isc:iec,jsc:jec)**2)

call comm%allreduce(tmp(1), rms_u_r, fckit_mpi_sum())
call comm%allreduce(tmp(2), rms_v_r, fckit_mpi_sum())
call comm%allreduce(tmp(3), rms_t_r, fckit_mpi_sum())
call comm%allreduce(tmp(4), rms_d_r, fckit_mpi_sum())
call comm%allreduce(tmp(5), rms_p_r, fckit_mpi_sum())
rms_u_r = sqrt(rms_u_r/gs3g)
rms_v_r = sqrt(rms_v_r/gs3g)
rms_t_r = sqrt(rms_t_r/gs3g)
rms_d_r = sqrt(rms_d_r/gs3g)
rms_p_r = sqrt(rms_p_r/gs3g)

! If root processor then perform checks
if (comm%rank() == 0) then

  ! Assert that min/max/rms are within tolerance (relative)
  ! ------------------------------------------------------
  if (abs(rms_u_r-rms_u_c)/abs(rms_u_c) > tol) then
      print *, 'u rms out of tolerance (run) (config) (rel diff) (tol)', rms_u_r, rms_u_c, &
               abs(rms_u_r-rms_u_c)/abs(rms_u_c), tol
      stop 1
  end if

  if (abs(rms_v_r-rms_v_c)/abs(rms_v_c) > tol) then
      print *, 'v rms out of tolerance (run) (config) (rel diff) (tol)', rms_v_r, rms_v_c, &
               abs(rms_v_r-rms_v_c)/abs(rms_v_c), tol
      stop 1
  end if

  if (abs(rms_t_r-rms_t_c)/abs(rms_t_c) > tol) then
      print *, 't rms out of tolerance (run) (config) (rel diff) (tol)', rms_t_r, rms_t_c, &
               abs(rms_t_r-rms_t_c)/abs(rms_t_c), tol
      stop 1
  end if

  if (abs(rms_d_r-rms_d_c)/abs(rms_d_c) > tol) then
      print *, 'delp rms out of tolerance (run) (config) (rel diff) (tol)', rms_d_r, rms_d_c, &
               abs(rms_d_r-rms_d_c)/abs(rms_d_c), tol
      stop 1
  end if

  if (abs(rms_p_r-rms_p_c)/abs(rms_p_c) > tol) then
      print *, 'phis rms out of tolerance (run) (config) (rel diff) (tol)', rms_p_r, rms_p_c, &
               abs(rms_p_r-rms_p_c)/abs(rms_p_c), tol
      stop 1
  end if

  print*, 'Variable rms values are within tolerance (PASSED)'

end if

end subroutine testing


! --------------------------------------------------------------------------------------------------

subroutine checker_tracers(i0, i1, j0, j1, ifirst, ilast, jfirst, jlast,  &
                           nq, km, q, lon, lat, nx, ny, rn)

! Arguments
integer,                       intent(in)  :: nq, km
integer,                       intent(in)  :: i0, i1, j0, j1
integer,                       intent(in)  :: ifirst, ilast, jfirst, jlast
real(kind=kind_real),           intent(in)  :: nx, ny
real(kind=kind_real), optional, intent(in)  :: rn
real(kind=kind_real),            intent(in)  :: lon(i0:i1,j0:j1), lat(i0:i1,j0:j1)
real(kind=kind_real),           intent(out) :: q(ifirst:ilast,jfirst:jlast,km,nq)

! Locals
real(kind=kind_real) :: qt(i0:i1,j0:j1)
real(kind=kind_real) :: qtmp, ftmp
integer :: i, j, k, iq

do j = j0, j1
  do i = i0, i1
    qtmp = sin(nx * lon(i,j)) * sin(ny * lat(i,j))
    if (qtmp < 0.) then
      qt(i,j) = 0.0
    else
      qt(i,j) = 1.0
    end if
  end do
end do

if (present(rn)) then
  do iq = 1, nq
    call random_seed()
    do k = 1, km
      do j = j0, j1
        do i = i0, i1
          call random_number(ftmp)
          q(i,j,k,iq) = qt(i,j) + rn * ftmp
        end do
      end do
    end do
  end do
else
  do iq = 1, nq
    do k = 1, km
      do j = j0, j1
        do i = i0, i1
          q(i,j,k,iq) = qt(i,j)
        end do
      end do
    end do
  end do
end if

end subroutine checker_tracers

! --------------------------------------------------------------------------------------------------

subroutine remap_dwinds(km, npz, ak0, bk0, psc, ud, vd, Atm)

  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: km, npz
  real(kind=kind_real),    intent(in):: ak0(km+1), bk0(km+1)
  real(kind=kind_real),    intent(in):: psc(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je)
  real(kind=kind_real),    intent(in)::  ud(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je+1,km)
  real(kind=kind_real),    intent(in)::  vd(Atm%bd%is:Atm%bd%ie+1,Atm%bd%js:Atm%bd%je,km)
! local:
  real(kind=kind_real), dimension(Atm%bd%isd:Atm%bd%ied,Atm%bd%jsd:Atm%bd%jed):: psd
  real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie+1, km+1):: pe0
  real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie+1,npz+1):: pe1
  real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie+1,npz):: qn1
  integer i,j,k
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed

!Not sure what this is for
  if (Atm%gridstruct%bounded_domain) then
     do j=jsd,jed
     do i=isd,ied
        psd(i,j) = Atm%ps(i,j)
     enddo
     enddo
  else
     do j=js,je
        do i=is,ie
           psd(i,j) = psc(i,j)
        enddo
     enddo
  endif
  call mpp_update_domains( psd,    Atm%domain, complete=.false. )
  call mpp_update_domains( Atm%ps, Atm%domain, complete=.true. )

  do 5000 j=js,je+1
!------
! map u
!------
     do k=1,km+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*0.5*(psd(i,j-1)+psd(i,j))
        enddo
     enddo
     do k=1,npz+1
        do i=is,ie
           pe1(i,k) = Atm%ak(k) + Atm%bk(k)*0.5*(Atm%ps(i,j-1)+Atm%ps(i,j))
        enddo
     enddo
     call mappm(km, pe0(is:ie,1:km+1), ud(is:ie,j,1:km), npz, pe1(is:ie,1:npz+1),   &
                qn1(is:ie,1:npz), is,ie, -1, 8, Atm%ptop)
     do k=1,npz
        do i=is,ie
           Atm%u(i,j,k) = qn1(i,k)
        enddo
     enddo
!------
! map v
!------
     if ( j/=(je+1) ) then

     do k=1,km+1
        do i=is,ie+1
           pe0(i,k) = ak0(k) + bk0(k)*0.5*(psd(i-1,j)+psd(i,j))
        enddo
     enddo
     do k=1,npz+1
        do i=is,ie+1
           pe1(i,k) = Atm%ak(k) + Atm%bk(k)*0.5*(Atm%ps(i-1,j)+Atm%ps(i,j))
        enddo
     enddo
     call mappm(km, pe0(is:ie+1,1:km+1), vd(is:ie+1,j,1:km), npz, pe1(is:ie+1,1:npz+1),  &
                qn1(is:ie+1,1:npz), is,ie+1, -1, 8, Atm%ptop)
     do k=1,npz
        do i=is,ie+1
           Atm%v(i,j,k) = qn1(i,k)
        enddo
     enddo

     endif

5000 continue

  if (mpp_pe()==mpp_root_pe()) write(*,*) 'done remap_dwinds'

 end subroutine remap_dwinds

! --------------------------------------------------------------------------------------------------

 subroutine remap_scalar(Atm, km, npz, ncnst, ak0, bk0, psc, qa, zh, omga, t_in)

   implicit none

   type(fv_atmos_type), intent(inout) :: Atm
   integer, intent(in):: km, npz, ncnst
   real(kind=kind_real),    intent(in):: ak0(km+1), bk0(km+1)
   real(kind=kind_real),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc
   real(kind=kind_real),    intent(in), optional, dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: omga, t_in
   real(kind=kind_real),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
   real(kind=kind_real),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km+1):: zh
 ! local:
   real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0
   real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1, dp2
   real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1
   real(kind=kind_real) qp(Atm%bd%is:Atm%bd%ie,km)
   real(kind=kind_real) wk(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je)
   real(kind=kind_real), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: z500
 !!! High-precision
   real(kind=R_GRID), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pn1
   real(kind=R_GRID):: gz_fv(npz+1)
   real(kind=R_GRID), dimension(2*km+1):: gz, pn
   real(kind=R_GRID), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pn0
   real(kind=R_GRID):: pst, pe0tmp
 !!! High-precision
   integer i,j,k,l,m, k2,iq
   integer  sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, hailwat, cld_amt, sgs_tke
   integer o3mr
   integer :: is,  ie,  js,  je

   is  = Atm%bd%is
   ie  = Atm%bd%ie
   js  = Atm%bd%js
   je  = Atm%bd%je

   sphum   = get_tracer_index(MODEL_ATMOS, 'sphum')
   liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
   ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
   rainwat = get_tracer_index(MODEL_ATMOS, 'rainwat')
   snowwat = get_tracer_index(MODEL_ATMOS, 'snowwat')
   graupel = get_tracer_index(MODEL_ATMOS, 'graupel')
   hailwat = get_tracer_index(MODEL_ATMOS, 'hailwat')
   cld_amt = get_tracer_index(MODEL_ATMOS, 'cld_amt')
   o3mr    = get_tracer_index(MODEL_ATMOS, 'o3mr')
   sgs_tke = get_tracer_index(MODEL_ATMOS, 'sgs_tke')

   if (mpp_pe()==1) then
     print *, 'In remap_scalar:'
     print *, 'ncnst = ', ncnst
     print *, 'nwat = ', Atm%flagstruct%nwat
     print *, 'sphum    = ', sphum
     print *, 'clwmr    = ', liq_wat
     print *, 'liq_wat  = ', liq_wat
     print *, ' o3mr    = ', o3mr
    if ( Atm%flagstruct%nwat .ge. 6 ) then
       print *, 'rainwat = ', rainwat
       print *, 'ice_wat = ', ice_wat
       print *, 'snowwat = ', snowwat
       print *, 'graupel = ', graupel
       IF ( Atm%flagstruct%nwat == 7 ) print *, 'hailwat = ', hailwat
     endif
     print *, 'sgs_tke = ', sgs_tke
     print *, 'cld_amt = ', cld_amt
   endif

   if ( sphum/=1 ) then
        call mpp_error(FATAL,'SPHUM must be 1st tracer')
   endif

   k2 = max(10, km/2)

   do 5000 j=js,je
      do k=1,km+1
         do i=is,ie
            pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
            pe0tmp = real(pe0(i,k), kind=R_GRID)
            pn0(i,k) = log(pe0tmp)
         enddo
      enddo

      do i=is,ie
         do k=1,km+1
            pn(k) = pn0(i,k)
            gz(k) = zh(i,j,k)*grav
         enddo
 ! Use log-p for interpolation/extrapolation
 ! mirror image method:
         do k=km+2, km+k2
                l = 2*(km+1) - k
            gz(k) = 2.*gz(km+1) - gz(l)
            pn(k) = 2.*pn(km+1) - pn(l)
         enddo

         do k=km+k2-1, 2, -1
           if( Atm%phis(i,j).le.gz(k) .and. Atm%phis(i,j).ge.gz(k+1) ) then
               pst = pn(k) + (pn(k+1)-pn(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
               go to 123
           endif
         enddo
 123     Atm%ps(i,j) = exp(pst)

 ! ------------------
 ! Find 500-mb height
 ! ------------------
         pst = log(500.e2)
         do k=km+k2-1, 2, -1
           if( pst.le.pn(k+1) .and. pst.ge.pn(k) ) then
               z500(i,j) = (gz(k+1) + (gz(k)-gz(k+1))*(pn(k+1)-pst)/(pn(k+1)-pn(k)))/grav
               go to 124
           endif
         enddo
 124     continue

      enddo   ! i-loop

      do i=is,ie
         pe1(i,1) = Atm%ak(1)
         pn1(i,1) = log(pe1(i,1))
      enddo
      do k=2,npz+1
        do i=is,ie
           pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
           pn1(i,k) = log(pe1(i,k))
        enddo
      enddo

 ! * Compute delp
      do k=1,npz
         do i=is,ie
            dp2(i,k) = pe1(i,k+1) - pe1(i,k)
            Atm%delp(i,j,k) = dp2(i,k)
         enddo
      enddo

 ! map tracers
       do iq=1,ncnst
         if (floor(qa(is,j,1,iq)) == -1000) cycle !skip missing scalars [floor(-999.99) is -1000]
          do k=1,km
             do i=is,ie
                qp(i,k) = qa(i,j,k,iq)
             enddo
          enddo
          call mappm(km, pe0, qp, npz, pe1,  qn1, is,ie, 0, 8, Atm%ptop)
          if ( iq==sphum ) then
             call fillq(ie-is+1, npz, 1, qn1, dp2)
          else
             call fillz(ie-is+1, npz, 1, qn1, dp2)
          endif
 ! The HiRam step of blending model sphum with NCEP data is obsolete because nggps is always cold starting...
          do k=1,npz
             do i=is,ie
                Atm%q(i,j,k,iq) = qn1(i,k)
             enddo
          enddo
       enddo

 !---------------------------------------------------
 ! Retrive temperature using  geopotential height from external data
 !---------------------------------------------------
    do i=is,ie
 ! Make sure FV3 top is lower than GFS; can not do extrapolation above the top at this point
       if ( pn1(i,1) .lt. pn0(i,1) ) then
            call mpp_error(FATAL,'FV3 top higher than external data')
       endif

       do k=1,km+1
          pn(k) = pn0(i,k)
          gz(k) = zh(i,j,k)*grav
       enddo
 !-------------------------------------------------
       do k=km+2, km+k2
          l = 2*(km+1) - k
          gz(k) = 2.*gz(km+1) - gz(l)
          pn(k) = 2.*pn(km+1) - pn(l)
       enddo
 !-------------------------------------------------

       gz_fv(npz+1) = Atm%phis(i,j)

       m = 1

       do k=1,npz
 ! Searching using FV3 log(pe): pn1
          do l=m,km+k2-1
             if ( (pn1(i,k).le.pn(l+1)) .and. (pn1(i,k).ge.pn(l)) ) then
                 gz_fv(k) = gz(l) + (gz(l+1)-gz(l))*(pn1(i,k)-pn(l))/(pn(l+1)-pn(l))
                 goto 555
             endif
          enddo
 555   m = l
       enddo

       do k=1,npz+1
          Atm%peln(i,k,j) = pn1(i,k)
       enddo

 !----------------------------------------------------
 ! Compute true temperature using hydrostatic balance
 !----------------------------------------------------
       if (.not. data_source_fv3gfs .or. .not. present(t_in)) then
         do k=1,npz
            Atm%pt(i,j,k) = (gz_fv(k)-gz_fv(k+1))/( rdgas*(pn1(i,k+1)-pn1(i,k))*(1.+zvir*Atm%q(i,j,k,sphum)) )
         enddo
 !------------------------------
 ! Remap input T logarithmically in p.
 !------------------------------
       else
         do k=1,km
             qp(i,k) = t_in(i,j,k)
         enddo

         call mappm(km, log(pe0), qp, npz, log(pe1), qn1, is,ie, 2, 4, Atm%ptop) ! pn0 and pn1 are higher-precision
                                                                                 ! and cannot be passed to mappm
         do k=1,npz
             Atm%pt(i,j,k) = qn1(i,k)
         enddo
       endif
       if ( .not. Atm%flagstruct%hydrostatic ) then
          do k=1,npz
             Atm%delz(i,j,k) = (gz_fv(k+1) - gz_fv(k)) / grav
          enddo
       endif

    enddo   ! i-loop

 !-----------------------------------------------------------------------
 ! seperate cloud water and cloud ice from Jan-Huey Chen's HiRAM code
 ! only use for NCEP IC and GFDL microphy
 !-----------------------------------------------------------------------
    if (.not. data_source_fv3gfs) then
       if ((Atm%flagstruct%nwat .eq. 3 .or. Atm%flagstruct%nwat .eq. 6) .and. &
            (Atm%flagstruct%ncep_ic .or. Atm%flagstruct%nggps_ic)) then
          do k=1,npz
             do i=is,ie

                qn1(i,k) = Atm%q(i,j,k,liq_wat)
                if (cld_amt .gt. 0) Atm%q(i,j,k,cld_amt) = 0.

                if ( Atm%pt(i,j,k) > 273.16 ) then       ! > 0C all liq_wat
                   Atm%q(i,j,k,liq_wat) = qn1(i,k)
                   Atm%q(i,j,k,ice_wat) = 0.
                else if ( Atm%pt(i,j,k) < 233.16 ) then  ! < -40C all ice_wat
                   Atm%q(i,j,k,liq_wat) = 0.
                   Atm%q(i,j,k,ice_wat) = qn1(i,k)
                else
                   if ( k.eq.1 ) then  ! between [-40,0]: linear interpolation
                      Atm%q(i,j,k,liq_wat) = qn1(i,k)*((Atm%pt(i,j,k)-233.16)/40.)
                      Atm%q(i,j,k,ice_wat) = qn1(i,k) - Atm%q(i,j,k,liq_wat)
                   else
                      if (Atm%pt(i,j,k)<258.16 .and. Atm%q(i,j,k-1,ice_wat)>1.e-5 ) then
                         Atm%q(i,j,k,liq_wat) = 0.
                         Atm%q(i,j,k,ice_wat) = qn1(i,k)
                      else  ! between [-40,0]: linear interpolation
                         Atm%q(i,j,k,liq_wat) = qn1(i,k)*((Atm%pt(i,j,k)-233.16)/40.)
                         Atm%q(i,j,k,ice_wat) = qn1(i,k) - Atm%q(i,j,k,liq_wat)
                      endif
                   endif
                endif

                if (Atm%flagstruct%nwat .eq. 6) then ! no need to check for nwat=7 (hail) since only nwat=3,6 treated here
                   Atm%q(i,j,k,rainwat) = 0.
                   Atm%q(i,j,k,snowwat) = 0.
                   Atm%q(i,j,k,graupel) = 0.
                   call mp_auto_conversion(Atm%q(i,j,k,liq_wat), Atm%q(i,j,k,rainwat),  &
                        Atm%q(i,j,k,ice_wat), Atm%q(i,j,k,snowwat) )
                endif
             enddo
          enddo

       endif
   endif ! data source /= FV3GFS GAUSSIAN NEMSIO/NETCDF and GRIB2 FILE

 ! For GFS spectral input, omega in pa/sec is stored as w in the input data so actual w(m/s) is calculated
 ! For GFS nemsio input, omega is 0, so best not to use for input since boundary data will not exist for w
 ! For FV3GFS NEMSIO input, w is already in m/s (but the code reads in as omga) and just needs to be remapped
 !-------------------------------------------------------------
 ! map omega or w
 !------- ------------------------------------------------------
    if ( (.not. Atm%flagstruct%hydrostatic) .and. (.not. Atm%flagstruct%ncep_ic) ) then
       do k=1,km
          do i=is,ie
             qp(i,k) = omga(i,j,k)
          enddo
       enddo
       call mappm(km, pe0, qp, npz, pe1, qn1, is,ie, -1, 4, Atm%ptop)
     if (data_source_fv3gfs) then
       do k=1,npz
          do i=is,ie
             atm%w(i,j,k) = qn1(i,k)
          enddo
       enddo
     else
       do k=1,npz
          do i=is,ie
             atm%w(i,j,k) = qn1(i,k)/atm%delp(i,j,k)*atm%delz(i,j,k)
          enddo
       enddo
      endif
    endif

 5000 continue

 ! Add some diagnostics:
   if (.not. Atm%flagstruct%hydrostatic) call p_maxmin('delz_model', Atm%delz, is, ie, js, je, npz, 1._kind_real)
   call p_maxmin('sphum_model', Atm%q(is:ie,js:je,1:npz,sphum), is, ie, js, je, npz, 1._kind_real)
   call p_maxmin('liq_wat_model', Atm%q(is:ie,js:je,1:npz,liq_wat), is, ie, js, je, npz, 1._kind_real)
   if (ice_wat .gt. 0) call p_maxmin('ice_wat_model', Atm%q(is:ie,js:je,1:npz,ice_wat), is, ie, js, je, npz, 1._kind_real)
   call p_maxmin('PS_model (mb)', Atm%ps(is:ie,js:je), is, ie, js, je, 1, 0.01_kind_real)
   call p_maxmin('PT_model', Atm%pt(is:ie,js:je,1:npz), is, ie, js, je, npz, 1._kind_real)
   call pmaxmn('ZS_model', Atm%phis(is:ie,js:je)/grav, is, ie, js, je, 1, 1._kind_real, Atm%gridstruct%area_64, Atm%domain)
   call pmaxmn('ZS_data', zh(is:ie,js:je,km+1), is, ie, js, je, 1, 1._kind_real, Atm%gridstruct%area_64, Atm%domain)
   do j=js,je
      do i=is,ie
         wk(i,j) = Atm%phis(i,j)/grav - zh(i,j,km+1)
   !      if ((wk(i,j) > 1800.).or.(wk(i,j)<-1600.)) then
   !         print *,'  '
   !         print *, 'Diff = ', wk(i,j), 'Atm%phis =', Atm%phis(i,j)/grav, 'zh = ', zh(i,j,km+1)
   !         print *, 'lat = ', Atm%gridstruct%agrid(i,j,2)/deg2rad, 'lon = ', Atm%gridstruct%agrid(i,j,1)/deg2rad
   !      endif
      enddo
   enddo
   call pmaxmn('ZS_diff (m)', wk, is, ie, js, je, 1, 1._kind_real, Atm%gridstruct%area_64, Atm%domain)

   !if (.not.Atm%gridstruct%bounded_domain) then
       !call prt_gb_nh_sh('DATA_IC Z500', is,ie, js,je, z500, Atm%gridstruct%area_64(is:ie,js:je), Atm%gridstruct%agrid_64(is:ie,js:je,2))
       !if ( .not. Atm%flagstruct%hydrostatic )  &
       !call prt_height('fv3_IC Z500', is,ie, js,je, 3, npz, 500.E2, Atm%phis, Atm%delz, Atm%peln,   &
       !                Atm%gridstruct%area_64(is:ie,js:je), Atm%gridstruct%agrid_64(is:ie,js:je,2))
   !endif

   do j=js,je
      do i=is,ie
         wk(i,j) = Atm%ps(i,j) - psc(i,j)
      enddo
   enddo
   call pmaxmn('PS_diff (mb)', wk, is, ie, js, je, 1, 0.01_kind_real, Atm%gridstruct%area_64, Atm%domain)

   if (mpp_pe()==mpp_root_pe()) write(*,*) 'done remap_scalar'

  end subroutine remap_scalar

! --------------------------------------------------------------------------------------------------

  subroutine pmaxmn(qname, q, is, ie, js, je, km, fac, area, domain)
   character(len=*), intent(in)::  qname
   integer, intent(in):: is, ie, js, je
   integer, intent(in):: km
   real(kind=kind_real), intent(in)::    q(is:ie, js:je, km)
   real(kind=kind_real), intent(in)::    fac
   real(kind=R_GRID), intent(IN)::  area(is-3:ie+3, js-3:je+3)
   type(domain2d), intent(INOUT) :: domain
!---local variables
   real(kind=kind_real) qmin, qmax, gmean
   integer i,j,k

   qmin = q(is,js,1)
   qmax = qmin
   gmean = 0.

   do k=1,km
   do j=js,je
      do i=is,ie
         if( q(i,j,k) < qmin ) then
             qmin = q(i,j,k)
         elseif( q(i,j,k) > qmax ) then
             qmax = q(i,j,k)
         endif
       enddo
   enddo
   enddo

   !call mp_reduce_min(qmin)
   !call mp_reduce_max(qmax)

   gmean = g_sum(domain, q(is,js,km), is, ie, js, je, 3, area, 1, reproduce=.true.)
   if(mpp_pe()==mpp_root_pe()) write(6,*) qname, qmax*fac, qmin*fac, gmean*fac

end subroutine pmaxmn

! --------------------------------------------------------------------------------------------------

subroutine p_maxmin(qname, q, is, ie, js, je, km, fac)
   character(len=*), intent(in)::  qname
   integer, intent(in):: is, ie, js, je, km
   real(kind=kind_real), intent(in)::    q(is:ie, js:je, km)
   real(kind=kind_real), intent(in)::    fac
   real(kind=kind_real) qmin, qmax
   integer i,j,k

   qmin = q(is,js,1)
   qmax = qmin
   do k=1,km
   do j=js,je
      do i=is,ie
         if( q(i,j,k) < qmin ) then
             qmin = q(i,j,k)
         elseif( q(i,j,k) > qmax ) then
             qmax = q(i,j,k)
         endif
       enddo
   enddo
   enddo
   !call mp_reduce_min(qmin)
   !call mp_reduce_max(qmax)
   if(mpp_pe()==mpp_root_pe()) write(6,*) qname, qmax*fac, qmin*fac

end subroutine p_maxmin

! --------------------------------------------------------------------------------------------------

subroutine mp_auto_conversion(ql, qr, qi, qs)
   real(kind=kind_real), intent(inout):: ql, qr, qi, qs
   real(kind=kind_real), parameter:: qi0_max = 2.0e-3
   real(kind=kind_real), parameter:: ql0_max = 2.5e-3

  ! Convert excess cloud water into rain:
    if ( ql > ql0_max ) then
         qr = ql - ql0_max
         ql = ql0_max
    endif
  ! Convert excess cloud ice into snow:
    if ( qi > qi0_max ) then
         qs = qi - qi0_max
         qi = qi0_max
    endif

   end subroutine mp_auto_conversion

   subroutine fillq(im, km, nq, q, dp)
      integer,  intent(in):: im            !< No. of longitudes
      integer,  intent(in):: km            !< No. of levels
      integer,  intent(in):: nq            !< Total number of tracers
      real(kind=kind_real) , intent(in)::  dp(im,km)       !< pressure thickness
      real(kind=kind_real) , intent(inout) :: q(im,km,nq)  !< tracer mixing ratio
   ! !LOCAL VARIABLES:
      integer i, k, ic, k1

      do ic=1,nq
   ! Bottom up:
         do k=km,2,-1
            k1 = k-1
            do i=1,im
              if( q(i,k,ic) < 0. ) then
                  q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                  q(i,k ,ic) = 0.
              endif
            enddo
         enddo
   ! Top down:
         do k=1,km-1
            k1 = k+1
            do i=1,im
               if( q(i,k,ic) < 0. ) then
                   q(i,k1,ic) = q(i,k1,ic) + q(i,k,ic)*dp(i,k)/dp(i,k1)
                   q(i,k ,ic) = 0.
               endif
            enddo
         enddo

      enddo

    end subroutine fillq

! --------------------------------------------------------------------------------------------------

    subroutine fillz(im, km, nq, q, dp)
      integer,  intent(in):: im                !< No. of longitudes
      integer,  intent(in):: km                !< No. of levels
      integer,  intent(in):: nq                !< Total number of tracers
      real(kind=kind_real) , intent(in)::  dp(im,km)           !< pressure thickness
      real(kind=kind_real) , intent(inout) :: q(im,km,nq)      !< tracer mixing ratio
   ! LOCAL VARIABLES:
      logical:: zfix(im)
      real(kind=kind_real) ::  dm(km)
      integer i, k, ic, k1
      real(kind=kind_real)  qup, qly, dup, dq, sum0, sum1, fac

      do ic=1,nq
   ! Top layer
         do i=1,im
            if( q(i,1,ic) < 0. ) then
                q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
                q(i,1,ic) = 0.
             endif
         enddo

   ! Interior
         zfix(:) = .false.
         do k=2,km-1
            do i=1,im
            if( q(i,k,ic) < 0. ) then
                zfix(i) = .true.
                if ( q(i,k-1,ic) > 0. ) then
   ! Borrow from above
                   dq = min ( q(i,k-1,ic)*dp(i,k-1), -q(i,k,ic)*dp(i,k) )
                   q(i,k-1,ic) = q(i,k-1,ic) - dq/dp(i,k-1)
                   q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
                endif
                if ( q(i,k,ic)<0.0 .and. q(i,k+1,ic)>0. ) then
   ! Borrow from below:
                   dq = min ( q(i,k+1,ic)*dp(i,k+1), -q(i,k,ic)*dp(i,k) )
                   q(i,k+1,ic) = q(i,k+1,ic) - dq/dp(i,k+1)
                   q(i,k  ,ic) = q(i,k  ,ic) + dq/dp(i,k  )
                endif
             endif
            enddo
         enddo

   ! Bottom layer
         k = km
         do i=1,im
            if( q(i,k,ic)<0. .and. q(i,k-1,ic)>0.) then
                zfix(i) = .true.
   ! Borrow from above
                qup =  q(i,k-1,ic)*dp(i,k-1)
                qly = -q(i,k  ,ic)*dp(i,k  )
                dup =  min(qly, qup)
                q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
                q(i,k,  ic) = q(i,k,  ic) + dup/dp(i,k  )
             endif
         enddo

   ! Perform final check and non-local fix if needed
         do i=1,im
            if ( zfix(i) ) then

              sum0 = 0.
              do k=2,km
                 dm(k) = q(i,k,ic)*dp(i,k)
                 sum0 = sum0 + dm(k)
              enddo

              if ( sum0 > 0. ) then
                sum1 = 0.
                do k=2,km
                   sum1 = sum1 + max(0., dm(k))
                enddo
                fac = sum0 / sum1
                do k=2,km
                   q(i,k,ic) = max(0., fac*dm(k)/dp(i,k))
                enddo
              endif

            endif
         enddo

      enddo
    end subroutine fillz

! --------------------------------------------------------------------------------------------------

    subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

      ! IV = 0: constituents
      ! IV = 1: potential temp
      ! IV =-1: winds

      ! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

      ! pe1: pressure at layer edges (from model top to bottom surface)
      !      in the original vertical coordinate
      ! pe2: pressure at layer edges (from model top to bottom surface)
      !      in the new vertical coordinate

       integer, intent(in):: i1, i2, km, kn, kord, iv
       real(kind=kind_real), intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1) !< pe1: pressure at layer edges from model top to bottom
                                                            !!      surface in the ORIGINAL vertical coordinate
                                                            !< pe2: pressure at layer edges from model top to bottom
                                                            !!      surface in the NEW vertical coordinate
      ! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
       real(kind=kind_real), intent(in )::  q1(i1:i2,km)
       real(kind=kind_real), intent(out)::  q2(i1:i2,kn)
       real(kind=kind_real), intent(IN) :: ptop
      ! local
            real(kind=kind_real)  qs(i1:i2)
            real(kind=kind_real) dp1(i1:i2,km)
            real(kind=kind_real) a4(4,i1:i2,km)
            integer i, k, l
            integer k0, k1
            real(kind=kind_real) pl, pr, tt, delp, qsum, dpsum, esl

            do k=1,km
               do i=i1,i2
                   dp1(i,k) = pe1(i,k+1) - pe1(i,k)
                  a4(1,i,k) = q1(i,k)
               enddo
            enddo

            if ( kord >7 ) then
                 call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
            else
                 call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
            endif

      !------------------------------------
      ! Lowest layer: constant distribution
      !------------------------------------

            do 5555 i=i1,i2
               k0 = 1
            do 555 k=1,kn

               if(pe2(i,k) .le. pe1(i,1)) then
      ! above old ptop
                  q2(i,k) = q1(i,1)
               elseif(pe2(i,k) .ge. pe1(i,km+1)) then
      ! Entire grid below old ps
                  q2(i,k) = q1(i,km)
               else

               do 45 L=k0,km
      ! locate the top edge at pe2(i,k)
               if( pe2(i,k) .ge. pe1(i,L) .and.        &
                   pe2(i,k) .le. pe1(i,L+1)    ) then
                   k0 = L
                   PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
                   if(pe2(i,k+1) .le. pe1(i,L+1)) then

      ! entire new grid is within the original grid
                     PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
                     TT = r3*(PR*(PR+PL)+PL**2)
                     q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                             - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
                    goto 555
                   else
      ! Fractional area...
                    delp = pe1(i,L+1) - pe2(i,k)
                    TT   = r3*(1.+PL*(1.+PL))
                    qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                           a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
                    dpsum = delp
                    k1 = L + 1
                   goto 111
                   endif
               endif
      45       continue

      111      continue
               do 55 L=k1,km
               if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

      ! Whole layer..

                  qsum  =  qsum + dp1(i,L)*q1(i,L)
                  dpsum = dpsum + dp1(i,L)
               else
                 delp = pe2(i,k+1)-pe1(i,L)
                 esl  = delp / dp1(i,L)
                 qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                       (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23*esl)) )
                dpsum = dpsum + delp
                 k0 = L
                 goto 123
               endif
      55       continue
              delp = pe2(i,k+1) - pe1(i,km+1)
              if(delp > 0.) then
      ! Extended below old ps
                 qsum = qsum + delp * q1(i,km)
                dpsum = dpsum + delp
              endif
      123     q2(i,k) = qsum / dpsum
            endif
      555   continue
      5555  continue

       end subroutine mappm

! --------------------------------------------------------------------------------------------------

       subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
         ! Optimized vertical profile reconstruction:
         ! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
          integer, intent(in):: i1, i2
          integer, intent(in):: km      !< vertical dimension
          integer, intent(in):: iv      !< iv =-1: winds
                                        !< iv = 0: positive definite scalars
                                        !< iv = 1: others
          integer, intent(in):: kord
          real(kind=kind_real), intent(in)   ::   qs(i1:i2)
          real(kind=kind_real), intent(in)   :: delp(i1:i2,km)     !< layer pressure thickness
          real(kind=kind_real), intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
         !-----------------------------------------------------------------------
          logical, dimension(i1:i2,km):: extm, ext5, ext6
          real(kind=kind_real)  gam(i1:i2,km)
          real(kind=kind_real)    q(i1:i2,km+1)
          real(kind=kind_real)   d4(i1:i2)
          real(kind=kind_real)   bet, a_bot, grat
          real(kind=kind_real)   pmp_1, lac_1, pmp_2, lac_2, x0, x1
          integer i, k, im

          if ( iv .eq. -2 ) then
               do i=i1,i2
                  gam(i,2) = 0.5
                    q(i,1) = 1.5*a4(1,i,1)
               enddo
               do k=2,km-1
                  do i=i1, i2
                           grat = delp(i,k-1) / delp(i,k)
                            bet =  2. + grat + grat - gam(i,k)
                         q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
                     gam(i,k+1) = grat / bet
                  enddo
               enddo
               do i=i1,i2
                     grat = delp(i,km-1) / delp(i,km)
                  q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                            (2. + grat + grat - gam(i,km))
                  q(i,km+1) = qs(i)
               enddo
               do k=km-1,1,-1
                 do i=i1,i2
                    q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
                 enddo
               enddo
          else
           do i=i1,i2
                  grat = delp(i,2) / delp(i,1)   ! grid ratio
                   bet = grat*(grat+0.5)
                q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
              gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
           enddo

           do k=2,km
              do i=i1,i2
                    d4(i) = delp(i,k-1) / delp(i,k)
                      bet =  2. + d4(i) + d4(i) - gam(i,k-1)
                   q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
                 gam(i,k) = d4(i) / bet
              enddo
           enddo

           do i=i1,i2
                  a_bot = 1. + d4(i)*(d4(i)+1.5)
              q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
                        / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
           enddo

           do k=km,1,-1
              do i=i1,i2
                 q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
              enddo
           enddo
          endif
         !----- Perfectly linear scheme --------------------------------
          if ( abs(kord) > 16 ) then
           do k=1,km
              do i=i1,i2
                 a4(2,i,k) = q(i,k  )
                 a4(3,i,k) = q(i,k+1)
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
           enddo
           return
          endif
         !----- Perfectly linear scheme --------------------------------

         !------------------
         ! Apply constraints
         !------------------
           im = i2 - i1 + 1

         ! Apply *large-scale* constraints
           do i=i1,i2
              q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
              q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
           enddo

           do k=2,km
              do i=i1,i2
                 gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
              enddo
           enddo

         ! Interior:
           do k=3,km-1
              do i=i1,i2
                 if ( gam(i,k-1)*gam(i,k+1)>0. ) then
         ! Apply large-scale constraint to ALL fields if not local max/min
                      q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
                      q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
                 else
                   if ( gam(i,k-1) > 0. ) then
         ! There exists a local max
                        q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
                   else
         ! There exists a local min
                          q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
                        if ( iv==0 ) q(i,k) = max(0., q(i,k))
                   endif
                 endif
              enddo
           enddo

         ! Bottom:
           do i=i1,i2
              q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
              q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
           enddo

           do k=1,km
              do i=i1,i2
                 a4(2,i,k) = q(i,k  )
                 a4(3,i,k) = q(i,k+1)
              enddo
           enddo

           do k=1,km
              if ( k==1 .or. k==km ) then
                do i=i1,i2
                   extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
                enddo
              else
                do i=i1,i2
                   extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
                enddo
              endif
              if ( abs(kord) > 9 ) then
                do i=i1,i2
                   x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
                   x1 = abs(a4(2,i,k)-a4(3,i,k))
                   a4(4,i,k) = 3.*x0
                   ext5(i,k) = abs(x0) > x1
                   ext6(i,k) = abs(a4(4,i,k)) > x1
                enddo
              endif
           enddo

         !---------------------------
         ! Apply subgrid constraints:
         !---------------------------
         ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
         ! Top 2 and bottom 2 layers always use monotonic mapping

           if ( iv==0 ) then
              do i=i1,i2
                 a4(2,i,1) = max(0., a4(2,i,1))
              enddo
           elseif ( iv==-1 ) then
               do i=i1,i2
                  if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
               enddo
           elseif ( iv==2 ) then
              do i=i1,i2
                 a4(2,i,1) = a4(1,i,1)
                 a4(3,i,1) = a4(1,i,1)
                 a4(4,i,1) = 0.
              enddo
           endif

           if ( iv/=2 ) then
              do i=i1,i2
                 a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
              enddo
              call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
           endif

         ! k=2
            do i=i1,i2
               a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
            enddo
            call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

         !-------------------------------------
         ! Huynh's 2nd constraint for interior:
         !-------------------------------------
           do k=3,km-2
              if ( abs(kord)<9 ) then
                do i=i1,i2
         ! Left  edges
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
         ! Right edges
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                                  max(a4(1,i,k), pmp_2, lac_2) )

                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo

              elseif ( abs(kord)==9 ) then
                do i=i1,i2
                   if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
         ! grid-scale 2-delta-z wave detected
                        a4(2,i,k) = a4(1,i,k)
                        a4(3,i,k) = a4(1,i,k)
                        a4(4,i,k) = 0.
                   else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
         ! grid-scale 2-delta-z wave detected
                        a4(2,i,k) = a4(1,i,k)
                        a4(3,i,k) = a4(1,i,k)
                        a4(4,i,k) = 0.
                   else
                     a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
         ! Check within the smooth region if subgrid profile is non-monotonic
                     if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                           pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                           lac_1 = pmp_1 + 1.5*gam(i,k+2)
                       a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                      max(a4(1,i,k), pmp_1, lac_1) )
                           pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                           lac_2 = pmp_2 - 1.5*gam(i,k-1)
                       a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                      max(a4(1,i,k), pmp_2, lac_2) )
                       a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
                     endif
                   endif
                enddo
              elseif ( abs(kord)==10 ) then
                do i=i1,i2
                   if( ext5(i,k) ) then
                       if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                            a4(2,i,k) = a4(1,i,k)
                            a4(3,i,k) = a4(1,i,k)
                       elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                            pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                            lac_1 = pmp_1 + 1.5*gam(i,k+2)
                            a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                           max(a4(1,i,k), pmp_1, lac_1) )
                            pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                            lac_2 = pmp_2 - 1.5*gam(i,k-1)
                            a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                           max(a4(1,i,k), pmp_2, lac_2) )
                       endif
                   elseif( ext6(i,k) ) then
                       if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                           pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                           lac_1 = pmp_1 + 1.5*gam(i,k+2)
                           a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                          max(a4(1,i,k), pmp_1, lac_1) )
                           pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                           lac_2 = pmp_2 - 1.5*gam(i,k-1)
                           a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                          max(a4(1,i,k), pmp_2, lac_2) )
                       endif
                   endif
                enddo
                do i=i1,i2
                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo
              elseif ( abs(kord)==12 ) then
                do i=i1,i2
                   if( extm(i,k) ) then
         ! grid-scale 2-delta-z wave detected
                       a4(2,i,k) = a4(1,i,k)
                       a4(3,i,k) = a4(1,i,k)
                       a4(4,i,k) = 0.
                   else        ! not a local extremum
                     a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
         ! Check within the smooth region if subgrid profile is non-monotonic
                     if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                           pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                           lac_1 = pmp_1 + 1.5*gam(i,k+2)
                       a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                      max(a4(1,i,k), pmp_1, lac_1) )
                           pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                           lac_2 = pmp_2 - 1.5*gam(i,k-1)
                       a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                      max(a4(1,i,k), pmp_2, lac_2) )
                       a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
                     endif
                   endif
                enddo
              elseif ( abs(kord)==13 ) then
                do i=i1,i2
                   if( ext6(i,k) ) then
                      if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
         ! grid-scale 2-delta-z wave detected
                          a4(2,i,k) = a4(1,i,k)
                          a4(3,i,k) = a4(1,i,k)
                      endif
                   endif
                enddo
                do i=i1,i2
                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo
              elseif ( abs(kord)==14 ) then

                do i=i1,i2
                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo

              elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
                do i=i1,i2
                   if ( ext5(i,k) ) then  ! c90_mp122
                      if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
         ! grid-scale 2-delta-z wave detected
                           a4(2,i,k) = a4(1,i,k)
                           a4(3,i,k) = a4(1,i,k)
                      endif
                   elseif( ext6(i,k) ) then
         ! Check within the smooth region if subgrid profile is non-monotonic
                           pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                           lac_1 = pmp_1 + 1.5*gam(i,k+2)
                       a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                      max(a4(1,i,k), pmp_1, lac_1) )
                           pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                           lac_2 = pmp_2 - 1.5*gam(i,k-1)
                       a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                      max(a4(1,i,k), pmp_2, lac_2) )
                   endif
                enddo
                do i=i1,i2
                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo
              elseif ( abs(kord)==16 ) then
                do i=i1,i2
                   if( ext5(i,k) ) then
                      if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                          a4(2,i,k) = a4(1,i,k)
                          a4(3,i,k) = a4(1,i,k)
                      elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                          ! Left  edges
                          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                          lac_1 = pmp_1 + 1.5*gam(i,k+2)
                          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                              max(a4(1,i,k), pmp_1, lac_1) )
                          ! Right edges
                          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                          lac_2 = pmp_2 - 1.5*gam(i,k-1)
                          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                              max(a4(1,i,k), pmp_2, lac_2) )
                      endif
                   endif
                enddo
                do i=i1,i2
                   a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                enddo
              else      ! kord = 11
                do i=i1,i2
                  if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
         ! Noisy region:
                       a4(2,i,k) = a4(1,i,k)
                       a4(3,i,k) = a4(1,i,k)
                       a4(4,i,k) = 0.
                  else
                       a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                  endif
                enddo
              endif

         ! Additional constraint to ensure positivity
              if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

           enddo      ! k-loop

         !----------------------------------
         ! Bottom layer subgrid constraints:
         !----------------------------------
           if ( iv==0 ) then
              do i=i1,i2
                 a4(3,i,km) = max(0., a4(3,i,km))
              enddo
           elseif ( iv .eq. -1 ) then
               do i=i1,i2
                  if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
               enddo
           endif

           do k=km-1,km
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
              if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
              if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
           enddo

          end subroutine cs_profile


          subroutine cs_limiters(im, extm, a4, iv)
          integer, intent(in) :: im
          integer, intent(in) :: iv
          logical, intent(in) :: extm(im)
          real(kind=kind_real) , intent(inout) :: a4(4,im)   !< PPM array
         ! LOCAL VARIABLES:
          real(kind=kind_real)  da1, da2, a6da
          integer i

          if ( iv==0 ) then
         ! Positive definite constraint
             do i=1,im
             if( a4(1,i)<=0.) then
                 a4(2,i) = a4(1,i)
                 a4(3,i) = a4(1,i)
                 a4(4,i) = 0.
             else
               if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
                  if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
         ! local minimum is negative
                      if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                          a4(3,i) = a4(1,i)
                          a4(2,i) = a4(1,i)
                          a4(4,i) = 0.
                      elseif( a4(3,i) > a4(2,i) ) then
                          a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                          a4(3,i) = a4(2,i) - a4(4,i)
                      else
                          a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                          a4(2,i) = a4(3,i) - a4(4,i)
                      endif
                  endif
               endif
             endif
             enddo
          elseif ( iv==1 ) then
             do i=1,im
               if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
                  a4(2,i) = a4(1,i)
                  a4(3,i) = a4(1,i)
                  a4(4,i) = 0.
               else
                  da1  = a4(3,i) - a4(2,i)
                  da2  = da1**2
                  a6da = a4(4,i)*da1
                  if(a6da < -da2) then
                     a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                     a4(3,i) = a4(2,i) - a4(4,i)
                  elseif(a6da > da2) then
                     a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                     a4(2,i) = a4(3,i) - a4(4,i)
                  endif
               endif
             enddo
          else
         ! Standard PPM constraint
             do i=1,im
               if( extm(i) ) then
                  a4(2,i) = a4(1,i)
                  a4(3,i) = a4(1,i)
                  a4(4,i) = 0.
               else
                  da1  = a4(3,i) - a4(2,i)
                  da2  = da1**2
                  a6da = a4(4,i)*da1
                  if(a6da < -da2) then
                     a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                     a4(3,i) = a4(2,i) - a4(4,i)
                  elseif(a6da > da2) then
                     a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                     a4(2,i) = a4(3,i) - a4(4,i)
                  endif
               endif
             enddo
          endif
          end subroutine cs_limiters



          subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

         ! !INPUT PARAMETERS:
          integer, intent(in):: iv      !< iv =-1: winds
                                        !! iv = 0: positive definite scalars
                                        !! iv = 1: others
                                        !! iv = 2: temp (if remap_t) and w (iv=-2)
          integer, intent(in):: i1      !< Starting longitude
          integer, intent(in):: i2      !< Finishing longitude
          integer, intent(in):: km      !< vertical dimension
          integer, intent(in):: kord    !< Order (or more accurately method no.):
                                        !!
          real(kind=kind_real) , intent(in):: delp(i1:i2,km)     !< layer pressure thickness

         ! !INPUT/OUTPUT PARAMETERS:
          real(kind=kind_real) , intent(inout):: a4(4,i1:i2,km)  !< Interpolated values

         ! DESCRIPTION:
         !
         !   Perform the piecewise parabolic reconstruction
         !
         ! !REVISION HISTORY:
         ! S.-J. Lin   revised at GFDL 2007
         !-----------------------------------------------------------------------
         ! local arrays:
               real(kind=kind_real)    dc(i1:i2,km)
               real(kind=kind_real)    h2(i1:i2,km)
               real(kind=kind_real)  delq(i1:i2,km)
               real(kind=kind_real)   df2(i1:i2,km)
               real(kind=kind_real)    d4(i1:i2,km)

         ! local scalars:
               integer i, k, km1, lmt, it
               real(kind=kind_real)  fac
               real(kind=kind_real)  a1, a2, c1, c2, c3, d1, d2
               real(kind=kind_real)  qm, dq, lac, qmp, pmp

               km1 = km - 1
                it = i2 - i1 + 1

               do k=2,km
                  do i=i1,i2
                     delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
                       d4(i,k  ) = delp(i,k-1) + delp(i,k)
                  enddo
               enddo

               do k=2,km1
                  do i=i1,i2
                          c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                          c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
                     df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                             (d4(i,k)+delp(i,k+1))
                     dc(i,k) = sign( min(abs(df2(i,k)),              &
                                     max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                           a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
                  enddo
               enddo

         !-----------------------------------------------------------
         ! 4th order interpolation of the provisional cell edge value
         !-----------------------------------------------------------

               do k=3,km1
                  do i=i1,i2
                     c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
                     a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
                     a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
                     a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                               ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                 delp(i,k-1)*a1*dc(i,k  ) )
                  enddo
               enddo

         !     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

         ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
         ! Top
               do i=i1,i2
                  d1 = delp(i,1)
                  d2 = delp(i,2)
                  qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
                  dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
                  c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
                  c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
                  a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         ! Top edge:
         !-------------------------------------------------------
                  a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
         !-------------------------------------------------------
         !        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
         !-------------------------------------------------------
         ! No over- and undershoot condition
                  a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
                  a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
                  dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
               enddo

         ! Enforce monotonicity  within the top layer

               if( iv==0 ) then
                  do i=i1,i2
                     a4(2,i,1) = max(0., a4(2,i,1))
                     a4(2,i,2) = max(0., a4(2,i,2))
                  enddo
               elseif( iv==-1 ) then
                  do i=i1,i2
                     if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
                  enddo
               elseif( abs(iv)==2 ) then
                  do i=i1,i2
                     a4(2,i,1) = a4(1,i,1)
                     a4(3,i,1) = a4(1,i,1)
                  enddo
               endif

         ! Bottom
         ! Area preserving cubic with 2nd deriv. = 0 at the surface
               do i=i1,i2
                  d1 = delp(i,km)
                  d2 = delp(i,km1)
                  qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
                  dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
                  c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
                  c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
                  a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         ! Bottom edge:
         !-----------------------------------------------------
                  a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         !        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
         !-----------------------------------------------------
         !        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
         ! No over- and under-shoot condition
                  a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
                  a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
                  dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
               enddo


         ! Enforce constraint on the "slope" at the surface

               if( iv==0 ) then
                   do i=i1,i2
                      a4(2,i,km) = max(0.,a4(2,i,km))
                      a4(3,i,km) = max(0.,a4(3,i,km))
                   enddo
               elseif( iv<0 ) then
                   do i=i1,i2
                      if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
                   enddo
               endif

            do k=1,km1
               do i=i1,i2
                  a4(3,i,k) = a4(2,i,k+1)
               enddo
            enddo

         !-----------------------------------------------------------
         ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
         !-----------------------------------------------------------
         ! Top 2 and bottom 2 layers always use monotonic mapping
               do k=1,2
                  do i=i1,i2
                     a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                  enddo
                  call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
               enddo

               if(kord >= 7) then
         !-----------------------
         ! Huynh's 2nd constraint
         !-----------------------
               do k=2,km1
                  do i=i1,i2
         ! Method#1
         !           h2(i,k) = delq(i,k) - delq(i,k-1)
         ! Method#2 - better
                     h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                              / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                              * delp(i,k)**2
         ! Method#3
         !!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
                  enddo
               enddo

               fac = 1.5           ! original quasi-monotone

               do k=3,km-2
                 do i=i1,i2
         ! Right edges
         !        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
         !        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
         !
                  pmp   = 2.*dc(i,k)
                  qmp   = a4(1,i,k) + pmp
                  lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                                 max(a4(1,i,k), qmp, lac) )
         ! Left  edges
         !        qmp   = a4(1,i,k) - 2.0*delq(i,k)
         !        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
         !
                  qmp   = a4(1,i,k) - pmp
                  lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
                  a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                              max(a4(1,i,k), qmp, lac))
         !-------------
         ! Recompute A6
         !-------------
                  a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                 enddo
         ! Additional constraint to ensure positivity when kord=7
                  if (iv == 0 .and. kord >= 6 )                      &
                      call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
               enddo

               else

                  lmt = kord - 3
                  lmt = max(0, lmt)
                  if (iv == 0) lmt = min(2, lmt)

                  do k=3,km-2
                     if( kord /= 4) then
                       do i=i1,i2
                          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                       enddo
                     endif
                     if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
                  enddo
               endif

               do k=km1,km
                  do i=i1,i2
                     a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
                  enddo
                  call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
               enddo

          end subroutine ppm_profile


          subroutine ppm_limiters(dm, a4, itot, lmt)

         ! INPUT PARAMETERS:
               real(kind=kind_real) , intent(in):: dm(*)     !< Linear slope
               integer, intent(in) :: itot      !< Total Longitudes
               integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
                                                !< (Lin) 2: Positive definite constraint
                                                !< 3: do nothing (return immediately)
         ! INPUT/OUTPUT PARAMETERS:
               real(kind=kind_real) , intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
         ! LOCAL VARIABLES:
               real(kind=kind_real)  qmp
               real(kind=kind_real)  da1, da2, a6da
               real(kind=kind_real)  fmin
               integer i

         ! Developer: S.-J. Lin

               if ( lmt == 3 ) return

               if(lmt == 0) then
         ! Standard PPM constraint
               do i=1,itot
               if(dm(i) == 0.) then
                  a4(2,i) = a4(1,i)
                  a4(3,i) = a4(1,i)
                  a4(4,i) = 0.
               else
                  da1  = a4(3,i) - a4(2,i)
                  da2  = da1**2
                  a6da = a4(4,i)*da1
                  if(a6da < -da2) then
                     a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                     a4(3,i) = a4(2,i) - a4(4,i)
                  elseif(a6da > da2) then
                     a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                     a4(2,i) = a4(3,i) - a4(4,i)
                  endif
               endif
               enddo

               elseif (lmt == 1) then

         ! Improved full monotonicity constraint (Lin 2004)
         ! Note: no need to provide first guess of A6 <-- a4(4,i)
               do i=1, itot
                    qmp = 2.*dm(i)
                  a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
                  a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
                  a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
               enddo

               elseif (lmt == 2) then

         ! Positive definite constraint
               do i=1,itot
               if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
               fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
                  if( fmin < 0. ) then
                  if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
                     a4(3,i) = a4(1,i)
                     a4(2,i) = a4(1,i)
                     a4(4,i) = 0.
                  elseif(a4(3,i) > a4(2,i)) then
                     a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                     a4(3,i) = a4(2,i) - a4(4,i)
                  else
                     a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                     a4(2,i) = a4(3,i) - a4(4,i)
                  endif
                  endif
               endif
               enddo

               endif

          end subroutine ppm_limiters


! --------------------------------------------------------------------------------------------------

end program cold_to_warm
