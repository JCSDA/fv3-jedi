program cold_to_warm

! fckit uses
use fckit_module
use fckit_mpi_module

! fms uses
use constants_mod,              only: grav
use fms_io_mod,                 only: nullify_domain
use fms_mod,                    only: fms_init
use mpp_mod,                    only: mpp_exit, mpp_pe, mpp_npes, mpp_error, FATAL, NOTE, &
                                      mpp_root_pe
use mpp_domains_mod,            only: domain2D, mpp_deallocate_domain, mpp_define_layout, &
                                      mpp_define_mosaic, mpp_define_io_domain, mpp_domains_exit, &
                                      mpp_domains_set_stack_size
use field_manager_mod,          only: fm_string_len, field_manager_init, MODEL_ATMOS
use tracer_manager_mod,         only: get_number_tracers, get_tracer_names, get_tracer_index, &
                                      NO_TRACER, set_tracer_profile
use fms_io_mod,                 only: restart_file_type, register_restart_field
use fms_io_mod,                 only: free_restart_type, restore_state, save_restart
use fms_io_mod,                 only: set_domain, nullify_domain
use mpp_domains_mod,            only: east, north, center

! fv3 uses
use fv_arrays_mod,              only: fv_atmos_type, deallocate_fv_atmos_type, R_GRID
use fv_grid_utils_mod,          only: mid_pt_sphere, get_unit_vect2, get_latlon_vector, inner_prod
use external_ic_mod,            only: remap_scalar, remap_dwinds
use test_cases_mod,             only: checker_tracers

! fv3jedi uses ***this code should not make use of state or increment type***
use fv3jedi_geom_mod,           only: initialize_fms => initialize, fv3jedi_geom
use fv_prec_mod,                only: kind_fv3
use fv_init_mod,                only: fv_init
use fv3jedi_fmsnamelist_mod,    only: fv3jedi_fmsnamelist

! Nothing implicit
implicit none

type :: state_cold_type
  integer :: isc, iec, jsc, jec, npz
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: u_w_cold, v_w_cold, u_s_cold, v_s_cold
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: ud_cold
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: vd_cold
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: t_cold
  real(kind=kind_fv3), dimension(:,:,:,:), allocatable :: q_cold
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: zh_cold
  real(kind=kind_fv3), dimension(:,:,:),   allocatable :: w_cold
  real(kind=kind_fv3), dimension(:,:),     allocatable :: ps_cold
  real(kind=kind_fv3), dimension(:,:),     allocatable :: orog_filt
end type state_cold_type

! Local variables
type(fv3jedi_geom) :: geom
integer :: arg_count, gtile, ptile=1, nlevs, ntracers, ntprog, ierr, isc, iec, jsc, jec, npz
character(len=256) :: yaml_file
type(fckit_configuration) :: config, config_fms, config_geom, config_input, config_remap, &
                             config_output, config_test
type(fckit_mpi_comm)      :: comm
type(fv_atmos_type), allocatable :: Atm(:)
logical, allocatable :: grids_on_this_pe(:)
type(fv3jedi_fmsnamelist) :: fmsnamelist

type(state_cold_type)       :: state_cold

real(kind=kind_fv3), dimension(:), allocatable :: phis_flat

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
call geom%create(config_geom, comm, nlevs)

! Create fv3 object
! -----------------
call fmsnamelist%replace_namelist(config_remap)
call fv_init(Atm, 300.0_kind_fv3, grids_on_this_pe, ptile, gtile, .true.)
call fmsnamelist%revert_namelist

Atm(1)%flagstruct%nggps_ic = .true.
Atm(1)%ak = real(geom%ak,kind_fv3)
Atm(1)%bk = real(geom%bk,kind_fv3)
Atm%ptop = real(geom%ak(1),kind_fv3)

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
allocate(state_cold%  ps_cold(isc:iec  , jsc:jec                  ))
allocate(state_cold%orog_filt(isc:iec  , jsc:jec                  ))

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
real(kind=kind_fv3), allocatable :: ak(:), bk(:)
real(kind=kind_fv3) :: wt, qt, m_fac


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
  ak = 0.0_kind_fv3
  bk = 0.0_kind_fv3
  ak(itoa:levp+1) = Atm(1)%ak(1:npz+1)
  bk(itoa:levp+1) = Atm(1)%bk(1:npz+1)
  ak(1) = max(1.e-9_kind_fv3, ak(1))

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
            qt = wt*(1.0_kind_fv3 + Atm(1)%q(i,j,k,liq_wat) + Atm(1)%q(i,j,k,ice_wat) + &
                                    Atm(1)%q(i,j,k,rainwat) + Atm(1)%q(i,j,k,snowwat) + &
                                    Atm(1)%q(i,j,k,graupel))
          else
            qt = wt*(1.0_kind_fv3 + sum(Atm(1)%q(i,j,k,2:Atm(1)%flagstruct%nwat)))
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
            qt = wt*(1.0_kind_fv3 + Atm(1)%q(i,j,k,liq_wat) + Atm(1)%q(i,j,k,ice_wat) + &
                                    Atm(1)%q(i,j,k,rainwat) + Atm(1)%q(i,j,k,snowwat) + &
                                    Atm(1)%q(i,j,k,graupel))
          else
             qt = wt*(1.0_kind_fv3 + sum(Atm(1)%q(i,j,k,2:Atm(1)%flagstruct%nwat)))
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
                       9.0_kind_fv3, 9.0_kind_fv3)
endif


end subroutine remap


! --------------------------------------------------------------------------------------------------


subroutine testing(comm, geom, conf, Atm)

type(fckit_mpi_comm),      intent(in) :: comm
type(fv3jedi_geom),        intent(in) :: geom
type(fckit_configuration), intent(in) :: conf
type(fv_atmos_type),       intent(in) :: Atm(:)

real(kind=kind_fv3) :: rms_u_c, rms_v_c, rms_t_c, rms_d_c, rms_p_c, tol
real(kind=kind_fv3) :: rms_u_r, rms_v_r, rms_t_r, rms_d_r, rms_p_r
integer :: isc, iec, jsc, jec, npz
real(kind=kind_fv3) :: tmp(5), gs3, gs3g

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
gs3 = real((geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)*geom%npz, kind_fv3)
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


end program cold_to_warm
