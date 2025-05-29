program ua_to_ud

! fckit uses
use fckit_module
use fckit_mpi_module

! fms uses
use fms_io_mod,                 only: register_restart_field, free_restart_type, restore_state, &
                                      save_restart, restart_file_type
use mpp_domains_mod,            only: east, north, center

! fv3jedi uses ***this code should not make use of state or increment type***
use fv3jedi_geom_mod,           only: initialize_fms => initialize, fv3jedi_geom
use fv3jedi_kinds_mod,          only: kind_real
use fv3jedi_fmsnamelist_mod,    only: fv3jedi_fmsnamelist
use wind_vt_mod,                only: a_to_d

! Nothing implicit
implicit none

! Local variables
type(fv3jedi_geom) :: geom
type(fckit_mpi_comm) :: comm
integer :: npx, npy, npz
character(len=256) :: yaml_file
type(fckit_configuration) :: config, config_fms, config_geom, config_input, config_output, &
                             config_test
type(fv3jedi_fmsnamelist) :: fmsnamelist
real(kind=kind_real), dimension(:,:,:),   allocatable :: u, v
real(kind=kind_real), dimension(:,:,:),   allocatable :: ua, va


! Initialize fckit and communicator
! ---------------------------------
call fckit_main%initialise()
comm = fckit_mpi_comm("world")

! Check if a YAML file was provided as an argument
! ------------------------------------------------
if (command_argument_count() < 1) then
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
call config%get_or_die("output",   config_output)

! Initialize fms
! --------------
call initialize_fms(config_fms, comm)

! Intialize fv3 geometry
! ---------------------
call fmsnamelist%replace_namelist(config_geom)
call geom%create(config_geom, comm, npx, npy, npz)

! Allocate state
! --------------
allocate( u(geom%isc:geom%iec  , geom%jsc:geom%jec+1, 1:geom%npz))
allocate( v(geom%isc:geom%iec+1, geom%jsc:geom%jec  , 1:geom%npz))
allocate(ua(geom%isc:geom%iec  , geom%jsc:geom%jec  , 1:geom%npz))
allocate(va(geom%isc:geom%iec  , geom%jsc:geom%jec  , 1:geom%npz))

! Read the restarts
! -----------------
call read_a(config_input, geom, ua, va)

! Convert the winds
! -----------------
call a_to_d(geom, ua, va, u, v)

! Write the restarts
! ------------------
call write_d(config_output, geom, u, v)

! Check output numbers
! --------------------
if (config%has("testing")) then

  ! Get testing config
  call config%get_or_die("testing", config_test)

  ! Check output
  call testing(comm, geom, config_test, u, v)

endif

! Delete geom
! -----------
call geom%delete()


! --------------------------------------------------------------------------------------------------


contains


! --------------------------------------------------------------------------------------------------


subroutine read_a(conf, geom, ua, va)

! Arguments
type(fckit_configuration), intent(in)    :: conf
type(fv3jedi_geom),        intent(in)    :: geom
real(kind=kind_real),      intent(inout) :: ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real),      intent(inout) :: va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

! Locals
character(len=1024) :: dpath, fcore
character(len=:), allocatable :: str
integer :: rc
type(restart_file_type) :: rst

! Read file path and names from config
! ------------------------------------
call conf%get_or_die("datapath", str)
dpath = str
deallocate(str)
call conf%get_or_die("filename_core", str)
fcore = str
deallocate(str)

! Register and read
! -----------------
rc = register_restart_field(rst, trim(fcore), 'ua', ua, domain=geom%domain, position=center  )
rc = register_restart_field(rst, trim(fcore), 'va', va, domain=geom%domain, position=center  )
call restore_state(rst, directory=trim(adjustl(dpath)))
call free_restart_type(rst)

end subroutine read_a


! --------------------------------------------------------------------------------------------------


subroutine write_d(conf, geom, u, v)

! Arguments
type(fckit_configuration), intent(in) :: conf
type(fv3jedi_geom),        intent(in) :: geom
real(kind=kind_real),      intent(in) :: u(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real),      intent(in) :: v(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)

! Locals
character(len=1024) :: dpath, fcore
character(len=:), allocatable :: str
integer :: rc
type(restart_file_type) :: rst

! Read file path and names from config
! ------------------------------------
call conf%get_or_die("datapath", str)
dpath = str
deallocate(str)
call conf%get_or_die("filename_core", str)
fcore = str
deallocate(str)

! Register and write
! -----------------
rc = register_restart_field( rst, trim(fcore), 'u', u, domain=geom%domain, position=north, &
                             longname = 'u_component_of_native_D_grid_wind', units = 'ms-1')
rc = register_restart_field( rst, trim(fcore), 'v', v, domain=geom%domain, position=east, &
                             longname = 'v_component_of_native_D_grid_wind', units = 'ms-1')

call save_restart(rst, directory=trim(adjustl(dpath)))
call free_restart_type(rst)

end subroutine write_d


! --------------------------------------------------------------------------------------------------


subroutine testing(comm, geom, conf, u, v)

type(fckit_mpi_comm),      intent(in) :: comm
type(fv3jedi_geom),        intent(in) :: geom
type(fckit_configuration), intent(in) :: conf
real(kind=kind_real),      intent(in) :: u(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real),      intent(in) :: v(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)

real(kind=kind_real) :: min_u_c, max_u_c, rms_u_c, min_v_c, max_v_c, rms_v_c, tol
real(kind=kind_real) :: min_u_r, max_u_r, rms_u_r, min_v_r, max_v_r, rms_v_r
integer :: isc, iec, jsc, jec, npz
real(kind=kind_real) :: tmp(3), gs3, gs3g

! Get min/max/rms from config
call conf%get_or_die("min_u", min_u_c)
call conf%get_or_die("max_u", max_u_c)
call conf%get_or_die("rms_u", rms_u_c)
call conf%get_or_die("min_v", min_v_c)
call conf%get_or_die("max_v", max_v_c)
call conf%get_or_die("rms_v", rms_v_c)
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
tmp(1) = minval(u(isc:iec,jsc:jec,1:npz), mask=u(isc:iec,jsc:jec,1:npz)<10e10)
tmp(2) = maxval(u(isc:iec,jsc:jec,1:npz), mask=u(isc:iec,jsc:jec,1:npz)<10e10)
tmp(3) =    sum(u(isc:iec,jsc:jec,1:npz)**2, mask=u(isc:iec,jsc:jec,1:npz)<10e10)
call comm%allreduce(tmp(1), min_u_r, fckit_mpi_min())
call comm%allreduce(tmp(2), max_u_r, fckit_mpi_max())
call comm%allreduce(tmp(3), rms_u_r, fckit_mpi_sum())
rms_u_r = sqrt(rms_u_r/gs3g)

! Min/Max/SumSquares (v)
tmp(1) = minval(v(isc:iec,jsc:jec,1:npz), mask=v(isc:iec,jsc:jec,1:npz)<10e10)
tmp(2) = maxval(v(isc:iec,jsc:jec,1:npz), mask=v(isc:iec,jsc:jec,1:npz)<10e10)
tmp(3) =    sum(v(isc:iec,jsc:jec,1:npz)**2, mask=v(isc:iec,jsc:jec,1:npz)<10e10)
call comm%allreduce(tmp(1), min_v_r, fckit_mpi_min())
call comm%allreduce(tmp(2), max_v_r, fckit_mpi_max())
call comm%allreduce(tmp(3), rms_v_r, fckit_mpi_sum())
rms_v_r = sqrt(rms_v_r/gs3g)

! If root processor then perform checks
if (comm%rank() == 0) then

  ! Assert that min/max/rms are within tolerance (relative)
  ! ------------------------------------------------------
  if (abs(min_u_r-min_u_c)/abs(min_u_c) > tol) then
      print *, 'u min out of tolerance (run) (config) (rel diff) (tol)', min_u_r, min_u_c, &
               abs(min_u_r-min_u_c)/abs(min_u_c), tol
      stop 1
  end if

  if (abs(max_u_r-max_u_c)/abs(max_u_c) > tol) then
      print *, 'u max out of tolerance (run) (config) (rel diff) (tol)', max_u_r, max_u_c, &
               abs(max_u_r-max_u_c)/abs(max_u_c), tol
      stop 1
  end if

  if (abs(rms_u_r-rms_u_c)/abs(rms_u_c) > tol) then
      print *, 'u rms out of tolerance (run) (config) (rel diff) (tol)', rms_u_r, rms_u_c, &
               abs(rms_u_r-rms_u_c)/abs(rms_u_c), tol
      stop 1
  end if

  if (abs(min_v_r-min_v_c)/abs(min_v_c) > tol) then
      print *, 'v min out of tolerance (run) (config) (rel diff) (tol)', min_v_r, min_v_c, &
               abs(min_v_r-min_v_c)/abs(min_v_c), tol
      stop 1
  end if

  if (abs(max_v_r-max_v_c)/abs(max_v_c) > tol) then
      print *, 'v max out of tolerance (run) (config) (rel diff) (tol)', max_v_r, max_v_c, &
               abs(max_v_r-max_v_c)/abs(max_v_c), tol
      stop 1
  end if

  if (abs(rms_v_r-rms_v_c)/abs(rms_v_c) > tol) then
      print *, 'v rms out of tolerance (run) (config) (rel diff) (tol)', rms_v_r, rms_v_c, &
               abs(rms_v_r-rms_v_c)/abs(rms_v_c), tol
      stop 1
  end if

  print*, 'u min/max/rms and v min/max/rms are within tolerance (PASSED)'

end if

end subroutine testing


! --------------------------------------------------------------------------------------------------


end program ua_to_ud
