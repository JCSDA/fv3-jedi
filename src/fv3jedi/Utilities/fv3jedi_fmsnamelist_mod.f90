module fv3jedi_fmsnamelist_mod

! fckit uses
use fckit_configuration_module, only: fckit_configuration

! fms uses
use fms_mod, only: check_nml_error
use mpp_mod, only: input_nml_file, get_ascii_file_num_lines_and_length, read_ascii_file

! fv3jedi uses
use fv3jedi_kinds_mod,           only: kind_real

! --------------------------------------------------------------------------------------------------

implicit none
private
public fv3jedi_fmsnamelist

type fv3jedi_fmsnamelist
  character(len=:), dimension(:), allocatable :: input_nml_file_fms
contains
  procedure, public :: replace_namelist
  procedure, public :: revert_namelist
  procedure :: file_to_namelist
  procedure :: conf_to_namelist
  final     :: dummy_final
end type fv3jedi_fmsnamelist

contains

! --------------------------------------------------------------------------------------------------

subroutine replace_namelist(self, conf)

! Arguments
class(fv3jedi_fmsnamelist), intent(inout) :: self
type(fckit_configuration),  intent(in)    :: conf

! Locals
character(len=1024) :: namelist_filename
character(len=:), allocatable :: str
logical :: use_int_namelist

! Optionally use the namelist that already exists in FMS
use_int_namelist = .false.
if (conf%has("use internal namelist")) then
  call conf%get_or_die("use internal namelist", use_int_namelist)
endif

! Replace the namelist with a file or individual yaml config
if (.not. use_int_namelist) then
  if (conf%has("namelist filename")) then
    call conf%get_or_die("namelist filename", str)
    namelist_filename = str
    call self%file_to_namelist(namelist_filename)
  else
    call self%conf_to_namelist(conf)
  endif
endif

end subroutine replace_namelist

! --------------------------------------------------------------------------------------------------

subroutine file_to_namelist(self, namelist_filename)

! Arguments
class(fv3jedi_fmsnamelist), intent(inout) :: self
character(len=*),           intent(in)    :: namelist_filename

! Locals
integer :: lines_and_length(2)
logical :: file_exists

! Check the file exists
! ---------------------
inquire( file=trim(namelist_filename), exist = file_exists )
if(.not. file_exists)then
  call abor1_ftn("fv3jedi_fmsnamelist_mod.file_to_namelist "//namelist_filename//" does not exist")
endif

! Deallocate local copy of namelist if needed
! -------------------------------------------
if (allocated(self%input_nml_file_fms)) deallocate(self%input_nml_file_fms)


! Move fms version of namelist to local
! -------------------------------------
if (allocated(input_nml_file)) call move_alloc(input_nml_file, self%input_nml_file_fms)


! Get number of lines and data length from new file
! -------------------------------------------------
lines_and_length = get_ascii_file_num_lines_and_length(namelist_filename)


! Allocate new fms namelist
! -------------------------
allocate(character(len=lines_and_length(2))::input_nml_file(lines_and_length(1)))


! Read file into fms namelist
! ---------------------------
call read_ascii_file(namelist_filename, lines_and_length(2), input_nml_file)


end subroutine file_to_namelist

! --------------------------------------------------------------------------------------------------

subroutine conf_to_namelist(self, conf)

! Arguments
class(fv3jedi_fmsnamelist), intent(inout) :: self
type(fckit_configuration),  intent(in)    :: conf

! Locals
integer :: nmls, ios, ierr
integer :: npx, npy, npz, ntiles, nwat
integer, allocatable :: layout(:), io_layout(:)
logical :: nested, regional, do_schmidt, hydrostatic
real(kind=kind_real) :: stretch_fac, target_lat, target_lon

namelist /fv_core_nml/ npx, npy, npz, ntiles, layout, io_layout, regional, nested, do_schmidt, &
                       target_lat, target_lon, stretch_fac, hydrostatic, nwat


! Replace FMS namelist file with geometry
! ---------------------------------------
allocate(layout(2))
allocate(io_layout(2))

nmls = 1

! Compulsory elements
nmls = nmls + 1
if (.not. conf%get('npx', npx)) call abor1_ftn("conf_to_namelist: did not find npx in config")
nmls = nmls + 1
if (.not. conf%get('npy', npy)) call abor1_ftn("conf_to_namelist: did not find npy in config")
nmls = nmls + 1
if (.not. conf%get('npz', npz)) call abor1_ftn("conf_to_namelist: did not find npz in config")

! Optional elements with defaults
nmls = nmls + 1; if (.not. conf%get('ntiles'     , ntiles     )) ntiles      = 6
nmls = nmls + 1; if (.not. conf%get('nwat'       , nwat       )) nwat        = 1
nmls = nmls + 1; if (.not. conf%get('layout'     , layout     )) layout      = (/1,1/)
nmls = nmls + 1; if (.not. conf%get('io_layout'  , io_layout  )) io_layout   = (/1,1/)
nmls = nmls + 1; if (.not. conf%get('nested'     , nested     )) nested      = .false.
nmls = nmls + 1; if (.not. conf%get('regional'   , regional   )) regional    = .false.
nmls = nmls + 1; if (.not. conf%get('do_schmidt' , do_schmidt )) do_schmidt  = .false.
nmls = nmls + 1; if (.not. conf%get('hydrostatic', hydrostatic)) hydrostatic = .true.
nmls = nmls + 1; if (.not. conf%get('stretch_fac', stretch_fac)) stretch_fac = 0.0_kind_real
nmls = nmls + 1; if (.not. conf%get('target_lat' , target_lat )) target_lat  = 0.0_kind_real
nmls = nmls + 1; if (.not. conf%get('target_lon' , target_lon )) target_lon  = 0.0_kind_real

! Deallocate existing FMS namelist array and reallocate
if (allocated(input_nml_file)) call move_alloc(input_nml_file, self%input_nml_file_fms)
allocate(character(len=255)::input_nml_file(nmls))

! Write new namelist to internal FMS namelist
write (input_nml_file, fv_core_nml, iostat=ios)
ierr = check_nml_error(ios,'fv3jedi geom writing input_nml_file')

end subroutine conf_to_namelist

! --------------------------------------------------------------------------------------------------

subroutine revert_namelist(self)

class(fv3jedi_fmsnamelist), target, intent(inout) :: self

! Replace fms namelist with the saved data
! ----------------------------------------
if (allocated(self%input_nml_file_fms) .and. allocated(input_nml_file)) then
  deallocate(input_nml_file)
  call move_alloc(self%input_nml_file_fms, input_nml_file)
else if (allocated(self%input_nml_file_fms)) then
  deallocate(self%input_nml_file_fms)
endif

end subroutine revert_namelist

! --------------------------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_fmsnamelist), intent(inout) :: self
call self%revert_namelist()
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fmsnamelist_mod
