module fv3jedi_io_gfs_mod

! iso
use iso_c_binding

! oops
use datetime_mod
use string_utils,                 only: swap_name_member

! fckit
use fckit_configuration_module,   only: fckit_configuration

! fms
use fms_io_mod,                   only: restart_file_type, register_restart_field
use fms_io_mod,                   only: free_restart_type, restore_state, save_restart
use mpp_domains_mod,              only: east, north, center
use mpp_mod,                      only: mpp_pe, mpp_root_pe

! fv3jedi
use fv3jedi_constants_mod,        only: rad2deg
use fv3jedi_field_mod,            only: fv3jedi_field, has_field
use fv3jedi_geom_mod,             only: fv3jedi_geom
use fv3jedi_io_utils_mod,         only: vdate_to_datestring, replace_text
use fv3jedi_kinds_mod,            only: kind_real

! --------------------------------------------------------------------------------------------------

implicit none
private
public fv3jedi_io_gfs

! If adding a new file it is added here and object and config in setup
integer, parameter :: numfiles = 9

type fv3jedi_io_gfs
 character(len=128) :: datapath
 character(len=128) :: filenames(numfiles)
 character(len=128) :: filenames_conf(numfiles)
 integer :: index_core = 1  ! Files like fv_core.res.tile<n>.nc
 integer :: index_trcr = 2  ! Files like fv_tracer.res.tile<n>.nc
 integer :: index_sfcd = 3  ! Files like sfc_data.tile<n>.nc
 integer :: index_sfcw = 4  ! Files like fv_srf_wnd.res.tile<n>.nc
 integer :: index_cplr = 5  ! Files like coupler.res
 integer :: index_spec = 6  ! Files like grid_spec.res.tile<n>.nc
 integer :: index_phys = 7  ! Files like phy_data.tile<n>.nc
 integer :: index_orog = 8  ! Files like C<npx-1>_oro_data.tile<n>.nc
 integer :: index_cold = 9  ! Files like gfs_data.tile<n>.nc
 logical :: ps_in_file
 logical :: skip_coupler
 logical :: prepend_date
 contains
  procedure :: setup_conf  ! Setup for when config is available, called from constructors
  procedure :: setup_date  ! Setup when datetime is available
  procedure :: read_meta
  procedure :: read_fields
  procedure :: write
  final     :: dummy_final
end type fv3jedi_io_gfs

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup_conf(self, f_conf)

class(fv3jedi_io_gfs),     intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_conf

integer :: n
character(len=:), allocatable :: str
character(len=13) :: fileconf(numfiles)

! Get path to files
! -----------------
call f_conf%get_or_die("datapath",str)
if (len(str) > 128) &
  call abor1_ftn('fv3jedi_io_gfs_mod.setup: datapath too long, max FMS char length= 128')

! For ensemble methods switch out member template
! -----------------------------------------------
call swap_name_member(f_conf, str)

self%datapath = str
deallocate(str)


!Set default filenames
!---------------------
self%filenames_conf(self%index_core) = 'fv_core.res.nc'
self%filenames_conf(self%index_trcr) = 'fv_tracer.res.nc'
self%filenames_conf(self%index_sfcd) = 'sfc_data.nc'
self%filenames_conf(self%index_sfcw) = 'srf_wnd.nc'
self%filenames_conf(self%index_cplr) = 'coupler.res'
self%filenames_conf(self%index_spec) = 'null'
self%filenames_conf(self%index_phys) = 'phy_data.nc'
self%filenames_conf(self%index_orog) = 'oro_data.nc'
self%filenames_conf(self%index_cold) = 'gfs_data.nc'

! Configuration to parse for the filenames
! ----------------------------------------
fileconf(self%index_core) = "filename_core"
fileconf(self%index_trcr) = "filename_trcr"
fileconf(self%index_sfcd) = "filename_sfcd"
fileconf(self%index_sfcw) = "filename_sfcw"
fileconf(self%index_cplr) = "filename_cplr"
fileconf(self%index_spec) = "filename_spec"
fileconf(self%index_phys) = "filename_phys"
fileconf(self%index_orog) = "filename_orog"
fileconf(self%index_cold) = "filename_cold"


! Set files names based on user input
! -----------------------------------
do n = 1, numfiles

  ! Retrieve user input filenames if available
  if (f_conf%has(fileconf(n))) then
    call f_conf%get_or_die(fileconf(n),str)
    if (len(str) > 128) call abor1_ftn("fv3jedi_io_gfs_mod.setup: "//fileconf(n)//&
                                        " too long, max FMS char length= 128")
    self%filenames_conf(n) = str
    deallocate(str)
  endif

  ! Config filenames to filenames
  self%filenames(n) = trim(self%filenames_conf(n))

enddo

! Option to retrieve Ps from delp
! -------------------------------
self%ps_in_file = .false.
if (f_conf%has("psinfile")) then
  call f_conf%get_or_die("psinfile",self%ps_in_file)
endif

! Option to skip read/write of coupler file
! -----------------------------------------
self%skip_coupler = .false.
if (f_conf%has("skip coupler file")) then
  call f_conf%get_or_die("skip coupler file",self%skip_coupler)
endif

! Option to turn off prepending file with date
! --------------------------------------------
if (.not.f_conf%get("prepend files with date", self%prepend_date)) then
  self%prepend_date = .true.
endif

end subroutine setup_conf

! --------------------------------------------------------------------------------------------------

subroutine setup_date(self, vdate)

class(fv3jedi_io_gfs),               intent(inout) :: self
type(datetime),                      intent(in)    :: vdate

integer :: n
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss

! Datetime to strings
! -------------------
call vdate_to_datestring(vdate, yyyy=yyyy, mm=mm, dd=dd, hh=hh, min=min, ss=ss)

do n = 1, numfiles

  ! Config filenames to filenames
  self%filenames(n) = trim(self%filenames_conf(n))

  ! Swap out datetime templates if needed
  if (index(self%filenames(n),"%yyyy") > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%yyyy',yyyy)
  if (index(self%filenames(n),"%mm"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%mm'  ,mm  )
  if (index(self%filenames(n),"%dd"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%dd'  ,dd  )
  if (index(self%filenames(n),"%hh"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%hh'  ,hh  )
  if (index(self%filenames(n),"%MM"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%MM'  ,min )
  if (index(self%filenames(n),"%ss"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%ss'  ,ss  )

enddo

end subroutine setup_date

! --------------------------------------------------------------------------------------------------

subroutine read_meta(self, geom, vdate, calendar_type, date_init)

class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom          !< Geometry
type(datetime),        intent(inout) :: vdate         !< DateTime
integer,               intent(inout) :: calendar_type !< GFS calendar type
integer,               intent(inout) :: date_init(6)  !< GFS date intialized

integer :: date(6)
integer(kind=c_int) :: idate, isecs

type(restart_file_type)  :: restart_spec
integer :: idrst
real(kind=kind_real), allocatable, dimension(:,:) :: grid_lat, grid_lon


! Read Lat-Lon and check consitency with geom
! -------------------------------------------
if (trim(self%filenames(self%index_spec)) .ne. "null" .and. trim(self%datapath) .ne. "null") then

  allocate(grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec))
  allocate(grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec))

  idrst = register_restart_field( restart_spec, trim(self%filenames(self%index_spec)), &
                                  "grid_latt", grid_lat, domain=geom%domain )
  idrst = register_restart_field( restart_spec, trim(self%filenames(self%index_spec)), &
                                  "grid_lont", grid_lon, domain=geom%domain )

  call restore_state(restart_spec, directory=trim(adjustl(self%datapath)))
  call free_restart_type(restart_spec)

  if ((maxval(abs(grid_lat-rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec)))>1.0e-4) &
  .or.(maxval(abs(grid_lon-rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec)))>1.0e-4))then
    call abor1_ftn("Grid in gridspec file does not match that in the geometry")
  endif
  deallocate(grid_lat)
  deallocate(grid_lon)
endif

! Get dates from coupler.res
!---------------------------
if (.not. self%skip_coupler) then
  open(101, file=trim(adjustl(self%datapath))//'/'//self%filenames(self%index_cplr), form='formatted')
  read(101, '(i6)')  calendar_type
  read(101, '(6i6)') date_init
  read(101, '(6i6)') date
  close(101)
  idate=date(1)*10000+date(2)*100+date(3)
  isecs=date(4)*3600+date(5)*60+date(6)
else
  idate = 20000101
  isecs = 0
endif

! Set datetime
call datetime_from_ifs(vdate, idate, isecs)

end subroutine read_meta

! --------------------------------------------------------------------------------------------------

subroutine read_fields(self, geom, fields)

implicit none
class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom
type(fv3jedi_field),   intent(inout) :: fields(:)

type(restart_file_type) :: restart(numfiles)
logical :: rstflag(numfiles)
integer :: n, indexrst, position, var, idrst

logical :: havedelp
integer :: indexof_ps, indexof_delp
real(kind=kind_real), allocatable :: delp(:,:,:)

! Register and read fields
! ------------------------
rstflag = .false.

! Check whether delp in fields
! ----------------------------
indexof_ps = -1
indexof_delp = -1
havedelp = has_field(fields, 'delp', indexof_delp)

! Loop over fields and register their restart file
! ------------------------------------------------
do var = 1,size(fields)

  ! If need ps and not in file will compute from delp so read delp in place of ps
  if (trim(fields(var)%fv3jedi_name) == 'ps' .and. .not.self%ps_in_file) then
    indexof_ps = var
    if (havedelp) cycle ! Do not register delp twice
    deallocate(fields(indexof_ps)%array)
    allocate(fields(indexof_ps)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
    fields(indexof_ps)%short_name = 'DELP'
  endif

  select case (trim(fields(var)%io_file))
    case("core")
      indexrst = self%index_core
    case("tracer")
      indexrst = self%index_trcr
    case("surface")
      indexrst = self%index_sfcd
    case("surface_wind")
      indexrst = self%index_sfcw
    case("physics")
      indexrst = self%index_phys
    case("orography")
      indexrst = self%index_orog
    case("cold")
      indexrst = self%index_cold
    case("default")
      call abor1_ftn("fv3jedi_io_gfs_mod: Abort, field "//trim(fields(var)%short_name)//&
                      " does not have IOFile specified in the FieldSets metadata or it"&
                      " does not match options in gfs IO module")
  end select

  ! Convert fv3jedi position to fms position
  position = center
  if (fields(var)%staggerloc == 'northsouth') then
    position = north
  elseif (fields(var)%staggerloc == 'eastwest') then
    position = east
  endif

  ! Flag to read this restart
  rstflag(indexrst) = .true.

  ! Register this restart
  idrst = register_restart_field( restart(indexrst), trim(self%filenames(indexrst)), &
                                  trim(fields(var)%short_name), fields(var)%array, &
                                  domain=geom%domain, position=position )

enddo

! Loop over files and read fields
! -------------------------------
do n = 1, numfiles
  if (rstflag(n)) then
    call restore_state(restart(n), directory=trim(adjustl(self%datapath)))
    call free_restart_type(restart(n))
  endif
enddo

! Compute ps from DELP
! --------------------
if (indexof_ps > 0) then
  allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  if (.not. havedelp) then
    delp = fields(indexof_ps)%array
    deallocate(fields(indexof_ps)%array)
    allocate(fields(indexof_ps)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1))
  else
    delp = fields(indexof_delp)%array
  endif
  fields(indexof_ps)%array(:,:,1) = sum(delp,3)
  fields(indexof_ps)%short_name = 'ps'
endif

end subroutine read_fields

! --------------------------------------------------------------------------------------------------

subroutine write(self, geom, fields, vdate, calendar_type, date_init)

implicit none
class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom          !< Geom
type(fv3jedi_field),   intent(in)    :: fields(:)     !< Fields to be written
type(datetime),        intent(in)    :: vdate         !< DateTime
integer,               intent(in)    :: calendar_type !< GFS calendar type
integer,               intent(in)    :: date_init(6)  !< GFS date intialized

logical :: rstflag(numfiles)
integer :: n, indexrst, position, var, idrst, date(6)
integer(kind=c_int) :: idate, isecs
type(restart_file_type) :: restart(numfiles)
character(len=64)  :: datefile

! Get datetime
! ------------
call datetime_to_ifs(vdate, idate, isecs)
date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

! Convert integer datetime into string and prepend file names
! -----------------------------------------------------------
write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2,A1)') date(1),date(2),date(3),".",&
                                                      date(4),date(5),date(6),"."

if (self%prepend_date) then
  do n = 1, numfiles
    self%filenames(n) = trim(datefile)//trim(self%filenames(n))
  enddo
endif

rstflag = .false.

! Loop over fields and register their restart file
! ------------------------------------------------
do var = 1,size(fields)

  select case (trim(fields(var)%io_file))
    case("core")
      indexrst = self%index_core
    case("tracer")
      indexrst = self%index_trcr
    case("surface")
      indexrst = self%index_sfcd
    case("surface_wind")
      indexrst = self%index_sfcw
    case("physics")
      indexrst = self%index_phys
    case("orography")
      indexrst = self%index_orog
    case("cold")
      indexrst = self%index_cold
    case("default")
      call abor1_ftn("fv3jedi_io_gfs_mod: Abort, field "//trim(fields(var)%short_name)//&
                      " does not have IOFile specified in the FieldSets metadata")
  end select

  ! Convert fv3jedi position to fms position
  position = center
  if (fields(var)%staggerloc == 'northsouth') then
    position = north
  elseif (fields(var)%staggerloc == 'eastwest') then
    position = east
  endif

  ! Flag to read this restart
  rstflag(indexrst) = .true.

  ! Register this restart
  idrst = register_restart_field( restart(indexrst), trim(self%filenames(indexrst)), &
                                  fields(var)%short_name, fields(var)%array, domain=geom%domain, &
                                  position=position, longname = trim(fields(var)%long_name), &
                                  units = trim(fields(var)%units) )

enddo


! Loop over files and write fields
! -------------------------------
do n = 1, numfiles
  if (rstflag(n)) then
    call save_restart(restart(n), directory=trim(adjustl(self%datapath)))
    call free_restart_type(restart(n))
  endif
enddo


!Write date/time info in coupler.res
!-----------------------------------
if (mpp_pe() == mpp_root_pe() .and. .not. self%skip_coupler) then
   open(101, file = trim(adjustl(self%datapath))//'/'// &
        trim(adjustl(self%filenames(self%index_cplr))), form='formatted')
   write( 101, '(i6,8x,a)' ) calendar_type, &
        '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
   write( 101, '(6i6,8x,a)') date_init, 'Model start time:   year, month, day, hour, minute, second'
   write( 101, '(6i6,8x,a)') date,      'Current model time: year, month, day, hour, minute, second'
   close(101)
endif

end subroutine write

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
type(fv3jedi_io_gfs), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_gfs_mod
