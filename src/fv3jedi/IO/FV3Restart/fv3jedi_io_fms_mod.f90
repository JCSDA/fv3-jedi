! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_fms_mod

! oops
use datetime_mod
use string_utils,                 only: swap_name_member

! fckit
use fckit_configuration_module,   only: fckit_configuration

! fms
use fms_io_mod,                   only: restart_file_type, register_restart_field
use fms_io_mod,                   only: free_restart_type, restore_state, save_restart
use fms_io_mod,                   only: set_domain, nullify_domain
use mpp_domains_mod,              only: east, north, center, domain2D
use mpp_mod,                      only: mpp_pe, mpp_root_pe

! fv3jedi
use fv3jedi_constants_mod,        only: rad2deg
use fv3jedi_field_mod,            only: fv3jedi_field, hasfield
use fv3jedi_io_utils_mod,         only: vdate_to_datestring, replace_text, add_iteration
use fv3jedi_kinds_mod,            only: kind_real

! --------------------------------------------------------------------------------------------------

implicit none
private
public fv3jedi_io_fms, read_fields

! If adding a new file it is added here and object and config in setup
integer, parameter :: numfiles = 9

type fv3jedi_io_fms
 logical :: input_is_date_templated
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
 logical :: has_prefix
 character(len=128) :: prefix
 integer :: calendar_type
 ! Geometry copies
 type(domain2D), pointer :: domain
 integer :: npz
 contains
   procedure :: create
   procedure :: delete
   procedure :: read
   procedure :: write
   final     :: dummy_final
end type fv3jedi_io_fms

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, domain, npz)

class(fv3jedi_io_fms),     intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(domain2D), target,    intent(in)    :: domain
integer,                   intent(in)    :: npz

integer :: n
character(len=:), allocatable :: str
character(len=13) :: fileconf(numfiles)

! Get path to files
! -----------------
call conf%get_or_die("datapath",str)
if (len(str) > 128) &
  call abor1_ftn('fv3jedi_io_fms_mod.setup: datapath too long, max FMS char length= 128')

! For ensemble methods switch out member template
! -----------------------------------------------
call swap_name_member(conf, str)

self%datapath = str
deallocate(str)


!Set default filenames
!---------------------
self%filenames_conf(self%index_core) = 'fv_core.res.nc'
self%filenames_conf(self%index_trcr) = 'fv_tracer.res.nc'
self%filenames_conf(self%index_sfcd) = 'sfc_data.nc'
self%filenames_conf(self%index_sfcw) = 'fv_srf_wnd.res.nc'
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
  if (conf%has(fileconf(n))) then
    call conf%get_or_die(fileconf(n),str)
    if (len(str) > 128) call abor1_ftn("fv3jedi_io_fms_mod.setup: "//fileconf(n)//&
                                        " too long, max FMS char length= 128")
    call add_iteration(conf,str)
    self%filenames_conf(n) = str
    deallocate(str)
  endif

  ! Config filenames to filenames
  self%filenames(n) = trim(self%filenames_conf(n))

enddo

! Option to retrieve Ps from delp
! -------------------------------
self%ps_in_file = .false.
if (conf%has("psinfile")) then
  call conf%get_or_die("psinfile",self%ps_in_file)
endif

! Option to skip read/write of coupler file
! -----------------------------------------
self%skip_coupler = .false.
if (conf%has("skip coupler file")) then
  call conf%get_or_die("skip coupler file",self%skip_coupler)
endif

! Option to turn off prepending file with date
! --------------------------------------------
if (.not.conf%get("prepend files with date", self%prepend_date)) then
  self%prepend_date = .true.
endif

! Optionally the file name to be read is datetime templated
! ---------------------------------------------------------
if (conf%has("filename is datetime templated")) then
  call conf%get_or_die("filename is datetime templated", self%input_is_date_templated)
else
  self%input_is_date_templated = .false.
endif

! Option to overwrite date etc...
! -------------------------------
self%has_prefix = conf%has("prefix")
if (self%has_prefix) then
  call conf%get_or_die("prefix",str)
  self%prefix = trim(str)
endif

! Calendar type
! -------------
self%calendar_type = 2
if (conf%has("calendar type")) then
  call conf%get_or_die("calendar type", self%calendar_type)
endif

! Geometry copies
! ---------------
self%domain => domain
self%npz = npz

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_io_fms), intent(inout) :: self

if (associated(self%domain)) nullify(self%domain)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine read(self, vdate, fields)

class(fv3jedi_io_fms), intent(inout) :: self
type(datetime),        intent(inout) :: vdate
type(fv3jedi_field),   intent(inout) :: fields(:)
integer :: n

! Overwrite any datetime templates in the file names
! --------------------------------------------------
if (self%input_is_date_templated) call setup_date(self, vdate)

! Use prefix if present
! ---------------------
if (self%has_prefix) then
  do n = 1, numfiles
    self%filenames(n) = trim(self%prefix)//"."//trim(self%filenames_conf(n))
  enddo
endif

! Read meta data
! --------------
if (.not. self%skip_coupler) call read_meta(self, vdate)

! Read fields
! -----------
call read_fields(self, fields)

end subroutine read

! --------------------------------------------------------------------------------------------------

subroutine write(self, vdate, fields)

class(fv3jedi_io_fms), intent(inout) :: self
type(datetime),        intent(in)    :: vdate
type(fv3jedi_field),   intent(in)    :: fields(:)

! Overwrite any datetime templates in the file names
! --------------------------------------------------
call setup_date(self, vdate)

! Write metadata and fields
! -------------------------
call write_all(self, fields, vdate)

end subroutine write

! --------------------------------------------------------------------------------------------------

subroutine setup_date(self, vdate)

type(fv3jedi_io_fms), intent(inout) :: self
type(datetime),       intent(in)    :: vdate

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

subroutine read_meta(self, vdate)

type(fv3jedi_io_fms), intent(inout) :: self
type(datetime),       intent(inout) :: vdate         !< DateTime

integer :: date(6)
integer :: idate, itime
character(len=8) :: cdate
character(len=6) :: ctime

integer :: calendar_type
integer :: date_init(6)
character(len=64) :: vdate_string_file, vdate_string

! Get datetime from coupler.res
open(101, file=trim(adjustl(self%datapath))//'/'//self%filenames(self%index_cplr), form='formatted')
read(101, '(i6)')  calendar_type
read(101, '(6i6)') date_init
read(101, '(6i6)') date
close(101)

! Pad and convert to string
idate=date(1)*10000+date(2)*100+date(3)
itime=date(4)*10000+date(5)*100+date(6)
write(cdate,"(I0.8)") idate  ! Looks like YYYYMMDD
write(ctime,"(I0.6)") itime  ! Looks like HHmmSS

! Compute string form of the datetime in the fields
call datetime_to_string(vdate, vdate_string)

! Convert to string that matches format returned by datetime_to_string YYYY-MM-DDTHH:mm:SS
vdate_string_file = cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//'T'// &
                    ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//'Z'

! Assert
if (trim(vdate_string_file) .ne. trim(vdate_string)) &
  call abor1_ftn("io_cube_sphere_history.read.check_datetime: Datetime set in config (" &
                 //trim(vdate_string)//") does not match that read from the file (" &
                 //trim(vdate_string_file)//").")

end subroutine read_meta

! --------------------------------------------------------------------------------------------------

subroutine read_fields(self, fields)

implicit none
type(fv3jedi_io_fms), intent(inout) :: self
type(fv3jedi_field),  intent(inout) :: fields(:)

type(restart_file_type) :: restart(numfiles)
logical :: rstflag(numfiles)
integer :: n, indexrst, position, var, idrst

logical :: havedelp
integer :: indexof_ps, indexof_delp
real(kind=kind_real), allocatable :: delp(:,:,:)

! Set FMS IO internal domain
! --------------------------
call set_domain(self%domain)

! Register and read fields
! ------------------------
rstflag = .false.

! Check whether delp in fields
! ----------------------------
indexof_ps = -1
indexof_delp = -1
havedelp = hasfield(fields, 'delp', indexof_delp)

! Loop over fields and register their restart file
! ------------------------------------------------
do var = 1,size(fields)

  ! If need ps and not in file will compute from delp so read delp in place of ps
  if (trim(fields(var)%fv3jedi_name) == 'ps' .and. .not.self%ps_in_file) then
    indexof_ps = var
    if (havedelp) cycle ! Do not register delp twice
    deallocate(fields(indexof_ps)%array)
    allocate(fields(indexof_ps)%array(fields(indexof_ps)%isc:fields(indexof_ps)%iec, &
                fields(indexof_ps)%jsc:fields(indexof_ps)%jec,1:self%npz))
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
      call abor1_ftn("fv3jedi_io_fms_mod: Abort, field "//trim(fields(var)%short_name)//&
                      " does not have IOFile specified in the FieldSets metadata or it"&
                      " does not match options in fms IO module")
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
                                  domain=self%domain, position=position )

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
  allocate(delp(fields(indexof_ps)%isc:fields(indexof_ps)%iec, &
                fields(indexof_ps)%jsc:fields(indexof_ps)%jec,1:self%npz))
  if (.not. havedelp) then
    delp = fields(indexof_ps)%array
    deallocate(fields(indexof_ps)%array)
    allocate(fields(indexof_ps)%array(fields(indexof_ps)%isc:fields(indexof_ps)%iec, &
                fields(indexof_ps)%jsc:fields(indexof_ps)%jec,1))
  else
    delp = fields(indexof_delp)%array
  endif
  fields(indexof_ps)%array(:,:,1) = sum(delp,3)
  fields(indexof_ps)%short_name = 'ps'
endif

! Nullify FMS IO internal domain
! ------------------------------
call nullify_domain()

end subroutine read_fields

! --------------------------------------------------------------------------------------------------

subroutine write_all(self, fields, vdate)

implicit none
type(fv3jedi_io_fms), intent(inout) :: self
type(fv3jedi_field),  intent(in)    :: fields(:)     !< Fields to be written
type(datetime),       intent(in)    :: vdate         !< DateTime

logical :: rstflag(numfiles)
integer :: n, indexrst, position, var, idrst, date(6)
integer :: idate, isecs
type(restart_file_type) :: restart(numfiles)
character(len=64)  :: datefile

! Set FMS IO internal domain
! --------------------------
call set_domain(self%domain)

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

! Use prefix if present
! ---------------------
if (self%has_prefix) then
  do n = 1, numfiles
    self%filenames(n) = trim(self%prefix)//"."//trim(self%filenames_conf(n))
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
      call abor1_ftn("fv3jedi_io_fms_mod: Abort, field "//trim(fields(var)%short_name)//&
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
                                  fields(var)%short_name, fields(var)%array, domain=self%domain, &
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
   write( 101, '(i6,8x,a)' ) self%calendar_type, &
        '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
   write( 101, '(6i6,8x,a)') date, 'Model start time:   year, month, day, hour, minute, second'
   write( 101, '(6i6,8x,a)') date, 'Current model time: year, month, day, hour, minute, second'
   close(101)
endif

! Nullify FMS IO internal domain
! ------------------------------
call nullify_domain()

end subroutine write_all

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
type(fv3jedi_io_fms), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_fms_mod
