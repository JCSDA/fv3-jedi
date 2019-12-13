module fv3jedi_io_gfs_mod

use datetime_mod
use iso_c_binding
use string_utils, only: swap_name_member

use fv3jedi_constants_mod,      only: rad2deg
use fv3jedi_geom_mod,           only: fv3jedi_geom
use fv3jedi_field_mod,          only: fv3jedi_field
use fv3jedi_kinds_mod,          only: kind_real
use mpp_mod,                    only: mpp_pe, mpp_root_pe
use fms_io_mod,                 only: restart_file_type, register_restart_field, &
                                      free_restart_type, restore_state, save_restart
use fckit_configuration_module, only: fckit_configuration

implicit none
private
public fv3jedi_io_gfs

type fv3jedi_io_gfs
 character(len=128) :: datapath_ti
 character(len=128) :: datapath_sp
 character(len=128) :: filename_spec
 character(len=128) :: filename_core
 character(len=128) :: filename_trcr
 character(len=128) :: filename_sfcd
 character(len=128) :: filename_sfcw
 character(len=128) :: filename_cplr
 logical :: ps_in_file
 contains
  procedure :: setup
  procedure :: read_meta
  procedure :: read_fields
  procedure :: write_all
  final     :: dummy_final
end type fv3jedi_io_gfs

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine setup(self,f_conf,psinfile)
use string_utils

class(fv3jedi_io_gfs),     intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_conf
integer, optional,         intent(in)    :: psinfile

character(len=:), allocatable :: str

!Set filenames
!--------------
self%filename_core = 'fv_core.res.nc'
self%filename_trcr = 'fv_tracer.res.nc'
self%filename_sfcd = 'sfc_data.nc'
self%filename_sfcw = 'srf_wnd.nc'
self%filename_cplr = 'coupler.res'

call f_conf%get_or_die("datapath_tile",str)
if (len(str) > 128) &
  call abor1_ftn('fv3jedi_io_gfs_mod.setup: datapath_tile too long, max FMS char length= 128')

call swap_name_member(f_conf, str)

self%datapath_ti = str
deallocate(str)

if (f_conf%has("filename_core")) then
   call f_conf%get_or_die("filename_core",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_core too long, max FMS char length= 128')
   self%filename_core = str
   deallocate(str)
endif
if (f_conf%has("filename_trcr")) then
   call f_conf%get_or_die("filename_trcr",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_trcr too long, max FMS char length= 128')
   self%filename_trcr = str
   deallocate(str)
endif
if (f_conf%has("filename_sfcd")) then
   call f_conf%get_or_die("filename_sfcd",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_sfcd too long, max FMS char length= 128')
   self%filename_sfcd = str
   deallocate(str)
endif
if (f_conf%has("filename_sfcw")) then
   call f_conf%get_or_die("filename_sfcw",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_sfcw too long, max FMS char length= 128')
   self%filename_sfcw = str
   deallocate(str)
endif
if (f_conf%has("filename_cplr")) then
   call f_conf%get_or_die("filename_cplr",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_cplr too long, max FMS char length= 128')
   self%filename_cplr = str
   deallocate(str)
endif

if (f_conf%has("filename_spec")) then
   call f_conf%get_or_die("filename_spec",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: filename_spec too long, max FMS char length= 128')
   self%filename_spec = str
   deallocate(str)
   call f_conf%get_or_die("datapath_spec",str)
   if (len(str) > 128) &
     call abor1_ftn('fv3jedi_io_gfs_mod.setup: datapath_spec too long, max FMS char length= 128')
   self%datapath_sp = str
   deallocate(str)
else
   self%filename_spec = "null"
   self%datapath_sp = "null"
endif

self%ps_in_file = .false.
if (present(psinfile)) then
  if (psinfile == 1) then
    self%ps_in_file = .true.
  endif
endif

end subroutine setup

! ------------------------------------------------------------------------------

subroutine read_meta(self, geom, vdate, calendar_type, date_init)

implicit none
class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom          !< Geometry
type(datetime),        intent(inout) :: vdate         !< DateTime
integer,               intent(inout) :: calendar_type !< GFS calendar type
integer,               intent(inout) :: date_init(6)  !< GFS date intialized

integer :: date(6)
integer(kind=c_int) :: idate, isecs

type(restart_file_type)  :: restart_spec
integer :: id_restart
real(kind=kind_real), allocatable, dimension(:,:) :: grid_lat, grid_lon


! Read Lat-Lon and check consitency with geom
! -------------------------------------------
if (trim(self%filename_spec) .ne. "null" .and. trim(self%datapath_sp) .ne. "null") then

  allocate(grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec))
  allocate(grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec))

  id_restart = register_restart_field( restart_spec, trim(self%filename_spec), "grid_latt", grid_lat, &
                                       domain=geom%domain )
  id_restart = register_restart_field( restart_spec, trim(self%filename_spec), "grid_lont", grid_lon, &
                                       domain=geom%domain )

  call restore_state(restart_spec, directory=trim(adjustl(self%datapath_sp)))
  call free_restart_type(restart_spec)

  if ( (maxval(abs(grid_lat-rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec))) > 1.0e-4) &
  .or. (maxval(abs(grid_lon-rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec))) > 1.0e-4) ) then
    call abor1_ftn("Grid in gridspec file does not match that in the geometry")
  endif

  deallocate(grid_lat)
  deallocate(grid_lon)
endif

! Get dates from coupler.res
!---------------------------
open(101, file=trim(adjustl(self%datapath_ti))//'/'//self%filename_cplr, form='formatted')
read(101, '(i6)')  calendar_type
read(101, '(6i6)') date_init
read(101, '(6i6)') date
close(101)
idate=date(1)*10000+date(2)*100+date(3)
isecs=date(4)*3600+date(5)*60+date(6)
call datetime_from_ifs(vdate, idate, isecs)

end subroutine read_meta

! ------------------------------------------------------------------------------

subroutine read_fields(self, geom, fields)

implicit none
class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom
type(fv3jedi_field),   intent(inout) :: fields(:)

type(restart_file_type), pointer :: restart
type(restart_file_type), target  :: restart_core
type(restart_file_type), target  :: restart_trcr
type(restart_file_type), target  :: restart_sfcd
type(restart_file_type), target  :: restart_sfcw
type(restart_file_type)  :: restart_spec
logical :: read_core, read_trcr, read_sfcd, read_sfcw
integer :: var, id_restart
character(len=255) :: filename
integer :: compute_ps, compute_ps_type
integer :: indexof_ps, indexof_delp
logical :: assocps, assocdelp
real(kind=kind_real), allocatable :: delp(:,:,:)

! Register and read fields
! ------------------------
read_core = .false.
read_trcr = .false.
read_sfcd = .false.
read_sfcw = .false.

assocdelp = .false.
assocps   = .false.
do var = 1,size(fields)
  if (trim(fields(var)%short_name) == 'DELP') assocdelp = .true.
  if (trim(fields(var)%short_name) == 'ps')   assocps   = .true.
enddo

if (assocps) then
  if (self%ps_in_file) then
    compute_ps_type = 0
  else
    compute_ps_type = 1
    if (assocdelp) compute_ps_type = 2
  endif
endif

do var = 1,size(fields)

  compute_ps = 0

  select case (trim(fields(var)%short_name))
  case("u","v","ud","vd","ua","va","phis","T","W","DZ","psi","chi","vort","divg","tv")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
  case("DELP","delp")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
    indexof_delp = var
  case("sphum","ice_wat","liq_wat","rainwat","snowwat","graupel","cld_amt","rh",&
       "o3mr","sulf","bc1","bc2","oc1","oc2",&
       "dust1","dust2","dust3","dust4","dust5","seas1","seas2","seas3","seas4")
    filename = self%filename_trcr
    restart => restart_trcr
    read_trcr = .true.
  case("slmsk","sheleg","tsea","vtype","stype","vfrac","stc","smc","snwdph","f10m")
    filename = self%filename_sfcd
    restart => restart_sfcd
    read_sfcd = .true.
  case("u_srf","v_srf")
    filename = self%filename_sfcw
    restart => restart_sfcw
    read_sfcw = .true.
  case("ps")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
    compute_ps = compute_ps_type
    if (compute_ps .ne. 0) then
      indexof_ps = var
      if (compute_ps == 1) allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
    endif
  case default
    call abor1_ftn("read_gfs: filename not set for "//trim(fields(var)%short_name))
  end select

  if (compute_ps == 0) then
    id_restart = register_restart_field( restart, trim(filename), trim(fields(var)%short_name), &
                                         fields(var)%array, domain=geom%domain, &
                                         position=fields(var)%staggerloc )
  elseif (compute_ps == 1) then
    id_restart = register_restart_field( restart, trim(filename), "DELP", &
                                         delp, domain=geom%domain, &
                                         position=fields(var)%staggerloc )
  endif

enddo

if (read_core) then
  call restore_state(restart_core, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_core)
endif
if (read_trcr) then
  call restore_state(restart_trcr, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_trcr)
endif
if (read_sfcd) then
  call restore_state(restart_sfcd, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_sfcd)
endif
if (read_sfcw) then
  call restore_state(restart_sfcw, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_sfcw)
endif

!Compute ps from DELP
if (compute_ps > 0) then
  if (compute_ps == 2) then
    fields(indexof_ps)%array(:,:,1) = sum(fields(indexof_delp)%array,3)
  elseif (compute_ps == 1) then
    fields(indexof_ps)%array(:,:,1) = sum(delp,3)
    deallocate(delp)
  endif
endif

end subroutine read_fields

! ------------------------------------------------------------------------------

subroutine write_all(self, geom, fields, vdate, calendar_type, date_init)

implicit none
class(fv3jedi_io_gfs), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom          !< Geom
type(fv3jedi_field),   intent(in)    :: fields(:)     !< Fields to be written
type(datetime),        intent(in)    :: vdate         !< DateTime
integer,               intent(in)    :: calendar_type !< GFS calendar type
integer,               intent(in)    :: date_init(6)  !< GFS date intialized

integer :: date(6)
integer(kind=c_int) :: idate, isecs
type(restart_file_type), pointer :: restart
type(restart_file_type), target  :: restart_core
type(restart_file_type), target  :: restart_trcr
type(restart_file_type), target  :: restart_sfcd
type(restart_file_type), target  :: restart_sfcw
integer :: var, id_restart, jlev
logical :: read_core, read_trcr, read_sfcd, read_sfcw, register
character(len=255) :: filename
character(len=64)  :: datefile

! Append file start with the datetime
! -----------------------------------
call datetime_to_ifs(vdate, idate, isecs)
date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2,A1)') date(1),date(2),date(3),".",date(4),date(5),date(6),"."
self%filename_core = trim(datefile)//trim(self%filename_core)
self%filename_trcr = trim(datefile)//trim(self%filename_trcr)
self%filename_sfcd = trim(datefile)//trim(self%filename_sfcd)
self%filename_sfcw = trim(datefile)//trim(self%filename_sfcw)
self%filename_cplr = trim(datefile)//trim(self%filename_cplr)

! Register and write fields
! -------------------------
read_core = .false.
read_trcr = .false.
read_sfcd = .false.
read_sfcw = .false.
do var = 1,size(fields)

  select case (trim(fields(var)%short_name))
  case("u","v","ud","vd","ua","va","phis","T","ps","DELP","delp","W","DZ","psi","chi","vort","divg","tv")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
  case("sphum","ice_wat","liq_wat","rainwat","snowwat","graupel","cld_amt","rh",&
       "o3mr","sulf","bc1","bc2","oc1","oc2",&
       "dust1","dust2","dust3","dust4","dust5","seas1","seas2","seas3","seas4")
    filename = self%filename_trcr
    restart => restart_trcr
    read_trcr = .true.
  case("slmsk","sheleg","tsea","vtype","stype","vfrac","stc","smc","snwdph","f10m")
    filename = self%filename_sfcd
    restart => restart_sfcd
    read_sfcd = .true.
  case("u_srf","v_srf")
    filename = self%filename_sfcw
    restart => restart_sfcw
    read_sfcw = .true.
  case default
    call abor1_ftn("write_gfs: filename not set for "//trim(fields(var)%short_name))
  end select

  id_restart = register_restart_field( restart, filename, fields(var)%short_name, fields(var)%array, &
                                         domain=geom%domain, position=fields(var)%staggerloc, &
                                         longname = trim(fields(var)%long_name), units = trim(fields(var)%units) )

enddo

if (read_core) then
  call save_restart(restart_core, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_core)
endif
if (read_trcr) then
  call save_restart(restart_trcr, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_trcr)
endif
if (read_sfcd) then
  call save_restart(restart_sfcd, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_sfcd)
endif
if (read_sfcw) then
  call save_restart(restart_sfcw, directory=trim(adjustl(self%datapath_ti)))
  call free_restart_type(restart_sfcw)
endif


!Write date/time info in coupler.res
!-----------------------------------
if (mpp_pe() == mpp_root_pe()) then
   open(101, file=trim(adjustl(self%datapath_ti))//'/'//trim(adjustl(self%filename_cplr)), form='formatted')
   write( 101, '(i6,8x,a)' ) calendar_type, &
        '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
   write( 101, '(6i6,8x,a)') date_init, 'Model start time:   year, month, day, hour, minute, second'
   write( 101, '(6i6,8x,a)') date,      'Current model time: year, month, day, hour, minute, second'
   close(101)
endif

end subroutine write_all

! ------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_io_gfs), intent(inout) :: self
end subroutine dummy_final

! ------------------------------------------------------------------------------

end module fv3jedi_io_gfs_mod
