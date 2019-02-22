module fv3jedi_io_gfs_mod

use config_mod
use datetime_mod
use iso_c_binding
use fv3jedi_constants_mod,   only: rad2deg
use fv3jedi_geom_mod,        only: fv3jedi_geom
use fv3jedi_field_mod,       only: fv3jedi_field
use fv3jedi_kinds_mod,       only: kind_real
use mpp_mod,                 only: mpp_pe, mpp_root_pe
use fms_io_mod,              only: restart_file_type, register_restart_field, &
                                   free_restart_type, restore_state, save_restart

implicit none
private
public fv3jedi_io_gfs

type fv3jedi_io_gfs
 character(len=255) :: datapath_ti
 character(len=255) :: datapath_sp
 character(len=255) :: filename_spec
 character(len=255) :: filename_core
 character(len=255) :: filename_trcr
 character(len=255) :: filename_sfcd
 character(len=255) :: filename_sfcw
 character(len=255) :: filename_cplr
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

subroutine setup(self,c_conf)

class(fv3jedi_io_gfs), intent(inout) :: self
type(c_ptr),           intent(in)    :: c_conf

!Set filenames
!--------------
self%filename_core = 'fv_core.res.nc'
self%filename_trcr = 'fv_tracer.res.nc'
self%filename_sfcd = 'sfc_data.nc'
self%filename_sfcw = 'srf_wnd.nc'
self%filename_cplr = 'coupler.res'

self%datapath_ti = config_get_string(c_conf,len(self%datapath_ti),"datapath_tile")

if (config_element_exists(c_conf,"filename_core")) then
   self%filename_core = config_get_string(c_conf,len(self%filename_core),"filename_core")
endif
if (config_element_exists(c_conf,"filename_trcr")) then
   self%filename_trcr = config_get_string(c_conf,len(self%filename_trcr),"filename_trcr")
endif
if (config_element_exists(c_conf,"filename_sfcd")) then
   self%filename_sfcd = config_get_string(c_conf,len(self%filename_sfcd),"filename_sfcd")
endif
if (config_element_exists(c_conf,"filename_sfcw")) then
   self%filename_sfcw = config_get_string(c_conf,len(self%filename_sfcw),"filename_sfcw")
endif
if (config_element_exists(c_conf,"filename_cplr")) then
   self%filename_cplr = config_get_string(c_conf,len(self%filename_cplr),"filename_cplr")
endif

if (config_element_exists(c_conf,"filename_spec")) then
   self%filename_spec = config_get_string(c_conf,len(self%filename_spec),"filename_spec")
   self%datapath_sp   = config_get_string(c_conf,len(self%datapath_sp),"datapath_spec")  
else
   self%filename_spec = "null"
   self%datapath_sp = "null"
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
open(101, file=trim(adjustl(self%datapath_ti))//self%filename_cplr, form='formatted')
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
logical :: compute_ps
integer :: indexof_ps
real(kind=kind_real), allocatable :: delp(:,:,:)

! Register and read fields
! ------------------------
read_core = .false.
read_trcr = .false.
read_sfcd = .false.
read_sfcw = .false.
do var = 1,size(fields)

  compute_ps = .false.

  select case (trim(fields(var)%short_name))
  case("u","v","ud","vd","ua","va","phis","T","DELP","W","DZ")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
  case("sphum","ice_wat","liq_wat","o3mr")
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
    compute_ps = .true.
    indexof_ps = var
    allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  case default
    call abor1_ftn("read_gfs: filename not set for "//trim(fields(var)%short_name))
  end select

  if (.not.compute_ps) then
    id_restart = register_restart_field( restart, trim(filename), trim(fields(var)%short_name), &
                                         fields(var)%array, domain=geom%domain, &
                                         position=fields(var)%staggerloc )
  else
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
if (allocated(delp)) then
  fields(indexof_ps)%array(:,:,1) = sum(delp,3)
  deallocate(delp)
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
integer :: var, id_restart
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

  register = .true.

  select case (trim(fields(var)%short_name))
  case("u","v","ud","vd","ua","va","phis","T","DELP","W","DZ")
    filename = self%filename_core
    restart => restart_core
    read_core = .true.
  case("sphum","ice_wat","liq_wat","o3mr")
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
  case("qls","qcn","cfcn","frocean","frland", &
       "varflt","ustar","bstar","zpbl","cm", &
       "ct","cq","kcbl","ts","khl","khu")
    register = .false. !Not currently available from GFS, do nothing
  case default
    call abor1_ftn("write_gfs: filename not set for "//trim(fields(var)%short_name))
  end select

  if (register) &
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
   open(101, file=trim(adjustl(self%datapath_ti))//trim(adjustl(self%filename_cplr)), form='formatted')
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
