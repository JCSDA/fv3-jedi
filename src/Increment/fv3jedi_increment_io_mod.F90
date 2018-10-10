module fv3jedi_increment_io_mod

use config_mod
use iso_c_binding
use datetime_mod
use fckit_log_module, only : log

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds, only: kind_real
use fv3jedi_increment_utils_mod, only: fv3jedi_increment

!For FMS like restarts
use mpp_domains_mod,   only: EAST, NORTH
use fms_io_mod,        only: restart_file_type, register_restart_field, &
                             free_restart_type, restore_state, save_restart

!For GEOS like restarts
use netcdf

implicit none
private
public read_fms_restart, write_fms_restart, &
       read_geos_restart, write_geos_restart

contains

! ------------------------------------------------------------------------------

subroutine read_fms_restart(geom, incr, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: incr      !< incr
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine read_fms_restart

! ------------------------------------------------------------------------------

subroutine write_fms_restart(geom, incr, c_conf, vdate)

use mpp_mod, only: mpp_pe, mpp_root_pe

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_increment), intent(in)    :: incr      !< incr
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime

!Locals
type(restart_file_type) :: Fv_restart
type(restart_file_type) :: Tr_restart

integer :: id_restart, iounit, io_status
integer :: date_init(6), date(6), calendar_type, layout(2)
integer(kind=c_int) :: idate, isecs

character(len=64) :: filename_core
character(len=64) :: filename_trcr
character(len=64) :: filename_cplr
character(len=64) :: datefile

character(len=20) :: sdate,validitydate

character(len=255) :: datapath_out


 ! Place to save restarts
 ! ----------------------
 datapath_out = "Data/"
 if (config_element_exists(c_conf,"datapath_write")) then
    datapath_out = config_get_string(c_conf,len(datapath_out),"datapath_write")
 endif


 ! Current date
 ! ------------
 call datetime_to_ifs(vdate, idate, isecs)
 date(1) = idate/10000
 date(2) = idate/100 - date(1)*100
 date(3) = idate - (date(1)*10000 + date(2)*100)
 date(4) = isecs/3600
 date(5) = (isecs - date(4)*3600)/60
 date(6) = isecs - (date(4)*3600 + date(5)*60)
 

 ! Naming convection for the file
 ! ------------------------------
 filename_core = 'fv_core.res.nc'
 if (config_element_exists(c_conf,"filename_core")) then
   filename_core = config_get_string(c_conf,len(filename_core),"filename_core")
 endif

 filename_trcr = 'fv_tracer.res.nc'
 if (config_element_exists(c_conf,"filename_trcr")) then
   filename_trcr = config_get_string(c_conf,len(filename_trcr),"filename_trcr")
 endif

 filename_cplr = 'coupler.res'
 if (config_element_exists(c_conf,"filename_cplr")) then
   filename_cplr = config_get_string(c_conf,len(filename_cplr),"filename_cplr")
 endif

 !Append with the date
 write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2,A1)') date(1),date(2),date(3),".",date(4),date(5),date(6),"."
 filename_core = trim(datefile)//trim(filename_core)
 filename_trcr = trim(datefile)//trim(filename_trcr)
 filename_cplr = trim(datefile)//trim(filename_cplr)


 ! Register the variables that should be written
 ! ---------------------------------------------

 !A-Grid winds, increment
 id_restart =  register_restart_field(Fv_restart, filename_core, 'ua', incr%ua, &
                                      domain=geom%domain )
 id_restart =  register_restart_field(Fv_restart, filename_core, 'va', incr%va, &
                                      domain=geom%domain )

 !Temperature
 id_restart = register_restart_field( Fv_restart, filename_core, 'T', incr%t, &
                                      domain=geom%domain )

 !Pressure thickness
 id_restart = register_restart_field( Fv_restart, filename_core, 'Ps', incr%ps, &
                                      domain=geom%domain )

 !Nonhydrostatic state
 if (.not. incr%hydrostatic) then
     id_restart =  register_restart_field( Fv_restart, filename_core, 'W', incr%w, &
                                           domain=geom%domain )
     id_restart =  register_restart_field( Fv_restart, filename_core, 'DZ', incr%delz, &
                                           domain=geom%domain )
 endif

 !Cell center lat/lon
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lat', geom%grid_lat, &
                                      domain=geom%domain )
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lon', geom%grid_lon, &
                                      domain=geom%domain )

 ! Write variables to file
 ! -----------------------`/
 call save_restart(Fv_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Fv_restart)


 !Write tracers to file
 !---------------------
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'sphum'  , incr%q, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'ice_wat', incr%qi, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'liq_wat', incr%ql, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'o3mr'   , incr%o3, &
                                       domain=geom%domain )

 call save_restart(Tr_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Tr_restart)


 !Write date/time info in coupler.res
 !-----------------------------------
 iounit = 101
 if (mpp_pe() == mpp_root_pe()) then
    print *,'write_file: date model init = ',incr%date_init
    print *,'write_file: date model now  = ',incr%date
    print *,'write_file: date vdate      = ',date
    open(iounit, file=trim(adjustl(datapath_out))//'RESTART/'//trim(adjustl(filename_cplr)), form='formatted')
    write( iounit, '(i6,8x,a)' ) incr%calendar_type, &
         '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
    write( iounit, '(6i6,8x,a)' )date, &
          'Model start time:   year, month, day, hour, minute, second'
    write( iounit, '(6i6,8x,a)' )date, &
          'Current model time: year, month, day, hour, minute, second'
    close(iounit)
 endif

 return

end subroutine write_fms_restart

! ------------------------------------------------------------------------------

subroutine read_geos_restart(geom, incr, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_increment), intent(inout) :: incr    !< incr
type(c_ptr),         intent(in)    :: c_conf   !< Configuration
type(datetime),      intent(inout) :: vdate    !< DateTime


end subroutine read_geos_restart

! ------------------------------------------------------------------------------

subroutine write_geos_restart(geom, incr, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_increment), intent(in)    :: incr      !< incr
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine write_geos_restart

! ------------------------------------------------------------------------------

end module fv3jedi_increment_io_mod
