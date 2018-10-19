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

!fckit mpi
use fckit_mpi_module, only : fckit_mpi_comm

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
type(fv3jedi_geom), intent(inout)   :: geom
type(fv3jedi_increment), intent(in) :: incr      !< incr
type(c_ptr), intent(in)             :: c_conf   !< Configuration
type(datetime), intent(inout)       :: vdate    !< DateTime


character(len=255) :: datapath
character(len=255) :: filename
character(len=64)  :: datefile


integer :: geostiledim = 0

integer, parameter :: nvar = 10
integer :: ncid, varid(nvar)

integer :: date(6)
integer(kind=c_int) :: idate, isecs

integer :: isc,iec,jsc,jec,im,jm,km
integer :: x_dimid, y_dimid, z_dimid
integer :: t_dimid, tile_dimid
integer, allocatable :: dimids2(:), dimids3(:)
integer, allocatable :: istart2(:), icount2(:)
integer, allocatable :: istart3(:), icount3(:)

integer :: writeprec
type(fckit_mpi_comm) :: f_comm

 !> Convenience
 !> -----------
 isc = incr%isc
 iec = incr%iec
 jsc = incr%jsc
 jec = incr%jec

 f_comm = fckit_mpi_comm()

 ! Place to save restarts
 ! ----------------------
 datapath = "Data/"
 if (config_element_exists(c_conf,"datapath")) then
    datapath = config_get_string(c_conf,len(datapath),"datapath")
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
 

 !Using tile as a dimension in the file?
 if (config_element_exists(c_conf,"geos_tile_dim")) then
    geostiledim = config_get_int(c_conf,"geos_tile_dim")
 endif

 ! Naming convection for the file
 ! ------------------------------
 filename = 'GEOS.eta.'

 if (config_element_exists(c_conf,"filename")) then
    filename = config_get_string(c_conf,len(filename),"filename")
 endif

 !Append with the date
 write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)
 filename = trim(datapath)//trim(filename)//trim(datefile)//trim("z.nc4")

 writeprec = NF90_FLOAT

 varid = 0

 if (f_comm%rank() == 0) then

   im = geom%npx - 1
   if (geostiledim == 0) then
     jm = 6*(geom%npy - 1)
   else
     jm = geom%npy - 1
   endif
   km = geom%npz
  
   !Create file to write to
   call nccheck ( nf90_create(trim(filename), NF90_CLOBBER, ncid), "nf90_create" )

   !Create dimensions
   call nccheck ( nf90_def_dim(ncid, "lon", im, x_dimid), "nf90_def_dim lon" )
   call nccheck ( nf90_def_dim(ncid, "lat", jm, y_dimid), "nf90_def_dim lat" )
   call nccheck ( nf90_def_dim(ncid, "lev", km, z_dimid), "nf90_def_dim lev" )
   call nccheck ( nf90_def_dim(ncid, "time", 1, t_dimid), "nf90_def_dim time" )
  
   !Add further dimension for the tile number if requested and create dimids
   if (geostiledim == 1) then
  
     call nccheck ( nf90_def_dim(ncid, "tile", geom%ntiles, tile_dimid), "nf90_def_dim tile"  )
  
     allocate(dimids2(4))
     dimids2 =  (/ x_dimid, y_dimid, t_dimid, tile_dimid /)
     allocate(dimids3(5))
     dimids3 =  (/ x_dimid, y_dimid, z_dimid, t_dimid, tile_dimid /)
  
   else
  
     allocate(dimids2(2))
     dimids2 =  (/ x_dimid, y_dimid /)
     allocate(dimids3(4))
     dimids3 =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)
  
   endif
  
   !Define variables to be written (geom)
print*, dimids2
   call nccheck( nf90_def_var(ncid, "lon", NF90_DOUBLE, dimids2, varid(1)), "nf90_def_var lon" )
   call nccheck( nf90_def_var(ncid, "lat", NF90_DOUBLE, dimids2, varid(2)), "nf90_def_var lat" )

   !Define variables to be written (increment)
   call nccheck( nf90_def_var(ncid, "ua" , writeprec  , dimids3, varid(3)), "nf90_def_var ua"  )

   !Close for this proc
   call nccheck( nf90_enddef(ncid), "nf90_enddef" )
   call nccheck( nf90_close(ncid), "nf90_close" )
  
   deallocate(dimids3,dimids2)

 endif !Root process


 !Broadcast the varids
 call f_comm%broadcast(varid,0)


 !Create local to this proc start/count
 if (geostiledim == 1) then

   allocate(istart3(5),istart2(4))
   allocate(icount3(5),icount2(4))

   istart3(1) = isc
   istart3(2) = jsc
   istart3(3) = geom%ntile
   istart3(4) = 1
   istart3(5) = 1
   icount3(1) = iec-isc+1
   icount3(2) = jec-jsc+1
   icount3(3) = 1
   icount3(4) = geom%npz
   icount3(5) = 1

   istart2(1) = isc
   istart2(2) = jsc
   istart2(3) = geom%ntile
   istart2(4) = 1
   icount2(1) = iec-isc+1
   icount2(2) = jec-jsc+1
   icount2(3) = 1
   icount2(4) = 1   

 else

   allocate(istart3(4),istart2(2))
   allocate(icount3(4),icount2(2))

   istart3(1) = isc
   istart3(2) = (geom%ntile-1)*(jm/geom%ntiles) + jsc
   istart3(3) = 1
   istart3(4) = 1
   icount3(1) = iec-isc+1
   icount3(2) = jec-jsc+1
   icount3(3) = geom%npz
   icount3(4) = 1

   istart2(1) = isc
   istart2(2) = (geom%ntile-1)*(jm/geom%ntiles) + jsc
!   istart2(3) = 1
   icount2(1) = iec-isc+1
   icount2(2) = jec-jsc+1
!   icount2(3) = 1

 endif

 call nccheck( nf90_open(trim(filename), NF90_WRITE, ncid), "nf90_open" )

print*, icount2
 call nccheck( nf90_put_var(ncid, varid(1), geom%grid_lon(isc:iec,jsc:jec), istart2, icount2), "nf90_put_var lon" )


! call nccheck( nf90_put_var(ncid, varid(2), geom%grid_lat(isc:iec,jsc:jec), istart2, icount2), "nf90_put_var lat"  )
! call nccheck( nf90_put_var(ncid, varid(3), incr%ua(isc:iec,jsc:jec,:)    , istart3, icount3), "nf90_put_var ua"  )

 call nccheck( nf90_close(ncid), "nf90_close" )

 deallocate(istart2,istart3)
 deallocate(icount2,icount3)

 call abor1_ftn("done")

end subroutine write_geos_restart

! ------------------------------------------------------------------------------

subroutine nccheck(status,iam)

implicit none
integer, intent ( in) :: status
character(len=*), optional :: iam 

character(len=255) :: error_descr

 if(status /= nf90_noerr) then

   error_descr = "fv3jedi_increment_io_mod: NetCDF error, aborting"

   if (present(iam)) then
     error_descr = trim(error_descr)//", "//trim(iam)
   endif

   error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

   call abor1_ftn(trim(error_descr))

 end if

end subroutine nccheck

! ------------------------------------------------------------------------------

end module fv3jedi_increment_io_mod
