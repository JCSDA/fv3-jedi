module fv3jedi_state_io_mod

use config_mod
use iso_c_binding
use datetime_mod
use fckit_log_module, only : log

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds, only: kind_real
use fv3jedi_state_utils_mod, only: fv3jedi_state

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

subroutine read_fms_restart(geom, state, c_conf, vdate)

use iso_c_binding
use datetime_mod

use mpp_domains_mod,     only: mpp_update_domains, DGRID_NE
use field_manager_mod,   only: MODEL_ATMOS
use tracer_manager_mod,  only: get_number_tracers, get_tracer_names, set_tracer_profile, get_tracer_index

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(inout) :: state      !< State
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime

!Locals
type(restart_file_type) :: Fv_restart
type(restart_file_type) :: Tr_restart
type(restart_file_type) :: Sf_restart

integer :: id_restart, iounit, io_status, layout(2)
integer :: date_init(6), date(6), calendar_type
integer(kind=c_int) :: idate, isecs

character(len=255) :: datapath_in, datapath_ti
character(len=255) :: filename_core
character(len=255) :: filename_trcr
character(len=255) :: filename_sfcd
character(len=255) :: filename_sfcw
character(len=255) :: filename_cplr

character(len=20) :: sdate,validitydate
character(len=1024)  :: buf

integer :: k

character(len=64):: tracer_name
integer :: ntracers, ntprog, nt, ierr

integer :: print_read_info = 0

integer :: read_crtm_surface

 !Set filenames
 !--------------
 filename_core = 'fv_core.res.nc'
 filename_trcr = 'fv_tracer.res.nc'
 filename_sfcd = 'sfc_data.nc'
 filename_sfcw = 'srf_wnd.nc'
 filename_cplr = 'coupler.res'

 if (config_element_exists(c_conf,"filename_core")) then
    filename_core = config_get_string(c_conf,len(filename_core),"filename_core")
 endif
 if (config_element_exists(c_conf,"filename_trcr")) then
    filename_trcr = config_get_string(c_conf,len(filename_trcr),"filename_trcr")
 endif
 if (config_element_exists(c_conf,"filename_sfcd")) then
    filename_sfcd = config_get_string(c_conf,len(filename_sfcd),"filename_sfcd")
 endif
 if (config_element_exists(c_conf,"filename_sfcw")) then
    filename_sfcw = config_get_string(c_conf,len(filename_sfcw),"filename_sfcw")
 endif
 if (config_element_exists(c_conf,"filename_cplr")) then
    filename_cplr = config_get_string(c_conf,len(filename_cplr),"filename_cplr")
 endif

 datapath_in = config_get_string(c_conf,len(datapath_in),"datapath_read")
 datapath_ti = config_get_string(c_conf,len(datapath_ti),"datapath_tile")

 ! Register the variables that should be read
 ! ------------------------------------------
 !D-Grid winds, nonlinear model only
 id_restart = register_restart_field(Fv_restart, filename_core, 'u', state%ud, &
              domain=geom%domain, position=NORTH)
 id_restart = register_restart_field(Fv_restart, filename_core, 'v', state%vd, &
              domain=geom%domain, position=EAST)

 !A-Grid winds, increment
 id_restart = register_restart_field(Fv_restart, filename_core, 'ua', state%ua, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Fv_restart, filename_core, 'va', state%va, &
                                     domain=geom%domain)

 !phis
 id_restart = register_restart_field(Fv_restart, filename_core, 'phis', state%phis, &
              domain=geom%domain)

 !Temperature
 id_restart = register_restart_field(Fv_restart, filename_core, 'T', state%t, &
              domain=geom%domain)

 !Pressure thickness
 id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', state%delp, &
              domain=geom%domain)

 !Nonhydrostatic variables
 if (.not. state%hydrostatic) then
    id_restart =  register_restart_field(Fv_restart, filename_core, 'W', state%w, &
                  domain=geom%domain)
    id_restart =  register_restart_field(Fv_restart, filename_core, 'DZ', state%delz, &
                  domain=geom%domain)
 endif

 ! Read file and fill variables
 ! ----------------------------
 call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Fv_restart)
 

 !Register and read tracers
 !-------------------------
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'sphum'  , state%q , &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'ice_wat', state%qi, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'liq_wat', state%ql, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'o3mr'   , state%o3, &
                                     domain=geom%domain)

 call restore_state(Tr_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Tr_restart)

 !Register and read surface state needed for crtm calculation
 !------------------------------------------------------------
 read_crtm_surface = 0
 if (config_element_exists(c_conf,"read_crtm_surface")) then
   read_crtm_surface = config_get_int(c_conf,"read_crtm_surface")
 endif

 if (read_crtm_surface == 1) then
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'slmsk' , state%slmsk , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'sheleg', state%sheleg, domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'tsea'  , state%tsea  , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vtype' , state%vtype , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stype' , state%stype , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vfrac' , state%vfrac , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stc'   , state%stc   , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'smc'   , state%smc   , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'snwdph', state%snwdph, domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'f10m'  , state%f10m  , domain=geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)

   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'u_srf' , state%u_srf , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'v_srf' , state%v_srf , domain=geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)
   state%havecrtmfields = .true.
 else
   state%havecrtmfields = .false.
   state%slmsk  = 0.0_kind_real
   state%sheleg = 0.0_kind_real
   state%tsea   = 0.0_kind_real
   state%vtype  = 0.0_kind_real
   state%stype  = 0.0_kind_real
   state%vfrac  = 0.0_kind_real
   state%stc    = 0.0_kind_real
   state%smc    = 0.0_kind_real
   state%u_srf  = 0.0_kind_real
   state%u_srf  = 0.0_kind_real
   state%v_srf  = 0.0_kind_real
   state%f10m   = 0.0_kind_real
 endif

 ! Get dates from file
 !--------------------

 ! read date from coupler.res text file.
 sdate = config_get_string(c_conf,len(sdate),"date")
 iounit = 101
 open(iounit, file=trim(adjustl(datapath_in))//filename_cplr, form='formatted')
 read(iounit, '(i6)')  calendar_type
 read(iounit, '(6i6)') date_init
 read(iounit, '(6i6)') date
 close(iounit)
 state%date = date
 state%date_init = date_init
 state%calendar_type = calendar_type
 idate=date(1)*10000+date(2)*100+date(3)
 isecs=date(4)*3600+date(5)*60+date(6)

 call datetime_from_ifs(vdate, idate, isecs)
 call datetime_to_string(vdate, validitydate)

 call log%info("read_file: validity date: "//trim(validitydate)) 
 call log%info("read_file: expected validity date: "//trim(sdate)) 

 return

end subroutine read_fms_restart

! ------------------------------------------------------------------------------

subroutine write_fms_restart(geom, state, c_conf, vdate)

use mpp_mod, only: mpp_pe, mpp_root_pe

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(in)    :: state      !< State
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
 !D-Grid winds, nonlinear model only
 id_restart = register_restart_field( Fv_restart, filename_core, 'u', state%ud, &
                                      domain=geom%domain,position=NORTH )
 id_restart = register_restart_field( Fv_restart, filename_core, 'v', state%vd, &
                                      domain=geom%domain,position=EAST )

 !A-Grid winds, increment
 id_restart =  register_restart_field(Fv_restart, filename_core, 'ua', state%ua, &
                                      domain=geom%domain )
 id_restart =  register_restart_field(Fv_restart, filename_core, 'va', state%va, &
                                      domain=geom%domain )

 !phis
 id_restart = register_restart_field( Fv_restart, filename_core, 'phis', state%phis, &
                                      domain=geom%domain )

 !Temperature
 id_restart = register_restart_field( Fv_restart, filename_core, 'T', state%t, &
                                      domain=geom%domain )

 !Pressure thickness
 id_restart = register_restart_field( Fv_restart, filename_core, 'DELP', state%delp, &
                                      domain=geom%domain )

 !Nonhydrostatic state
 if (.not. state%hydrostatic) then
     id_restart =  register_restart_field( Fv_restart, filename_core, 'W', state%w, &
                                           domain=geom%domain )
     id_restart =  register_restart_field( Fv_restart, filename_core, 'DZ', state%delz, &
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
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'sphum'  , state%q, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'ice_wat', state%qi, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'liq_wat', state%ql, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'o3mr'   , state%o3, &
                                       domain=geom%domain )

 call save_restart(Tr_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Tr_restart)


 !Write date/time info in coupler.res
 !-----------------------------------
 iounit = 101
 if (mpp_pe() == mpp_root_pe()) then
    print *,'write_file: date model init = ',state%date_init
    print *,'write_file: date model now  = ',state%date
    print *,'write_file: date vdate      = ',date
    open(iounit, file=trim(adjustl(datapath_out))//'RESTART/'//trim(adjustl(filename_cplr)), form='formatted')
    write( iounit, '(i6,8x,a)' ) state%calendar_type, &
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

subroutine read_geos_restart(state, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_state), intent(inout) :: state    !< State
type(c_ptr),         intent(in)    :: c_conf   !< Configuration
type(datetime),      intent(inout) :: vdate    !< DateTime

character(len=255) :: datapath
character(len=255) :: filename_eta

integer :: ncid, ncstat, dimid, varid

integer :: im, jm, lm, nm, l

integer :: date(6)
integer :: intdate, inttime
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: idate, isecs
character(len=20) :: sdate, validitydate

integer, allocatable :: istart(:), icount(:)

integer :: tileoff
logical :: tiledimension = .false.

integer :: isc,iec,jsc,jec

character(len=20)  :: var

 !> Convenience
 !> -----------
 isc = state%isc
 iec = state%iec
 jsc = state%jsc
 jec = state%jec


 !> Set filenames
 !> -------------
 filename_eta = 'GEOS.bkg.eta.nc4'

 if (config_element_exists(c_conf,"filename_eta")) then
    filename_eta = config_get_string(c_conf,len(filename_eta),"filename_eta")
 endif

 datapath = config_get_string(c_conf,len(datapath),"datapath_read")

 filename_eta  = trim(datapath)//trim("/")//trim(filename_eta )

 !> Open the file
 ncstat = nf90_open(filename_eta, NF90_NOWRITE, ncid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

 !> Get dimensions, lon,lat,lev,time
 ncstat = nf90_inq_dimid(ncid, "lon", dimid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
 ncstat = nf90_inquire_dimension(ncid, dimid, len = im)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

 ncstat = nf90_inq_dimid(ncid, "lat", dimid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
 ncstat = nf90_inquire_dimension(ncid, dimid, len = jm)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

 ncstat = nf90_inq_dimid(ncid, "lev", dimid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
 ncstat = nf90_inquire_dimension(ncid, dimid, len = lm)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

 ncstat = nf90_inq_dimid(ncid, "time", dimid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
 ncstat = nf90_inquire_dimension(ncid, dimid, len = nm)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))


 !> Read the time
 !> -------------
 allocate(istart(1))
 allocate(icount(1))
 istart = 1
 icount = 1

 !> Get time attributes
 ncstat = nf90_inq_varid(ncid, "time", varid)
 if(ncstat /= nf90_noerr) print *, "time: "//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_att(ncid, varid, "begin_date", intdate)
 if(ncstat /= nf90_noerr) print *, "time: "//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_att(ncid, varid, "begin_time", inttime)
 if(ncstat /= nf90_noerr) print *, "time: "//trim(nf90_strerror(ncstat))

 !> Pad with leading zeros if need be
 write(cdate,"(I0.8)") intdate
 write(ctime,"(I0.6)") inttime

 !> Back to integer
 read(cdate(1:4),*) date(1)
 read(cdate(5:6),*) date(2)
 read(cdate(7:8),*) date(3)
 read(ctime(1:2),*) date(4)
 read(ctime(3:4),*) date(5)
 read(ctime(5:6),*) date(6)

 !> To idate/isecs for Jedi
 idate = date(1)*10000 + date(2)*100 + date(3)
 isecs = date(4)*3600  + date(5)*60  + date(6)

 call datetime_from_ifs(vdate, idate, isecs)
 call datetime_to_string(vdate, validitydate)

 !> Print info to user
 sdate = config_get_string(c_conf,len(sdate),"date")
 call log%info("read_file: validity date: "//trim(validitydate)) 
 call log%info("read_file: expected validity date: "//trim(sdate)) 

 !> Make sure file dimensions equal to geometry
 if ( im /= state%npx-1 .or. lm /= state%npz) then
   call abor1_ftn("GEOS restarts: restart dimension not compatible with geometry")
 endif

 !> GEOS can use concatenated tiles or tile as a dimension
 if ( (im == state%npx-1) .and. (jm == 6*(state%npy-1) ) ) then
   tiledimension = .false.
   tileoff = (state%ntile-1)*(jm/state%ntiles)
 else
   tiledimension = .true.
   tileoff = 0
   call abor1_ftn("GEOS restarts: tile dimension in file not done yet")
 endif


 !Rank three variables
 !--------------------
 deallocate(istart,icount)

 if (.not. tiledimension) then
   allocate(istart(4))
   allocate(icount(4))
   istart(1) = isc
   istart(2) = tileoff + jsc
   istart(3) = 1
   istart(4) = 1

   icount(1) = iec-isc+1
   icount(2) = jec-jsc+1
   icount(3) = state%npz
   icount(4) = 1
 else
   allocate(istart(5))
   allocate(icount(5))
   istart(1) = isc
   istart(2) = jsc
   istart(3) = state%ntile
   istart(4) = 1
   istart(5) = 1

   icount(1) = iec-isc+1
   icount(2) = jec-jsc+1
   icount(3) = 1
   icount(4) = state%npz
   icount(5) = 1
 endif

 var = 'ud'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ud(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'vd'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%vd(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'ua'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ua(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'va'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%va(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 't'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%t(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'delp'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%delp(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'q'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%q(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'qi'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%qi(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'ql'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ql(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'o3mr'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%o3(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'qls'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%qls(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'qcn'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%qcn(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'cfcn'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%cfcn(isc:iec,jsc:jec,:), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 !Rank two variables
 !------------------
 deallocate(istart,icount)

 if (.not. tiledimension) then
   allocate(istart(3))
   allocate(icount(3))
   istart(1) = isc
   istart(2) = tileoff + jsc
   istart(3) = 1

   icount(1) = iec-isc+1
   icount(2) = jec-jsc+1
   icount(3) = 1
 else
   allocate(istart(4))
   allocate(icount(4))
   istart(1) = isc
   istart(2) = jsc
   istart(3) = state%ntile
   istart(4) = 1

   icount(1) = iec-isc+1
   icount(2) = jec-jsc+1
   icount(3) = 1
   icount(4) = 1
 endif

 var = 'phis'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%phis(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'frland'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%frland(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'frocean'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%frocean(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'kcbl'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%kcbl(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'ts'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ts(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'khl'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%khl(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'khu'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%khu(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'varflt'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%varflt(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'ustar'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ustar(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'bstar'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%bstar(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'zpbl'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%zpbl(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'cm'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%cm(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'ct'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%ct(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

 var = 'cq'
 ncstat = nf90_inq_varid (ncid, trim(var), varid)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state%cq(isc:iec,jsc:jec), istart, icount)
 if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))


 !Close this file
 ncstat = nf90_close(ncid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

 deallocate(istart,icount)


end subroutine read_geos_restart

! ------------------------------------------------------------------------------

subroutine write_geos_restart(geom, state, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(in)    :: state      !< State
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine write_geos_restart

! ------------------------------------------------------------------------------

end module fv3jedi_state_io_mod
