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

subroutine read_fms_restart(geom, fld, c_conf, vdate)

use iso_c_binding
use datetime_mod

use mpp_domains_mod,     only: mpp_update_domains, DGRID_NE
use field_manager_mod,   only: MODEL_ATMOS
use tracer_manager_mod,  only: get_number_tracers, get_tracer_names, set_tracer_profile, get_tracer_index

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(inout) :: fld      !< State
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
 id_restart = register_restart_field(Fv_restart, filename_core, 'u', fld%ud, &
              domain=geom%domain, position=NORTH)
 id_restart = register_restart_field(Fv_restart, filename_core, 'v', fld%vd, &
              domain=geom%domain, position=EAST)

 !A-Grid winds, increment
 id_restart = register_restart_field(Fv_restart, filename_core, 'ua', fld%ua, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Fv_restart, filename_core, 'va', fld%va, &
                                     domain=geom%domain)

 !phis
 id_restart = register_restart_field(Fv_restart, filename_core, 'phis', fld%phis, &
              domain=geom%domain)

 !Temperature
 id_restart = register_restart_field(Fv_restart, filename_core, 'T', fld%t, &
              domain=geom%domain)

 !Pressure thickness
 id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', fld%delp, &
              domain=geom%domain)

 !Nonhydrostatic variables
 if (.not. fld%hydrostatic) then
    id_restart =  register_restart_field(Fv_restart, filename_core, 'W', fld%w, &
                  domain=geom%domain)
    id_restart =  register_restart_field(Fv_restart, filename_core, 'DZ', fld%delz, &
                  domain=geom%domain)
 endif

 ! Read file and fill variables
 ! ----------------------------
 call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Fv_restart)
 

 !Register and read tracers
 !-------------------------
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'sphum'  , fld%q , &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'ice_wat', fld%qi, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'liq_wat', fld%ql, &
                                     domain=geom%domain)
 id_restart = register_restart_field(Tr_restart, filename_trcr, 'o3mr'   , fld%o3, &
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
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'slmsk' , fld%slmsk , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'sheleg', fld%sheleg, domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'tsea'  , fld%tsea  , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vtype' , fld%vtype , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stype' , fld%stype , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vfrac' , fld%vfrac , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stc'   , fld%stc   , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'smc'   , fld%smc   , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'snwdph', fld%snwdph, domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'f10m'  , fld%f10m  , domain=geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)

   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'u_srf' , fld%u_srf , domain=geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'v_srf' , fld%v_srf , domain=geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)
   fld%havecrtmfields = .true.
 else
   fld%havecrtmfields = .false.
   fld%slmsk  = 0.0_kind_real
   fld%sheleg = 0.0_kind_real
   fld%tsea   = 0.0_kind_real
   fld%vtype  = 0.0_kind_real
   fld%stype  = 0.0_kind_real
   fld%vfrac  = 0.0_kind_real
   fld%stc    = 0.0_kind_real
   fld%smc    = 0.0_kind_real
   fld%u_srf  = 0.0_kind_real
   fld%u_srf  = 0.0_kind_real
   fld%v_srf  = 0.0_kind_real
   fld%f10m   = 0.0_kind_real
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
 fld%date = date
 fld%date_init = date_init
 fld%calendar_type = calendar_type
 idate=date(1)*10000+date(2)*100+date(3)
 isecs=date(4)*3600+date(5)*60+date(6)

 call datetime_from_ifs(vdate, idate, isecs)
 call datetime_to_string(vdate, validitydate)

 call log%info("read_file: validity date: "//trim(validitydate)) 
 call log%info("read_file: expected validity date: "//trim(sdate)) 

 return

end subroutine read_fms_restart

! ------------------------------------------------------------------------------

subroutine write_fms_restart(geom, fld, c_conf, vdate)

use mpp_mod, only: mpp_pe, mpp_root_pe

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(in)    :: fld      !< State
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
 id_restart = register_restart_field( Fv_restart, filename_core, 'u', fld%ud, &
                                      domain=geom%domain,position=NORTH )
 id_restart = register_restart_field( Fv_restart, filename_core, 'v', fld%vd, &
                                      domain=geom%domain,position=EAST )

 !A-Grid winds, increment
 id_restart =  register_restart_field(Fv_restart, filename_core, 'ua', fld%ua, &
                                      domain=geom%domain )
 id_restart =  register_restart_field(Fv_restart, filename_core, 'va', fld%va, &
                                      domain=geom%domain )

 !phis
 id_restart = register_restart_field( Fv_restart, filename_core, 'phis', fld%phis, &
                                      domain=geom%domain )

 !Temperature
 id_restart = register_restart_field( Fv_restart, filename_core, 'T', fld%t, &
                                      domain=geom%domain )

 !Pressure thickness
 id_restart = register_restart_field( Fv_restart, filename_core, 'DELP', fld%delp, &
                                      domain=geom%domain )

 !Nonhydrostatic state
 if (.not. fld%hydrostatic) then
     id_restart =  register_restart_field( Fv_restart, filename_core, 'W', fld%w, &
                                           domain=geom%domain )
     id_restart =  register_restart_field( Fv_restart, filename_core, 'DZ', fld%delz, &
                                           domain=geom%domain )
 endif

 !Cell center lat/lon
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lat', geom%grid_lat, &
                                      domain=geom%domain )
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lon', geom%grid_lon, &
                                      domain=geom%domain )

 id_restart =  register_restart_field( Fv_restart, filename_core, 'ua', fld%ua, &
                                       domain=geom%domain )
 id_restart =  register_restart_field( Fv_restart, filename_core, 'va', fld%va, &
                                       domain=geom%domain )

 ! Write variables to file
 ! -----------------------`/
 call save_restart(Fv_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Fv_restart)


 !Write tracers to file
 !---------------------
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'sphum'  , fld%q, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'ice_wat', fld%qi, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'liq_wat', fld%ql, &
                                       domain=geom%domain )
  id_restart =  register_restart_field( Tr_restart, filename_trcr, 'o3mr'   , fld%o3, &
                                       domain=geom%domain )

 call save_restart(Tr_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Tr_restart)


 !Write date/time info in coupler.res
 !-----------------------------------
 iounit = 101
 if (mpp_pe() == mpp_root_pe()) then
    print *,'write_file: date model init = ',fld%date_init
    print *,'write_file: date model now  = ',fld%date
    print *,'write_file: date vdate      = ',date
    open(iounit, file=trim(adjustl(datapath_out))//'RESTART/'//trim(adjustl(filename_cplr)), form='formatted')
    write( iounit, '(i6,8x,a)' ) fld%calendar_type, &
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

subroutine read_geos_restart(geom, fld, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(inout) :: fld      !< State
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime

character(len=255) :: datapath
character(len=255) :: filename_eta
character(len=255) :: filename_sfc

integer :: ncid, ncstat, dimid, varid

integer :: im, jm, lm, nm, l

integer :: date(6)
integer :: intdate, inttime
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: idate, isecs
character(len=20) :: sdate, validitydate

real(kind=kind_real), allocatable :: state2d_global(:,:)

integer, allocatable :: istart(:), icount(:)

integer :: tileoff
logical :: tiledimension = .false.

integer :: isc,iec,jsc,jec


 !> Convenience
 !> -----------
 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec


 !> Set filenames
 !> -------------
 filename_eta = 'GEOS.bkg.eta.nc4'
 filename_sfc = 'GEOS.sfc.sfc.nc4'

 if (config_element_exists(c_conf,"filename_eta")) then
    filename_eta = config_get_string(c_conf,len(filename_eta),"filename_eta")
 endif
 if (config_element_exists(c_conf,"filename_sfc")) then
    filename_sfc = config_get_string(c_conf,len(filename_sfc),"filename_sfc")
 endif

 datapath = config_get_string(c_conf,len(datapath),"datapath_read")

 filename_eta  = trim(datapath)//trim("/")//trim(filename_eta )
 filename_sfc = trim(datapath)//trim("/")//trim(filename_sfc)

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
 if ( im /= geom%npx-1 .or. lm /= geom%npz) then
   call abor1_ftn("GEOS restarts: restart dimension not compatible with geometry")
 endif

 !> GEOS can use concatenated tiles or tile as a dimension
 if ( (im == geom%npx-1) .and. (jm == 6*(geom%npy-1) ) ) then
   tiledimension = .false.
   tileoff = (geom%ntile-1)*(jm/geom%ntiles)
 else
   tiledimension = .true.
   tileoff = 0
   call abor1_ftn("GEOS restarts: tile dimension in file not done yet")
 endif


 !> Read the state level by level
 !> -----------------------------
 allocate(state2d_global(im,jm))

 !> starts and counts
 deallocate(istart,icount)
 allocate(istart(4))
 allocate(icount(4))
 istart = 1
 icount(1) = im
 icount(2) = jm
 icount(3) = 1
 icount(4) = 1

 !> Loop over levels
 do l = 1,lm

   istart(3) = l

   !fld%ud (D-Grid winds, nonlinear model only)
   ncstat = nf90_inq_varid (ncid, 'u', varid)
   if(ncstat /= nf90_noerr) print *, "u: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "u: "//trim(nf90_strerror(ncstat))
   fld%ud(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%vd (D-Grid winds, nonlinear model only)
   ncstat = nf90_inq_varid (ncid, 'v', varid)
   if(ncstat /= nf90_noerr) print *, "v: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "v: "//trim(nf90_strerror(ncstat))
   fld%vd(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%ua (A-Grid winds, increment)
   ncstat = nf90_inq_varid (ncid, 'ua', varid)
   if(ncstat /= nf90_noerr) print *, "ua: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "ua: "//trim(nf90_strerror(ncstat))
   fld%ua(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%va (A-Grid winds, increment)
   ncstat = nf90_inq_varid (ncid, 'va', varid)
   if(ncstat /= nf90_noerr) print *, "va: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "va: "//trim(nf90_strerror(ncstat))
   fld%va(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%t
   ncstat = nf90_inq_varid (ncid, 'T', varid)
   if(ncstat /= nf90_noerr) print *, "T: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "T: "//trim(nf90_strerror(ncstat))
   fld%t(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%delp
   ncstat = nf90_inq_varid (ncid, 'delp', varid)
   if(ncstat /= nf90_noerr) print *, "delp: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "delp: "//trim(nf90_strerror(ncstat))
   fld%delp(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%q (sphum)
   ncstat = nf90_inq_varid (ncid, 'sphum', varid)
   if(ncstat /= nf90_noerr) print *, "sphum: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "sphum: "//trim(nf90_strerror(ncstat))
   fld%q(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%qi (liq_wat)
   ncstat = nf90_inq_varid (ncid, 'liq_wat', varid)
   if(ncstat /= nf90_noerr) print *, "liq_wat: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "liq_wat: "//trim(nf90_strerror(ncstat))
   fld%qi(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%ql (ice_wat)
   ncstat = nf90_inq_varid (ncid, 'ice_wat', varid)
   if(ncstat /= nf90_noerr) print *, "ice_wat: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "ice_wat: "//trim(nf90_strerror(ncstat))
   fld%ql(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

   !fld%o3 (o3mr)
   ncstat = nf90_inq_varid (ncid, 'o3mr', varid)
   if(ncstat /= nf90_noerr) print *, "o3mr: "//trim(nf90_strerror(ncstat))
   ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
   if(ncstat /= nf90_noerr) print *, "o3mr: "//trim(nf90_strerror(ncstat))
   fld%o3(isc:iec,jsc:jec,l) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

 enddo

 deallocate(istart,icount)

 !Two dimensional variables
 
 allocate(istart(3))
 allocate(icount(3))
 
 istart = 1
 icount(1) = im
 icount(2) = jm
 icount(3) = 1
 
 !fld%phis
 ncstat = nf90_inq_varid (ncid, 'phis', varid)
 if(ncstat /= nf90_noerr) print *, "phis: "//trim(nf90_strerror(ncstat))
 ncstat = nf90_get_var(ncid, varid, state2d_global, istart, icount)
 if(ncstat /= nf90_noerr) print *, "phis: "//trim(nf90_strerror(ncstat))
 fld%phis(isc:iec,jsc:jec) = state2d_global(isc:iec,tileoff+jsc:tileoff+jec)

 !Close this file
 ncstat = nf90_close(ncid)
 if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))


 !TODO: ADD THE SURFACE VARIABLES AND TRANSFORM THEM TO FV3FGFS STYLE


 deallocate(state2d_global)
 deallocate(istart,icount)


end subroutine read_geos_restart

! ------------------------------------------------------------------------------

subroutine write_geos_restart(geom, fld, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(in)    :: fld      !< State
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine write_geos_restart

! ------------------------------------------------------------------------------

end module fv3jedi_state_io_mod
