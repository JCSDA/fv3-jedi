module fv3jedi_state_io_mod

use config_mod
use iso_c_binding
use datetime_mod
use fckit_log_module, only : log

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_state_utils_mod, only: fv3jedi_state

!For FMS like restarts
use mpp_domains_mod,   only: EAST, NORTH
use fms_io_mod,        only: restart_file_type, register_restart_field, &
                             free_restart_type, restore_state, save_restart

!For GEOS like restarts
use netcdf
use mpi
use fckit_mpi_module, only: fckit_mpi_comm

implicit none
private
public read_fms_state, write_fms_state, &
       read_geos_state, write_geos_state

contains

! ------------------------------------------------------------------------------

subroutine read_fms_state(geom, state, c_conf, vdate)

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

 !u on D-Grid
 if (allocated(state%ud)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'u', state%ud, &
                domain=geom%domain, position=NORTH)
 endif

 !v on D-Grid
 if (allocated(state%vd)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'v', state%vd, &
                domain=geom%domain, position=EAST)
 endif

 !u on A-Grid
 if (allocated(state%ua)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'ua', state%ua, &
                                       domain=geom%domain)
 endif

 !v on A-Grid
 if (allocated(state%va)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'va', state%va, &
                                       domain=geom%domain)
 endif

 !phis
 if (allocated(state%phis)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'phis', state%phis, &
                domain=geom%domain)
 endif

 !Temperature
 if (allocated(state%t)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'T', state%t, &
                domain=geom%domain)
 endif

 !Pressure thickness
 if (allocated(state%delp)) then
   id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', state%delp, &
                domain=geom%domain)
 endif

 !Nonhydrostatic variables - w
 if (allocated(state%w)) then
   id_restart =  register_restart_field(Fv_restart, filename_core, 'W', state%w, &
                 domain=geom%domain)
 endif

 !Nonhydrostatic variables - delz
 if (allocated(state%delz)) then
   id_restart =  register_restart_field(Fv_restart, filename_core, 'DZ', state%delz, &
                 domain=geom%domain)
 endif

 ! Read file and fill variables
 ! ----------------------------
 call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Fv_restart)
 

 !Register and read tracers
 !-------------------------
 if (allocated(state%q)) then
   id_restart = register_restart_field(Tr_restart, filename_trcr, 'sphum'  , state%q , &
                                       domain=geom%domain)
 endif
 if (allocated(state%qi)) then
   id_restart = register_restart_field(Tr_restart, filename_trcr, 'ice_wat', state%qi, &
                                       domain=geom%domain)
 endif
 if (allocated(state%ql)) then
   id_restart = register_restart_field(Tr_restart, filename_trcr, 'liq_wat', state%ql, &
                                       domain=geom%domain)
 endif
 if (allocated(state%o3)) then
   id_restart = register_restart_field(Tr_restart, filename_trcr, 'o3mr'   , state%o3, &
                                       domain=geom%domain)
 endif

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

! call log%info("read_file: validity date: "//trim(validitydate)) 
! call log%info("read_file: expected validity date: "//trim(sdate)) 

 return

end subroutine read_fms_state

! ------------------------------------------------------------------------------

subroutine write_fms_state(geom, state, c_conf, vdate)

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

end subroutine write_fms_state

! ------------------------------------------------------------------------------

subroutine read_geos_state(geom, state, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom),  intent(in)    :: geom
type(fv3jedi_state), intent(inout) :: state    !< State
type(c_ptr),         intent(in)    :: c_conf   !< Configuration
type(datetime),      intent(inout) :: vdate    !< DateTime

character(len=255) :: datapath
character(len=255) :: filename
integer :: ncid, ncstat, dimid, varid
integer :: im, jm, lm, nm, l
integer :: date(6)
integer :: intdate, inttime
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: idate, isecs
character(len=20) :: sdate, validitydate

integer, allocatable :: istart3(:), icount3(:)
integer, allocatable :: istart2(:), icount2(:)

integer :: read_tlad_traj
integer :: geostiledim, tileoff
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
 filename = 'GEOS.bkg.eta.nc4'

 if (config_element_exists(c_conf,"filename")) then
    filename = config_get_string(c_conf,len(filename),"filename")
 endif

 read_tlad_traj = 0
 if (config_element_exists(c_conf,"read_tlad_traj")) then
   read_tlad_traj = config_get_int(c_conf,"read_tlad_traj")
 endif

 !> Open the file
 call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, ncid), "nf90_open"//trim(filename) )

 !> Get dimensions, XDim,YDim,lev,time
 call nccheck ( nf90_inq_dimid(ncid, "Xdim", dimid), "nf90_inq_dimid Xdim" )
 call nccheck ( nf90_inquire_dimension(ncid, dimid, len = im), "nf90_inquire_dimension Xdim" )

 call nccheck ( nf90_inq_dimid(ncid, "Ydim", dimid), "nf90_inq_dimid YDim" )
 call nccheck ( nf90_inquire_dimension(ncid, dimid, len = jm), "nf90_inquire_dimension YDim" )

 call nccheck ( nf90_inq_dimid(ncid, "lev", dimid), "nf90_inq_dimid lev" )
 call nccheck ( nf90_inquire_dimension(ncid, dimid, len = lm), "nf90_inquire_dimension lev" )

 call nccheck ( nf90_inq_dimid(ncid, "time", dimid), "nf90_inq_dimid time" )
 call nccheck ( nf90_inquire_dimension(ncid, dimid, len = nm), "nf90_inquire_dimension time" )


 !> Get time attributes
 call nccheck ( nf90_inq_varid(ncid, "time", varid), "nf90_inq_varid time" )
 call nccheck ( nf90_get_att(ncid, varid, "begin_date", intdate), "nf90_get_att begin_date" )
 call nccheck ( nf90_get_att(ncid, varid, "begin_time", inttime), "nf90_get_att begin_time" )

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
! call log%info("read_file: validity date: "//trim(validitydate)) 
! call log%info("read_file: expected validity date: "//trim(sdate)) 

 !> Make sure file dimensions equal to geometry
 if ( im /= state%npx-1 .or. lm /= state%npz) then
   call abor1_ftn("GEOS restarts: restart dimension not compatible with geometry")
 endif

 !> GEOS can use concatenated tiles or tile as a dimension
 if ( (im == state%npx-1) .and. (jm == 6*(state%npy-1) ) ) then
   geostiledim = 0
   tileoff = (state%ntile-1)*(jm/state%ntiles)
 else
   geostiledim = 1
   tileoff = 0
 endif

 ! Create local to this proc start/count
 if (geostiledim == 1) then
 
   allocate(istart3(5),icount3(5))
   allocate(istart2(4),icount2(4))
 
   istart3(1) = isc;          icount3(1) = iec-isc+1
   istart3(2) = jsc;          icount3(2) = jec-jsc+1
   istart3(3) = geom%ntile;   icount3(3) = 1
   istart3(4) = 1;            icount3(4) = geom%npz
   istart3(5) = 1;            icount3(5) = 1
   
   istart2(1) = isc;          icount2(1) = iec-isc+1
   istart2(2) = jsc;          icount2(2) = jec-jsc+1
   istart2(3) = geom%ntile;   icount2(3) = 1
   istart2(4) = 1;            icount2(4) = 1
 
 else
 
   allocate(istart3(4),icount3(4))
   allocate(istart2(3),icount2(3))
 
   istart3(1) = isc;          icount3(1) = iec-isc+1
   istart3(2) = tileoff+jsc;  icount3(2) = jec-jsc+1
   istart3(3) = 1;            icount3(3) = geom%npz
   istart3(4) = 1;            icount3(4) = 1
   
   istart2(1) = isc;          icount2(1) = iec-isc+1
   istart2(2) = tileoff+jsc;  icount2(2) = jec-jsc+1
   istart2(3) = 1;            icount2(3) = 1
 
 endif

 var = 'ud'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%ud(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'vd'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%vd(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'ua'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%ua(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'va'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%va(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 't'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%t(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'delp'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%delp(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'q'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%q(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'qi'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%qi(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'ql'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%ql(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'o3mr'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%o3(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )

 var = 'phis'
 call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
 call nccheck ( nf90_get_var(ncid, varid, state%phis(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )

 if (read_tlad_traj == 1) then

   !Optional TLM/ADM trajectory variables
  
   var = 'qls'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%qls(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )
  
   var = 'qcn'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%qcn(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )
  
   var = 'cfcn'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%cfcn(isc:iec,jsc:jec,:), istart3, icount3), "nf90_get_var"//trim(var) )
   
   var = 'frland'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%frland(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'frocean'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%frocean(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'kcbl'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%kcbl(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'ts'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%ts(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'khl'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%khl(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'khu'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%khu(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'varflt'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%varflt(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'ustar'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%ustar(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'bstar'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%bstar(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'zpbl'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%zpbl(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'cm'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%cm(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'ct'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%ct(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )
  
   var = 'cq'
   call nccheck ( nf90_inq_varid (ncid, trim(var), varid), "nf90_inq_varid"//trim(var) )
   call nccheck ( nf90_get_var(ncid, varid, state%cq(isc:iec,jsc:jec), istart2, icount2), "nf90_get_var"//trim(var) )

 endif

 !Close this file
 call nccheck ( nf90_close(ncid), "nf90_close" )

 ! Deallocate
 deallocate ( istart2, icount2 )
 deallocate ( istart3, icount3 )

end subroutine read_geos_state

! ------------------------------------------------------------------------------

subroutine write_geos_state(geom, state, c_conf, vdate)

implicit none

! Arguments
type(fv3jedi_geom), intent(inout)  :: geom
type(fv3jedi_state), intent(in)    :: state      !< State
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime

! Locals
character(len=255) :: datapath
character(len=255) :: filename
character(len=64)  :: datefile
integer :: geostiledim, tileoff
integer :: ncid, varid(1000), vc = 0
integer :: date(6)
integer(kind=c_int) :: idate, isecs
integer :: isc,iec,jsc,jec,im,jm,km
integer :: x_dimid, y_dimid, z_dimid
integer :: t_dimid, tile_dimid
integer, allocatable :: dimidsv(:), dimidsg(:), dimids2(:), dimids3(:)
integer, allocatable :: istart2(:), icount2(:)
integer, allocatable :: istart3(:), icount3(:)
type(fckit_mpi_comm) :: f_comm
integer :: k, levs(geom%npz)

! Convenience
isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec

f_comm = fckit_mpi_comm()

do k = 1,geom%npz
  levs(k) = k
enddo

! Place to save restarts
datapath = "Data/"
if (config_element_exists(c_conf,"datapath")) then
   datapath = config_get_string(c_conf,len(datapath),"datapath")
endif

! Current date
call datetime_to_ifs(vdate, idate, isecs)

date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

! Using tile as a dimension in the file?
geostiledim = 0
if (config_element_exists(c_conf,"geos_tile_dim")) then
   geostiledim = config_get_int(c_conf,"geos_tile_dim")
endif

! File total dims
im = geom%npx-1
jm = geom%npy-1
km = state%npz
if (geostiledim == 0) then
  jm = 6*jm
  tileoff = (state%ntile-1)*(jm/state%ntiles)
endif

! Naming convention for the file
filename = 'GEOS.eta.'

if (config_element_exists(c_conf,"filename")) then
   filename = config_get_string(c_conf,len(filename),"filename")
endif

! Append with the date
write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)
filename = trim(datapath)//trim(filename)//trim(datefile)//trim("z.nc4")

call nccheck( nf90_create( filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                           comm = f_comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

! Create dimensions
call nccheck ( nf90_def_dim(ncid, "Xdim",  im, x_dimid), "nf90_def_dim Xdim" )
call nccheck ( nf90_def_dim(ncid, "Ydim",  jm, y_dimid), "nf90_def_dim Ydim" )

! Add dimension for the tile number
if (geostiledim == 1) then
  call nccheck ( nf90_def_dim(ncid, "nf", geom%ntiles, tile_dimid), "nf90_def_dim nf"  )
endif

call nccheck ( nf90_def_dim(ncid, "lev",  km, z_dimid), "nf90_def_dim lev" )
call nccheck ( nf90_def_dim(ncid, "time", 1,  t_dimid), "nf90_def_dim time" )

! DimId arrays
if (geostiledim == 1) then

  allocate(dimidsv(1))
  dimidsv =  (/ z_dimid /)
  allocate(dimidsg(3))
  dimidsg =  (/ x_dimid, y_dimid, tile_dimid /)
  allocate(dimids2(4))
  dimids2 =  (/ x_dimid, y_dimid, tile_dimid, t_dimid /)
  allocate(dimids3(5))
  dimids3 =  (/ x_dimid, y_dimid, tile_dimid, z_dimid, t_dimid /)

else

  allocate(dimidsv(1))
  dimidsv =  (/ z_dimid /)
  allocate(dimidsg(2))
  dimidsg =  (/ x_dimid, y_dimid /)
  allocate(dimids2(3))
  dimids2 =  (/ x_dimid, y_dimid, t_dimid /)
  allocate(dimids3(4))
  dimids3 =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)

endif

! Define fields to be written (geom)
vc=vc+1;call nccheck( nf90_def_var(ncid, "lons", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lons" )
vc=vc+1;call nccheck( nf90_def_var(ncid, "lats", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lats" )

! Define fields to be written (state)
vc=vc+1;call nccheck( nf90_def_var(ncid, "ud",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var ud"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "vd",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var vd"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "ua",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var ua"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "va",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var va"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "t",    NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var t"    )
vc=vc+1;call nccheck( nf90_def_var(ncid, "delp", NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var delp" )
vc=vc+1;call nccheck( nf90_def_var(ncid, "q",    NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var q"    )
vc=vc+1;call nccheck( nf90_def_var(ncid, "qi",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var qi"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "ql",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var ql"   )
vc=vc+1;call nccheck( nf90_def_var(ncid, "o3",   NF90_DOUBLE, dimids3, varid(vc)), "nf90_def_var o3"   )

! End define mode
vc=0
call nccheck( nf90_enddef(ncid), "nf90_enddef" )

! Create local to this proc start/count
if (geostiledim == 1) then

  allocate(istart3(5),icount3(5))
  allocate(istart2(4),icount2(4))

  istart3(1) = isc;          icount3(1) = iec-isc+1
  istart3(2) = jsc;          icount3(2) = jec-jsc+1
  istart3(3) = geom%ntile;   icount3(3) = 1
  istart3(4) = 1;            icount3(4) = geom%npz
  istart3(5) = 1;            icount3(5) = 1
  
  istart2(1) = isc;          icount2(1) = iec-isc+1
  istart2(2) = jsc;          icount2(2) = jec-jsc+1
  istart2(3) = geom%ntile;   icount2(3) = 1
  istart2(4) = 1;            icount2(4) = 1

else

  allocate(istart3(4),icount3(4))
  allocate(istart2(3),icount2(3))

  istart3(1) = isc;          icount3(1) = iec-isc+1
  istart3(2) = tileoff+jsc;  icount3(2) = jec-jsc+1
  istart3(3) = 1;            icount3(3) = geom%npz
  istart3(4) = 1;            icount3(4) = 1
  
  istart2(1) = isc;          icount2(1) = iec-isc+1
  istart2(2) = tileoff+jsc;  icount2(2) = jec-jsc+1
  istart2(3) = 1;            icount2(3) = 1

endif

! Write fields (geom)
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), geom%grid_lon(isc:iec,jsc:jec), &
                                    start = istart2(1:3), count = icount2(1:3) ), "nf90_put_var lons" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), geom%grid_lat(isc:iec,jsc:jec), &
                                    start = istart2(1:3), count = icount2(1:3) ), "nf90_put_var lats" )

! Write fields (state)
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%ud  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ud" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%vd  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var vd" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%ua  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%va  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%t   (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%delp(isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%q   (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%qi  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%ql  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), state%o3  (isc:iec,jsc:jec,1:geom%npz), &
                                    start = istart3, count = icount3 ), "nf90_put_var ua" )

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

! Deallocate
deallocate ( dimidsv, dimidsg, dimids2, dimids3 )
deallocate ( istart2, icount2 )
deallocate ( istart3, icount3 )

end subroutine write_geos_state

! ------------------------------------------------------------------------------

subroutine nccheck(status,iam)

implicit none
integer, intent ( in) :: status
character(len=*), optional :: iam 

character(len=1024) :: error_descr

 if(status /= nf90_noerr) then

   error_descr = "fv3jedi_state_io_mod: NetCDF error, aborting"

   if (present(iam)) then
     error_descr = trim(error_descr)//", "//trim(iam)
   endif

   error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

   call abor1_ftn(trim(error_descr))

 end if

end subroutine nccheck

! ------------------------------------------------------------------------------

end module fv3jedi_state_io_mod
