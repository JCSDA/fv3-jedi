module fv3jedi_fields_io_mod

use config_mod
use iso_c_binding
use datetime_mod

use fv3jedi_kinds, only: kind_real
use fv3jedi_fields_utils_mod, only: fv3jedi_field

!For FMS like restarts
use mpp_domains_mod,   only: EAST, NORTH
use fms_io_mod,        only: restart_file_type, register_restart_field, &
                             free_restart_type, restore_state, save_restart

implicit none
private
public read_fms_restart, write_fms_restart, &
       read_geos_restart, write_geos_restart

contains

! ------------------------------------------------------------------------------

subroutine read_fms_restart(fld, c_conf, vdate)

use iso_c_binding
use datetime_mod

use mpp_domains_mod,     only: mpp_update_domains, DGRID_NE
use field_manager_mod,   only: MODEL_ATMOS
use tracer_manager_mod,  only: get_number_tracers, get_tracer_names, set_tracer_profile, get_tracer_index

implicit none

!Arguments
type(fv3jedi_field), intent(inout) :: fld      !< Fields
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
 !Winds
 if (trim(fld%geom%wind_type) == 'D-grid') then
    id_restart = register_restart_field(Fv_restart, filename_core, 'u', fld%Atm%u, &
                 domain=fld%geom%domain, position=NORTH)
    id_restart = register_restart_field(Fv_restart, filename_core, 'v', fld%Atm%v, &
                 domain=fld%geom%domain, position=EAST)
 elseif (trim(fld%geom%wind_type) == 'A-grid') then
    id_restart = register_restart_field(Fv_restart, filename_core, 'ua', fld%Atm%u, &
                                        domain=fld%geom%domain)
    id_restart = register_restart_field(Fv_restart, filename_core, 'va', fld%Atm%v, &
                                        domain=fld%geom%domain)
 endif

 !phis
 id_restart = register_restart_field(Fv_restart, filename_core, 'phis', fld%Atm%phis, &
              domain=fld%geom%domain)

 !Temperature
 id_restart = register_restart_field(Fv_restart, filename_core, 'T', fld%Atm%pt, &
              domain=fld%geom%domain)

 !Pressure thickness
 id_restart = register_restart_field(Fv_restart, filename_core, 'DELP', fld%Atm%delp, &
              domain=fld%geom%domain)

 !Nonhydrostatic variables
 if (.not. fld%Atm%hydrostatic) then
    id_restart =  register_restart_field(Fv_restart, filename_core, 'W', fld%Atm%w, &
                  domain=fld%geom%domain)
    id_restart =  register_restart_field(Fv_restart, filename_core, 'DZ', fld%Atm%delz, &
                  domain=fld%geom%domain)
 endif

 ! Read file and fill variables
 ! ----------------------------
 call restore_state(Fv_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Fv_restart)
 

 !Register and read tracers
 !-------------------------
 call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)

 if (ntracers /= fld%geom%ntracers) then
   call abor1_ftn("fv3jedi_fields: tracer table does not match geomtry")
 endif
 
 do nt = 1, ntprog
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    call set_tracer_profile (MODEL_ATMOS, nt, fld%Atm%q(:,:,:,nt) )
    id_restart = register_restart_field(Tr_restart, filename_trcr, tracer_name, fld%Atm%q(:,:,:,nt), &
                                        domain=fld%geom%domain)
 enddo

 call restore_state(Tr_restart, directory=trim(adjustl(datapath_ti)))
 call free_restart_type(Tr_restart)

 !Register and read surface fields needed for crtm calculation
 !------------------------------------------------------------
 read_crtm_surface = 0
 if (config_element_exists(c_conf,"read_crtm_surface")) then
   read_crtm_surface = config_get_int(c_conf,"read_crtm_surface")
 endif

 if (read_crtm_surface == 1) then
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'slmsk' , fld%Atm%slmsk , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'sheleg', fld%Atm%sheleg, domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'tsea'  , fld%Atm%tsea  , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vtype' , fld%Atm%vtype , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stype' , fld%Atm%stype , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'vfrac' , fld%Atm%vfrac , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'stc'   , fld%Atm%stc   , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'smc'   , fld%Atm%smc   , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'snwdph', fld%Atm%snwdph, domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcd, 'f10m'  , fld%Atm%f10m  , domain=fld%geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)

   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'u_srf' , fld%Atm%u_srf , domain=fld%geom%domain)
   id_restart = register_restart_field( Sf_restart, filename_sfcw, 'v_srf' , fld%Atm%v_srf , domain=fld%geom%domain)

   call restore_state(Sf_restart, directory=trim(adjustl(datapath_ti)))
   call free_restart_type(Sf_restart)
   fld%havecrtmfields = .true.
 else
   fld%havecrtmfields = .false.
   fld%Atm%slmsk  = 0.0_kind_real
   fld%Atm%sheleg = 0.0_kind_real
   fld%Atm%tsea   = 0.0_kind_real
   fld%Atm%vtype  = 0.0_kind_real
   fld%Atm%stype  = 0.0_kind_real
   fld%Atm%vfrac  = 0.0_kind_real
   fld%Atm%stc    = 0.0_kind_real
   fld%Atm%smc    = 0.0_kind_real
   fld%Atm%u_srf  = 0.0_kind_real
   fld%Atm%u_srf  = 0.0_kind_real
   fld%Atm%v_srf  = 0.0_kind_real
   fld%Atm%f10m   = 0.0_kind_real
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
 fld%Atm%date = date
 fld%Atm%date_init = date_init
 fld%Atm%calendar_type = calendar_type
 idate=date(1)*10000+date(2)*100+date(3)
 isecs=date(4)*3600+date(5)*60+date(6)

 if (fld%root_pe == 1 .and. print_read_info == 1 ) then
    print *,'read_file: integer time from coupler.res: ',date,idate,isecs
 endif

 call datetime_from_ifs(vdate, idate, isecs)
 call datetime_to_string(vdate, validitydate)

 if (fld%root_pe == 1 .and. print_read_info == 1 ) then
    print *,'read_file: validity date: ',trim(validitydate)
    print *,'read_file: expected validity date: ',trim(sdate)
 endif

 return

end subroutine read_fms_restart

! ------------------------------------------------------------------------------

subroutine write_fms_restart(fld, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_field), intent(in)    :: fld      !< Fields
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
 !Winds
 if (trim(fld%geom%wind_type) == 'D-grid') then
    id_restart = register_restart_field( Fv_restart, filename_core, 'u', fld%Atm%u, &
                                         domain=fld%geom%domain,position=NORTH )
    id_restart = register_restart_field( Fv_restart, filename_core, 'v', fld%Atm%v, &
                                         domain=fld%geom%domain,position=EAST )
 elseif (trim(fld%geom%wind_type) == 'A-grid') then
     id_restart =  register_restart_field(Fv_restart, filename_core, 'ua', fld%Atm%u, &
                                           domain=fld%geom%domain )
     id_restart =  register_restart_field(Fv_restart, filename_core, 'va', fld%Atm%v, &
                                           domain=fld%geom%domain )
 endif

 !phis
 id_restart = register_restart_field( Fv_restart, filename_core, 'phis', fld%Atm%phis, &
                                      domain=fld%geom%domain )

 !Temperature
 id_restart = register_restart_field( Fv_restart, filename_core, 'T', fld%Atm%pt, &
                                      domain=fld%geom%domain )

 !Pressure thickness
 id_restart = register_restart_field( Fv_restart, filename_core, 'DELP', fld%Atm%delp, &
                                      domain=fld%geom%domain )

 !Nonhydrostatic fields
 if (.not. fld%Atm%hydrostatic) then
     id_restart =  register_restart_field( Fv_restart, filename_core, 'W', fld%Atm%w, &
                                           domain=fld%geom%domain )
     id_restart =  register_restart_field( Fv_restart, filename_core, 'DZ', fld%Atm%delz, &
                                           domain=fld%geom%domain )
 endif

 !Cell center lat/lon
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lat', fld%geom%grid_lat, &
                                      domain=fld%geom%domain )
 id_restart = register_restart_field( Fv_restart, filename_core, 'grid_lon', fld%geom%grid_lon, &
                                      domain=fld%geom%domain )

 id_restart =  register_restart_field( Fv_restart, filename_core, 'ua', fld%Atm%ua, &
                                       domain=fld%geom%domain )
 id_restart =  register_restart_field( Fv_restart, filename_core, 'va', fld%Atm%va, &
                                       domain=fld%geom%domain )

 ! Write variables to file
 ! -----------------------
 call save_restart(Fv_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Fv_restart)


 !Write tracers to file
 !---------------------
 id_restart = register_restart_field( Tr_restart, filename_trcr, 'sphum', fld%Atm%q(:,:,:,1), &
                                      domain=fld%geom%domain )

 call save_restart(Tr_restart, directory=trim(adjustl(datapath_out))//'RESTART')
 call free_restart_type(Tr_restart)


 !Write date/time info in coupler.res
 !-----------------------------------
 iounit = 101
 if (fld%root_pe == 1) then
    print *,'write_file: date model init = ',fld%Atm%date_init
    print *,'write_file: date model now  = ',fld%Atm%date
    print *,'write_file: date vdate      = ',date
    open(iounit, file=trim(adjustl(datapath_out))//'RESTART/'//trim(adjustl(filename_cplr)), form='formatted')
    write( iounit, '(i6,8x,a)' ) fld%Atm%calendar_type, &
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

subroutine read_geos_restart(fld, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_field), intent(inout) :: fld      !< Fields
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine read_geos_restart

! ------------------------------------------------------------------------------

subroutine write_geos_restart(fld, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_field), intent(in)    :: fld      !< Fields
type(c_ptr), intent(in)            :: c_conf   !< Configuration
type(datetime), intent(inout)      :: vdate    !< DateTime


end subroutine write_geos_restart

! ------------------------------------------------------------------------------

end module fv3jedi_fields_io_mod
