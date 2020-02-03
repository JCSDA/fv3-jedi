! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_pseudo_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use duration_mod
use netcdf
use oops_variables_mod

use fckit_mpi_module, only: fckit_mpi_comm

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_io_gfs_mod, only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos

implicit none
private

public :: pseudo_model
public :: pseudo_create
public :: pseudo_delete
public :: pseudo_initialize
public :: pseudo_step
public :: pseudo_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: pseudo_model
  character(len=255)  :: pseudo_type !geos of gfs
  character(len=255)  :: pseudo_path
  character(len=255)  :: filename_bkgd
  character(len=255)  :: filename_crtm
  logical :: tiledim_bkgd, tiledim_crtm
  logical :: need_crtm_fields
end type pseudo_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine pseudo_create(self, geom, c_conf, vars)

implicit none
type(pseudo_model),   intent(inout) :: self
type(c_ptr),          intent(in)    :: c_conf
type(fv3jedi_geom),   intent(in)    :: geom
type(oops_variables), intent(in)    :: vars

type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str
integer :: var

! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)

!File types, paths and names
call f_conf%get_or_die("pseudo_type",str)
self%pseudo_type = str
deallocate(str)
call f_conf%get_or_die("pseudo_path",str)
self%pseudo_path = str
deallocate(str)

! Geos specific details
if (trim(self%pseudo_type) == "geos") then

  self%need_crtm_fields = .false.
  do var = 1, vars%nvars()
    if ( trim(vars%variable(var)) == 'vtype' .or. &
         trim(vars%variable(var)) == 'stype' .or. &
         trim(vars%variable(var)) == 'vfrac' ) then
      self%need_crtm_fields = .true.
    endif
  enddo

  call f_conf%get_or_die("filename_bkgd",str)
  self%filename_bkgd = str
  deallocate(str)
  if (self%need_crtm_fields) then
    call f_conf%get_or_die("filename_crtm",str)
    self%filename_crtm = str
    deallocate(str)
  endif
  if (.not. f_conf%get('tiledim_bkgd',self%tiledim_bkgd)) self%tiledim_bkgd = .true.
  if (.not. f_conf%get('tiledim_crtm',self%tiledim_crtm)) self%tiledim_crtm = .true.

endif



end subroutine pseudo_create

! ------------------------------------------------------------------------------

subroutine pseudo_delete(self)

implicit none
type(pseudo_model), intent(inout) :: self

end subroutine pseudo_delete

! ------------------------------------------------------------------------------

subroutine pseudo_initialize(self, state)

implicit none
type(pseudo_model),  intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

end subroutine pseudo_initialize

! ------------------------------------------------------------------------------

subroutine pseudo_step(self, state, geom, vdate)

implicit none
type(pseudo_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_geom),  intent(inout) :: geom
type(datetime),      intent(inout) :: vdate !< Valid datetime after step

character(len=20)  :: vdatec
character(len=255) :: date
character(len=4)   :: yyyy
character(len=2)   :: mm,dd,hh,mn,ss

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos

! Convert datetime to string
call datetime_to_string(vdate, vdatec)

! Write character form of date
write(date,*) vdatec(1:4),vdatec(6:7),vdatec(9:10),'_',vdatec(12:13),vdatec(15:16),'z.nc4'

! Date part of the filename
yyyy = vdatec(1 :4 )
mm   = vdatec(6 :7 )
dd   = vdatec(9 :10)
hh   = vdatec(12:13)
mn   = vdatec(15:16)
ss   = vdatec(18:19)

! File path/filename
if (trim(self%pseudo_type) == "gfs") then

  gfs%datapath_ti = trim(self%pseudo_path)
  gfs%filename_core = yyyy//mm//dd//"."//hh//mn//ss//'.fv_core.res.nc'
  gfs%filename_trcr = yyyy//mm//dd//"."//hh//mn//ss//'.fv_tracer.res.nc'
  gfs%filename_sfcd = yyyy//mm//dd//"."//hh//mn//ss//'.sfc_data.nc'
  gfs%filename_sfcw = yyyy//mm//dd//"."//hh//mn//ss//'.fv_srf_wnd.res.nc'
  gfs%filename_cplr = yyyy//mm//dd//"."//hh//mn//ss//'.coupler.res'
  gfs%datapath_sp = 'null'
  gfs%datapath_sp = 'null'

  call print_filename(self,gfs%filename_core,geom)
  call gfs%read_meta  ( geom, vdate, state%calendar_type, state%date_init )
  call gfs%read_fields( geom, state%fields )

elseif (trim(self%pseudo_type) == "geos") then

  geos%datapath = trim(self%pseudo_path)
  geos%filenames(1) = trim(self%filename_bkgd) ! Use background file for trajectory
  geos%filenames(2) = trim(self%filename_crtm) ! In case CRTM is being called

  call geos%setup(geom, state%fields, vdate, 'read')

  geos%tiledim(1) = self%tiledim_bkgd
  if (self%need_crtm_fields) geos%tiledim(2) = self%tiledim_crtm

  call geos%read_meta(geom, vdate, state%calendar_type, state%date_init)
  call geos%read_fields(geom, state%fields)
  call geos%delete()

else

  call abor1_ftn("fv3jedi_pseudo_mod: pseudo_model, model choice must be geos or gfs")

endif

end subroutine pseudo_step

! ------------------------------------------------------------------------------

subroutine pseudo_finalize(self, state)

implicit none
type(pseudo_model)  :: self
type(fv3jedi_state) :: state

end subroutine pseudo_finalize

! ------------------------------------------------------------------------------

subroutine print_filename(self,filename,geom)

implicit none
type(pseudo_model), intent(in) :: self
character(len=*),   intent(in) :: filename
type(fv3jedi_geom), intent(in) :: geom

! Print filename to the user
if (geom%f_comm%rank()==0) then
  print*, ' '
  print*, 'Pseudo model from file: ', trim(filename)
  print*, ' '
endif

end subroutine print_filename

! ------------------------------------------------------------------------------

end module fv3jedi_pseudo_mod
