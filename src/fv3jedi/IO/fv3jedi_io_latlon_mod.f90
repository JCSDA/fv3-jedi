! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_latlon_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use datetime_mod
use fckit_log_module, only : log
use fckit_mpi_module
use fv3jedi_constants_mod, only: rad2deg

use type_bump, only: bump_type

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

use netcdf
use mpi

implicit none

private
public fv3jedi_llgeom, create_latlon, delete_latlon, &
       write_latlon_metadata, write_latlon_field

interface write_latlon_field
 module procedure write_latlon_field_r3
 module procedure write_latlon_field_r2
end interface

!Module level type
type fv3jedi_llgeom
 type(fckit_mpi_comm) :: f_comm
 integer :: llcomm
 integer :: layout(2)
 logical :: thispe = .false.
 integer :: nx, ny, nxg, nyg
 integer :: npes
 real(kind=kind_real), allocatable :: lons(:)
 real(kind=kind_real), allocatable :: lats(:)
 type(bump_type) :: bump
 character(len=1024) :: filename
 integer, allocatable :: istart2(:), icount2(:)
 integer, allocatable :: istart3(:), icount3(:)
 contains
  final :: dummy_final
end type fv3jedi_llgeom

contains

! ------------------------------------------------------------------------------

subroutine create_latlon(geom, llgeom)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)    :: geom     !< Geometry
type(fv3jedi_llgeom), intent(inout) :: llgeom   !< LatLon Geometry

integer :: color

real(kind=kind_real) :: dx, dy
integer :: i, j, ji, jj, ii, locs_nlocs, ierr
real(kind=kind_real), allocatable :: locs_lat(:), locs_lon(:)

! Create lat lon grid and interpolation object for going from cube to lat-lon
! --------------------------------------------------------------------------

llgeom%f_comm = geom%f_comm

!Maximum of 12 IO processors
if (llgeom%f_comm%size() >= 12) then
  llgeom%layout(1) = 12
  llgeom%layout(2) = 1
elseif (llgeom%f_comm%size() >= 6) then
  llgeom%layout(1) = 6
  llgeom%layout(2) = 1
else
  call abor1_ftn("creat_latlon error: fewer than 6 npes not anticipated")
endif
llgeom%npes = llgeom%layout(1) * llgeom%layout(2)

!Since the lat lon grid is only for IO purposes it is only
!defined on a subset of PEs - those that will do the writing.
!This is generally more efficient than having many PEs trying
!to write to the same file.

llgeom%nxg = 4*(geom%npx - 1)
llgeom%nyg = 2*(geom%npy - 1) + 1
llgeom%nx = 0
llgeom%ny = 0

color = MPI_UNDEFINED

if (llgeom%f_comm%rank() <= llgeom%npes-1) then

  !Split communicator
  color = int(llgeom%f_comm%communicator() / 2)

  llgeom%thispe = .true.

  !Resolution
  dx = 360.0_kind_real / (real(llgeom%nxg,kind_real) - 1.0_kind_real)
  dy = 180.0_kind_real / (real(llgeom%nyg,kind_real) - 1.0_kind_real)

  llgeom%nx = llgeom%nxg / llgeom%layout(1)
  llgeom%ny = llgeom%nyg / llgeom%layout(2)

  allocate(llgeom%lons(llgeom%nx))
  allocate(llgeom%lats(llgeom%ny))

  !Each processor has subset of lons
  llgeom%lons(1) = dx * llgeom%nx * llgeom%f_comm%rank()
  do i = 2,llgeom%nx
    llgeom%lons(i) = llgeom%lons(i-1) + dx
  enddo

  !Each processor has all lats
  llgeom%lats(1) = -90.0_kind_real
  do i = 2,llgeom%ny
    llgeom%lats(i) = llgeom%lats(i-1) + dy
  enddo

  locs_nlocs = llgeom%nx*llgeom%ny
  allocate(locs_lon(locs_nlocs))
  allocate(locs_lat(locs_nlocs))

  ii = 0
  do jj = 1,llgeom%ny
    do ji = 1,llgeom%nx
      ii = ii + 1
      locs_lon(ii) = llgeom%lons(ji)
      locs_lat(ii) = llgeom%lats(jj)
    enddo
  enddo

else

  locs_nlocs = 0
  allocate(locs_lon(0))
  allocate(locs_lat(0))

endif

call initialize_bump(geom, locs_nlocs, locs_lat, locs_lon, llgeom%bump)

deallocate(locs_lon)
deallocate(locs_lat)

!IO communicator
!call MPI_Comm_split(llgeom%f_comm%communicator(), color, llgeom%f_comm%rank(), llgeom%llcomm, ierr)

! NC arrays
allocate(llgeom%istart3(4),llgeom%icount3(4))
allocate(llgeom%istart2(3),llgeom%icount2(3))
llgeom%istart3(1) = llgeom%nx * llgeom%f_comm%rank() + 1;  llgeom%icount3(1) = llgeom%nx
llgeom%istart3(2) = 1;                                     llgeom%icount3(2) = llgeom%ny
llgeom%istart3(3) = 1;                                     llgeom%icount3(3) = geom%npz
llgeom%istart3(4) = 1;                                     llgeom%icount3(4) = 1
llgeom%istart2(1) = llgeom%istart3(1);                     llgeom%icount2(1) = llgeom%icount3(1)
llgeom%istart2(2) = llgeom%istart3(2);                     llgeom%icount2(2) = llgeom%icount3(2)
llgeom%istart2(3) = llgeom%istart3(4);                     llgeom%icount2(3) = llgeom%icount3(4)

end subroutine create_latlon

! ------------------------------------------------------------------------------

subroutine delete_latlon(llgeom)

implicit none

!Arguments
type(fv3jedi_llgeom), intent(inout) :: llgeom  !< LatLon Geometry

if (llgeom%thispe) then
  deallocate(llgeom%lons)
  deallocate(llgeom%lats)

  deallocate ( llgeom%istart2, llgeom%icount2 )
  deallocate ( llgeom%istart3, llgeom%icount3 )
endif

call llgeom%bump%dealloc


end subroutine delete_latlon

! ------------------------------------------------------------------------------

subroutine write_latlon_metadata(geom, llgeom, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)    :: geom     !< Geometry
type(fv3jedi_llgeom), intent(inout) :: llgeom   !< LatLon Geometry
type(c_ptr),          intent(in)    :: c_conf   !< Configuration
type(datetime),       intent(in)    :: vdate    !< DateTime

integer :: date(6)
integer(kind=c_int) :: idate, isecs
character(len=64)   :: datefile

integer :: ncid, varid(2)
integer :: x_dimid, y_dimid, z_dimid, t_dimid
type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)


! Current date
call datetime_to_ifs(vdate, idate, isecs)
date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

! Naming convention for the file
llgeom%filename = 'Data/fv3jedi.latlon.'
if (f_conf%has("filename")) then
   call f_conf%get_or_die("filename",str)
   llgeom%filename = str
   deallocate(str)
endif

! Append filename with the curent datetime from fv3jedi
write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)
llgeom%filename = trim(llgeom%filename)//trim(datefile)//trim("z.nc4")

call nccheck( nf90_create( llgeom%filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                           comm = llgeom%f_comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

! Create dimensions
call nccheck ( nf90_def_dim(ncid, "lon" , llgeom%nxg , x_dimid), "nf90_def_dim lon" )
call nccheck ( nf90_def_dim(ncid, "lat" , llgeom%nyg , y_dimid), "nf90_def_dim lat" )
call nccheck ( nf90_def_dim(ncid, "lev" , geom%npz   , z_dimid), "nf90_def_dim lev" )
call nccheck ( nf90_def_dim(ncid, "time", 1          , t_dimid), "nf90_def_dim time" )

! Define fields to be written (geom)
call nccheck( nf90_def_var(ncid, "lons", NF90_DOUBLE, x_dimid, varid(1)), "nf90_def_var lons" )
call nccheck( nf90_def_var(ncid, "lats", NF90_DOUBLE, y_dimid, varid(2)), "nf90_def_var lats" )

call nccheck( nf90_enddef(ncid), "nf90_enddef" )

if (llgeom%thispe) then

  ! Write fields (geom)
  call nccheck( nf90_put_var( ncid, varid(1), llgeom%lons, &
                              start = llgeom%istart2(1:1), count = llgeom%icount2(1:1) ), "nf90_put_var lons" )
  call nccheck( nf90_put_var( ncid, varid(2), llgeom%lats, &
                              start = llgeom%istart2(2:2), count = llgeom%icount2(2:2) ), "nf90_put_var lats" )

endif

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

!Let LatLon PEs catch up
call llgeom%f_comm%barrier()

end subroutine write_latlon_metadata

! ------------------------------------------------------------------------------

subroutine write_latlon_field_r3(geom, llgeom, csfield, fieldname, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)    :: geom           !< Geometry
type(fv3jedi_llgeom), intent(inout) :: llgeom         !< LatLon Geometry
real(kind=kind_real), intent(in)    :: csfield(:,:,:) !< Field to write
character(len=*),     intent(in)    :: fieldname      !< Name for field
type(c_ptr),          intent(in)    :: c_conf         !< Configuration
type(datetime),       intent(in)    :: vdate          !< DateTime

integer :: ji, jj, jk, csngrid, llngrid, ii, i, j, k, n
real(kind=kind_real), allocatable :: csfield_bump(:,:), llfield_bump(:,:)
real(kind=kind_real), allocatable :: llfield(:,:,:)

integer :: ncid, varid
integer :: x_dimid, y_dimid, z_dimid, t_dimid, dimids(4)

! Interpolate the field to the lat-lon grid
! -----------------------------------------

if (llgeom%thispe) then
  allocate(llfield(1:llgeom%nx,1:llgeom%ny,1:geom%npz))
else
  allocate(llfield(0,0,0))
endif
llfield = 0.0_kind_real

csngrid = (geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)
llngrid = llgeom%nx*llgeom%ny

!Allocate fields with BUMP arrangement
allocate(csfield_bump(csngrid,1))
allocate(llfield_bump(llngrid,1))

do jk = 1, geom%npz

  !Pack Cube sphere field as BUMP wants
  ii = 0
  do jj = geom%jsc, geom%jec
    do ji = geom%isc, geom%iec
      ii = ii + 1
      csfield_bump(ii, 1) = csfield(ji, jj, jk)
    enddo
  enddo

  !Bilinear interpolation to latlon grid
  call llgeom%bump%apply_obsop(csfield_bump,llfield_bump)

  !Unpack BUMP latlon field
  if (llgeom%thispe) then
    ii = 0
    do jj = 1,llgeom%ny
      do ji = 1,llgeom%nx
        ii = ii + 1
        llfield(ji,jj,jk) = llfield_bump(ii,1)
      enddo
    enddo
  endif

enddo

!Open file
call nccheck( nf90_open( llgeom%filename, ior(NF90_WRITE, NF90_MPIIO), ncid, &
                         comm = llgeom%f_comm%communicator(), info = MPI_INFO_NULL), "nf90_open" )

!Dimension ID
call nccheck( nf90_inq_dimid(ncid, "lon" , x_dimid), "nf90_inq_dimid lon" )
call nccheck( nf90_inq_dimid(ncid, "lat" , y_dimid), "nf90_inq_dimid lat" )
call nccheck( nf90_inq_dimid(ncid, "lev" , z_dimid), "nf90_inq_dimid lev" )
call nccheck( nf90_inq_dimid(ncid, "time", t_dimid), "nf90_inq_dimid time" )
dimids =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)

! Write field to the file
call nccheck( nf90_def_var( ncid, trim(fieldname), NF90_DOUBLE, dimids, varid), "nf90_def_var"//trim(fieldname)   )
call nccheck( nf90_enddef(ncid), "nf90_enddef" )

if (llgeom%thispe) then
  call nccheck( nf90_put_var( ncid, varid, llfield, start = llgeom%istart3, count = llgeom%icount3 ), "nf90_put_var"//trim(fieldname) )
endif

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

!Let LatLon PEs catch up
call llgeom%f_comm%barrier()

deallocate(llfield)
deallocate(csfield_bump,llfield_bump)

end subroutine write_latlon_field_r3

! ------------------------------------------------------------------------------

subroutine write_latlon_field_r2(geom, llgeom, field, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)    :: geom         !< Geometry
type(fv3jedi_llgeom), intent(inout) :: llgeom       !< LatLon Geometry
real(kind=kind_real), intent(in)    :: field(:,:)   !< Field to write
type(c_ptr),          intent(in)    :: c_conf       !< Configuration
type(datetime),       intent(in)    :: vdate        !< DateTime

if (llgeom%thispe) then
endif

end subroutine write_latlon_field_r2

! ------------------------------------------------------------------------------

subroutine initialize_bump(geom, locs_nlocs, locs_lat, locs_lon, bump)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: locs_nlocs
real(kind=kind_real), intent(in)    :: locs_lat(locs_nlocs)
real(kind=kind_real), intent(in)    :: locs_lon(locs_nlocs)
type(bump_type),      intent(inout) :: bump

!Locals
integer :: mod_num
real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:)
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

character(len=5)   :: cbumpcount
character(len=255) :: bump_nam_prefix

type(fckit_mpi_comm) :: f_comm


! Communicator from geometry
! --------------------------
f_comm = geom%f_comm

! Each bump%nam%prefix must be distinct
! -------------------------------------
write(cbumpcount,"(I0.5)") 1
bump_nam_prefix = 'fv3jedi_bump_c2latllo_'//cbumpcount

!Get the Solution dimensions
!---------------------------
mod_num = (geom%iec - geom%isc + 1) * (geom%jec - geom%jsc + 1)


!Calculate interpolation weight using BUMP
!-----------------------------------------
allocate(mod_lat(mod_num))
allocate(mod_lon(mod_num))
mod_lat = reshape( rad2deg*geom%grid_lat(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] )
mod_lon = reshape( rad2deg*geom%grid_lon(geom%isc:geom%iec,      &
                                         geom%jsc:geom%jec),     &
                                        [mod_num] )


! Namelist options
! ----------------

!Important namelist options
call bump%nam%init

!Less important namelist options (should not be changed)
bump%nam%prefix = trim(bump_nam_prefix)   ! Prefix for files output
bump%nam%default_seed = .true.
bump%nam%new_obsop = .true.

bump%nam%write_obsop = .false.
bump%nam%verbosity = "none"

! Initialize geometry
! -------------------
allocate(area(mod_num))
allocate(vunit(mod_num,1))
allocate(lmask(mod_num,1))
area = 1.0           ! Dummy area
vunit = 1.0          ! Dummy vertical unit
lmask = .true.       ! Mask

! Initialize BUMP
! ---------------
call bump%setup_online( f_comm,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                        nobs=locs_nlocs,lonobs=locs_lon(:),latobs=locs_lat(:) )

!Run BUMP drivers
call bump%run_drivers

!Partial deallocate option
call bump%partial_dealloc

! Release memory
! --------------
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate(mod_lat)
deallocate(mod_lon)

end subroutine initialize_bump

! ------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_llgeom), intent(inout) :: self
end subroutine dummy_final

! ------------------------------------------------------------------------------

subroutine nccheck(status,iam)

implicit none
integer, intent ( in) :: status
character(len=*), optional :: iam

character(len=1024) :: error_descr

 if(status /= nf90_noerr) then

   error_descr = "fv3jedi_io_latlon_mod: NetCDF error, aborting"

   if (present(iam)) then
     error_descr = trim(error_descr)//", "//trim(iam)
   endif

   error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

   call abor1_ftn(trim(error_descr))

 end if

end subroutine nccheck

! ------------------------------------------------------------------------------

end module fv3jedi_io_latlon_mod
