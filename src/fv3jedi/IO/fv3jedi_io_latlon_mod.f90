! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_latlon_mod

use iso_c_binding

use netcdf
use mpi

use datetime_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : log
use fckit_mpi_module
use fv3jedi_constants_mod, only: rad2deg

use interpolatorbump_mod,         only: bump_interpolator

use fv3jedi_constants_mod,    only: rad2deg
use fv3jedi_field_mod,        only: fv3jedi_field
use fv3jedi_geom_mod,         only: fv3jedi_geom
use fv3jedi_kinds_mod,        only: kind_int, kind_real
use fv3jedi_io_utils_mod,     only: vdate_to_datestring, replace_text
use fv3jedi_netcdf_utils_mod, only: nccheck

implicit none

private
public fv3jedi_llgeom!, create_latlon, delete_latlon, write_latlon_metadata, write_latlon_fields, &


!Module level type
type fv3jedi_llgeom
 type(fckit_mpi_comm) :: f_comm
 type(bump_interpolator) :: bumpinterp
 integer :: llcomm
 integer :: layout(2)
 logical :: thispe = .false.
 integer :: nx, ny, nxg, nyg
 integer :: npes
 real(kind=kind_real), allocatable :: lons(:)
 real(kind=kind_real), allocatable :: lats(:)
 character(len=1024) :: filename
 integer, allocatable :: istart2(:), icount2(:)
 integer, allocatable :: istart3(:), icount3(:)
 integer, allocatable :: istarte(:), icounte(:)
 contains
  procedure :: setup_conf
  procedure :: setup_date
  procedure :: write
  procedure :: delete
  final :: dummy_final
end type fv3jedi_llgeom

contains

! --------------------------------------------------------------------------------------------------
subroutine setup_conf(self, geom)

implicit none

!Arguments
class(fv3jedi_llgeom),     intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom          !< Geometry

call log%debug("fv3jedi_io_latlon_mod%setup_conf starting")

call create_latlon(self, geom)

call log%debug("fv3jedi_io_latlon_mod%setup_conf done")

end subroutine setup_conf

! --------------------------------------------------------------------------------------------------

subroutine setup_date(self, vdate)

implicit none

!Arguments
class(fv3jedi_llgeom),     intent(inout) :: self
type(datetime),            intent(in)    :: vdate         !< DateTime

integer :: n
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss

call log%debug("fv3jedi_io_latlon_mod%setup_date starting")

call log%debug("fv3jedi_io_latlon_mod%setup_date done")

end subroutine setup_date

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3jedi_llgeom), intent(inout) :: self

call log%debug("fv3jedi_io_latlon_mod%delete starting")

call delete_latlon(self)

call log%debug("fv3jedi_io_latlon_mod%delete done")

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine write(self, geom, fields, f_conf, vdate)

implicit none
class(fv3jedi_llgeom),     intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom          !< Geom
type(fv3jedi_field),       intent(in)    :: fields(:)     !< Fields to be written
type(fckit_configuration), intent(in)    :: f_conf        !< Configuration
type(datetime),            intent(in)    :: vdate         !< DateTime

integer :: var

call log%debug("fv3jedi_io_latlon_mod%write starting")

call write_latlon_metadata(self, geom, f_conf, vdate)

call write_latlon_fields(self, geom, fields)

call log%debug("fv3jedi_io_latlon_mod%write done")

end subroutine write

! --------------------------------------------------------------------------------------------------

subroutine create_latlon(llgeom, geom)

implicit none

!Arguments
type(fv3jedi_llgeom), intent(inout) :: llgeom   !< LatLon Geometry
type(fv3jedi_geom),   intent(in)    :: geom     !< Geometry

integer :: color

real(kind=kind_real) :: dx, dy
integer :: i, ierr

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
  call abor1_ftn("create_latlon error: fewer than 6 npes not anticipated")
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

else

  llgeom%nx = 0
  llgeom%ny = 0
  allocate(llgeom%lons(llgeom%nx))
  allocate(llgeom%lats(llgeom%ny))

endif

! Initialize bump interpolator
call llgeom%bumpinterp%init(geom%f_comm, afunctionspace_in=geom%afunctionspace, lon1d_out=llgeom%lons, lat1d_out=llgeom%lats, &
   & nl=geom%npz)

!IO communicator
!call MPI_Comm_split(llgeom%f_comm%communicator(), color, llgeom%f_comm%rank(), llgeom%llcomm, ierr)

! NC arrays
allocate(llgeom%istart3(4),llgeom%icount3(4))
allocate(llgeom%istarte(4),llgeom%icounte(4))
allocate(llgeom%istart2(3),llgeom%icount2(3))
llgeom%istart3(1) = llgeom%nx * llgeom%f_comm%rank() + 1;  llgeom%icount3(1) = llgeom%nx
llgeom%istart3(2) = 1;                                     llgeom%icount3(2) = llgeom%ny
llgeom%istart3(3) = 1;                                     llgeom%icount3(3) = geom%npz
llgeom%istart3(4) = 1;                                     llgeom%icount3(4) = 1
llgeom%istart2(1) = llgeom%istart3(1);                     llgeom%icount2(1) = llgeom%icount3(1)
llgeom%istart2(2) = llgeom%istart3(2);                     llgeom%icount2(2) = llgeom%icount3(2)
llgeom%istart2(3) = llgeom%istart3(4);                     llgeom%icount2(3) = llgeom%icount3(4)

llgeom%istarte = llgeom%istarte
llgeom%icounte = llgeom%icount3
llgeom%istarte(3) = geom%npz + 1

end subroutine create_latlon

! --------------------------------------------------------------------------------------------------

subroutine delete_latlon(llgeom)

implicit none

!Arguments
type(fv3jedi_llgeom), intent(inout) :: llgeom  !< LatLon Geometry

if (llgeom%thispe) then
  deallocate(llgeom%lons)
  deallocate(llgeom%lats)

  deallocate ( llgeom%istart2, llgeom%icount2 )
  deallocate ( llgeom%istart3, llgeom%icount3 )
  deallocate ( llgeom%istarte, llgeom%icounte )
endif

call llgeom%bumpinterp%delete()


end subroutine delete_latlon

! --------------------------------------------------------------------------------------------------

subroutine write_latlon_metadata(llgeom, geom, f_conf, vdate)

implicit none

!Arguments
type(fv3jedi_llgeom),      intent(inout) :: llgeom   !< LatLon Geometry
type(fv3jedi_geom),        intent(in)    :: geom     !< Geometry
type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
type(datetime),            intent(in)    :: vdate    !< DateTime

integer :: date(6)
integer(kind=c_int) :: idate, isecs
character(len=64)   :: datefile

integer :: ncid, varid(2)
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid
character(len=:), allocatable :: str

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
call nccheck ( nf90_def_dim(ncid, "edge", geom%npz+1 , e_dimid), "nf90_def_dim lev" )
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

! --------------------------------------------------------------------------------------------------

subroutine write_latlon_fields(llgeom, geom, fields)

implicit none

!Arguments
type(fv3jedi_llgeom), target, intent(inout)   :: llgeom         !< LatLon Geometry
type(fv3jedi_geom),           intent(in)      :: geom           !< Geometry
type(fv3jedi_field),          intent(in)      :: fields(:)      !< Field to write

integer :: var, ji, jj, jk, csngrid, llngrid, ii, i, j, k, n
real(kind=kind_real), allocatable :: llfield(:,:,:)

integer :: ncid, varid
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid
integer, target  :: dimids3(4), dimids2(3), dimidse(4)
integer, pointer :: dimids(:), istart(:), icount(:)

! Loop over fields
! ----------------
do var = 1, size(fields)

  ! Only certain fields can be written for now
  ! ------------------------------------------
  if ( trim(fields(var)%space) == 'magnitude' .and. &
       trim(fields(var)%staggerloc) == 'center' .and. &
       .not. fields(var)%integerfield ) then

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

    if (llgeom%thispe) then
      ! Interpolate
      call llgeom%bumpinterp%apply(fields(var)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1:fields(var)%npz), &
         & llfield(:,:,1:fields(var)%npz))
    endif

    !Open file
    call nccheck( nf90_open( llgeom%filename, ior(NF90_WRITE, NF90_MPIIO), ncid, &
                            comm = llgeom%f_comm%communicator(), info = MPI_INFO_NULL), "nf90_open" )

    !Dimension ID
    call nccheck( nf90_inq_dimid(ncid, "lon" , x_dimid), "nf90_inq_dimid lon" )
    call nccheck( nf90_inq_dimid(ncid, "lat" , y_dimid), "nf90_inq_dimid lat" )
    call nccheck( nf90_inq_dimid(ncid, "lev" , z_dimid), "nf90_inq_dimid lev" )
    call nccheck( nf90_inq_dimid(ncid, "edge", e_dimid), "nf90_inq_dimid edge" )
    call nccheck( nf90_inq_dimid(ncid, "time", t_dimid), "nf90_inq_dimid time" )

    dimids3 =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)
    dimidse =  (/ x_dimid, y_dimid, e_dimid, t_dimid /)
    dimids2 =  (/ x_dimid, y_dimid, t_dimid /)

    ! Pointer to dimensions based on number of levels
    if (associated(dimids)) nullify(dimids)
    if (associated(istart)) nullify(istart)
    if (associated(icount)) nullify(icount)
    if (fields(var)%npz == geom%npz) then
      dimids => dimids3
      istart => llgeom%istart3
      icount => llgeom%icount3
    elseif (fields(var)%npz == 1) then
      dimids => dimids2
      istart => llgeom%istart2
      icount => llgeom%icount2
    elseif (fields(var)%npz == geom%npz+1) then
      dimids => dimidse
      istart => llgeom%istarte
      icount => llgeom%icounte
    endif

    if (associated(dimids)) then

      ! Write field to the file
      call nccheck( nf90_def_var( ncid, trim(fields(var)%short_name), NF90_DOUBLE, dimids, varid), &
                    "nf90_def_var"//trim(fields(var)%short_name) )
      call nccheck( nf90_put_att(ncid, varid, "long_name", trim(fields(var)%long_name) ), "nf90_put_att" )
      call nccheck( nf90_put_att(ncid, varid, "units"    , trim(fields(var)%units)     ), "nf90_put_att" )
      call nccheck( nf90_enddef(ncid), "nf90_enddef" )

      if (llgeom%thispe) then
       call nccheck( nf90_put_var( ncid, varid, llfield, start = istart, count = icount), &
                     "nf90_put_var"//trim(fields(var)%short_name) )
      endif

    endif

    ! Close file
    call nccheck( nf90_close(ncid), "nf90_close" )

    !Let LatLon PEs catch up
    call llgeom%f_comm%barrier()

    deallocate(llfield)

  endif

enddo

end subroutine write_latlon_fields

! --------------------------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_llgeom), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_latlon_mod
