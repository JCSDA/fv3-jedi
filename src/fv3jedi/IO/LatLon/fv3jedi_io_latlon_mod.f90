! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_latlon_mod

! libs
use netcdf
use mpi

! fckit
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm

! oops
use datetime_mod
use string_utils,               only: swap_name_member
use interpolatorbump_mod,       only: bump_interpolator

! fv3jedi
use fv3jedi_field_mod,          only: fv3jedi_field
use fv3jedi_geom_mod,           only: fv3jedi_geom
use fv3jedi_kinds_mod,          only: kind_int, kind_real
use fv3jedi_io_utils_mod,       only: add_iteration
use fv3jedi_netcdf_utils_mod,   only: nccheck

implicit none
private
public fv3jedi_io_latlon

type fv3jedi_io_latlon
 type(fckit_mpi_comm) :: comm
 type(bump_interpolator) :: bumpinterp
 integer :: llcomm
 integer :: layout(2)
 logical :: thispe = .false.
 integer :: nx, ny, nxg, nyg
 integer :: npes
 integer :: isc, iec, jsc, jec, npx, npy, npz
 integer :: isd, ied, jsd, jed
 real(kind=kind_real), allocatable :: lons(:)
 real(kind=kind_real), allocatable :: lats(:)
 character(len=1024) :: filename
 integer, allocatable :: istart2(:), icount2(:)
 integer, allocatable :: istart3(:), icount3(:)
 integer, allocatable :: istarte(:), icounte(:)
 contains
  procedure :: create
  procedure :: delete
  procedure :: write
  final :: dummy_final
end type fv3jedi_io_latlon

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

! Arguments
class(fv3jedi_io_latlon),  intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

! Locals
character(len=:), allocatable :: str

! Local copy of geometry things
self%isc  = geom%isc
self%iec  = geom%iec
self%jsc  = geom%jsc
self%jec  = geom%jec
self%isd  = geom%isd
self%ied  = geom%ied
self%jsd  = geom%jsd
self%jed  = geom%jed
self%npx  = geom%npx
self%npy  = geom%npy
self%npz  = geom%npz
self%comm = geom%f_comm

call create_latlon(self, geom)

! Naming convention for the file
! For ensemble methods switch out member template
self%filename = 'Data/fv3jedi.latlon.'
call conf%get_or_die("filename",str)
call swap_name_member(conf, str)
call add_iteration(conf, str)
self%filename = str
deallocate(str)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(fv3jedi_io_latlon), intent(inout) :: self

call delete_latlon(self)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine write(self, vdate, fields)

! Arguments
class(fv3jedi_io_latlon),  intent(inout) :: self
type(datetime),            intent(in)    :: vdate
type(fv3jedi_field),       intent(in)    :: fields(:)

! Write metadata
! --------------
call write_latlon_metadata(self, vdate)

! Write fields
! ------------
call write_latlon_fields(self, fields)

end subroutine write

! --------------------------------------------------------------------------------------------------

subroutine create_latlon(self, geom)

! Arguments
type(fv3jedi_io_latlon), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom

! Locals
integer :: color
real(kind=kind_real) :: dx, dy
integer :: i, ierr


! Create lat lon grid and interpolation object for going from cube to lat-lon
! --------------------------------------------------------------------------
!Maximum of 12 IO processors
if (self%comm%size() >= 12) then
  self%layout(1) = 12
  self%layout(2) = 1
elseif (self%comm%size() >= 6) then
  self%layout(1) = 6
  self%layout(2) = 1
else
  call abor1_ftn("create_latlon error: fewer than 6 npes not anticipated")
endif
self%npes = self%layout(1) * self%layout(2)

!Since the lat lon grid is only for IO purposes it is only
!defined on a subset of PEs - those that will do the writing.
!This is generally more efficient than having many PEs trying
!to write to the same file.

self%nxg = 4*(self%npx - 1)
self%nyg = 2*(self%npy - 1) + 1
self%nx = 0
self%ny = 0

color = MPI_UNDEFINED

if (self%comm%rank() <= self%npes-1) then

  !Split communicator
  color = int(self%comm%communicator() / 2)

  self%thispe = .true.

  !Resolution
  dx = 360.0_kind_real / (real(self%nxg,kind_real) - 1.0_kind_real)
  dy = 180.0_kind_real / (real(self%nyg,kind_real) - 1.0_kind_real)

  self%nx = self%nxg / self%layout(1)
  self%ny = self%nyg / self%layout(2)

  allocate(self%lons(self%nx))
  allocate(self%lats(self%ny))

  !Each processor has subset of lons
  self%lons(1) = dx * self%nx * self%comm%rank()
  do i = 2,self%nx
    self%lons(i) = self%lons(i-1) + dx
  enddo

  !Each processor has all lats
  self%lats(1) = -90.0_kind_real
  do i = 2,self%ny
    self%lats(i) = self%lats(i-1) + dy
  enddo

else

  self%nx = 0
  self%ny = 0
  allocate(self%lons(self%nx))
  allocate(self%lats(self%ny))

endif

! Initialize bump interpolator
call self%bumpinterp%init(self%comm, afunctionspace_in=geom%afunctionspace, lon1d_out=self%lons, &
                          lat1d_out=self%lats, nl0=self%npz)

!IO communicator
!call MPI_Comm_split(self%comm%communicator(), color, self%comm%rank(), self%llcomm, ierr)

! NC arrays
allocate(self%istart3(4),self%icount3(4))
allocate(self%istarte(4),self%icounte(4))
allocate(self%istart2(3),self%icount2(3))
self%istart3(1) = self%nx * self%comm%rank() + 1;  self%icount3(1) = self%nx
self%istart3(2) = 1;                                     self%icount3(2) = self%ny
self%istart3(3) = 1;                                     self%icount3(3) = self%npz
self%istart3(4) = 1;                                     self%icount3(4) = 1
self%istart2(1) = self%istart3(1);                     self%icount2(1) = self%icount3(1)
self%istart2(2) = self%istart3(2);                     self%icount2(2) = self%icount3(2)
self%istart2(3) = self%istart3(4);                     self%icount2(3) = self%icount3(4)

self%istarte = self%istarte
self%icounte = self%icount3
self%istarte(3) = self%npz + 1

end subroutine create_latlon

! --------------------------------------------------------------------------------------------------

subroutine delete_latlon(self)

!Arguments
type(fv3jedi_io_latlon), intent(inout) :: self  !< LatLon Geometry

if (self%thispe) then
  deallocate(self%lons)
  deallocate(self%lats)

  deallocate ( self%istart2, self%icount2 )
  deallocate ( self%istart3, self%icount3 )
  deallocate ( self%istarte, self%icounte )
endif

call self%bumpinterp%delete()

end subroutine delete_latlon

! --------------------------------------------------------------------------------------------------

subroutine write_latlon_metadata(self, vdate)

!Arguments
type(fv3jedi_io_latlon), intent(inout) :: self
type(datetime),          intent(in)    :: vdate

integer :: date(6)
integer :: idate, isecs
character(len=64)   :: datefile

integer :: ncid, varid(2)
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid

! Current date
call datetime_to_ifs(vdate, idate, isecs)
date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

! Append filename with the curent datetime from fv3jedi
write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)
self%filename = trim(self%filename)//trim(datefile)//trim("z.nc4")

call nccheck( nf90_create( self%filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                           comm = self%comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

! Create dimensions
call nccheck ( nf90_def_dim(ncid, "lon" , self%nxg , x_dimid), "nf90_def_dim lon" )
call nccheck ( nf90_def_dim(ncid, "lat" , self%nyg , y_dimid), "nf90_def_dim lat" )
call nccheck ( nf90_def_dim(ncid, "lev" , self%npz   , z_dimid), "nf90_def_dim lev" )
call nccheck ( nf90_def_dim(ncid, "edge", self%npz+1 , e_dimid), "nf90_def_dim edge")
call nccheck ( nf90_def_dim(ncid, "time", 1          , t_dimid), "nf90_def_dim time")

! Define fields to be written
call nccheck( nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, varid(1)), "nf90_def_var lon" )
call nccheck( nf90_put_att(ncid, varid(1), "long_name", "longitude") )
call nccheck( nf90_put_att(ncid, varid(1), "units", "degrees_east") )

call nccheck( nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, varid(2)), "nf90_def_var lat" )
call nccheck( nf90_put_att(ncid, varid(2), "long_name", "latitude") )
call nccheck( nf90_put_att(ncid, varid(2), "units", "degrees_north") )

call nccheck( nf90_enddef(ncid), "nf90_enddef" )

if (self%thispe) then

  ! Write fields
  call nccheck( nf90_put_var( ncid, varid(1), self%lons, &
                              start = self%istart2(1:1), count = self%icount2(1:1) ), "nf90_put_var lons" )
  call nccheck( nf90_put_var( ncid, varid(2), self%lats, &
                              start = self%istart2(2:2), count = self%icount2(2:2) ), "nf90_put_var lats" )

endif

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

!Let LatLon PEs catch up
call self%comm%barrier()

end subroutine write_latlon_metadata

! --------------------------------------------------------------------------------------------------

subroutine write_latlon_fields(self, fields)

!Arguments
type(fv3jedi_io_latlon), target, intent(inout)   :: self
type(fv3jedi_field),             intent(in)      :: fields(:)

integer :: var, ji, jj, jk, llngrid, ii, i, j, k, n
real(kind=kind_real), allocatable :: llfield(:,:,:)

integer :: ncid, varid
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid
integer, target  :: dimids3(4), dimids2(3), dimidse(4)
integer, pointer :: dimids(:), istart(:), icount(:)
real(kind_real),allocatable :: array_with_halo(:,:,:)

! Loop over fields
! ----------------
do var = 1, size(fields)

  ! Only certain fields can be written for now
  ! ------------------------------------------
  if ( trim(fields(var)%space) == 'magnitude' .and. &
       trim(fields(var)%horizontal_stagger_location) == 'center' .and. &
       .not. fields(var)%kind == 'integer' ) then

    ! Interpolate the field to the lat-lon grid
    ! -----------------------------------------
    allocate(llfield(1:self%nx,1:self%ny,1:self%npz))
    llfield = 0.0_kind_real
    llngrid = self%nx*self%ny

    ! Copy field in array with halo points
    allocate(array_with_halo(self%isd:self%ied,self%jsd:self%jed,1:fields(var)%npz))
    array_with_halo = 0.0_kind_real
    array_with_halo(self%isc:self%iec,self%jsc:self%jec,1:fields(var)%npz) = &
       fields(var)%array(self%isc:self%iec,self%jsc:self%jec,1:fields(var)%npz)

    ! Interpolate
    call self%bumpinterp%apply(array_with_halo,llfield(:,:,1:fields(var)%npz))

    !Open file
    call nccheck( nf90_open( self%filename, ior(NF90_WRITE, NF90_MPIIO), ncid, &
                            comm = self%comm%communicator(), info = MPI_INFO_NULL), "nf90_open" )

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
    if (fields(var)%npz == self%npz) then
      dimids => dimids3
      istart => self%istart3
      icount => self%icount3
    elseif (fields(var)%npz == 1) then
      dimids => dimids2
      istart => self%istart2
      icount => self%icount2
    elseif (fields(var)%npz == self%npz+1) then
      dimids => dimidse
      istart => self%istarte
      icount => self%icounte
    endif

    if (associated(dimids)) then

      ! Write field to the file
      call nccheck( nf90_def_var( ncid, trim(fields(var)%io_name), NF90_DOUBLE, dimids, varid), &
                    "nf90_def_var"//trim(fields(var)%io_name) )
      call nccheck( nf90_put_att(ncid, varid, "long_name", trim(fields(var)%long_name) ), "nf90_put_att" )
      call nccheck( nf90_put_att(ncid, varid, "units"    , trim(fields(var)%units)     ), "nf90_put_att" )
      call nccheck( nf90_enddef(ncid), "nf90_enddef" )

      if (self%thispe) then
       call nccheck( nf90_put_var( ncid, varid, llfield, start = istart, count = icount), &
                     "nf90_put_var"//trim(fields(var)%io_name) )
      endif

    endif

    ! Close file
    call nccheck( nf90_close(ncid), "nf90_close" )

    !Let LatLon PEs catch up
    call self%comm%barrier()

    deallocate(llfield)
    deallocate(array_with_halo)
  endif

enddo

end subroutine write_latlon_fields

! --------------------------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_io_latlon), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_latlon_mod
