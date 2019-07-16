! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_geos_mod

use config_mod
use datetime_mod
use iso_c_binding

use fckit_mpi_module

use fv3jedi_constants_mod,    only: rad2deg
use fv3jedi_geom_mod,         only: fv3jedi_geom
use fv3jedi_field_mod,        only: fv3jedi_field
use fv3jedi_kinds_mod,        only: kind_real
use fv3jedi_netcdf_utils_mod, only: nccheck

use mpi
use netcdf

implicit none
private
public fv3jedi_io_geos

type fv3jedi_io_geos
 logical :: iam_io_proc        !Flag for procs doing IO
 type(fckit_mpi_comm) :: ccomm !Component communicator
 integer :: tcomm, ocomm       !Communicator for each tile and for output
 integer :: trank, tsize       !Tile come info
 integer :: crank, csize       !Component comm info
 integer :: orank, osize       !Output comm info
 integer :: vindex
 integer, allocatable :: ncid(:)
 logical :: tiledim
 character(len=10) :: filetype
 character(len=4) :: XdimVar, YdimVar
 integer :: ncdim3, ncdim2
 integer, allocatable :: istart3(:), icount3(:)
 integer, allocatable :: istart2(:), icount2(:)
 integer :: numfiles
  contains
  procedure :: create
  procedure :: delete
  procedure :: read_time
  procedure :: read_fields
  procedure :: write_all
  final     :: dummy_final
end type fv3jedi_io_geos


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, accesstype, filetype, filename, tiledim)

implicit none

class(fv3jedi_io_geos),     intent(inout) :: self
type(fv3jedi_geom),         intent(in)    :: geom
character(len=*),           intent(in)    :: accesstype
character(len=10),          intent(in)    :: filetype
character(len=*), optional, intent(in)    :: filename(:)
integer, optional,          intent(in)    :: tiledim

integer :: ierr, ncstat, nf
integer :: im,jm,lm,nm,dimid
integer :: tileoffset

! Tile dimension in file
self%tiledim = .true.
if (present(tiledim)) then
  self%tiledim = tiledim == 1
endif

! Component communicator
self%ccomm = fckit_mpi_comm()
self%csize = self%ccomm%size()
self%crank = self%ccomm%rank()

! Split comm and set flag for IO procs
self%iam_io_proc = .true.

if (self%csize > 6) then

  ! Tile communicator
  call mpi_comm_split(self%ccomm%communicator(), geom%ntile, self%ccomm%rank(), self%tcomm, ierr)
  call mpi_comm_rank(self%tcomm, self%trank, ierr)
  call mpi_comm_size(self%tcomm, self%tsize, ierr)

  ! Write communicator
  call mpi_comm_split(self%ccomm%communicator(), self%trank, geom%ntile, self%ocomm, ierr)
  call mpi_comm_rank(self%ocomm, self%orank, ierr)
  call mpi_comm_size(self%ocomm, self%osize, ierr)

  if (self%trank .ne. 0) self%iam_io_proc = .false.

else

  ! For 6 cores output comm is componenent comm
  call mpi_comm_dup(self%ccomm%communicator(), self%ocomm, ierr)

endif

! Allocatables based on type of file
self%filetype = filetype
if (self%filetype == 'geos') then
  self%numfiles = 1
elseif (self%filetype == 'geos-rst') then
  self%numfiles = 3
else
  call abor1_ftn("fv3jedi_io_geos_mod.create support filetpe geos or geos-rst only")
endif
allocate(self%ncid(self%numfiles))

! Prepare the file handle for later access
if (self%iam_io_proc) then

  if (trim(accesstype) == 'read') then

    ! Open the file for reading
    do nf = 1,self%numfiles
      call nccheck ( nf90_open(trim(filename(nf)), NF90_NOWRITE, self%ncid(nf)), "nf90_open"//trim(filename(nf)) )
    enddo

    ! Get dimensions, XDim,YDim,lev,time
    ncstat = nf90_inq_dimid(self%ncid(1), "Xdim", dimid)
    if(ncstat /= nf90_noerr) &
    ncstat = nf90_inq_dimid(self%ncid(1), "lon", dimid)
    if(ncstat /= nf90_noerr) &
    call abor1_ftn("Failed to find Xdim or lon in GEOS read")

    call nccheck ( nf90_inquire_dimension(self%ncid(1), dimid, len = im), "nf90_inquire_dimension Xdim/lon" )

    ncstat = nf90_inq_dimid(self%ncid(1), "Ydim", dimid)
    if(ncstat /= nf90_noerr) &
    ncstat = nf90_inq_dimid(self%ncid(1), "lat", dimid)
    if(ncstat /= nf90_noerr) &
    call abor1_ftn("Failed to find Ydim or lat in GEOS read")

    call nccheck ( nf90_inquire_dimension(self%ncid(1), dimid, len = jm), "nf90_inquire_dimension YDim/lat" )

    call nccheck ( nf90_inq_dimid(self%ncid(1), "lev",  dimid), "nf90_inq_dimid lev" )
    call nccheck ( nf90_inquire_dimension(self%ncid(1), dimid, len = lm), "nf90_inquire_dimension lev" )

    call nccheck ( nf90_inq_dimid(self%ncid(1), "time", dimid), "nf90_inq_dimid time" )
    call nccheck ( nf90_inquire_dimension(self%ncid(1), dimid, len = nm), "nf90_inquire_dimension time" )

    ! GEOS can use concatenated tiles or tile as a dimension
    if ( (im == geom%npx-1) .and. (jm == 6*(geom%npy-1) ) ) then
      tileoffset = (geom%ntile-1)*(jm/geom%ntiles)
      self%ncdim3 = 4
      self%ncdim2 = 3
    elseif ( (im == geom%npx-1) .and. (jm == geom%npy-1 ) ) then
      tileoffset = 0
      self%ncdim3 = 5
      self%ncdim2 = 4
    else
      call abor1_ftn("fv3jedi_io_geos_mod.create: dimension mismatch between geometry and file")
    endif

    if (geom%npz .ne. lm) &
    call abor1_ftn("fv3jedi_io_geos_mod.create: level mismatch between geometry and file")

  elseif (trim(accesstype) == 'write') then

    ! Choose dimension based on tiledim
    if ( .not. self%tiledim ) then
      tileoffset = (geom%ntile-1)*(6*(geom%npy-1)/geom%ntiles)
      self%ncdim3 = 4
      self%ncdim2 = 3
      self%XdimVar = 'lon'
      self%YdimVar = 'lat'
    else
      tileoffset = 0
      self%ncdim3 = 5
      self%ncdim2 = 4
      self%XdimVar = 'Xdim'
      self%YdimVar = 'Ydim'
    endif

  endif

  allocate(self%istart3(self%ncdim3),self%icount3(self%ncdim3))
  allocate(self%istart2(self%ncdim2),self%icount2(self%ncdim2))

  ! Create local to this proc start/count
  if (self%ncdim3 == 5) then
    self%istart3(1) = 1;           self%icount3(1) = geom%npx-1  !X
    self%istart3(2) = 1;           self%icount3(2) = geom%npy-1  !Y
    self%istart3(3) = geom%ntile;  self%icount3(3) = 1           !Tile
    self%istart3(4) = 1;           self%icount3(4) = 1           !Lev
    self%istart3(5) = 1;           self%icount3(5) = 1           !Time
    self%istart2(1) = 1;           self%icount2(1) = geom%npx-1
    self%istart2(2) = 1;           self%icount2(2) = geom%npy-1
    self%istart2(3) = geom%ntile;  self%icount2(3) = 1
    self%istart2(4) = 1;           self%icount2(4) = 1
    self%vindex = 4
  elseif (self%ncdim3 == 4) then
    self%istart3(1) = 1;              self%icount3(1) = geom%npx-1
    self%istart3(2) = tileoffset+1;   self%icount3(2) = geom%npy-1
    self%istart3(3) = 1;              self%icount3(3) = 1
    self%istart3(4) = 1;              self%icount3(4) = 1
    self%istart2(1) = 1;              self%icount2(1) = geom%npx-1
    self%istart2(2) = tileoffset+1;   self%icount2(2) = geom%npy-1
    self%istart2(3) = 1;              self%icount2(3) = 1
    self%vindex = 3
  else
    call abor1_ftn("fv3jedi_io_geos_mod.create: ncdim3 set incorrectly")
  endif

endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none

class(fv3jedi_io_geos), intent(inout) :: self

integer :: ierr, nf

if (self%iam_io_proc) then
  !Close the file
  do nf = 1,self%numfiles
    call nccheck ( nf90_close(self%ncid(nf)), "nf90_close" )
  enddo
endif
deallocate(self%ncid)

! Release split comms
if (self%csize > 6) call MPI_Comm_free(self%tcomm, ierr)
call MPI_Comm_free(self%ocomm, ierr)

! Deallocate start/count
if (self%iam_io_proc) then
  deallocate ( self%istart2, self%icount2 )
  deallocate ( self%istart3, self%icount3 )
endif

end subroutine delete

! ------------------------------------------------------------------------------

subroutine read_time(self, vdate)

implicit none

class(fv3jedi_io_geos), intent(in)    :: self
type(datetime),         intent(inout) :: vdate

integer :: varid, date(6), intdate, inttime, idate, isecs
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: cidate, cisecs

idate = 0
isecs = 0

if (self%iam_io_proc) then

  ! Get time attributes
  call nccheck ( nf90_inq_varid(self%ncid(1), "time", varid), "nf90_inq_varid time" )
  call nccheck ( nf90_get_att(self%ncid(1), varid, "begin_date", intdate), "nf90_get_att begin_date" )
  call nccheck ( nf90_get_att(self%ncid(1), varid, "begin_time", inttime), "nf90_get_att begin_time" )

  ! Pad with leading zeros if need be
  write(cdate,"(I0.8)") intdate
  write(ctime,"(I0.6)") inttime

  ! Back to integer
  read(cdate(1:4),*) date(1)
  read(cdate(5:6),*) date(2)
  read(cdate(7:8),*) date(3)
  read(ctime(1:2),*) date(4)
  read(ctime(3:4),*) date(5)
  read(ctime(5:6),*) date(6)

  ! To idate/isecs for Jedi
  idate = date(1)*10000 + date(2)*100 + date(3)
  isecs = date(4)*3600  + date(5)*60  + date(6)

endif

! Set the object date from the date of the file
call self%ccomm%broadcast(idate,0)
call self%ccomm%broadcast(isecs,0)
cidate = idate
cisecs = isecs
call datetime_from_ifs(vdate, cidate, cisecs)

end subroutine read_time

! ------------------------------------------------------------------------------

subroutine read_fields(self, geom, fields)

implicit none

class(fv3jedi_io_geos), target, intent(in)    :: self
type(fv3jedi_geom),             intent(in)    :: geom
type(fv3jedi_field),            intent(inout) :: fields(:)

integer :: varid, n, lev, fieldncid
integer, pointer :: istart(:), icount(:)
integer, allocatable, target :: istart3(:)
real(kind=kind_real), allocatable :: arrayg(:,:), delp(:,:,:)


! Array for level of whole tile
allocate(arrayg(1:geom%npx-1,1:geom%npy-1))

do n = 1,size(fields)

  ! For GEOS restarts specify file for variables
  fieldncid = 1
  if (self%filetype == 'geos-rst') then
    select case (trim(fields(n)%short_name))
    case("U","V","W","PT","PKZ","PE","DZ")
      fieldncid = 1
    case("Q","QILS","QICN","QLLS","QLCN","CLLS","CLCN")
      fieldncid = 2
    case("PHIS")
      fieldncid = 3
    case default
      call abor1_ftn("fv3jedi_io_geos_mod.read_fields: no geos restart for "//trim(fields(n)%short_name))
    end select
  endif


  if (self%iam_io_proc) then

    ! Local counts for 3 to allow changing start point
    if (.not.allocated(istart3)) then
      allocate(istart3(size(self%istart3)))
      istart3 = self%istart3
    endif

    if (fields(n)%npz == 1) then
      istart => self%istart2; icount => self%icount2
    elseif (fields(n)%npz > 1) then
      istart => istart3; icount => self%icount3
    endif

  endif

  !If ps then switch to delp
  if (trim(fields(n)%fv3jedi_name) == 'ps') then
    fields(n)%short_name = 'delp'
    fields(n)%npz = geom%npz
    if (self%iam_io_proc) then
      istart => istart3; icount => self%icount3
    endif
    allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  endif


  do lev = 1,fields(n)%npz

    arrayg = 0.0_kind_real

    if (self%iam_io_proc) then

      !Set counter to current level
      istart3(self%vindex) = lev

      call nccheck ( nf90_inq_varid (self%ncid(fieldncid), trim(fields(n)%short_name), varid), &
                    "nf90_inq_varid "//trim(fields(n)%short_name) )
      call nccheck ( nf90_get_var( self%ncid(fieldncid), varid, arrayg, istart, icount), &
                    "nf90_get_var "//trim(fields(n)%short_name) )

    endif

    if (trim(fields(n)%fv3jedi_name) .ne. 'ps') then
      if (self%csize > 6) then
        call scatter_tile(geom, self%tcomm, 1, arrayg, fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev))
      else
        fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev) = arrayg(geom%isc:geom%iec,geom%jsc:geom%jec)
      endif
    else
      if (self%csize > 6) then
        call scatter_tile(geom, self%tcomm, 1, arrayg, delp(geom%isc:geom%iec,geom%jsc:geom%jec,lev))
      else
        delp(geom%isc:geom%iec,geom%jsc:geom%jec,lev) = arrayg(geom%isc:geom%iec,geom%jsc:geom%jec)
      endif
    endif

  enddo

  if (self%iam_io_proc) then
    nullify(istart,icount)
    if (n == size(fields)) deallocate(istart3)
  endif

  if (trim(fields(n)%fv3jedi_name) == 'ps') then
    fields(n)%short_name = 'ps'
    fields(n)%npz = 1
    fields(n)%array(:,:,1) = sum(delp,3)
    deallocate(delp)
  endif

enddo

deallocate(arrayg)

end subroutine read_fields

! ------------------------------------------------------------------------------

subroutine write_all(self, geom, fields, c_conf, vdate)

implicit none

class(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_geom),             intent(in)    :: geom
type(fv3jedi_field),            intent(in)    :: fields(:)
type(c_ptr),                    intent(in)    :: c_conf
type(datetime),                 intent(in)    :: vdate

! Locals
character(len=255), allocatable :: filename(:)
character(len=255) :: datapath
character(len=15)  :: datefile
character(len=8)   :: date8s, cubesize
character(len=6)   :: time6s
integer, allocatable :: vc(:)
integer :: n, nf, lev, date8, time6, ymult, dfend, fieldncid
integer :: varid(1000), date(6)
integer(kind=c_int) :: idate, isecs
integer :: x_dimid, y_dimid, n_dimid, z_dimid, ze_dimid, t_dimid, z4_dimid, c_dimid, o_dimid
integer, allocatable :: dimidsv(:,:), dimidsg(:,:), dimids2(:,:), dimids3(:,:), dimids3e(:,:), dimids3_4(:,:)
integer, target, allocatable :: istart3(:)
integer, allocatable :: dimids(:), intarray(:)
real(kind=kind_real), allocatable :: arrayg(:,:), realarray(:)
integer, pointer :: istart(:), icount(:)
logical :: write_lev4 = .false.


! Whole level of tile array
allocate(arrayg(1:geom%npx-1,1:geom%npy-1))

! Define all the variables and meta data to be written
! ----------------------------------------------------
if (self%iam_io_proc) then

  ! Place to save restarts
  datapath = "Data/"
  if (config_element_exists(c_conf,"datapath")) then
     datapath = config_get_string(c_conf,len(datapath),"datapath")
  endif

  ! Base filename
  allocate(filename(self%numfiles))
  if (self%filetype == 'geos') then
    filename(1) = 'geos.'
    if (config_element_exists(c_conf,"filename")) then
       filename(1) = config_get_string(c_conf,len(filename(1)),"filename")
    endif
    dfend = 15
  elseif (self%filetype == 'geos-rst') then
    filename(1) = 'fvcore_internal_rst.'
    filename(2) = 'moist_internal_rst.'
    filename(3) = 'surf_import_rst.'
    if (config_element_exists(c_conf,"filename-fvcore")) then
       filename(1) = config_get_string(c_conf,len(filename(1)),"filename-fvcore")
    endif
    if (config_element_exists(c_conf,"filename-moist")) then
       filename(2) = config_get_string(c_conf,len(filename(2)),"filename-moist")
    endif
    if (config_element_exists(c_conf,"filename-surf")) then
       filename(3) = config_get_string(c_conf,len(filename(3)),"filename-surf")
    endif
    dfend = 11
  endif

  ! Append with the date
  call datetime_to_ifs(vdate, idate, isecs)
  date(1) = idate/10000
  date(2) = idate/100 - date(1)*100
  date(3) = idate - (date(1)*10000 + date(2)*100)
  date(4) = isecs/3600
  date(5) = (isecs - date(4)*3600)/60
  date(6) = isecs - (date(4)*3600 + date(5)*60)
  write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)

  write(date8s,'(I4,I0.2,I0.2)')   date(1),date(2),date(3)
  write(time6s,'(I0.2,I0.2,I0.2)') date(4),date(5),date(6)
  read(date8s,*)  date8
  read(time6s,*)  time6

  write(cubesize,'(I8)') geom%npx-1

  ! DimId arrays
  allocate(dimidsv(self%numfiles,1))
  allocate(dimidsg(self%numfiles,self%ncdim2-1))
  allocate(dimids2(self%numfiles,self%ncdim2))
  allocate(dimids3(self%numfiles,self%ncdim3))
  allocate(dimids3e(self%numfiles,self%ncdim3))
  allocate(dimids3_4(self%numfiles,self%ncdim3))

  allocate(vc(self%numfiles))

  do nf = 1,self%numfiles

    filename(nf) = trim(datapath)//trim(filename(nf))//trim(datefile(1:dfend))//trim("z.nc4")

    ! Create and open the file for parallel write
    call nccheck( nf90_create( trim(filename(nf)), ior(NF90_NETCDF4, NF90_MPIIO), self%ncid(nf), &
                               comm = self%ocomm, info = MPI_INFO_NULL), "nf90_create" )

    ! Create dimensions
    ymult = 1
    if (.not. self%tiledim) ymult = 6

    call nccheck ( nf90_def_dim(self%ncid(nf), trim(self%XdimVar), geom%npx-1,  x_dimid), "nf90_def_dim "//trim(self%XdimVar) )
    call nccheck ( nf90_def_dim(self%ncid(nf), trim(self%YdimVar), ymult*(geom%npy-1),  y_dimid), "nf90_def_dim "//trim(self%YdimVar) )
    if (self%tiledim) call nccheck ( nf90_def_dim(self%ncid(nf), "nf",   geom%ntiles, n_dimid), "nf90_def_dim nf"   )
    call nccheck ( nf90_def_dim(self%ncid(nf), "lev",  geom%npz,    z_dimid), "nf90_def_dim lev"  )
    call nccheck ( nf90_def_dim(self%ncid(nf), "edge",  geom%npz+1, ze_dimid), "nf90_def_dim edge"  )
    call nccheck ( nf90_def_dim(self%ncid(nf), "time", 1,           t_dimid), "nf90_def_dim time" )

    do n = 1,size(fields)
      if (fields(n)%npz == 4) then
        write_lev4 = .true.
      endif
    enddo
    if (write_lev4) call nccheck ( nf90_def_dim(self%ncid(nf), "lev4", 4, z4_dimid), "nf90_def_dim lev"  )

    !Needed by GEOS for ingesting cube sphere field
    call nccheck ( nf90_def_dim(self%ncid(nf), "ncontact", 4, c_dimid), "nf90_def_dim ncontact" )
    call nccheck ( nf90_def_dim(self%ncid(nf), "orientationStrLen", 5, o_dimid), "nf90_def_dim orientationStrLend" )

    dimidsv(nf,:)   =  (/ z_dimid /)

    if ( self%tiledim ) then
      dimidsg(nf,:)   =  (/ x_dimid, y_dimid, n_dimid /)
      dimids2(nf,:)   =  (/ x_dimid, y_dimid, n_dimid, t_dimid /)
      dimids3(nf,:)   =  (/ x_dimid, y_dimid, n_dimid, z_dimid, t_dimid /)
      dimids3e(nf,:)  =  (/ x_dimid, y_dimid, n_dimid, ze_dimid, t_dimid /)
      dimids3_4(nf,:) =  (/ x_dimid, y_dimid, n_dimid, z4_dimid, t_dimid /)
    else
      dimidsg(nf,:)   =  (/ x_dimid, y_dimid /)
      dimids2(nf,:)   =  (/ x_dimid, y_dimid, t_dimid /)
      dimids3(nf,:)   =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)
      dimids3e(nf,:)  =  (/ x_dimid, y_dimid, ze_dimid, t_dimid /)
      dimids3_4(nf,:) =  (/ x_dimid, y_dimid, z4_dimid, t_dimid /)
    endif

    ! Define fields to be written (geom)
    vc(nf)=0;

    if (self%tiledim) then
      vc(nf)=vc(nf)+1;
      call nccheck( nf90_def_var(self%ncid(nf), "nf", NF90_INT, n_dimid, varid(vc(nf))), "nf90_def_var nf" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "cubed-sphere face") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "axis", "e") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "grads_dim", "e") )
    endif

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), trim(self%XdimVar), NF90_DOUBLE, x_dimid, varid(vc(nf))), "nf90_def_var "//trim(self%XdimVar) )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "Fake Longitude for GrADS Compatibility") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "degrees_east") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), trim(self%YdimVar), NF90_DOUBLE, y_dimid, varid(vc(nf))), "nf90_def_var "//trim(self%YdimVar) )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "Fake Latitude for GrADS Compatibility") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "degrees_north") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), "lons", NF90_DOUBLE, dimidsg(nf,:), varid(vc(nf))), "nf90_def_var lons" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "longitude") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "degrees_east") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), "lats", NF90_DOUBLE, dimidsg(nf,:), varid(vc(nf))), "nf90_def_var lats" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "latitude") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "degrees_north") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), "lev", NF90_DOUBLE, z_dimid, varid(vc(nf))), "nf90_def_var lev" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "vertical level") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "layer") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "positive", "down") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "coordinate", "eta") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "standard_name", "model_layers") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), "edge", NF90_DOUBLE, ze_dimid, varid(vc(nf))), "nf90_def_var edge" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "vertical level edges") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "units", "layer") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "positive", "down") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "coordinate", "eta") )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "standard_name", "model_layers") )

    vc(nf)=vc(nf)+1;
    call nccheck( nf90_def_var(self%ncid(nf), "time", NF90_INT, t_dimid, varid(vc(nf))), "nf90_def_var time" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "long_name", "time"), "nf90_def_var time long_name" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "begin_date", date8), "nf90_def_var time begin_date" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "begin_time", time6), "nf90_def_var time begin_time" )

    vc(nf)=vc(nf)+1; !(Needed by GEOS to ingest cube sphere analysis)
    call nccheck( nf90_def_var(self%ncid(nf), "cubed_sphere", NF90_CHAR, varid(vc(nf))), "nf90_def_var cubed_sphere" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "grid_mapping_name", "gnomonic cubed-sphere"), "nf90_def_var time grid_mapping_name" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "file_format_version", "2.90"), "nf90_def_var time file_format_version" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "additional_vars", "contacts,orientation,anchor"), "nf90_def_var time additional_vars" )
    call nccheck( nf90_put_att(self%ncid(nf), varid(vc(nf)), "gridspec_file", "C"//trim(cubesize)//"_gridspec.nc4"), "nf90_def_var gridspec_file" )

    !vc(nf)=vc(nf)+1; !(Needed by GEOS to ingest cube sphere analysis)
    !call nccheck( nf90_def_var(self%ncid(nf), "ncontact", NF90_INT, varid(vc(nf))), "nf90_def_var ncontact" )

  enddo

  ! Define fields to be written
  do n = 1,size(fields)

    ! For GEOS restarts specify file for variables
    fieldncid = 1
    if (self%filetype == 'geos-rst') then
      select case (trim(fields(n)%short_name))
      case("U","V","W","PT","PKZ","PE","DZ")
        fieldncid = 1
      case("Q","QILS","QICN","QLLS","QLCN","CLLS","CLCN")
        fieldncid = 2
      case("PHIS")
        fieldncid = 3
      case default
        call abor1_ftn("fv3jedi_io_geos_mod.write_all: no geos restart for "//trim(fields(n)%short_name))
      end select
    endif

    if (fields(n)%npz == 1) then
      allocate(dimids(self%ncdim2))
      dimids = dimids2(fieldncid,:)
    elseif (fields(n)%npz == geom%npz) then
      allocate(dimids(self%ncdim3))
      dimids = dimids3(fieldncid,:)
    elseif (fields(n)%npz == geom%npz+1) then
      allocate(dimids(self%ncdim3))
      dimids = dimids3e(fieldncid,:)
    elseif (fields(n)%npz == 4) then
      allocate(dimids(self%ncdim3))
      dimids = dimids3_4(fieldncid,:)
    else
      call abor1_ftn("write_geos: vertical dimension not supported")
    endif

    vc(fieldncid)=vc(fieldncid)+1
    call nccheck( nf90_def_var(self%ncid(fieldncid), trim(fields(n)%short_name), NF90_DOUBLE, dimids, varid(vc(fieldncid))), &
                  "nf90_def_var"//trim(fields(n)%short_name)   )
    call nccheck( nf90_put_att(self%ncid(fieldncid), varid(vc(fieldncid)), "long_name"    , trim(fields(n)%long_name) ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid(fieldncid), varid(vc(fieldncid)), "units"        , trim(fields(n)%units)     ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid(fieldncid), varid(vc(fieldncid)), "standard_name", trim(fields(n)%long_name) ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid(fieldncid), varid(vc(fieldncid)), "coordinates"  , "lons lats"               ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid(fieldncid), varid(vc(fieldncid)), "grid_mapping" , "cubed_sphere"            ), "nf90_put_att" )

    deallocate(dimids)

  enddo

  do nf = 1,self%numfiles

    ! End define mode
    call nccheck( nf90_enddef(self%ncid(nf)), "nf90_enddef" )

  enddo

endif


! Write Xdim/lon, YDim/lat, and nf (tile)
! ---------------------------------------
if (self%iam_io_proc) then

  do nf = 1,self%numfiles

    vc(nf)=0

    if (self%tiledim) then
      allocate(intarray(6))
      do n = 1,6
        intarray(n) = n
      enddo
      vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), intarray ), "nf90_put_var nf" )
      deallocate(intarray)
    endif

    allocate(realarray(geom%npx-1))
    do n = 1,geom%npx-1
      realarray(n) = real(n,kind_real)
    enddo

    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), realarray ), "nf90_put_var "//trim(self%XdimVar) )
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), realarray ), "nf90_put_var "//trim(self%YdimVar) )

    deallocate(realarray)

  enddo

endif


! Gather longitudes to write
! --------------------------
if (self%csize > 6) then
  call gather_tile(geom, self%tcomm, 1, rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec), arrayg)
else
  arrayg = rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec)
endif

if (self%iam_io_proc) then
  do nf = 1,self%numfiles
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), arrayg, &
                                        start = self%istart2(1:self%ncdim2-1), count = self%icount2(1:self%ncdim2-1) ), "nf90_put_var lons" )
  enddo
endif


! Gather latitudes and write
! --------------------------
if (self%csize > 6) then
  call gather_tile(geom, self%tcomm, 1, rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec), arrayg)
else
  arrayg = rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec)
endif

if (self%iam_io_proc) then
  do nf = 1,self%numfiles
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), arrayg, &
                                        start = self%istart2(1:self%ncdim2-1), count = self%icount2(1:self%ncdim2-1) ), "nf90_put_var lats" )
  enddo
endif


! Write model levels & time
! -------------------------
if (self%iam_io_proc) then
  allocate(intarray(geom%npz))
  do n = 1,geom%npz
    intarray(n) = n
  enddo
  do nf = 1,self%numfiles
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), intarray ), "nf90_put_var lev" )
  enddo
  deallocate(intarray)

  allocate(intarray(geom%npz+1))
  do n = 1,geom%npz+1
    intarray(n) = n
  enddo
  do nf = 1,self%numfiles
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), intarray ), "nf90_put_var edge" )
  enddo
  deallocate(intarray)

  do nf = 1,self%numfiles
    vc(nf)=vc(nf)+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc(nf)), 0 ), "nf90_put_var time" )
  enddo
endif

! Cube sphere thing for GEOS
vc=vc+1

! Loop over fields and levels and write fields to file
! ----------------------------------------------------
do n = 1,size(fields)

  ! For GEOS restarts specify file for variables
  fieldncid = 1
  if (self%filetype == 'geos-rst') then
    select case (trim(fields(n)%short_name))
    case("U","V","W","PT","PKZ","PE","DZ")
      fieldncid = 1
    case("Q","QILS","QICN","QLLS","QLCN","CLLS","CLCN")
      fieldncid = 2
    case("PHIS")
      fieldncid = 3
    case default
      call abor1_ftn("fv3jedi_io_geos_mod.write_all: no geos restart for "//trim(fields(n)%short_name))
    end select
  endif


  if (self%iam_io_proc) then

    ! Local counts for 3 to allow changing start point
    if (.not.allocated(istart3)) then
      allocate(istart3(size(self%istart3)))
      istart3 = self%istart3
    endif

    if (fields(n)%npz == 1) then
      istart => self%istart2; icount => self%icount2
    elseif (fields(n)%npz > 1) then
      istart => istart3; icount => self%icount3
    endif

    !Up field counter outside of level loop
    vc(fieldncid) = vc(fieldncid) + 1

  endif

  do lev = 1,fields(n)%npz

    if (self%csize > 6) then
      call gather_tile(geom, self%tcomm, 1, fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev), arrayg)
    else
      arrayg = fields(n)%array(:,:,lev)
    endif

    if (self%iam_io_proc) then

      istart3(self%vindex) = lev

      call nccheck( nf90_put_var( self%ncid(fieldncid), varid(vc(fieldncid)), arrayg, start = istart, count = icount ), &
                                  "nf90_put_var "//trim(fields(n)%short_name) )

    endif

  enddo

  if (self%iam_io_proc .and. n==size(fields)) deallocate(istart3)

enddo

if (self%iam_io_proc) &
deallocate ( dimidsv, dimidsg, dimids2, dimids3, dimids3e )

deallocate(arrayg, vc)

end subroutine write_all

! ------------------------------------------------------------------------------

subroutine gather_tile(geom, comm, nlev, array_l, array_g)

implicit none

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec,1:nlev)  ! Local array
real(kind=kind_real), intent(inout) :: array_g(1:geom%npx-1,1:geom%npy-1,1:nlev)            ! Gathered array (only valid on root)

real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)
integer :: comm_size, ierr
integer :: ji, jj, jk, jc, n
integer :: npx_g, npy_g, npx_l, npy_l
integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:)
integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)

!Get comm size
call mpi_comm_size(comm, comm_size, ierr)

!Array of counts and displacement
allocate(counts(comm_size), displs(comm_size))
do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

!Horizontal size for global and local
npx_g = geom%npx-1
npy_g = geom%npy-1
npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(geom%isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
deallocate(counts,displs)

! Pack whole tile array into vector
allocate(vectorcounts(comm_size), vectordispls(comm_size))

!Gather counts and displacement
n = 0
do jc = 1,comm_size
  vectordispls(jc) = n
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
      enddo
    enddo
  enddo
  vectorcounts(jc) = n - vectordispls(jc)
enddo

! Pack local array into vector
allocate(vector_l(npx_l*npy_l*nlev))
n = 0
do jk = 1,nlev
  do jj = geom%jsc,geom%jec
    do ji = geom%isc,geom%iec
      n = n+1
      vector_l(n) = array_l(ji,jj,jk)
    enddo
  enddo
enddo

! Gather the full field
allocate(vector_g(npx_g*npy_g*nlev))
call mpi_gatherv( vector_l, npx_l*npy_l, mpi_double_precision, &
                  vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                  0, comm, ierr)
deallocate(vector_l,vectorcounts,vectordispls)

!Unpack global vector into array
n = 0
do jc = 1,comm_size
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
        array_g(ji,jj,jk) = vector_g(n)
      enddo
    enddo
  enddo
enddo
deallocate(isc_l, iec_l, jsc_l, jec_l)

deallocate(vector_g)

end subroutine gather_tile

! ------------------------------------------------------------------------------

subroutine scatter_tile(geom, comm, nlev, array_g, array_l)

implicit none

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_g(1:geom%npx-1,1:geom%npy-1,nlev)            ! Gathered array (only valid on root)
real(kind=kind_real), intent(inout) :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec,nlev)  ! Local array

real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)
integer :: comm_size, ierr
integer :: ji, jj, jk, jc, n
integer :: npx_g, npy_g, npx_l, npy_l
integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:)
integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)

!Get comm size
call mpi_comm_size(comm, comm_size, ierr)

!Array of counts and displacement
allocate(counts(comm_size), displs(comm_size))
do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

!Horizontal size for global and local
npx_g = geom%npx-1
npy_g = geom%npy-1
npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(geom%isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
deallocate(counts,displs)

! Pack whole tile array into vector
allocate(vector_g(npx_g*npy_g*nlev))
allocate(vectorcounts(comm_size), vectordispls(comm_size))

n = 0
do jc = 1,comm_size
  vectordispls(jc) = n
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
        vector_g(n) = array_g(ji,jj,jk)
      enddo
    enddo
  enddo
  vectorcounts(jc) = n - vectordispls(jc)
enddo
deallocate(isc_l, iec_l, jsc_l, jec_l)

! Scatter tile array to processors
allocate(vector_l(npx_l*npy_l*nlev))

call mpi_scatterv( vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                   vector_l, npx_l*npy_l, mpi_double_precision, &
                   0, comm, ierr )

deallocate(vector_g,vectorcounts,vectordispls)

! Unpack local vector into array
n = 0
do jk = 1,nlev
  do jj = geom%jsc,geom%jec
    do ji = geom%isc,geom%iec
      n = n+1
      array_l(ji,jj,jk) = vector_l(n)
    enddo
  enddo
enddo
deallocate(vector_l)

end subroutine scatter_tile

! ------------------------------------------------------------------------------

subroutine dummy_final(self)
implicit none
type(fv3jedi_io_geos), intent(inout) :: self
end subroutine dummy_final

! ------------------------------------------------------------------------------

end module fv3jedi_io_geos_mod
