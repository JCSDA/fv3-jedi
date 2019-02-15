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
 integer :: ncid, vindex
 integer, allocatable :: istart3(:), icount3(:)
 integer, allocatable :: istart2(:), icount2(:)
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

subroutine create(self, geom, accesstype, filename)

class(fv3jedi_io_geos),     intent(inout) :: self
type(fv3jedi_geom),         intent(in)    :: geom
character(len=*),           intent(in)    :: accesstype
character(len=*), optional, intent(in)    :: filename

integer :: ierr
integer :: im,jm,lm,nm,dimid
integer :: tileoffset, ncdim3, ncdim2

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


! Prepare the file handle for later access
if (self%iam_io_proc) then

  if (trim(accesstype) == 'read') then

    ! Open the file for reading
    call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, self%ncid), "nf90_open"//trim(filename) )

    ! Get dimensions, XDim,YDim,lev,time
    call nccheck ( nf90_inq_dimid(self%ncid, "Xdim", dimid), "nf90_inq_dimid Xdim" )
    call nccheck ( nf90_inquire_dimension(self%ncid, dimid, len = im), "nf90_inquire_dimension Xdim" )
    
    call nccheck ( nf90_inq_dimid(self%ncid, "Ydim", dimid), "nf90_inq_dimid YDim" )
    call nccheck ( nf90_inquire_dimension(self%ncid, dimid, len = jm), "nf90_inquire_dimension YDim" )
    
    call nccheck ( nf90_inq_dimid(self%ncid, "lev",  dimid), "nf90_inq_dimid lev" )
    call nccheck ( nf90_inquire_dimension(self%ncid, dimid, len = lm), "nf90_inquire_dimension lev" )
    
    call nccheck ( nf90_inq_dimid(self%ncid, "time", dimid), "nf90_inq_dimid time" )
    call nccheck ( nf90_inquire_dimension(self%ncid, dimid, len = nm), "nf90_inquire_dimension time" )

    ! GEOS can use concatenated tiles or tile as a dimension
    if ( (im == geom%npx-1) .and. (jm == 6*(geom%npy-1) ) ) then
      tileoffset = (geom%ntile-1)*(jm/geom%ntiles)
      ncdim3 = 4
      ncdim2 = 3
    elseif ( (im == geom%npx-1) .and. (jm == geom%npy-1 ) ) then
      tileoffset = 0
      ncdim3 = 5
      ncdim2 = 4
    else
      call abor1_ftn("fv3jedi_io_geos_mod.create: dimension mismatch between geometry and file")
    endif

    if (geom%npz .ne. lm) &
    call abor1_ftn("fv3jedi_io_geos_mod.create: level mismatch between geometry and file")

  elseif (trim(accesstype) == 'write') then

    ! Always write with tile dimension
    tileoffset = 0
    ncdim3 = 5
    ncdim2 = 4

  endif

  allocate(self%istart3(ncdim3),self%icount3(ncdim3))
  allocate(self%istart2(ncdim2),self%icount2(ncdim2))
  
  ! Create local to this proc start/count
  if (tileoffset == 0) then
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
  else
    self%istart3(1) = 1;              self%icount3(1) = geom%npx-1
    self%istart3(2) = tileoffset+1;   self%icount3(2) = geom%npy-1
    self%istart3(3) = 1;              self%icount3(3) = 1
    self%istart3(4) = 1;              self%icount3(4) = 1
    self%istart2(1) = 1;              self%icount2(1) = geom%npx-1
    self%istart2(2) = tileoffset+1;   self%icount2(2) = geom%npy-1
    self%istart2(3) = 1;              self%icount2(3) = 1
    self%vindex = 3
  endif
  
endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_io_geos), intent(inout) :: self

integer :: ierr

if (self%iam_io_proc) then
  !Close the file
  call nccheck ( nf90_close(self%ncid), "nf90_close" )
endif

! Release split comms
if (self%csize > 6) call MPI_Comm_free(self%tcomm, ierr)
call MPI_Comm_free(self%ocomm, ierr)

! Deallocate start/count
deallocate ( self%istart2, self%icount2 )
deallocate ( self%istart3, self%icount3 )

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
  call nccheck ( nf90_inq_varid(self%ncid, "time", varid), "nf90_inq_varid time" )
  call nccheck ( nf90_get_att(self%ncid, varid, "begin_date", intdate), "nf90_get_att begin_date" )
  call nccheck ( nf90_get_att(self%ncid, varid, "begin_time", inttime), "nf90_get_att begin_time" )
  
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

integer :: varid, n, lev
integer, pointer :: istart(:), icount(:)
integer, allocatable, target :: istart3(:)
real(kind=kind_real), allocatable :: arrayg(:,:)

! Array for level of whole tile
allocate(arrayg(1:geom%npx-1,1:geom%npy-1))

do n = 1,size(fields)

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

  do lev = 1,fields(n)%npz

    arrayg = 0.0_kind_real
 
    if (self%iam_io_proc) then

      !Set counter to current level
      istart3(self%vindex) = lev

      call nccheck ( nf90_inq_varid (self%ncid, trim(fields(n)%short_name), varid), &
                    "nf90_inq_varid "//trim(fields(n)%short_name) )
      call nccheck ( nf90_get_var( self%ncid, varid, arrayg, istart, icount), &
                    "nf90_get_var "//trim(fields(n)%short_name) )

    endif

    if (self%csize > 6) then
      call scatter_tile(geom, self%tcomm, arrayg, fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev))
    else
      fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev) = arrayg(geom%isc:geom%iec,geom%jsc:geom%jec)
    endif

  enddo

  if (self%iam_io_proc) then
    nullify(istart,icount)
    if (n == size(fields)) deallocate(istart3)    
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
character(len=255) :: datapath, filename
character(len=64)  :: datefile
character(len=8)   :: date8s
character(len=6)   :: time6s
integer :: n, lev, date8, time6
integer :: varid(1000), date(6), vc
integer(kind=c_int) :: idate, isecs
integer :: x_dimid, y_dimid, n_dimid, z_dimid, t_dimid
integer :: ndimidsv, ndimidsg, ndimids2, ndimids3
integer, allocatable :: dimidsv(:), dimidsg(:), dimids2(:), dimids3(:)
integer, target, allocatable :: istart3(:)
integer, allocatable :: dimids(:), intarray(:)
real(kind=kind_real), allocatable :: arrayg(:,:), realarray(:)
integer, pointer :: istart(:), icount(:)

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
  filename = 'GEOS.eta.'
  if (config_element_exists(c_conf,"filename")) then
     filename = config_get_string(c_conf,len(filename),"filename")
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
  filename = trim(datapath)//trim(filename)//trim(datefile)//trim("z.nc4")

  write(datefile,'(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)') date(1),"-",date(2),"-",date(3)," "&
                                                    ,date(4),":",date(5),":",date(6)
  
  write(date8s,'(I4,I0.2,I0.2)')   date(1),date(2),date(3)
  write(time6s,'(I0.2,I0.2,I0.2)') date(4),date(5),date(6)
  read(date8s,*)  date8
  read(time6s,*)  time6
   
  ! Create and open the file for parallel write
  call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_MPIIO), self%ncid, &
                             comm = self%ocomm, info = MPI_INFO_NULL), "nf90_create" )

  ! Create dimensions
  call nccheck ( nf90_def_dim(self%ncid, "Xdim", geom%npx-1,  x_dimid), "nf90_def_dim Xdim" )
  call nccheck ( nf90_def_dim(self%ncid, "Ydim", geom%npy-1,  y_dimid), "nf90_def_dim Ydim" )
  call nccheck ( nf90_def_dim(self%ncid, "nf",   geom%ntiles, n_dimid), "nf90_def_dim nf"   )
  call nccheck ( nf90_def_dim(self%ncid, "lev",  geom%npz,    z_dimid), "nf90_def_dim lev"  )
  call nccheck ( nf90_def_dim(self%ncid, "time", 1,           t_dimid), "nf90_def_dim time" )
  
  ! DimId arrays
  ndimidsv = 1
  ndimidsg = 3
  ndimids2 = 4
  ndimids3 = 5

  allocate(dimidsv(ndimidsv))
  dimidsv =  (/ z_dimid /)
  allocate(dimidsg(ndimidsg))
  dimidsg =  (/ x_dimid, y_dimid, n_dimid /)
  allocate(dimids2(ndimids2))
  dimids2 =  (/ x_dimid, y_dimid, n_dimid, t_dimid /)
  allocate(dimids3(ndimids3))
  dimids3 =  (/ x_dimid, y_dimid, n_dimid, z_dimid, t_dimid /)
  
  ! Define fields to be written (geom) 
  vc=1;
  call nccheck( nf90_def_var(self%ncid, "nf", NF90_INT, n_dimid, varid(vc)), "nf90_def_var nf" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "cubed-sphere face") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "axis", "e") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "grads_dim", "e") )

  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "Xdim", NF90_DOUBLE, x_dimid, varid(vc)), "nf90_def_var Xdim" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "Fake Longitude for GrADS Compatibility") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "units", "degrees_east") )

  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "Ydim", NF90_DOUBLE, y_dimid, varid(vc)), "nf90_def_var Ydim" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "Fake Latitude for GrADS Compatibility") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "units", "degrees_north") )

  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "lons", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lons" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "longitude") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "units", "degrees_east") )
  
  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "lats", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lats" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "latitude") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "units", "degrees_north") )
  
  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "lev", NF90_DOUBLE, z_dimid, varid(vc)), "nf90_def_var lev" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "vertical level") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "units", "layer") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "positive", "down") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "coordinate", "eta") )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "standard_name", "model_layers") )

  vc=vc+1;
  call nccheck( nf90_def_var(self%ncid, "time", NF90_INT, t_dimid, varid(vc)), "nf90_def_var time" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name", "time"), "nf90_def_var time long_name" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "begin_date", date8), "nf90_def_var time begin_date" )
  call nccheck( nf90_put_att(self%ncid, varid(vc), "begin_time", time6), "nf90_def_var time begin_time" )
  
  ! Define fields to be written
  do n = 1,size(fields)
  
    if (fields(n)%npz == 1) then
      allocate(dimids(ndimids2))
      dimids = dimids2
    elseif (fields(n)%npz == geom%npz) then
      allocate(dimids(ndimids3))
      dimids = dimids3
    else
      call abor1_ftn("read_geos: vertical dimension not supported")
    endif
  
    vc=vc+1
    call nccheck( nf90_def_var(self%ncid, trim(fields(n)%short_name), NF90_DOUBLE, dimids, varid(vc)), &
                  "nf90_def_var"//trim(fields(n)%short_name)   )
    call nccheck( nf90_put_att(self%ncid, varid(vc), "long_name"    , trim(fields(n)%long_name) ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid, varid(vc), "units"        , trim(fields(n)%units)     ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid, varid(vc), "standard_name", trim(fields(n)%long_name) ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid, varid(vc), "coordinates"  , "lons lats"               ), "nf90_put_att" )
    call nccheck( nf90_put_att(self%ncid, varid(vc), "grid_mapping" , "cubed_sphere"            ), "nf90_put_att" )
  
    deallocate(dimids)
  
  enddo
  
  ! End define mode
  call nccheck( nf90_enddef(self%ncid), "nf90_enddef" )
   
endif


! Write Xdim, YDim, and nf (tile) 
! -------------------------------
if (self%iam_io_proc) then

  vc=0

  allocate(intarray(6))
  do n = 1,6
    intarray(n) = n
  enddo
  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), intarray ), "nf90_put_var nf" )  
  deallocate(intarray)

  allocate(realarray(geom%npx-1))
  do n = 1,geom%npx-1
    realarray(n) = real(n,kind_real)
  enddo

  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), realarray ), "nf90_put_var Xdim" )  
  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), realarray ), "nf90_put_var Ydim" )  

  deallocate(realarray)

endif


! Gather longitudes to write
! --------------------------
if (self%csize > 6) then
  call gather_tile(geom, self%tcomm, rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec), arrayg)
else
  arrayg = rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec)
endif

if (self%iam_io_proc) then
  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), arrayg, &
                                      start = self%istart2(1:3), count = self%icount2(1:3) ), "nf90_put_var lons" )
endif


! Gather latitudes and write
! --------------------------
if (self%csize > 6) then
  call gather_tile(geom, self%tcomm, rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec), arrayg)
else
  arrayg = rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec)
endif

if (self%iam_io_proc) then
  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), arrayg, &
                                      start = self%istart2(1:3), count = self%icount2(1:3) ), "nf90_put_var lats" )
endif


! Write model levels
! ------------------
if (self%iam_io_proc) then
  allocate(intarray(geom%npz))
  do n = 1,geom%npz
    intarray(n) = n
  enddo
  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), intarray ), "nf90_put_var lev" )  
  deallocate(intarray)

  vc=vc+1;call nccheck( nf90_put_var( self%ncid, varid(vc), 0 ), "nf90_put_var time" )  
endif


! Loop over fields and levels and write fields to file
! ----------------------------------------------------
do n = 1,size(fields)
  
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
    vc = vc + 1

  endif 
 
  do lev = 1,fields(n)%npz

    if (self%csize > 6) then
      call gather_tile(geom, self%tcomm, fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev), arrayg)
    else
      arrayg = fields(n)%array(:,:,lev)
    endif

    if (self%iam_io_proc) then

      istart3(self%vindex) = lev

      call nccheck( nf90_put_var( self%ncid, varid(vc), arrayg, start = istart, count = icount ), &
                                  "nf90_put_var "//trim(fields(n)%short_name) )
 
    endif

  enddo

  if (self%iam_io_proc .and. n==size(fields)) deallocate(istart3)

enddo

if (self%iam_io_proc) &
deallocate ( dimidsv, dimidsg, dimids2, dimids3 )

deallocate(arrayg)

end subroutine write_all

! ------------------------------------------------------------------------------

subroutine gather_tile(geom, comm, array_l, array_g)

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
real(kind=kind_real), intent(in)    :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec)  ! Local array
real(kind=kind_real), intent(inout) :: array_g(1:geom%npx-1,1:geom%npy-1)            ! Gathered array (only valid on root)

integer :: comm_size
integer :: n, ierr, npx_l, npy_l, subarray, resized_subarray
integer :: sizes_g(2), sizes_l(2), start_l(2), arraydispls_me
integer, allocatable :: counts(:), displs(:), arraydispls(:)
integer(kind=MPI_ADDRESS_KIND) :: extent, lb
real(kind=kind_real) :: forsize

call mpi_comm_size(comm, comm_size, ierr)

npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

sizes_g = [geom%npx-1, geom%npy-1]
sizes_l = [npx_l, npy_l]
start_l = [geom%isc-1, geom%jsc-1]

! Create recieving array
call mpi_type_create_subarray(2, sizes_g, sizes_l, start_l, mpi_order_fortran, mpi_double_precision, subarray, ierr)
call mpi_type_commit(subarray, ierr)

! Perform resizing
extent = sizeof(forsize)
lb = 0
call mpi_type_create_resized(subarray, lb, extent, resized_subarray, ierr)
call mpi_type_commit(resized_subarray,ierr)

! Set counts and displacement and gather
allocate(counts(comm_size), arraydispls(comm_size), displs(comm_size))

do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

arraydispls = 0
arraydispls_me = (geom%isc - 1) + (geom%jsc - 1) * (geom%npy-1)
call mpi_allgatherv(arraydispls_me, 1, mpi_int, arraydispls, counts, displs, mpi_int, comm, ierr)

! Gather the full field
call mpi_gatherv( array_l, npx_l*npy_l, mpi_double_precision, &
                  array_g, counts, arraydispls, resized_subarray, &
                  0, comm, ierr)

! Deallocate
deallocate(counts,displs,arraydispls)

end subroutine gather_tile

! ------------------------------------------------------------------------------

subroutine scatter_tile(geom, comm, array_g, array_l)

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
real(kind=kind_real), intent(in)    :: array_g(1:geom%npx-1,1:geom%npy-1)            ! Gathered array (only valid on root)
real(kind=kind_real), intent(inout) :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec)  ! Local array

integer :: comm_size
integer :: n, ierr, npx_l, npy_l, subarray, resized_subarray
integer :: sizes_g(2), sizes_l(2), start_l(2), arraydispls_me
integer, allocatable :: counts(:), displs(:), arraydispls(:)
integer(kind=MPI_ADDRESS_KIND) :: extent, lb
real(kind=kind_real) :: forsize

call mpi_comm_size(comm, comm_size, ierr)

npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

sizes_g = [geom%npx-1, geom%npy-1]
sizes_l = [npx_l, npy_l]
start_l = [geom%isc-1, geom%jsc-1]

! Create recieving array
call mpi_type_create_subarray(2, sizes_g, sizes_l, start_l, mpi_order_fortran, mpi_double_precision, subarray, ierr)
call mpi_type_commit(subarray, ierr)

! Perform resizing
extent = sizeof(forsize)
lb = 0
call mpi_type_create_resized(subarray, lb, extent, resized_subarray, ierr)
call mpi_type_commit(resized_subarray,ierr)

! Set counts and displacement and gather
allocate(counts(comm_size), arraydispls(comm_size), displs(comm_size))

do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

arraydispls = 0
arraydispls_me = (geom%isc - 1) + (geom%jsc - 1) * (geom%npy-1)
call mpi_allgatherv(arraydispls_me, 1, mpi_int, arraydispls, counts, displs, mpi_int, comm, ierr)

! Scatter the full field
call mpi_scatterv( array_g, counts, arraydispls, resized_subarray, &
                   array_l, npx_l*npy_l, mpi_double_precision, &
                   0, comm, ierr )

! Deallocate
deallocate(counts,displs,arraydispls)

end subroutine scatter_tile

! ------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_io_geos), intent(inout) :: self
end subroutine dummy_final

! ------------------------------------------------------------------------------

end module fv3jedi_io_geos_mod
