! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_io_geos_mod

use config_mod
use datetime_mod
use iso_c_binding

use fckit_mpi_module

use fv3jedi_geom_mod,         only: fv3jedi_geom
use fv3jedi_field_mod,        only: fv3jedi_field
use fv3jedi_kinds_mod,        only: kind_real
use fv3jedi_netcdf_utils_mod, only: nccheck

use mpi
use netcdf

implicit none
private
public read_geos, write_geos

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine read_geos(geom, fields, c_conf, vdate)

implicit none

!Arguments
type(fv3jedi_geom),  intent(in)    :: geom       !< Geom
type(fv3jedi_field), intent(inout) :: fields(:)  !< Fields
type(c_ptr),         intent(in)    :: c_conf     !< Configuration
type(datetime),      intent(inout) :: vdate      !< DateTime

logical :: iread !Am I an IO processor?
character(len=255) :: filename
integer :: ncid, dimid, varid
integer :: im, jm, lm, nm
integer :: date(6), intdate, inttime
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: idate, isecs
integer :: ncdim3, ncdim2
integer, target, allocatable :: istart3(:), icount3(:)
integer, target, allocatable :: istart2(:), icount2(:)
integer, pointer :: istart(:), icount(:)
integer :: n, geostiledim, tileoff

! Six processors with corner of tile will handle IO
! -------------------------------------------------
iread = .false.
if (geom%isc == 1 .and. geom%jsc == 1) then
  iread = .true.
endif

!> Set filenames
!> -------------
filename = 'Data/GEOS.bkg.eta.nc4'

if (config_element_exists(c_conf,"filename")) then
   filename = config_get_string(c_conf,len(filename),"filename")
endif

!> Open the file
call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, ncid), "nf90_open"//trim(filename) )

!> Get dimensions, XDim,YDim,lev,time
call nccheck ( nf90_inq_dimid(ncid, "Xdim", dimid), "nf90_inq_dimid Xdim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = im), "nf90_inquire_dimension Xdim" )

call nccheck ( nf90_inq_dimid(ncid, "Ydim", dimid), "nf90_inq_dimid YDim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = jm), "nf90_inquire_dimension YDim" )

call nccheck ( nf90_inq_dimid(ncid, "lev",  dimid), "nf90_inq_dimid lev" )
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

!> Set the object date from the date of the file
call datetime_from_ifs(vdate, idate, isecs)

!> Make sure file dimensions equal to geometry
if ( im /= geom%npx-1 .or. lm /= geom%npz) then
  call abor1_ftn("GEOS restarts: restart dimension not compatible with geometry")
endif

!> GEOS can use concatenated tiles or tile as a dimension
if ( (im == geom%npx-1) .and. (jm == 6*(geom%npy-1) ) ) then
  geostiledim = 0
  tileoff = (geom%ntile-1)*(jm/geom%ntiles)
  ncdim3 = 4
  ncdim2 = 3
else
  geostiledim = 1
  tileoff = 0
  ncdim3 = 5
  ncdim2 = 4
endif

allocate(istart3(ncdim3),icount3(ncdim3))
allocate(istart2(ncdim2),icount2(ncdim2))

! Create local to this proc start/count
if (geostiledim == 1) then
  istart3(1) = geom%isc;          icount3(1) = geom%iec-geom%isc+1
  istart3(2) = geom%jsc;          icount3(2) = geom%jec-geom%jsc+1
  istart3(3) = geom%ntile;        icount3(3) = 1
  istart3(4) = 1;                 icount3(4) = geom%npz
  istart3(5) = 1;                 icount3(5) = 1
  
  istart2(1) = geom%isc;          icount2(1) = geom%iec-geom%isc+1
  istart2(2) = geom%jsc;          icount2(2) = geom%jec-geom%jsc+1
  istart2(3) = geom%ntile;        icount2(3) = 1
  istart2(4) = 1;                 icount2(4) = 1
else
  istart3(1) = geom%isc;          icount3(1) = geom%iec-geom%isc+1
  istart3(2) = tileoff+geom%jsc;  icount3(2) = geom%jec-geom%jsc+1
  istart3(3) = 1;                 icount3(3) = geom%npz
  istart3(4) = 1;                 icount3(4) = 1
  
  istart2(1) = geom%isc;          icount2(1) = geom%iec-geom%isc+1
  istart2(2) = tileoff+geom%jsc;  icount2(2) = geom%jec-geom%jsc+1
  istart2(3) = 1;                 icount2(3) = 1
endif

do n = 1,size(fields)

  if (fields(n)%npz == 1) then
    istart => istart2; icount => icount2
  elseif (fields(n)%npz == geom%npz) then
    istart => istart3; icount => icount3
  else
    call abor1_ftn("read_geos: vertical dimension not supported")
  endif

  call nccheck ( nf90_inq_varid (ncid, trim(fields(n)%short_name), varid), "nf90_inq_varid"//trim(fields(n)%short_name) )
  call nccheck ( nf90_get_var( ncid, varid, &
                               fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1:fields(n)%npz), &
                               istart, icount), "nf90_get_var "//trim(fields(n)%short_name) )

  nullify(istart,icount)

enddo

!Close this file
call nccheck ( nf90_close(ncid), "nf90_close" )

! Deallocate
deallocate ( istart2, icount2 )
deallocate ( istart3, icount3 )

end subroutine read_geos

! ------------------------------------------------------------------------------

subroutine write_geos(geom, fields, c_conf, vdate)

implicit none

! Arguments
type(fv3jedi_geom),  intent(in)    :: geom       !< Geom
type(fv3jedi_field), intent(in)    :: fields(:)  !< Fields
type(c_ptr),         intent(in)    :: c_conf     !< Configuration
type(datetime),      intent(in)    :: vdate      !< DateTime

! Locals
character(len=255) :: datapath
character(len=255) :: filename
character(len=64)  :: datefile
integer :: geostiledim, tileoff
integer :: ncid, varid(1000), vc
integer :: date(6)
integer(kind=c_int) :: idate, isecs
integer :: im,jm,km,n
integer :: x_dimid, y_dimid, z_dimid
integer :: t_dimid, tile_dimid
integer :: ncdim3, ncdim2
integer, target, allocatable :: dimidsv(:), dimidsg(:), dimids2(:), dimids3(:)
integer, target, allocatable :: istart2(:), icount2(:)
integer, target, allocatable :: istart3(:), icount3(:)
integer, pointer :: dimids(:), istart(:), icount(:)
type(fckit_mpi_comm) :: f_comm
integer :: k, levs(geom%npz)

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
geostiledim = 1
if (config_element_exists(c_conf,"geos_tile_dim")) then
   geostiledim = config_get_int(c_conf,"geos_tile_dim")
endif

! File total dims
im = geom%npx-1
jm = geom%npy-1
km = geom%npz
if (geostiledim == 0) then
  jm = 6*jm
  tileoff = (geom%ntile-1)*(jm/geom%ntiles)
  ncdim3 = 4
  ncdim2 = 3
else
  ncdim3 = 5
  ncdim2 = 4
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
vc = 0

vc=vc+1;
call nccheck( nf90_def_var(ncid, "lons", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lons" )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "longitude") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_east") )

vc=vc+1;
call nccheck( nf90_def_var(ncid, "lats", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lats" )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "latitude") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_north") )

! Define fields to be written
do n = 1,size(fields)

  if (fields(n)%npz == 1) then
    dimids => dimids2
  elseif (fields(n)%npz == geom%npz) then
    dimids => dimids3
  else
    call abor1_ftn("read_geos: vertical dimension not supported")
  endif

  vc=vc+1
  call nccheck( nf90_def_var(ncid, trim(fields(n)%short_name), NF90_DOUBLE, dimids, varid(vc)), &
                "nf90_def_var"//trim(fields(n)%short_name)   )
  call nccheck( nf90_put_att(ncid, varid(vc), "long_name"   , trim(fields(n)%long_name) ), "nf90_put_att" )
  call nccheck( nf90_put_att(ncid, varid(vc), "units"       , trim(fields(n)%units)     ), "nf90_put_att" )
  call nccheck( nf90_put_att(ncid, varid(vc), "coordinates" , "lons lats"               ), "nf90_put_att" )
  call nccheck( nf90_put_att(ncid, varid(vc), "grid_mapping", "cubed_sphere"            ), "nf90_put_att" )

  nullify(dimids)

enddo

! End define mode
call nccheck( nf90_enddef(ncid), "nf90_enddef" )

allocate(istart3(ncdim3),icount3(ncdim3))
allocate(istart2(ncdim2),icount2(ncdim2))

! Create local to this proc start/count
if (geostiledim == 1) then
  istart3(1) = geom%isc;          icount3(1) = geom%iec-geom%isc+1
  istart3(2) = geom%jsc;          icount3(2) = geom%jec-geom%jsc+1
  istart3(3) = geom%ntile;        icount3(3) = 1
  istart3(4) = 1;                 icount3(4) = geom%npz
  istart3(5) = 1;                 icount3(5) = 1
  
  istart2(1) = geom%isc;          icount2(1) = geom%iec-geom%isc+1
  istart2(2) = geom%jsc;          icount2(2) = geom%jec-geom%jsc+1
  istart2(3) = geom%ntile;        icount2(3) = 1
  istart2(4) = 1;                 icount2(4) = 1
else
  istart3(1) = geom%isc;          icount3(1) = geom%iec-geom%isc+1
  istart3(2) = tileoff+geom%jsc;  icount3(2) = geom%jec-geom%jsc+1
  istart3(3) = 1;                 icount3(3) = geom%npz
  istart3(4) = 1;                 icount3(4) = 1

  istart2(1) = geom%isc;          icount2(1) = geom%iec-geom%isc+1
  istart2(2) = tileoff+geom%jsc;  icount2(2) = geom%jec-geom%jsc+1
  istart2(3) = 1;                 icount2(3) = 1
endif

vc=0

! Write fields (geom)
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec), &
                                    start = istart2(1:3), count = icount2(1:3) ), "nf90_put_var lons" )
vc=vc+1;call nccheck( nf90_put_var( ncid, varid(vc), geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec), &
                                    start = istart2(1:3), count = icount2(1:3) ), "nf90_put_var lats" )

! Write fields 
do n = 1,size(fields)

  if (fields(n)%npz == 1) then
    istart => istart2; icount => icount2
  elseif (fields(n)%npz == geom%npz) then
    istart => istart3; icount => icount3
  else
    call abor1_ftn("read_geos: vertical dimension not supported")
  endif

  vc=vc+1;
  call nccheck( nf90_put_var( ncid, varid(vc), &
                              fields(n)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1:fields(n)%npz), &
                              start = istart, count = icount ), "nf90_put_var "//trim(fields(n)%short_name) )

  nullify(istart,icount)

enddo

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

! Deallocate
deallocate ( dimidsv, dimidsg, dimids2, dimids3 )
deallocate ( istart2, icount2 )
deallocate ( istart3, icount3 )

end subroutine write_geos

! ------------------------------------------------------------------------------

end module fv3jedi_io_geos_mod
