! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_fv3_mod

use netcdf

use femps_kinds_mod, only: kind_real
use femps_grid_mod, only: fempsgrid
use femps_utils_mod
use femps_const_mod

implicit none

private
public fv3grid_to_ugrid, fv3field_to_ufield, ufield_to_fv3field, write_field_fv3format, read_field_fv3format

interface fv3grid_to_ugrid
   module procedure fv3grid_to_ugrid_read
   module procedure fv3grid_to_ugrid_
end interface

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_to_ugrid_read(grid,path)

implicit none

type(fempsgrid),  intent(inout) :: grid
character(len=*), intent(in)    :: path

integer :: iv, fxdim, vxdim
character(len=2055) :: filename

real(kind=kind_real), allocatable, dimension(:,:,:) :: flong
real(kind=kind_real), allocatable, dimension(:,:,:) :: flat
real(kind=kind_real), allocatable, dimension(:,:,:) :: vlong
real(kind=kind_real), allocatable, dimension(:,:,:) :: vlat

do iv = 1,grid%ngrids

  write(filename,"(A,A10,I0.4,A4)") trim(path),'/fv3grid_c',grid%ncube(iv),'.nc4'

  call message('Reading fv3 grid from '//trim(filename), trace)

  call fv3grid_read(filename,fxdim,vxdim,flong,flat,vlong,vlat)

  call fv3grid_to_ugrid_(grid,iv,fxdim,vxdim,flat,flong,vlat,vlong)

  deallocate(flong,flat,vlong,vlat)

enddo

end subroutine fv3grid_to_ugrid_read

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_read(filename,fxdim,vxdim,flong,flat,vlong,vlat)

implicit none

character(len=*),                                    intent(in)    :: filename
integer,                                             intent(out)   :: fxdim
integer,                                             intent(out)   :: vxdim
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: flong
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: flat
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: vlong
real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout) :: vlat

integer :: ncid, dimid, varid

! Read the fv3 grid from file
! ---------------------------

call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, ncid), "nf90_open"//trim(filename) )

! Get size of cube
call nccheck ( nf90_inq_dimid(ncid, "fxdim",  dimid), "nf90_inq_dimid fxdim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = fxdim), "nf90_inquire_dimension fxdim" )

call nccheck ( nf90_inq_dimid(ncid, "vxdim",  dimid), "nf90_inq_dimid vxdim" )
call nccheck ( nf90_inquire_dimension(ncid, dimid, len = vxdim), "nf90_inquire_dimension vxdim" )

! Allocate arrays
allocate(flong(fxdim,fxdim,6))
allocate(flat (fxdim,fxdim,6))
allocate(vlong(vxdim,vxdim,6))
allocate(vlat (vxdim,vxdim,6))

! Get flong and flat
call nccheck ( nf90_inq_varid (ncid, 'flons', varid), "nf90_inq_varid flons" )
call nccheck ( nf90_get_var   (ncid, varid, flong), "nf90_get_var flons" )
call nccheck ( nf90_inq_varid (ncid, 'flats', varid), "nf90_inq_varid flats" )
call nccheck ( nf90_get_var   (ncid, varid, flat), "nf90_get_var flats" )

! Get vlong and vlat
call nccheck ( nf90_inq_varid (ncid, 'vlons', varid), "nf90_inq_varid vlons" )
call nccheck ( nf90_get_var   (ncid, varid, vlong), "nf90_get_var vlons" )
call nccheck ( nf90_inq_varid (ncid, 'vlats', varid), "nf90_inq_varid vlats" )
call nccheck ( nf90_get_var   (ncid, varid, vlat), "nf90_get_var vlats" )

call nccheck ( nf90_close   (ncid), "nf90_close" )

end subroutine fv3grid_read

! --------------------------------------------------------------------------------------------------

subroutine fv3grid_to_ugrid_(grid,iv,fxdim,vxdim,flat,flong,vlat,vlong)

implicit none

type(fempsgrid),      intent(inout) :: grid
integer,              intent(in)    :: iv
integer,              intent(in)    :: fxdim
integer,              intent(in)    :: vxdim
real(kind=kind_real), intent(in)    :: flong(fxdim,fxdim,6)
real(kind=kind_real), intent(in)    :: flat (fxdim,fxdim,6)
real(kind=kind_real), intent(in)    :: vlong(vxdim,vxdim,6)
real(kind=kind_real), intent(in)    :: vlat (vxdim,vxdim,6)

integer :: ji,jj,jt,n

! Subroutine for converting the fv3 grid into the
! unstructured grid used by femps

! Grid faces
n = 0
do jt = 1,6
  do jj = 1,fxdim
    do ji = 1,fxdim
      n = n+1
      grid%flong(n,iv) = flong(ji,jj,jt)
      grid%flat(n,iv)  = flat (ji,jj,jt)
    enddo
  enddo
enddo

! Grid vertices
n = 0
do jt = 1,6
  do jj = 1,fxdim
    do ji = 1,fxdim
      n = n+1
      grid%vlong(n,iv) = vlong(ji,jj,jt)
      grid%vlat(n,iv)  = vlat (ji,jj,jt)
    enddo
  enddo
enddo

! Two vertices are neglected in above sweep
grid%vlong(grid%nvert(iv)-1,iv) = vlong(1,fxdim+1,1)
grid%vlat (grid%nvert(iv)-1,iv) = vlat (1,fxdim+1,1)
grid%vlong(grid%nvert(iv)  ,iv) = vlong(fxdim+1,1,2)
grid%vlat (grid%nvert(iv)  ,iv) = vlat (fxdim+1,1,2)

end subroutine fv3grid_to_ugrid_

! --------------------------------------------------------------------------------------------------

subroutine fv3field_to_ufield(grid,npx,field_fv3,field_femps)

implicit none

type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: npx
real(kind=kind_real), intent(in)  :: field_fv3  (npx,npx,6)
real(kind=kind_real), intent(out) :: field_femps(grid%nface(grid%ngrids))

integer :: ji,jj,jt,n

! Pack fv3 field into unstructured field (Always A-Grid)
n = 0
do jt = 1,6
  do jj = 1,npx
    do ji = 1,npx
      n = n+1
      field_femps(n) = field_fv3(ji,jj,jt)
    enddo
  enddo
enddo

end subroutine fv3field_to_ufield

! --------------------------------------------------------------------------------------------------

subroutine ufield_to_fv3field(grid,npx,field_femps,field_fv3)

implicit none

type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: npx
real(kind=kind_real), intent(in)  :: field_femps(grid%nface(grid%ngrids))
real(kind=kind_real), intent(out) :: field_fv3  (npx,npx,6)

integer :: ji,jj,jt,n

! Pack fv3 field into unstructured field (Always A-Grid)
n = 0
do jt = 1,6
  do jj = 1,npx
    do ji = 1,npx
      n = n+1
      field_fv3(ji,jj,jt) = field_femps(n)
    enddo
  enddo
enddo

end subroutine ufield_to_fv3field

! --------------------------------------------------------------------------------------------------

subroutine write_field_fv3format(grid,igrid,filename,fieldname,field)

implicit none
type(fempsgrid),      intent(in) :: grid
integer,              intent(in) :: igrid
character(len=*),     intent(in) :: filename
character(len=*),     intent(in) :: fieldname
real(kind=kind_real), intent(in) :: field(grid%nface(igrid))

real(kind=kind_real), allocatable :: lon_fv3(:,:,:)
real(kind=kind_real), allocatable :: lat_fv3(:,:,:)
real(kind=kind_real), allocatable :: field_fv3(:,:,:)

integer :: ncid, vc, varid(1000), n, tiles(6)
integer :: xdim_dimid, ydim_dimid, tile_dimid
real(kind=kind_real), allocatable :: xdims(:), ydims(:)

allocate(lon_fv3(grid%ncube(igrid),grid%ncube(igrid),6))
allocate(lat_fv3(grid%ncube(igrid),grid%ncube(igrid),6))
allocate(field_fv3(grid%ncube(igrid),grid%ncube(igrid),6))

call ufield_to_fv3field(grid,grid%ncube(igrid),grid%flong(:,igrid),lon_fv3)
call ufield_to_fv3field(grid,grid%ncube(igrid),grid%flat (:,igrid),lat_fv3)
call ufield_to_fv3field(grid,grid%ncube(igrid),field,field_fv3)

allocate(xdims(grid%ncube(igrid)))
allocate(ydims(grid%ncube(igrid)))

do n = 1, grid%ncube(igrid)
  xdims(n) = real(n,kind_real)
  ydims(n) = real(n,kind_real)
enddo

do n = 1,6
  tiles(n) = n
enddo

! Create file
! -----------
call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), ncid), "nf90_create" )

! Define dimensions
! -----------------
call nccheck( nf90_def_dim(ncid, "Xdim", grid%ncube(igrid), xdim_dimid), "nf90_def_dim Xdim" )
call nccheck( nf90_def_dim(ncid, "Ydim", grid%ncube(igrid), ydim_dimid), "nf90_def_dim Ydim" )
call nccheck( nf90_def_dim(ncid, "nf",   6,                 tile_dimid), "nf90_def_dim nf"   )


! Define variables
! ----------------
vc = 1

call nccheck( nf90_def_var(ncid, "Xdim",  NF90_DOUBLE, (/ xdim_dimid /), varid(vc)), "nf90_def_var Xdim"  )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "Fake Longitude for GrADS Compatibility") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_east") )

vc = vc + 1
call nccheck( nf90_def_var(ncid, "Ydim",  NF90_DOUBLE, (/ ydim_dimid /), varid(vc)), "nf90_def_var Ydim"  )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "Fake Longitude for GrADS Compatibility") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_north") )

vc = vc + 1
call nccheck( nf90_def_var(ncid, "lons",  NF90_DOUBLE, (/ xdim_dimid, ydim_dimid, tile_dimid /), varid(vc)), "nf90_def_var lon"  )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "longitude") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_east") )

vc = vc + 1
call nccheck( nf90_def_var(ncid, "lats",  NF90_DOUBLE, (/ xdim_dimid, ydim_dimid, tile_dimid /), varid(vc)), "nf90_def_var lat"  )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "latitude") )
call nccheck( nf90_put_att(ncid, varid(vc), "units", "degrees_north") )

vc = vc + 1
call nccheck( nf90_def_var(ncid, "nf", NF90_INT, (/ tile_dimid /), varid(vc)), "nf90_def_var nf" )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name", "cubed-sphere face") )
call nccheck( nf90_put_att(ncid, varid(vc), "axis", "e") )
call nccheck( nf90_put_att(ncid, varid(vc), "grads_dim", "e") )

vc = vc + 1
call nccheck( nf90_def_var(ncid, fieldname, NF90_DOUBLE, &
    (/ xdim_dimid, ydim_dimid, tile_dimid /), varid(vc)), "nf90_def_var field" )
call nccheck( nf90_put_att(ncid, varid(vc), "long_name"    , fieldname      ), "nf90_put_att" )
call nccheck( nf90_put_att(ncid, varid(vc), "units"        , '1'            ), "nf90_put_att" )
call nccheck( nf90_put_att(ncid, varid(vc), "standard_name", fieldname      ), "nf90_put_att" )
call nccheck( nf90_put_att(ncid, varid(vc), "coordinates"  , "lons lats"    ), "nf90_put_att" )
call nccheck( nf90_put_att(ncid, varid(vc), "grid_mapping" , "cubed_sphere" ), "nf90_put_att" )

! End define
! ----------
call nccheck( nf90_enddef(ncid), "nf90_enddef" )


! Write variables
! ---------------
vc = 1
call nccheck( nf90_put_var( ncid, varid(vc), xdims           ), "nf90_put_var xdims"   ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), ydims           ), "nf90_put_var ydims"   ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), rad2deg*lon_fv3 ), "nf90_put_var lon_fv3" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), rad2deg*lat_fv3 ), "nf90_put_var lat_fv3" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), tiles           ), "nf90_put_var nf"      ); vc = vc + 1


call nccheck( nf90_put_var( ncid, varid(vc), field_fv3 ), "nf90_put_var field"   ); vc = vc + 1

! Close file
! ----------
call nccheck ( nf90_close(ncid), "nf90_close" )


deallocate(lon_fv3,lat_fv3,field_fv3,xdims,ydims)

end subroutine write_field_fv3format

! --------------------------------------------------------------------------------------------------

subroutine read_field_fv3format(grid,igrid,filename,fieldname,field)

implicit none
type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: igrid
character(len=*),     intent(in)  :: filename
character(len=*),     intent(in)  :: fieldname
real(kind=kind_real), intent(out) :: field(grid%nface(igrid))

integer :: ncid, varid
real(kind=kind_real), allocatable :: field_fv3(:,:,:)


allocate(field_fv3(grid%ncube(igrid),grid%ncube(igrid),6))

call nccheck ( nf90_open(trim(filename), NF90_NOWRITE, ncid), "nf90_open"//trim(filename) )

call nccheck ( nf90_inq_varid (ncid, trim(fieldname), varid), "nf90_inq_varid "//trim(fieldname) )

call nccheck ( nf90_get_var( ncid, varid, field_fv3), "nf90_get_var "//trim(fieldname) )

call nccheck ( nf90_close(ncid), "nf90_close" )

call fv3field_to_ufield(grid,grid%ncube(igrid),field_fv3,field)

deallocate(field_fv3)

end subroutine read_field_fv3format

! --------------------------------------------------------------------------------------------------

end module femps_fv3_mod
