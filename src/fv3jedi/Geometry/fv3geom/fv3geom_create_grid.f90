module fv3geom_create_grid_mod

use fv3geom_constants_mod,     only: kind_real, pi
use fv3geom_gnomonicgrids_mod, only: gnomonic_grids
use fv3geom_utils_mod,         only: mirror_grid, direct_transform

implicit none

private
public create_fv3_grid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create_fv3_grid(grid_type, npx, npy, ng, isc, iec, jsc, jec, tile, lat, lon, &
                           shift_fac_in, do_schmidt_in, stretch_fac, target_lat, target_lon)

integer,                        intent(in)  :: grid_type
integer,                        intent(in)  :: npx
integer,                        intent(in)  :: npy
integer,                        intent(in)  :: ng
integer,                        intent(in)  :: isc, iec, jsc, jec, tile
real(kind=kind_real),           intent(out) :: lat(isc:iec+1,jsc:jec+1)
real(kind=kind_real),           intent(out) :: lon(isc:iec+1,jsc:jec+1)
real(kind=kind_real), optional, intent(in)  :: shift_fac_in
logical,              optional, intent(in)  :: do_schmidt_in
real(kind=kind_real), optional, intent(in)  :: stretch_fac
real(kind=kind_real), optional, intent(in)  :: target_lat
real(kind=kind_real), optional, intent(in)  :: target_lon

! Locals
integer :: i, j, n, nregions = 6
real(kind=kind_real), allocatable :: xs(:,:)
real(kind=kind_real), allocatable :: ys(:,:)
real(kind=kind_real), allocatable :: grid_global(:,:,:,:)
real(kind=kind_real) :: shift_fac
logical :: do_schmidt


! Create grid for a face of the cube centered on (0,0)
! ----------------------------------------------------
allocate(xs(npx,npy))
allocate(ys(npx,npy))

call gnomonic_grids(grid_type, npx-1, xs, ys)


! Set up remaining tiles
! ----------------------
allocate(grid_global(1-ng:npx+ng, 1-ng:npy+ng, 2, nregions))

do j = 1,npy
  do i = 1,npx
    grid_global(i,j,1,1) = xs(i,j)
    grid_global(i,j,2,1) = ys(i,j)
  enddo
enddo

call mirror_grid(grid_global, ng, npx, npy, 2, nregions)


! Adjustments
! -----------
shift_fac = 0.0_kind_real
do_schmidt = .false.
if (present(shift_fac_in)) shift_fac = shift_fac_in
if (present(do_schmidt_in)) do_schmidt = do_schmidt_in

do n=1,nregions

  ! Shift the grid, e.g. to move corners away from certain places
  if ( .not.do_schmidt .and. shift_fac > 1.0e-4_kind_real ) then
    grid_global(1:npx,1:npy,1,n) = grid_global(1:npx,1:npy,1,n) - pi/shift_fac
  endif

  do j=1,npy
    do i=1,npx

      if ( grid_global(i,j,1,n) < 0. ) grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
      if (abs(grid_global(i,j,1,1)) < 1.0e-10_kind_real) grid_global(i,j,1,1) = 0.0_kind_real
      if (abs(grid_global(i,j,2,1)) < 1.0e-10_kind_real) grid_global(i,j,2,1) = 0.0_kind_real

    enddo
  enddo
enddo

grid_global(  1,1:npy,:,2) = grid_global(npx,1:npy,:,1)
grid_global(  1,1:npy,:,3) = grid_global(npx:1:-1,npy,:,1)
grid_global(1:npx,npy,:,5) = grid_global(1,npy:1:-1,:,1)
grid_global(1:npx,npy,:,6) = grid_global(1:npx,1,:,1)

grid_global(1:npx,  1,:,3) = grid_global(1:npx,npy,:,2)
grid_global(1:npx,  1,:,4) = grid_global(npx,npy:1:-1,:,2)
grid_global(npx,1:npy,:,6) = grid_global(npx:1:-1,1,:,2)

grid_global(  1,1:npy,:,4) = grid_global(npx,1:npy,:,3)
grid_global(  1,1:npy,:,5) = grid_global(npx:1:-1,npy,:,3)

grid_global(npx,1:npy,:,3) = grid_global(1,1:npy,:,4)
grid_global(1:npx,  1,:,5) = grid_global(1:npx,npy,:,4)
grid_global(1:npx,  1,:,6) = grid_global(npx,npy:1:-1,:,4)

grid_global(  1,1:npy,:,6) = grid_global(npx,1:npy,:,5)


! Schmidt transform
! -----------------
if ( do_schmidt ) then
  do n = 1,nregions
    call direct_transform( stretch_fac, 1, npx, 1, npy, target_lon, target_lat, &
                           n, grid_global(1:npx,1:npy,1,n), grid_global(1:npx,1:npy,2,n) )
  enddo
endif


! Fill up outputs
! ---------------
do j = jsc,jec+1
   do i = isc,iec+1
      lon(i,j) = grid_global(i,j,1,tile)
      lat(i,j) = grid_global(i,j,2,tile)
   enddo
enddo


end subroutine create_fv3_grid

! --------------------------------------------------------------------------------------------------

end module fv3geom_create_grid_mod
