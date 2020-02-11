module fv3geom_gnomonicgrids_mod

use fv3geom_constants_mod
use fv3geom_utils_mod

implicit none

private
public gnomonic_grids

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine gnomonic_grids(grid_type, im, lon, lat)

integer,              intent(in)  :: im, grid_type
real(kind=kind_real), intent(out) :: lon(im+1,im+1)
real(kind=kind_real), intent(out) :: lat(im+1,im+1)

! Locals
integer :: i, j

if(grid_type==0) call gnomonic_ed  (im, lon, lat)
if(grid_type==1) call gnomonic_dist(im, lon, lat)
if(grid_type==2) call gnomonic_angl(im, lon, lat)

call symm_ed(im, lon, lat)
do j = 1, im+1
  do i = 1, im+1
    lon(i,j) = lon(i,j) - pi
  enddo
enddo

end subroutine gnomonic_grids

! --------------------------------------------------------------------------------------------------

subroutine gnomonic_ed(im, lamda, theta)

!-----------------------------------------------------
! Equal distance along the 4 edges of the cubed sphere
!-----------------------------------------------------
! Properties:
!            * defined by intersections of great circles
!            * max(dx,dy; global) / min(dx,dy; global) = sqrt(2) = 1.4142
!            * Max(aspect ratio) = 1.06089
!            * the N-S coordinate curves are const longitude on the 4 faces with equator
! For C2000: (dx_min, dx_max) = (3.921, 5.545)    in km unit
! This is the grid of choice for global cloud resolving

integer,              intent(in)  :: im
real(kind=kind_real), intent(out) :: lamda(im+1,im+1)
real(kind=kind_real), intent(out) :: theta(im+1,im+1)

! Local:
real(kind=kind_real) :: pp(3,im+1,im+1)
real(kind=kind_real) :: rsq3, alpha, dely
integer :: i, j, k

rsq3  = 1.0_8/sqrt(3.0_8)
alpha = asin( rsq3 )

! Ranges:
! lamda = [0.75*pi, 1.25*pi]
! theta = [-alpha, alpha]

dely = 2.0_8*alpha / real(im, 8)

! Define East-West edges:
do j=1,im+1
  lamda(1,   j) = 0.75_8*pi                  ! West edge
  lamda(im+1,j) = 1.25_8*pi                  ! East edge
  theta(1,   j) = -alpha + dely*real(j-1,8)  ! West edge
  theta(im+1,j) = theta(1,j)                 ! East edge
enddo

! Get North-South edges by symmetry:

do i=2,im
  call mirror_latlon(lamda(1,1), theta(1,1), lamda(im+1,im+1), theta(im+1,im+1), &
                     lamda(1,i), theta(1,i), lamda(i,1),       theta(i,      1) )
   lamda(i,im+1) =  lamda(i,1)
   theta(i,im+1) = -theta(i,1)
enddo

! Set 4 corners:
call latlon2xyz2(lamda(1    ,  1), theta(1,      1), pp(1,   1,   1))
call latlon2xyz2(lamda(im+1,   1), theta(im+1,   1), pp(1,im+1,   1))
call latlon2xyz2(lamda(1,   im+1), theta(1,   im+1), pp(1,   1,im+1))
call latlon2xyz2(lamda(im+1,im+1), theta(im+1,im+1), pp(1,im+1,im+1))

! Map edges on the sphere back to cube:
! Intersections at x=-rsq3

i = 1
do j=2,im
  call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,j))
  pp(2,i,j) = -pp(2,i,j)*rsq3/pp(1,i,j)
  pp(3,i,j) = -pp(3,i,j)*rsq3/pp(1,i,j)
enddo

j=1
do i=2,im
  call latlon2xyz2(lamda(i,j), theta(i,j), pp(1,i,1))
  pp(2,i,1) = -pp(2,i,1)*rsq3/pp(1,i,1)
  pp(3,i,1) = -pp(3,i,1)*rsq3/pp(1,i,1)
enddo

do j=1,im+1
  do i=1,im+1
    pp(1,i,j) = -rsq3
  enddo
enddo

do j=2,im+1
  do i=2,im+1
    !Copy y-z face of the cube along j=1
    pp(2,i,j) = pp(2,i,1)
    !Copy along i=1
    pp(3,i,j) = pp(3,1,j)
  enddo
enddo

call cart_to_latlon( (im+1)*(im+1), pp, lamda, theta)

end subroutine gnomonic_ed

! --------------------------------------------------------------------------------------------------

subroutine gnomonic_dist(im, lamda, theta)

! This is the commonly known equi-distance grid

integer,              intent(in)  :: im
real(kind=kind_real), intent(out) :: lamda(im+1,im+1)
real(kind=kind_real), intent(out) :: theta(im+1,im+1)

! Locals
real(kind=kind_real) p(3,im+1,im+1)
real(kind=kind_real) rsq3, xf, y0, z0, y, x, z, ds
real(kind=kind_real) dy, dz
integer j,k

! Face-2

rsq3 = 1.0_kind_real/sqrt(3.0_kind_real)
xf = -rsq3
y0 =  rsq3
dy = -2.0_kind_real*rsq3/im
z0 = -rsq3
dz =  2.0_kind_real*rsq3/im

do k=1,im+1
   do j=1,im+1
      p(1,j,k) = xf
      p(2,j,k) = y0 + (j-1)*dy
      p(3,j,k) = z0 + (k-1)*dz
   enddo
enddo
call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

end subroutine gnomonic_dist

! --------------------------------------------------------------------------------------------------

subroutine gnomonic_angl(im, lamda, theta)

! This is the commonly known equi-angular grid

integer,              intent(in)  :: im
real(kind=kind_real), intent(out) :: lamda(im+1,im+1)
real(kind=kind_real), intent(out) :: theta(im+1,im+1)

! Locals
integer :: j,k
real(kind=kind_real) :: p(3,im+1,im+1), rsq3, dp

dp = 0.5_kind_real*pi/real(im,kind=kind_real)

rsq3 = 1.0_kind_real/sqrt(3.0_kind_real)

do k=1,im+1
  do j=1,im+1
    p(1,j,k) =-rsq3               ! constant
    p(2,j,k) =-rsq3*tan(-0.25_kind_real*pi+(j-1)*dp)
    p(3,j,k) = rsq3*tan(-0.25_kind_real*pi+(k-1)*dp)
  enddo
enddo

call cart_to_latlon( (im+1)*(im+1), p, lamda, theta)

end subroutine gnomonic_angl

! --------------------------------------------------------------------------------------------------

end module fv3geom_gnomonicgrids_mod
