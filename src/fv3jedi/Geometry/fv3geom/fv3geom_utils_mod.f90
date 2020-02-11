module fv3geom_utils_mod

use fv3geom_constants_mod

implicit none

public

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine mirror_latlon(lon1, lat1, lon2, lat2, lon0, lat0, lon3, lat3)

! Given the "mirror" as defined by (lon1, lat1), (lon2, lat2), and center
! of the sphere, compute the mirror image of (lon0, lat0) as  (lon3, lat3)

real(kind=kind_real), intent(in)  :: lon1, lat1, lon2, lat2, lon0, lat0
real(kind=kind_real), intent(out) :: lon3, lat3

! Locals
real(kind=kind_real) :: p0(3), p1(3), p2(3), nb(3), pp(3), sp(2)
real(kind=kind_real) :: pdot
integer :: k

call latlon2xyz2(lon0, lat0, p0)
call latlon2xyz2(lon1, lat1, p1)
call latlon2xyz2(lon2, lat2, p2)
call vect_cross(nb, p1, p2)

pdot = sqrt(nb(1)**2+nb(2)**2+nb(3)**2)
do k=1,3
  nb(k) = nb(k) / pdot
enddo

pdot = p0(1)*nb(1) + p0(2)*nb(2) + p0(3)*nb(3)
do k=1,3
  pp(k) = p0(k) - 2.0_kind_real*pdot*nb(k)
enddo

call cart_to_latlon(1, pp, sp(1), sp(2))
lon3 = sp(1)
lat3 = sp(2)

end subroutine  mirror_latlon

! --------------------------------------------------------------------------------------------------

subroutine cart_to_latlon(np, q, xs, ys)

integer,              intent(in)    :: np
real(kind=kind_real), intent(inout) :: q(3,np)
real(kind=kind_real), intent(inout) :: xs(np), ys(np)

! Locals
real(kind=kind_real), parameter :: esl=1.0e-10_kind_real
real(kind=kind_real) :: p(3)
real(kind=kind_real) :: dist, lat, lon
integer i,k

do i=1,np

  do k=1,3
    p(k) = q(k,i)
  enddo

  dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)

  do k=1,3
    p(k) = p(k) / dist
  enddo

  if ( (abs(p(1))+abs(p(2)))  < esl ) then
    lon = 0.0_kind_real
  else
    lon = atan2( p(2), p(1) )   ! range [-pi,pi]
  endif

  if ( lon < 0.) lon = 2.0_kind_real*pi + lon

  ! Right hand system:
  lat = asin(p(3))

  xs(i) = lon
  ys(i) = lat

  ! q Normalized:
  do k=1,3
    q(k,i) = p(k)
  enddo

enddo

end subroutine cart_to_latlon

! --------------------------------------------------------------------------------------------------

subroutine vect_cross(e, p1, p2)

real(kind=kind_real), intent(in)  :: p1(3), p2(3)
real(kind=kind_real), intent(out) :: e(3)

! Perform cross products of 3D vectors: e = P1 X P2
e(1) = p1(2)*p2(3) - p1(3)*p2(2)
e(2) = p1(3)*p2(1) - p1(1)*p2(3)
e(3) = p1(1)*p2(2) - p1(2)*p2(1)

end subroutine vect_cross

! --------------------------------------------------------------------------------------------------

subroutine latlon2xyz2(lon, lat, p3)

real(kind=kind_real), intent(in)  :: lon, lat
real(kind=kind_real), intent(out) :: p3(3)

! Locals
real(kind=kind_real) :: e(2)

e(1) = lon
e(2) = lat
call latlon2xyz(e, p3)

end subroutine latlon2xyz2

! --------------------------------------------------------------------------------------------------

subroutine latlon2xyz(p, e, id)

! Routine to map (lon, lat) to (x,y,z)

real(kind=kind_real),  intent(in)  :: p(2)
real(kind=kind_real),  intent(out) :: e(3)
integer, optional,     intent(in)  :: id   ! id=0 do nothing; id=1, right_hand

! Locals
integer :: n
real(kind=kind_real) :: q(2)
real(kind=kind_real) :: e1, e2, e3

do n=1,2
  q(n) = p(n)
enddo

e1 = cos(q(2)) * cos(q(1))
e2 = cos(q(2)) * sin(q(1))
e3 = sin(q(2))

! Truncate to the desired precision:
e(1) = e1
e(2) = e2
e(3) = e3

end subroutine latlon2xyz

! --------------------------------------------------------------------------------------------------

subroutine symm_ed(im, lamda, theta)

! Make grid symmetrical to i=im/2+1

integer,              intent(in)  :: im
real(kind=kind_real), intent(out) :: lamda(im+1,im+1)
real(kind=kind_real), intent(out) :: theta(im+1,im+1)

! Locals
integer :: i,j,ip,jp
real(kind=kind_real) :: avg

do j=2,im+1
  do i=2,im
    lamda(i,j) = lamda(i,1)
  enddo
enddo

do j=1,im+1
  do i=1,im/2
    ip = im + 2 - i
    avg = 0.5_kind_real*(lamda(i,j)-lamda(ip,j))
    lamda(i, j) = avg + pi
    lamda(ip,j) = pi - avg
    avg = 0.5_kind_real*(theta(i,j)+theta(ip,j))
    theta(i, j) = avg
    theta(ip,j) = avg
  enddo
enddo

! Make grid symmetrical to j=im/2+1
do j=1,im/2
  jp = im + 2 - j
  do i=2,im
    avg = 0.5_kind_real*(lamda(i,j)+lamda(i,jp))
    lamda(i, j) = avg
    lamda(i,jp) = avg
    avg = 0.5_kind_real*(theta(i,j)-theta(i,jp))
    theta(i, j) =  avg
    theta(i,jp) = -avg
  enddo
enddo

end subroutine symm_ed

! --------------------------------------------------------------------------------------------------

subroutine mirror_grid(grid_global, ng, npx, npy, ndims, nregions)

integer,              intent(in)    :: ng, npx, npy, ndims, nregions
real(kind=kind_real), intent(inout) :: grid_global(1-ng:npx  +ng,1-ng:npy  +ng,ndims,1:nregions)

! Locals
integer :: i,j,n,n1,n2,nreg
real(kind=kind_real) :: x1, y1, z1, x2, y2, z2, ang

! Mirror Across the 0-longitude

nreg = 1

do j = 1,ceiling(npy/2.0_kind_real)
  do i = 1,ceiling(npx/2.0_kind_real)

    x1 = 0.25_kind_real * (abs(grid_global(i        ,j        ,1,nreg)) + &
                           abs(grid_global(npx-(i-1),j        ,1,nreg)) + &
                           abs(grid_global(i        ,npy-(j-1),1,nreg)) + &
                           abs(grid_global(npx-(i-1),npy-(j-1),1,nreg)))
    grid_global(i        ,j        ,1,nreg) = sign(x1,grid_global(i        ,j        ,1,nreg))
    grid_global(npx-(i-1),j        ,1,nreg) = sign(x1,grid_global(npx-(i-1),j        ,1,nreg))
    grid_global(i        ,npy-(j-1),1,nreg) = sign(x1,grid_global(i        ,npy-(j-1),1,nreg))
    grid_global(npx-(i-1),npy-(j-1),1,nreg) = sign(x1,grid_global(npx-(i-1),npy-(j-1),1,nreg))

    y1 = 0.25_kind_real * (abs(grid_global(i        ,j        ,2,nreg)) + &
                           abs(grid_global(npx-(i-1),j        ,2,nreg)) + &
                           abs(grid_global(i        ,npy-(j-1),2,nreg)) + &
                           abs(grid_global(npx-(i-1),npy-(j-1),2,nreg)))
    grid_global(i        ,j        ,2,nreg) = sign(y1,grid_global(i        ,j        ,2,nreg))
    grid_global(npx-(i-1),j        ,2,nreg) = sign(y1,grid_global(npx-(i-1),j        ,2,nreg))
    grid_global(i        ,npy-(j-1),2,nreg) = sign(y1,grid_global(i        ,npy-(j-1),2,nreg))
    grid_global(npx-(i-1),npy-(j-1),2,nreg) = sign(y1,grid_global(npx-(i-1),npy-(j-1),2,nreg))

    ! force dateline/greenwich-meridion consitency
    if (mod(npx,2) /= 0) then
      if ( (i==1+(npx-1)/2.0_kind_real) ) then
        grid_global(i,j        ,1,nreg) = 0.0_kind_real
        grid_global(i,npy-(j-1),1,nreg) = 0.0_kind_real
      endif
    endif

  enddo
enddo

do nreg=2,nregions
  do j=1,npy
    do i=1,npx

      x1 = grid_global(i,j,1,1)
      y1 = grid_global(i,j,2,1)
      z1 = radius

      if (nreg == 2) then

        ang = -90.0_kind_real
        call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis

      elseif (nreg == 3) then

        ang = -90.0_kind_real
        call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
        ang = 90.0_kind_real
        call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
        x2=x1
        y2=y1
        z2=z1

        ! Force North Pole and dateline/greenwich-meridion consitency
        if (mod(npx,2) /= 0) then
          if ( (i==1+(npx-1)/2.0_kind_real) .and. (i==j) ) then
            x2 = 0.0_kind_real
            y2 = pi/2.0_kind_real
          endif
          if ( (j==1+(npy-1)/2.0_kind_real) .and. (i < 1+(npx-1)/2.0_kind_real) ) then
            x2 = 0.0_kind_real
          endif
          if ( (j==1+(npy-1)/2.0_kind_real) .and. (i > 1+(npx-1)/2.0_kind_real) ) then
            x2 = pi
          endif
        endif

      elseif (nreg == 4) then

        ang = -180.0_kind_real
        call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
        ang = 90.0_kind_real
        call rot_3d( 1, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the x-axis
        x2=x1
        y2=y1
        z2=z1

        ! Force dateline/greenwich-meridion consitency
        if (mod(npx,2) /= 0) then
          if ( (j==1+(npy-1)/2.0_kind_real) ) then
            x2 = pi
          endif
        endif

      elseif (nreg == 5) then

        ang = 90.0_kind_real
        call rot_3d( 3, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the z-axis
        ang = 90.0_kind_real
        call rot_3d( 2, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the y-axis
        x2=x1
        y2=y1
        z2=z1

      elseif (nreg == 6) then

        ang = 90.0_kind_real
        call rot_3d( 2, x1, y1, z1, ang, x2, y2, z2, 1, 1)  ! rotate about the y-axis
        ang = 0.0_kind_real
        call rot_3d( 3, x2, y2, z2, ang, x1, y1, z1, 1, 1)  ! rotate about the z-axis
        x2=x1
        y2=y1
        z2=z1

        ! Force South Pole and dateline/greenwich-meridion consitency
        if (mod(npx,2) /= 0) then
          if ( (i==1+(npx-1)/2.0_kind_real) .and. (i==j) ) then
            x2 = 0.0_kind_real
            y2 = -pi/2.0_kind_real
          endif
          if ( (i==1+(npx-1)/2.0_kind_real) .and. (j > 1+(npy-1)/2.0_kind_real) ) then
            x2 = 0.0_kind_real
          endif
          if ( (i==1+(npx-1)/2.0_kind_real) .and. (j < 1+(npy-1)/2.0_kind_real) ) then
            x2 = pi
          endif
        endif

      endif

      grid_global(i,j,1,nreg) = x2
      grid_global(i,j,2,nreg) = y2

    enddo
  enddo
enddo

end subroutine mirror_grid

! --------------------------------------------------------------------------------------------------

subroutine rot_3d(axis, x1in, y1in, z1in, angle, x2out, y2out, z2out, degrees, convert)

integer,              intent(in)    :: axis    ! axis of rotation 1=x, 2=y, 3=z
real(kind=kind_real), intent(in)    :: x1in
real(kind=kind_real), intent(in)    :: y1in
real(kind=kind_real), intent(in)    :: z1in
real(kind=kind_real), intent(inout) :: angle   ! angle to rotate in radians
real(kind=kind_real), intent(out)   :: x2out
real(kind=kind_real), intent(out)   :: y2out
real(kind=kind_real), intent(out)   :: z2out
integer, optional,    intent(in)    :: degrees ! if present convert angle from degrees to radians
integer, optional,    intent(in)    :: convert ! if present convert input point from spherical to
                                               ! cartesian, rotate, and convert back

real(kind=kind_real)  :: c, s
real(kind=kind_real)  :: x1,y1,z1, x2,y2,z2

if ( present(convert) ) then

  call spherical_to_cartesian(x1in, y1in, z1in, x1, y1, z1)

else

  x1=x1in
  y1=y1in
  z1=z1in

endif

if ( present(degrees) ) then
  angle = angle*torad
endif

c = cos(angle)
s = sin(angle)

select case(axis)

case(1)
  x2 =  x1
  y2 =  c*y1 + s*z1
  z2 = -s*y1 + c*z1
case(2)
  x2 = c*x1 - s*z1
  y2 = y1
  z2 = s*x1 + c*z1
case(3)
  x2 =  c*x1 + s*y1
  y2 = -s*x1 + c*y1
  z2 = z1
case default
  write(*,*) "invalid axis: must be 1 for x, 2 for y, 3 for z."
end select

if ( present(convert) ) then
  call cartesian_to_spherical(x2, y2, z2, x2out, y2out, z2out)
else
  x2out=x2
  y2out=y2
  z2out=z2
endif

end subroutine rot_3d

! --------------------------------------------------------------------------------------------------

subroutine cartesian_to_spherical(x, y, z, lon, lat, r)

real(kind=kind_real) , intent(in)  :: x, y, z
real(kind=kind_real) , intent(out) :: lon, lat, r

r = sqrt(x*x + y*y + z*z)
if ( (abs(x) + abs(y)) < 1.e-10 ) then       ! poles:
  lon = 0.0_kind_real
else
  lon = atan2(y,x)    ! range: [-pi,pi]
endif

#ifdef RIGHT_HAND
  lat = asin(z/r)
#else
  lat = acos(z/r) - pi/2.
#endif

end subroutine cartesian_to_spherical

! --------------------------------------------------------------------------------------------------

subroutine spherical_to_cartesian(lon, lat, r, x, y, z)

real(kind=kind_real) , intent(in)  :: lon, lat, r
real(kind=kind_real) , intent(out) :: x, y, z

x = r * cos(lon) * cos(lat)
y = r * sin(lon) * cos(lat)

#ifdef RIGHT_HAND
  z =  r * sin(lat)
#else
  z = -r * sin(lat)
#endif

end subroutine spherical_to_cartesian

! --------------------------------------------------------------------------------------------------

subroutine direct_transform(c, i1, i2, j1, j2, lon_p, lat_p, n, lon, lat)

! This is a direct transformation of the standard (symmetrical) cubic grid
! to a locally enhanced high-res grid on the sphere; it is an application
! of the Schmidt transformation at the south pole followed by a
! pole_shift_to_target (rotation) operation

real(kind=kind_real),                         intent(in) :: c                ! Stretching factor
real(kind=kind_real),                         intent(in) :: lon_p, lat_p     ! center location of the target face, radian
integer,                                   intent(in) :: n                ! grid face number
integer,                                   intent(in) :: i1, i2, j1, j2   !  0 <= lon <= 2*pi ;    -pi/2 <= lat <= pi/2
real(kind=kind_real), dimension(i1:i2,j1:j2), intent(inout):: lon, lat

! Locals
real(kind_real) :: lat_t, sin_p, cos_p, sin_lat, cos_lat, sin_o, p2, two_pi
real(kind_real) :: c2p1, c2m1
integer :: i, j

p2 = 0.5_kind_real*pi
two_pi = 2.0_kind_real*pi

c2p1 = 1.0_kind_real + c*c
c2m1 = 1.0_kind_real - c*c

sin_p = sin(lat_p)
cos_p = cos(lat_p)

do j=j1,j2
  do i=i1,i2
    if ( abs(c2m1) > 1.0e-7_kind_real ) then
      sin_lat = sin(lat(i,j))
      lat_t = asin( (c2m1+c2p1*sin_lat)/(c2p1+c2m1*sin_lat) )
    else         ! no stretching
      lat_t = lat(i,j)
    endif
    sin_lat = sin(lat_t)
    cos_lat = cos(lat_t)
    sin_o = -(sin_p*sin_lat + cos_p*cos_lat*cos(lon(i,j)))
    if ( (1.-abs(sin_o)) < 1.0e-7_kind_real ) then    ! poles
      lon(i,j) = 0.0_kind_real
      lat(i,j) = sign( p2, sin_o )
    else
      lat(i,j) = asin( sin_o )
      lon(i,j) = lon_p + atan2( -cos_lat*sin(lon(i,j)), -sin_lat*cos_p+cos_lat*sin_p*cos(lon(i,j)))
      if ( lon(i,j) < 0.0_kind_real ) then
        lon(i,j) = lon(i,j) + two_pi
      elseif( lon(i,j) >= two_pi ) then
        lon(i,j) = lon(i,j) - two_pi
      endif
    endif
  enddo
enddo

end subroutine direct_transform

! --------------------------------------------------------------------------------------------------

end module fv3geom_utils_mod
