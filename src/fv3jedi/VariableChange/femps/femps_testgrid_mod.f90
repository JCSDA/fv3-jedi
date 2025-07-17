! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_testgrid_mod

use femps_kinds_mod
use femps_grid_mod
use femps_const_mod
use femps_utils_mod

implicit none
private
public cstestgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine cstestgrid(grid)

implicit none
type(fempsgrid),   intent(inout) :: grid

integer :: igrid, i, j, ixv, iv, n, n2
real(kind=kind_real) :: dlambda, lambda1, lambda2, t1, t2, long, lat, &
                        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, &
                        lambdaf1, lambdaf2, tf1, tf2, xf1, yf1, zf1, xf2, yf2, zf2, &
                        xf3, yf3, zf3, xf4, yf4, zf4, xf5, yf5, zf5, xf6, yf6, zf6

! Generate the grid
! -----------------
do igrid = 1, grid%ngrids

  ! Size of panels on this grid
  n = grid%ncube(igrid)
  n2 = n*n
  dlambda = 0.5_kind_real*pi/n

  !
  ! Loop over vertices/faces of one panel
  do j = 1, n
    lambda2 = (j-1)*dlambda - piby4
    lambdaf2 = (j-0.5_kind_real)*dlambda - piby4
    t2 = tan(lambda2)
    tf2 = tan(lambdaf2)
    do i = 1, n
      lambda1 = (i-1)*dlambda - piby4
      lambdaf1 = (i-0.5_kind_real)*dlambda - piby4
      t1 = tan(lambda1)
      tf1 = tan(lambdaf1)

  !   Set up coordinates of vertices and faces
  !   Panel 1
  !   Index of vertex
      ixv = (j-1)*n + i
  !   Cartesian coordinates of vertex
      x1 = 1.0_kind_real/sqrt(1.0_kind_real + t1*t1 + t2*t2)
      y1 = x1*t1
      z1 = x1*t2
  !   Lat long coordinates of vertex
      call xyz2ll(x1,y1,z1,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf1 = 1.0_kind_real/sqrt(1.0_kind_real + tf1*tf1 + tf2*tf2)
      yf1 = xf1*tf1
      zf1 = xf1*tf2
  !   Lat long coordinates of face
      call xyz2ll(xf1,yf1,zf1,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

  !   Panel 2
  !   Index of vertex
      ixv = ixv + n2
  !   Cartesian coordinates of vertex
      x2 = -y1
      y2 = x1
      z2 = z1
  !   Lat long coordinates of vertex
      call xyz2ll(x2,y2,z2,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf2 = -yf1
      yf2 = xf1
      zf2 = zf1
  !   Lat long coordinates of face
      call xyz2ll(xf2,yf2,zf2,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

  !   Panel 3
  !   Index of vertex
      ixv = ixv + n2
  !   Cartesian coordinates of vertex
      x3 = x2
      y3 = -z2
      z3 = y2
  !   Lat long coordinates of vertex
      call xyz2ll(x3,y3,z3,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf3 = xf2
      yf3 = -zf2
      zf3 = yf2
  !   Lat long coordinates of face
      call xyz2ll(xf3,yf3,zf3,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

  !   Panel 4
  !   Index of vertex
      ixv = ixv + n2
  !   Cartesian coordinates of vertex
      x4 = -z3
      y4 = y3
      z4 = x3
  !   Lat long coordinates of vertex
      call xyz2ll(x4,y4,z4,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf4 = -zf3
      yf4 = yf3
      zf4 = xf3
  !   Lat long coordinates of face
      call xyz2ll(xf4,yf4,zf4,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

  !   Panel 5
  !   Index of vertex
      ixv = ixv + n2
  !   Cartesian coordinates of vertex
      x5 = -y4
      y5 = x4
      z5 = z4
  !   Lat long coordinates of vertex
      call xyz2ll(x5,y5,z5,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf5 = -yf4
      yf5 = xf4
      zf5 = zf4
  !   Lat long coordinates of face
      call xyz2ll(xf5,yf5,zf5,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

  !   Panel 6
  !   Index of vertex
      ixv = ixv + n2
  !   Cartesian coordinates of vertex
      x6 = x5
      y6 = -z5
      z6 = y5
  !   Lat long coordinates of vertex
      call xyz2ll(x6,y6,z6,long,lat)
      grid%vlong(ixv,igrid) = long
      grid%vlat(ixv,igrid) = lat
  !   Cartesian coordinates of face
      xf6 = xf5
      yf6 = -zf5
      zf6 = yf5
  !   Lat long coordinates of face
      call xyz2ll(xf6,yf6,zf6,long,lat)
      grid%flong(ixv,igrid) = long
      grid%flat(ixv,igrid) = lat

    enddo
  enddo

  ! Vertex 6*n2 + 1 is at top left of panels 1, 3, and 5
  iv = 6*n2 + 1
  lambda2 = piby4
  t2 = tan(lambda2)
  lambda1 = - piby4
  t1 = tan(lambda1)
  ! Cartesian coordinates of vertex
  x1 = 1.0_kind_real/sqrt(1.0_kind_real + t1*t1 + t2*t2)
  y1 = x1*t1
  z1 = x1*t2
  ! Lat long coordinates of vertex
  call xyz2ll(x1,y1,z1,long,lat)
  grid%vlong(iv,igrid) = long
  grid%vlat(iv,igrid) = lat

  ! Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
  iv = 6*n2 + 2
  x1 = -x1
  y1 = -y1
  z1 = -z1
  ! Lat long coordinates of vertex
  call xyz2ll(x1,y1,z1,long,lat)
  grid%vlong(iv,igrid) = long
  grid%vlat(iv,igrid) = lat

enddo ! End of main loop over grids

end subroutine cstestgrid

! --------------------------------------------------------------------------------------------------

end module femps_testgrid_mod
