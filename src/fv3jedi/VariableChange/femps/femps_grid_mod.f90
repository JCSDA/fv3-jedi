! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_grid_mod

use netcdf
use femps_kinds_mod
use femps_utils_mod
use femps_const_mod

implicit none
private
public fempsgrid

! Class containing the grid
! -------------------------
type fempsgrid

  ! Grid type
  ! ---------
  character(len=2) :: gtype  ! cs (cubed-sphere), ih (icosahedral hexagons), di (diamonds)

  ! Mpi information
  ! ---------------
  integer :: comm, rank, csize

  ! Check convergence
  ! -----------------
  logical :: check_convergence

  ! Dimensions of each grid and hierarchy
  ! -------------------------------------
  integer :: ngrids      ! Number of grids in multigrid hierarchy
  integer :: ncubex      ! Number of faces on highest resolution grid
  integer :: nfacex      ! Number of faces on highest resolution grid
  integer :: nedgex      ! Number of edges on highest resolution each grid
  integer :: nvertx      ! Number of vertices on highest resolution edge

  integer, allocatable, dimension(:) :: ncube ! typical cube value (e.g. c96)
  integer, allocatable, dimension(:) :: nface ! faces on each grid
  integer, allocatable, dimension(:) :: nedge ! edges on each grid
  integer, allocatable, dimension(:) :: nvert ! vertices on each edge

  ! User provided lat/lon (radians)
  ! -------------------------------
  real(kind=kind_real), allocatable, dimension(:,:) :: flong ! longitude of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: flat  ! latitude of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: vlong ! longitude of vertices on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: vlat  ! latitude of vertices on each grid

  ! Coordinates and geometrical information
  ! ---------------------------------------
  real(kind=kind_real), allocatable, dimension(:,:) :: farea ! area of faces on each grid
  real(kind=kind_real), allocatable, dimension(:,:) :: ldist ! primal edge length, i.e. distance between neighbouring vertices
  real(kind=kind_real), allocatable, dimension(:,:) :: ddist ! dual edge length, i.e. distance between neighbouring face centres

  ! Connectivity
  ! ------------
  integer :: nefmx, nevmx                       ! maximum in neoff and neofv
  integer, allocatable, dimension(:,:) :: neoff ! number of edges and vertices of each face on each grid
  integer, allocatable, dimension(:,:) :: neofv ! number of edges of each vertex on each grid

  integer, private :: dimfnxtf, dimeoff, dimvoff, dimfnxte, dimvofe, dimfofv, dimeofv

  integer, allocatable, dimension(:,:,:) :: fnxtf   ! faces next to each face on each grid
  integer, allocatable, dimension(:,:,:) :: eoff    ! edges of each face on each grid
  integer, allocatable, dimension(:,:,:) :: voff    ! vertices of each face on each grid
  integer, allocatable, dimension(:,:,:) :: fnxte   ! faces either side of each edge on each grid
  integer, allocatable, dimension(:,:,:) :: vofe    ! vertices at the ends of each edge on each grid
  integer, allocatable, dimension(:,:,:) :: fofv    ! faces around each vertex on each grid
  integer, allocatable, dimension(:,:,:) :: eofv    ! edges incident on each vertex on each grid

  integer :: niter  ! Number of iterations for the solver

  contains

   procedure, public :: setup
   procedure, public :: delete
   procedure, public :: build_cs
   procedure, public :: writegrid
   procedure, public :: centroid
   procedure, public :: dual_centroid

end type fempsgrid

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup(self,gridtype,ngrids,cube,niter,comm,rank,csize,check_convergence)

implicit none
class(fempsgrid),  intent(inout) :: self
character(len=2),  intent(in)    :: gridtype
integer, optional, intent(in)    :: ngrids
integer, optional, intent(in)    :: cube
integer, optional, intent(in)    :: niter
integer, optional, intent(in)    :: comm
integer, optional, intent(in)    :: rank
integer, optional, intent(in)    :: csize
logical, optional, intent(in)    :: check_convergence

integer :: igrid, ncube

self%gtype = gridtype

self%niter = 10
if (present(niter)) self%niter = niter

! MPI information
self%comm = 0
self%rank = 0
self%csize = 1
if (present(comm)) self%comm = comm
if (present(rank)) self%rank = rank
if (present(csize)) self%csize = csize

! Check convergence
self%check_convergence = .false.
if (present(check_convergence)) self%check_convergence = check_convergence


if (self%gtype=='cs') then

  ! Set fempsgrid for a cubed sphere geometry
  ! -----------------------------------------
  if (.not.present(ngrids)) then
    self%ngrids = 6
  else
    self%ngrids = ngrids
  endif

  ! Allocate arrays based on ngrids only
  ! ------------------------------------
  allocate(self%ncube(self%ngrids))
  allocate(self%nface(self%ngrids))
  allocate(self%nedge(self%ngrids))
  allocate(self%nvert(self%ngrids))

  ! Set cube value for highest resolution grid
  ! ------------------------------------------
  if (.not.present(cube)) then
    self%ncubex = 96
  else
    self%ncubex = cube
  endif

  ! Other indicies for faces, verts, connectivity etc
  ! -------------------------------------------------
  self%nfacex = 6*self%ncubex*self%ncubex
  self%nedgex = 2*self%nfacex
  self%nvertx = self%nfacex + 2

  ! Fill arrays and report error if bad division of grid found
  ! ----------------------------------------------------------
  ncube = self%ncubex
  do igrid = self%ngrids,1,-1
    if (mod(ncube,2) == 0 .or. igrid == 1) then
      self%ncube(igrid) = ncube
      self%nface(igrid) = 6*ncube*ncube
      self%nedge(igrid) = 2*self%nface(igrid)
      self%nvert(igrid) = self%nface(igrid) + 2
    else
      call message('femps_grid_mod.setup: this number of grids is not possible', fatal)
    endif
    ncube = ncube / 2
  enddo

  ! Dimension for connectivity arrays
  ! ---------------------------------
  self%dimfnxtf = 4
  self%dimeoff  = 4
  self%dimvoff  = 4
  self%dimfnxte = 2
  self%dimvofe  = 2
  self%dimfofv  = 4
  self%dimeofv  = 4

elseif (self%gtype=='ih') then

  ! Set fempsgrid for a icosahedral hexagons geometry
  ! -------------------------------------------------
  self%ngrids = 7

  self%nfacex=5*2**(2*self%ngrids-1)+2
  self%nedgex=15*2**(2*self%ngrids-1)
  self%nvertx=5*2**(2*self%ngrids)

  self%dimfnxtf = 6
  self%dimeoff  = 6
  self%dimvoff  = 6
  self%dimfnxte = 2
  self%dimvofe  = 2
  self%dimfofv  = 3
  self%dimeofv  = 3

elseif (self%gtype=='di') then

  call message('Diamond grid cell not implemented here', fatal)

else

  call message('Grid type should be cs (cubed-sphere), ih (icosahedral hexagons) or di (diamonds)', fatal)

endif

! Allocate grid variables
allocate(self%neoff   (self%nfacex,self%ngrids))
allocate(self%neofv   (self%nvertx,self%ngrids))
allocate(self%fnxtf   (self%nfacex,self%dimfnxtf,self%ngrids))
allocate(self%eoff    (self%nfacex,self%dimeoff ,self%ngrids))
allocate(self%voff    (self%nfacex,self%dimvoff ,self%ngrids))
allocate(self%fnxte   (self%nedgex,self%dimfnxte,self%ngrids))
allocate(self%vofe    (self%nedgex,self%dimvofe ,self%ngrids))
allocate(self%fofv    (self%nvertx,self%dimfofv ,self%ngrids))
allocate(self%eofv    (self%nvertx,self%dimeofv ,self%ngrids))
allocate(self%flong   (self%nfacex,self%ngrids))
allocate(self%flat    (self%nfacex,self%ngrids))
allocate(self%vlong   (self%nvertx,self%ngrids))
allocate(self%vlat    (self%nvertx,self%ngrids))
allocate(self%farea   (self%nfacex,self%ngrids))
allocate(self%ldist   (self%nedgex,self%ngrids))
allocate(self%ddist   (self%nedgex,self%ngrids))

self%neoff = 0
self%neofv = 0
self%fnxtf = 0
self%eoff  = 0
self%voff  = 0
self%fnxte = 0
self%vofe  = 0
self%fofv  = 0
self%eofv  = 0
self%flong = 0.0_kind_real
self%flat  = 0.0_kind_real
self%vlong = 0.0_kind_real
self%vlat  = 0.0_kind_real
self%farea = 0.0_kind_real
self%ldist = 0.0_kind_real
self%ddist = 0.0_kind_real

end subroutine setup

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fempsgrid), intent(inout) :: self

! Deallocate grid variables
if(allocated(self%ncube)) deallocate(self%ncube)
if(allocated(self%nface)) deallocate(self%nface)
if(allocated(self%nedge)) deallocate(self%nedge)
if(allocated(self%nvert)) deallocate(self%nvert)
if(allocated(self%neoff)) deallocate(self%neoff)
if(allocated(self%neofv)) deallocate(self%neofv)
if(allocated(self%fnxtf)) deallocate(self%fnxtf)
if(allocated(self%eoff )) deallocate(self%eoff )
if(allocated(self%voff )) deallocate(self%voff )
if(allocated(self%fnxte)) deallocate(self%fnxte)
if(allocated(self%vofe )) deallocate(self%vofe )
if(allocated(self%fofv )) deallocate(self%fofv )
if(allocated(self%eofv )) deallocate(self%eofv )
if(allocated(self%flong)) deallocate(self%flong)
if(allocated(self%flat )) deallocate(self%flat )
if(allocated(self%vlong)) deallocate(self%vlong)
if(allocated(self%vlat )) deallocate(self%vlat )
if(allocated(self%farea)) deallocate(self%farea)
if(allocated(self%ldist)) deallocate(self%ldist)
if(allocated(self%ddist)) deallocate(self%ddist)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine build_cs(self,flavour_in,nsmooth_in)

implicit none
class(fempsgrid),  intent(inout) :: self
integer, optional, intent(in)    :: flavour_in
integer, optional, intent(in)    :: nsmooth_in

! flavour = cube sphere grid flavour
! 1 - equiangular cube (this is used with FV3)
! 2 - barycentric as in TCD14
! 3 - centroidal

! nsmooth = number of smoothing iterations.
! For flavour 1, no iterations are taken. nsmooth is ignored.
! For flavour 2, nsmooth = 1 is recommended. It must be at least 1 for a consistent FV H operator.
! For flavour 3, nsmooth = 3 is recommended.

! Module to generate a cubed sphere grid, including cross-
! reference tables of adjacent faces, edges and vertices,
! coordinates of faces and vertices, and lengths of edges
! and areas of faces.
!         .....
!        :  3  |
!        :     |
!         -----
!  -----  .....  .....  -----
! |  5  :|  1  :|  2  :|  4  :
! |     :|     :|     :|     :
!  .....  -----  -----  .....
!         .....
!        |  6  :
!        |     :
!         -----
! Solid lines: left and bottom edges of panel
! Dotted lines: right and top edges of panel
! John Thuburn Nov 2011
! Updated to include equiangular and centroidal variants 25/6/15

integer :: igrid, i, j, ixv, p1, p2, pp, jr, iv, iv0, n, n2, &
           ie0, ie1, ie2, if0, if1, if2, if3, iv1, iv2, ix1, ix2, &
           ixmin, ifmin, if21, if22, iv11, iv12, iv21, iv22, &
           ismooth, flavour, nsmooth

real(kind=kind_real) :: dlambda, lambda1, lambda2, t1, t2, long, lat, &
                        x1, y1, z1, x2, y2, z2, xc, yc, zc, x0, y0, z0, &
                        rmag, atri, aface, lmn, lmx, dmn, dmx, dav, cs, s, &
                        theta, thetamin, sn, d1x, d1y, d1z, d2x, d2y, d2z, &
                        lambdaf1, lambdaf2, tf1, tf2, &
                        amin, amax

logical :: lfound
character(len=2056) :: gridmessage

! Default setup
! -------------
flavour = 1
if (present(flavour_in)) flavour = flavour_in

nsmooth = 3
if (present(nsmooth_in)) nsmooth = nsmooth_in


! Generate the grid
! -----------------
do igrid = 1, self%ngrids

  ! Size of panels on this grid
  n = self%ncube(igrid)
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

  !   Set up incidence tables ignoring complications at
  !   panel edges
      do p1 = 1, 6
        ixv = (p1 - 1)*n2 + (j-1)*n + i

  !     Edges of the face
        self%eoff(ixv,1,igrid) = 2*ixv - 1
        self%eoff(ixv,2,igrid) = 2*ixv
        self%eoff(ixv,3,igrid) = 2*ixv + 1
        self%eoff(ixv,4,igrid) = 2*ixv + 2*n
  !     Vertices of the face
        self%voff(ixv,1,igrid) = ixv
        self%voff(ixv,2,igrid) = ixv + 1
        self%voff(ixv,3,igrid) = ixv + n + 1
        self%voff(ixv,4,igrid) = ixv + n
  !     Faces neighboring this face
        self%fnxtf(ixv,1,igrid) = ixv - 1
        self%fnxtf(ixv,2,igrid) = ixv - n
        self%fnxtf(ixv,3,igrid) = ixv + 1
        self%fnxtf(ixv,4,igrid) = ixv + n
  !     Edges incident on the vertex
        self%eofv(ixv,1,igrid) = 2*ixv - 2
        self%eofv(ixv,2,igrid) = 2*(ixv-n) - 1
        self%eofv(ixv,3,igrid) = 2*ixv
        self%eofv(ixv,4,igrid) = 2*ixv - 1

      enddo

    enddo
  enddo


  ! Now sort out complications at panel edges
  do j = 1, n
    jr = n + 1 - j

    do pp = 1, 3

      ! Odd numbered panels
      p1 = 2*pp - 1

      ! Left edge of panel p1 joins to top edge of panel p1 - 2
      ! Reverse order
      p2 = modulo(p1 + 3, 6) + 1
      ixv = (p1 - 1)*n2 + n*(j - 1) + 1
      self%fnxtf(ixv,1,igrid) = p2*n2 - n + jr
      self%eofv(ixv,1,igrid) = 2*p2*n2 - 2*n - 1 + 2*(jr + 1)

      ! Bottom edge of panel p1 joins to top edge of panel p1 - 1
      p2 = modulo(p1 + 4, 6) + 1
      ixv = (p1 - 1)*n2 + j
      self%fnxtf(ixv,2,igrid) = p2*n2 - n + j
      self%eofv(ixv,2,igrid) = 2*p2*n2 - 2*n - 1 + 2*j

      ! Right edge of panel p1 joins to left edge of panel p1 + 1
      p2 = modulo(p1, 6) + 1
      ixv = (p1 - 1)*n2 + n*j
      self%eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*(j - 1)*n + 1
      self%voff(ixv,2,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1
      self%voff(ixv,3,igrid) = (p2 - 1)*n2 + j*n + 1
      self%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + (j - 1)*n + 1

      ! Top edge of panel p1 joins to left edge of panel p1 + 2
      ! Reverse order
      p2 = modulo(p1 + 1, 6) + 1
      ixv = p1*n2 - n + j
      self%eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*(jr - 1)*n + 1
      self%voff(ixv,3,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1
      self%voff(ixv,4,igrid) = (p2 - 1)*n2 + jr*n + 1
      self%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + (jr - 1)*n + 1

      ! Even numbered panels
      p1 = 2*pp

      ! Left edge of panel p1 joins to right edge of panel p1 - 1
      p2 = modulo(p1 + 4, 6) + 1
      ixv = (p1 - 1)*n2 + n*(j - 1) + 1
      self%fnxtf(ixv,1,igrid) = (p2 - 1)*n2 + n*j
      self%eofv(ixv,1,igrid) = 2*(p2 - 1)*n2 + 2*n*j

      ! Bottom edge of panel p1 joins to right edge of panel p1 - 2
      ! Reverse order
      p2 = modulo(p1 + 3, 6) + 1
      ixv = (p1 - 1)*n2 + j
      self%fnxtf(ixv,2,igrid) = (p2 - 1)*n2 + n*jr
      self%eofv(ixv,2,igrid) = 2*(p2 - 1)*n2 + 2*n*(jr + 1)

      ! Top edge of panel p1 joins to bottom edge of panel p1 + 1
      p2 = modulo(p1, 6) + 1
      ixv = p1*n2 - n + j
      self%eoff(ixv,4,igrid) = 2*(p2 - 1)*n2 + 2*j
      self%voff(ixv,3,igrid) = (p2 - 1)*n2 + j + 1
      self%voff(ixv,4,igrid) = (p2 - 1)*n2 + j
      self%fnxtf(ixv,4,igrid) = (p2 - 1)*n2 + j

      ! Right edge of panel p1 joins to bottom edge of panel p1 + 2
      ! Reverse order
      p2 = modulo(p1 + 1, 6) + 1
      ixv = (p1-1)*n2 + n*j
      self%eoff(ixv,3,igrid) = 2*(p2 - 1)*n2 + 2*jr
      self%voff(ixv,2,igrid) = (p2 - 1)*n2 + jr + 1
      self%voff(ixv,3,igrid) = (p2 - 1)*n2 + jr
      self%fnxtf(ixv,3,igrid) = (p2 - 1)*n2 + jr

    enddo

  enddo

  ! All faces have 4 edges and vertices
  self%neoff(1:self%nface(igrid),igrid) = 4

  ! Almost all vertices have 4 edges (exceptions dealt with below)
  self%neofv(1:self%nvert(igrid),igrid) = 4

  ! Vertices not correctly captured by the above
  ! Corner vertices only have 3 edges
  do pp = 1, 3

    ! Bottom left of odd numbered panels
    p1 = 2*pp - 1
    ixv = (p1 - 1)*n2 + 1
    ! First edge needs to be deleted
    self%eofv(ixv,1,igrid) = self%eofv(ixv,2,igrid)
    self%eofv(ixv,2,igrid) = self%eofv(ixv,3,igrid)
    self%eofv(ixv,3,igrid) = self%eofv(ixv,4,igrid)
    self%eofv(ixv,4,igrid) = 0
    self%neofv(ixv,igrid) = 3

    ! Bottom left of even numbered panels
    p1 = 2*pp
    ixv = (p1 - 1)*n2 + 1
    ! Second edge needs to be deleted
    self%eofv(ixv,2,igrid) = self%eofv(ixv,3,igrid)
    self%eofv(ixv,3,igrid) = self%eofv(ixv,4,igrid)
    self%eofv(ixv,4,igrid) = 0
    self%neofv(ixv,igrid) = 3

  enddo

  ! Vertex 6*n2 + 1 is at top left of panels 1, 3, and 5
  iv = 6*n2 + 1
  do pp = 1, 3
    p1 = 2*pp - 1
    ixv = p1*n2 - n + 1
    self%voff(ixv,4,igrid) = iv
    self%eofv(iv,pp,igrid) = 2*p1*n2 - 2*n + 1
  enddo
  self%neofv(iv,igrid) = 3

  ! Vertex 6*n2 + 2 is at bottom right of panels 2, 4, and 6
  iv = 6*n2 + 2
  do pp = 1, 3
    p1 = 2*pp
    ixv = (p1 - 1)*n2 + n
    self%voff(ixv,2,igrid) = iv
    self%eofv(iv,pp,igrid) = 2*(p1 - 1)*n2 + 2*n
  enddo
  self%neofv(iv,igrid) = 3


enddo ! End of main loop over grids


! Now construct inverse tables
! First initialize entries to zero
self%fnxte = 0
self%vofe = 0
self%fofv = 0

do igrid = 1, self%ngrids

  do j = 1, self%nface(igrid)
    do i = 1, 4
      ie1 = self%eoff(j,i,igrid)
      call addtab(self%nedgex,2,self%fnxte(1,1,igrid),ie1,j)
      ixv = self%voff(j,i,igrid)
      call addtab(self%nvertx,4,self%fofv(1,1,igrid),ixv,j)
    enddo
  enddo

  do j = 1, self%nvert(igrid)
    do i = 1, 4
      ixv = self%eofv(j,i,igrid)
      if (ixv > 0) THEN
        call addtab(self%nedgex,2,self%vofe(1,1,igrid),ixv,j)
      endif
    enddo
  enddo

enddo

! Calculate geometrical quantities
do igrid = 1, self%ngrids

  ! Smoothing iterations
  if (flavour == 1) THEN

    ! No smoothing for equiangular cube

  elseif (flavour == 2) THEN

    ! TCD barycentric cube
    do ismooth = 1, nsmooth

      ! First locate face centres at barycentres of
      ! surrounding vertices
      do if1 = 1, self%nface(igrid)
        xc = 0.0_kind_real
        yc = 0.0_kind_real
        zc = 0.0_kind_real
        do i = 1, 4
          ixv = self%voff(if1,i,igrid)
          long = self%vlong(ixv,igrid)
          lat = self%vlat(ixv,igrid)
          call ll2xyz(long,lat,x1,y1,z1)
          xc = xc + x1
          yc = yc + y1
          zc = zc + z1
        enddo
        rmag = 1.0_kind_real/sqrt(xc*xc + yc*yc + zc*zc)
        xc = xc*rmag
        yc = yc*rmag
        zc = zc*rmag
        call xyz2ll(xc,yc,zc,long,lat)
        self%flong(if1,igrid) = long
        self%flat(if1,igrid) = lat
      enddo

      ! Next relocate vertices at barycentres of
      ! surrounding face centres - needed for H operator
      do iv1 = 1, self%nvert(igrid)
        xc = 0.0_kind_real
        yc = 0.0_kind_real
        zc = 0.0_kind_real
        do i = 1, self%neofv(iv1,igrid)
          if1 = self%fofv(iv1,i,igrid)
          long = self%flong(if1,igrid)
          lat = self%flat(if1,igrid)
          call ll2xyz(long,lat,x1,y1,z1)
          xc = xc + x1
          yc = yc + y1
          zc = zc + z1
        enddo
        rmag = 1.0_kind_real/sqrt(xc*xc + yc*yc + zc*zc)
        xc = xc*rmag
        yc = yc*rmag
        zc = zc*rmag
        call xyz2ll(xc,yc,zc,long,lat)
        self%vlong(iv1,igrid) = long
        self%vlat(iv1,igrid) = lat
      enddo

    enddo

  elseif (flavour == 3) THEN

    ! Centroidal cube
    do ismooth = 1, nsmooth
      ! Move faces to centroids
      do if0 = 1, self%nface(igrid)
        call self%centroid(if0,long,lat,igrid)
        self%flong(if0,igrid) = long
        self%flat(if0,igrid) = lat
      enddo
      ! Move vertices to centroids
      do iv0 = 1, self%nvert(igrid)
        call self%dual_centroid(iv0,long,lat,igrid)
        self%vlong(iv0,igrid) = long
        self%vlat(iv0,igrid) = lat
      enddo
    enddo

  else

    call message('flavour not recognized. Please pick another.', fatal)

  endif

! Tabulate areas
  do if1 = 1, self%nface(igrid)
    long = self%flong(if1,igrid)
    lat = self%flat(if1,igrid)
    call ll2xyz(long,lat,x0,y0,z0)
!   Compute face area
    aface = 0.0_kind_real
    do i = 1, 4
      ie1 = self%eoff(if1,i,igrid)
      ixv = self%vofe(ie1,1,igrid)
      long = self%vlong(ixv,igrid)
      lat = self%vlat(ixv,igrid)
      call ll2xyz(long,lat,x1,y1,z1)
      ixv = self%vofe(ie1,2,igrid)
      long = self%vlong(ixv,igrid)
      lat = self%vlat(ixv,igrid)
      call ll2xyz(long,lat,x2,y2,z2)
      call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,atri)
      aface = aface + atri
    enddo
    self%farea(if1,igrid) = aface
  enddo
  amin = MINVAL(self%farea(1:self%nface(igrid),igrid))
  amax = MAXVAL(self%farea(1:self%nface(igrid),igrid))
  write(gridmessage,*) 'Grid ',igrid,' min and max face area ',amin,amax,' ratio ',amin/amax
  call message(gridmessage, trace)

! Tabulate lengths of edges and distances between face centres
! across each edge
  lmn=5.0_kind_real
  lmx=0.0_kind_real
  dmn=5.0_kind_real
  dmx=0.0_kind_real
  dav=0.0_kind_real
  do ie0 = 1, self%nedge(igrid)
!   Vertices at ends of this edge
    iv1 = self%vofe(ie0,1,igrid)
    iv2 = self%vofe(ie0,2,igrid)
    long = self%vlong(iv1,igrid)
    lat = self%vlat(iv1,igrid)
    call ll2xyz(long,lat,x1,y1,z1)
    long = self%vlong(iv2,igrid)
    lat = self%vlat(iv2,igrid)
    call ll2xyz(long,lat,x2,y2,z2)
    call spdist(x1,y1,z1,x2,y2,z2,s)
    self%ldist(ie0,igrid) = s
    lmn = MIN(lmn,self%ldist(ie0,igrid))
    lmx = MIN(lmx,self%ldist(ie0,igrid))
!   Faces either side of this edge
    if1 = self%fnxte(ie0,1,igrid)
    if2 = self%fnxte(ie0,2,igrid)
    long = self%flong(if1,igrid)
    lat = self%flat(if1,igrid)
    call ll2xyz(long,lat,x1,y1,z1)
    long = self%flong(if2,igrid)
    lat = self%flat(if2,igrid)
    call ll2xyz(long,lat,x2,y2,z2)
    call spdist(x1,y1,z1,x2,y2,z2,s)
    self%ddist(ie0,igrid) = s
    dmn = MIN(dmn,self%ddist(ie0,igrid))
    dmx = MIN(dmx,self%ddist(ie0,igrid))
    dav = dav + self%ddist(ie0,igrid)/self%nedge(igrid)
  enddo

enddo


! Sort FNXTF into anticlockwise order on each grid
! and sort EOFF to correspond to FNXTF
! Also sort fofv into anticlockwise order
do igrid = 1, self%ngrids

  do if0 = 1, self%nface(igrid)

!   Coordinates of face if0
    long = self%flong(if0,igrid)
    lat = self%flat(if0,igrid)
    call ll2xyz(long,lat,x0,y0,z0)
    do ix1 = 1, 2
!     Coordinates of IX1'th neighbour
      if1 = self%fnxtf(if0,ix1,igrid)
      long = self%flong(if1,igrid)
      lat = self%flat(if1,igrid)
      call ll2xyz(long,lat,x1,y1,z1)
      d1x = x1 - x0
      d1y = y1 - y0
      d1z = z1 - z0
!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0
      do ix2 = ix1 + 1, 4
!       Coordinates of IX2'th neighbour
        if2 = self%fnxtf(if0,ix2,igrid)
        long = self%flong(if2,igrid)
        lat = self%flat(if2,igrid)
        call ll2xyz(long,lat,x2,y2,z2)
        d2x=x2 - x0
        d2y=y2 - y0
        d2z=z2 - z0
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0*(d1y*d2z - d1z*d2y) &
           + y0*(d1z*d2x - d1x*d2z) &
           + z0*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        if ((theta < thetamin) .AND. (theta > 0.0_kind_real)) THEN
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        endif
      enddo
!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = self%fnxtf(if0,ix1+1,igrid)
      self%fnxtf(if0,ix1+1,igrid) = ifmin
      self%fnxtf(if0,ixmin,igrid) = if3
    enddo

    do ix1 = 1, 4
      if1 = self%fnxtf(if0,ix1,igrid)
      ix2 = ix1 - 1
      lfound = .FALSE.
      do WHILE (.NOT. lfound)
         ix2 = ix2 + 1
         ie1 = self%eoff(if0,ix2,igrid)
         if21 = self%fnxte(ie1,1,igrid)
         if22 = self%fnxte(ie1,2,igrid)
         if ((if21 + if22) == (if0 + if1)) lfound = .TRUE.
      enddo
!     Edge IE2 corresponds to face IF1
      self%eoff(if0,ix2,igrid) = self%eoff(if0,ix1,igrid)
      self%eoff(if0,ix1,igrid) = ie1
    enddo

  enddo

  do iv0 = 1, self%nvert(igrid)
!   Coordinates of vertex iv0
    long = self%vlong(iv0,igrid)
    lat = self%vlat(iv0,igrid)
    call ll2xyz(long,lat,x0,y0,z0)
    do ix1 = 1, self%neofv(iv0,igrid) - 2
!     Coordinates of IX1'th face
      if1 = self%fofv(iv0,ix1,igrid)
      long = self%flong(if1,igrid)
      lat = self%flat(if1,igrid)
      call ll2xyz(long,lat,x1,y1,z1)
      d1x = x1 - x0
      d1y = y1 - y0
      d1z = z1 - z0
!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0
      do ix2 = ix1 + 1, self%neofv(iv0,igrid)
!       Coordinates of IX2'th neighbour
        if2 = self%fofv(iv0,ix2,igrid)
        long = self%flong(if2,igrid)
        lat = self%flat(if2,igrid)
        call ll2xyz(long,lat,x2,y2,z2)
        d2x=x2 - x0
        d2y=y2 - y0
        d2z=z2 - z0
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0*(d1y*d2z - d1z*d2y) &
           + y0*(d1z*d2x - d1x*d2z) &
           + z0*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        if ((theta < thetamin) .AND. (theta > 0.0_kind_real)) THEN
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        endif
      enddo
!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = self%fofv(iv0,ix1+1,igrid)
      self%fofv(iv0,ix1+1,igrid) = ifmin
      self%fofv(iv0,ixmin,igrid) = if3
    enddo
  enddo

enddo

!
! Order VOFF so that the k'th vertex is between the
! k'th and (k+1)'th edges in EOFF
do igrid = 1, self%ngrids
  do if0 = 1, self%nface(igrid)
    do ix1 = 1, self%neoff(if0,igrid)
      ix2 = ix1 + 1
      if (ix2 > self%neoff(if0,igrid)) ix2 = 1
      ie1 = self%eoff(if0,ix1,igrid)
      ie2 = self%eoff(if0,ix2,igrid)
      ! Find the common vertex of IE1 and IE2
      iv11 = self%vofe(ie1,1,igrid)
      iv12 = self%vofe(ie1,2,igrid)
      iv21 = self%vofe(ie2,1,igrid)
      iv22 = self%vofe(ie2,2,igrid)
      if ((iv11 == iv21) .OR. (iv11 == iv22)) THEN
        iv0 = iv11
      elseif ((iv12 == iv21) .OR. (iv12 == iv22)) THEN
        iv0 = iv12
      else
        call message('Common vertex not found',fatal)
      endif
      self%voff(if0,ix1,igrid) = iv0
    enddo
  enddo
enddo


! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
do igrid = 1, self%ngrids
  do ie0 = 1, self%nedge(igrid)
    if1 = self%fnxte(ie0,1,igrid)
    if2 = self%fnxte(ie0,2,igrid)
    long = self%flong(if1,igrid)
    lat = self%flat(if1,igrid)
    call ll2xyz(long,lat,x0,y0,z0)
    long = self%flong(if2,igrid)
    lat = self%flat(if2,igrid)
    call ll2xyz(long,lat,x1,y1,z1)
    d1x = x1 - x0
    d1y = y1 - y0
    d1z = z1 - z0
    iv1 = self%vofe(ie0,1,igrid)
    iv2 = self%vofe(ie0,2,igrid)
    long = self%vlong(iv1,igrid)
    lat = self%vlat(iv1,igrid)
    call ll2xyz(long,lat,x0,y0,z0)
    long = self%vlong(iv2,igrid)
    lat = self%vlat(iv2,igrid)
    call ll2xyz(long,lat,x1,y1,z1)
    d2x = x1 - x0
    d2y = y1 - y0
    d2z = z1 - z0
    sn = x0*(d1y*d2z - d1z*d2y) &
       + y0*(d1z*d2x - d1x*d2z) &
       + z0*(d1x*d2y - d1y*d2x)
    if (sn < 0.0_kind_real) THEN
      ! Swap the two vertices
      self%vofe(ie0,1,igrid) = iv2
      self%vofe(ie0,2,igrid) = iv1
    endif
  enddo
enddo

self%nefmx = maxval(self%neoff)
self%nevmx = maxval(self%neofv)

end subroutine build_cs

! --------------------------------------------------------------------------------------------------

subroutine writegrid(self,filename)

implicit none
class(fempsgrid), intent(in) :: self
character(len=*), intent(in) :: filename

integer :: ncid, vc, varid(1000)
integer ::   ngrids_dimid,  nfacex_dimid,  nedgex_dimid,   nvertx_dimid, &
           dimfnxtf_dimid, dimeoff_dimid, dimvoff_dimid, dimfnxte_dimid, &
            dimvofe_dimid, dimfofv_dimid, dimeofv_dimid


! Create file
! -----------
call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), ncid), "nf90_create" )


! Define dimensions
! -----------------
call nccheck( nf90_def_dim(ncid, "ngrids",   self%ngrids,     ngrids_dimid), "nf90_def_dim ngrids"  )
call nccheck( nf90_def_dim(ncid, "nfacex",   self%nfacex,     nfacex_dimid), "nf90_def_dim nfacex"  )
call nccheck( nf90_def_dim(ncid, "nvertx",   self%nvertx,     nvertx_dimid), "nf90_def_dim nvertx"  )
call nccheck( nf90_def_dim(ncid, "nedgex",   self%nedgex,     nedgex_dimid), "nf90_def_dim nedgex"  )
call nccheck( nf90_def_dim(ncid, "dimfnxtf", self%dimfnxtf, dimfnxtf_dimid), "nf90_def_dim dimfnxtf")
call nccheck( nf90_def_dim(ncid, "dimeoff",  self%dimeoff,   dimeoff_dimid), "nf90_def_dim dimeoff" )
call nccheck( nf90_def_dim(ncid, "dimvoff",  self%dimvoff,   dimvoff_dimid), "nf90_def_dim dimvoff" )
call nccheck( nf90_def_dim(ncid, "dimfnxte", self%dimfnxte, dimfnxte_dimid), "nf90_def_dim dimfnxte")
call nccheck( nf90_def_dim(ncid, "dimvofe",  self%dimvofe,   dimvofe_dimid), "nf90_def_dim dimvofe" )
call nccheck( nf90_def_dim(ncid, "dimfofv",  self%dimfofv,   dimfofv_dimid), "nf90_def_dim dimfofv" )
call nccheck( nf90_def_dim(ncid, "dimeofv",  self%dimeofv,   dimeofv_dimid), "nf90_def_dim dimeofv" )


! Define variables
! ----------------
vc = 1

call nccheck( nf90_def_var(ncid, "nface", NF90_INT,   (/ ngrids_dimid /), varid(vc)), "nf90_def_var nface" ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "nedge", NF90_INT,   (/ ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "nvert", NF90_INT,   (/ ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "neoff", NF90_INT,   (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "neofv", NF90_INT,   (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fnxtf", NF90_INT,   &
    (/ nfacex_dimid, dimfnxtf_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "eoff" , NF90_INT,   &
    (/ nfacex_dimid, dimeoff_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "voff" , NF90_INT,   &
    (/ nfacex_dimid, dimvoff_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fnxte", NF90_INT,   &
    (/ nedgex_dimid, dimfnxte_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vofe" , NF90_INT,   &
    (/ nedgex_dimid, dimvofe_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "fofv" , NF90_INT,   &
    (/ nvertx_dimid, dimfofv_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "eofv" , NF90_INT,   &
    (/ nvertx_dimid, dimeofv_dimid , ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "flong", NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "flat" , NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vlong", NF90_FLOAT, (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "vlat" , NF90_FLOAT, (/ nvertx_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "farea", NF90_FLOAT, (/ nfacex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "ldist", NF90_FLOAT, (/ nedgex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1
call nccheck( nf90_def_var(ncid, "ddist", NF90_FLOAT, (/ nedgex_dimid, ngrids_dimid /), varid(vc)), "nf90_def_var " ); vc = vc + 1


! End define
! ----------
call nccheck( nf90_enddef(ncid), "nf90_enddef" )


! Write variables
! ---------------
vc = 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nface ), "nf90_put_var nface" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nedge ), "nf90_put_var nedge" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%nvert ), "nf90_put_var nvert" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%neoff ), "nf90_put_var neoff" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%neofv ), "nf90_put_var neofv" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fnxtf ), "nf90_put_var fnxtf" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%eoff  ), "nf90_put_var eoff " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%voff  ), "nf90_put_var voff " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fnxte ), "nf90_put_var fnxte" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vofe  ), "nf90_put_var vofe " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%fofv  ), "nf90_put_var fofv " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%eofv  ), "nf90_put_var eofv " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%flong ), "nf90_put_var flong" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%flat  ), "nf90_put_var flat " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vlong ), "nf90_put_var vlong" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%vlat  ), "nf90_put_var vlat " ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%farea ), "nf90_put_var farea" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%ldist ), "nf90_put_var ldist" ); vc = vc + 1
call nccheck( nf90_put_var( ncid, varid(vc), self%ddist ), "nf90_put_var ddist" ); vc = vc + 1


! Close file
! ----------
call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine writegrid

! --------------------------------------------------------------------------------------------------

subroutine centroid(self,if0,long,lat,igrid)

implicit none
class(fempsgrid),     intent(in)  :: self
integer,              intent(in)  :: if0
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(out) :: long
real(kind=kind_real), intent(out) :: lat

integer :: ixe, ie1, iv1, iv2
real(kind=kind_real) :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
                        xc, yc, zc, a, aby3, mag

! Find the centroid of cell if0 on grid igrid
! -------------------------------------------

! Coordinates of `centre' of face (i.e. dual vertex)
long1 = self%flong(if0,igrid)
lat1 = self%flat(if0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the face
! Hence find area of face and centroid
xc = 0.0_kind_real
yc = 0.0_kind_real
zc = 0.0_kind_real
do ixe = 1, self%neoff(if0,igrid)
  ie1 = self%eoff(if0,ixe,igrid)
  iv1 = self%vofe(ie1,1,igrid)
  iv2 = self%vofe(ie1,2,igrid)
  long1 = self%vlong(iv1,igrid)
  lat1 = self%vlat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = self%vlong(iv2,igrid)
  lat1 = self%vlat(iv2,igrid)
  call ll2xyz(long1,lat1,x2,y2,z2)
  call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0_kind_real
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
enddo
mag = sqrt(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag
call xyz2ll(xc,yc,zc,long,lat)

end subroutine centroid

! --------------------------------------------------------------------------------------------------

subroutine dual_centroid(self,iv0,long,lat,igrid)

implicit none
class(fempsgrid),     intent(in)  :: self
integer,              intent(in)  :: iv0
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(out) :: long
real(kind=kind_real), intent(out) :: lat

integer :: ixe, ie1, iv1, iv2
real(kind=kind_real) :: long1, lat1, x0, y0, z0, x1, y1, z1, x2, y2, z2, &
                        xc, yc, zc, a, aby3, mag

! Find the centroid of dual cell iv0 on grid igrid
! ------------------------------------------------

! Coordinates of `centre' of dual cell (i.e. vertex)
long1 = self%vlong(iv0,igrid)
lat1 = self%vlat(iv0,igrid)
call ll2xyz(long1,lat1,x0,y0,z0)

! Loop over edges in turn and calculate area of triangle
! formed by the edge and the centre of the dual cell
! Hence find area of dual cell and centroid
xc = 0.0_kind_real
yc = 0.0_kind_real
zc = 0.0_kind_real
do ixe = 1, self%neofv(iv0,igrid)
  ie1 = self%eofv(iv0,ixe,igrid)
  iv1 = self%fnxte(ie1,1,igrid)
  iv2 = self%fnxte(ie1,2,igrid)
  long1 = self%flong(iv1,igrid)
  lat1 = self%flat(iv1,igrid)
  call ll2xyz(long1,lat1,x1,y1,z1)
  long1 = self%flong(iv2,igrid)
  lat1 = self%flat(iv2,igrid)
  call ll2xyz(long1,lat1,x2,y2,z2)
  call starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,a)
  aby3 = a/3.0_kind_real
  xc = xc + (x0 + x1 + x2)*aby3
  yc = yc + (y0 + y1 + y2)*aby3
  zc = zc + (z0 + z1 + z2)*aby3
enddo

mag = sqrt(xc*xc + yc*yc + zc*zc)
xc = xc/mag
yc = yc/mag
zc = zc/mag

call xyz2ll(xc,yc,zc,long,lat)

end subroutine dual_centroid

! --------------------------------------------------------------------------------------------------

end module femps_grid_mod

! Conventions
!
! Latitude and longitude in radians.
! Area and lengths are for the unit sphere.
!
! 1. fnxtf and eoff are ordered anticlockwise and such that
! the i'th edge lies between the face in question and its i'th
! neighbour.
!
! 2. voff are ordered anticlockwise such that the k'th vertex
! is the common vertex of the k'th and (k+1)'th edge.
!
! 3. eofv are ordered anticlockwise.
!
! 4. fofv are ordered anticlockwise such that the k'th face lies
! between the k'th and (k+1)'th edge.
!
! 5. The positive normal direction n points from
! fnxte(:,1,:) -> fnxte(:,2,:)
! and the positive tangential direction t points from
! vofe(:,1,:) -> vofe(:,2,:)
! such that t = k x n (where k is vertical).
