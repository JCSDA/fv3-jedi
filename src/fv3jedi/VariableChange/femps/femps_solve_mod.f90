! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_solve_mod

use femps_const_mod
use femps_utils_mod
use femps_grid_mod
use femps_operators_mod
use femps_kinds_mod

implicit none

private
public preliminary, laplace, inverselaplace

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine preliminary(grid,oprs)

! Preliminary calculations and setting up

implicit none
type(fempsgrid), intent(inout) :: grid
type(fempsoprs), intent(inout) :: oprs

! Build the operators
! -------------------
call oprs%setup(grid)
call oprs%build(grid)

! Dimensionalize
! --------------
grid%farea = grid%farea*rearth*rearth
grid%ldist = grid%ldist*rearth
grid%ddist = grid%ddist*rearth
oprs%varea = oprs%varea*rearth*rearth
oprs%lmass = oprs%lmass/(rearth*rearth)

! Build a lumped version of the J, M and H matrices
! -------------------------------------------------
call buildjlump(grid,oprs)
call buildmlump(grid,oprs)
call buildhlump(grid,oprs)

! And build a local operator to approximate the inverse of M
! ----------------------------------------------------------
call buildxminv(grid,oprs)

! Build coefficients used in Laplacian operator on all grids
! ----------------------------------------------------------
call buildlap(grid,oprs)

end subroutine preliminary

! --------------------------------------------------------------------------------------------------

subroutine buildlap(grid,oprs)

! Extract the diagonal coefficients of the Laplacian operator needed
! for the multigrid Poisson solver at all resolutions

implicit none
type(fempsgrid), intent(in)    :: grid
type(fempsoprs), intent(inout) :: oprs

integer :: igrid, if1, if2, if3, ie1, ie2, &
           ixd2, ixm, ixd1, ixl
real(kind=kind_real) :: temp, cd2, cm, cd1, cl


! Under-relaxation parameter. ( *** There is scope to optimize here *** )
! -----------------------------------------------------------------------
do igrid = 1, grid%ngrids
  if (grid%nefmx == 6) then
    oprs%underrel(igrid) = 0.8_kind_real
  elseif (grid%nefmx == 4) then
    oprs%underrel(igrid) = 0.8_kind_real
  else
    call message('Choose a sensible value for underrel in buildlap',fatal)
  endif
enddo

! Extract diagonal coefficient of Laplacian operator
do igrid = 1, grid%ngrids
  ! Loop over cells
  do if1 = 1, grid%nface(igrid)

    temp = 0.0_kind_real
    ! Loop over edges of if1 involved in Dprimal2 operator
    do ixd2 = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ixd2,igrid)
      cd2 = -oprs%eoffin(if1,ixd2,igrid)
      ! Loop over edges involved in approximate M-inverse operator
      do ixm = 1, oprs%nxminvsten(ie1,igrid)
        ie2 = oprs%xminvsten(ie1,ixm,igrid)
        cm = oprs%xminv(ie1,ixm,igrid)
        ! Loop over cells involved in Ddual1 operator
        cd1 = 1.0_kind_real
        do ixd1 = 1, 2
          if2 = grid%fnxte(ie2,ixd1,igrid)
          cd1 = -cd1
          ! Loop over cells in L operator
          do ixl = 1, oprs%nlsten(if2,igrid)
            if3 = oprs%lsten(if2,ixl,igrid)
            cl = oprs%lmass(if2,ixl,igrid)
            IF (if3 == if1) temp = temp + cd2*cm*cd1*cl
          enddo
        enddo
      enddo
    enddo
    oprs%lapdiag(if1,igrid) = temp

  enddo
enddo

end subroutine buildlap

! --------------------------------------------------------------------------------------------------

subroutine Dprimal1(grid,f,df,igrid,nv,ne)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises pointwise values
! at vertices; df comprises integrals of the derivative
! of f along primal cell edges.

implicit none
type(fempsgrid),      intent(in)  :: grid
integer,              intent(in)  :: igrid, nv, ne
real(kind=kind_real), intent(in)  :: f(nv)
real(kind=kind_real), intent(out) :: df(ne)

integer :: ie1, iv1, iv2

do ie1 = 1, ne
  iv1 = grid%vofe(ie1,1,igrid)
  iv2 = grid%vofe(ie1,2,igrid)
  df(ie1) = f(iv2) - f(iv1)
enddo

end subroutine Dprimal1

! --------------------------------------------------------------------------------------------------

subroutine Dprimal2(grid,oprs,f,df,igrid,ne,nf)

! To compute the exterior derivative df of the field f
! on primal grid number igrid. f comprises integrals along
! primal edges; df comprises integrals of the derivative
! over primal cells.

implicit none
type(fempsgrid), intent(in) :: grid
type(fempsoprs), intent(in) :: oprs
integer,         intent(in) :: igrid, ne, nf

real(kind=kind_real),       intent(in)  :: f(ne)
real(kind=kind_real),       intent(out) :: df(nf)

integer :: if1, ix, ie1
real(kind=kind_real) :: temp

do if1 = 1, nf
  temp = 0.0_kind_real
  do ix = 1, grid%neoff(if1,igrid)
    ie1 = grid%eoff(if1,ix,igrid)
    temp = temp - f(ie1)*oprs%eoffin(if1,ix,igrid)
  enddo
  df(if1) = temp
enddo

end subroutine Dprimal2

! --------------------------------------------------------------------------------------------------

subroutine Ddual1(grid,f,df,igrid,nf,ne)

! To compute the exterior derivative df of the field f
! on dual grid number igrid. f comprises pointwise values
! at face centres; df comprises integrals of the derivative
! of f along dual cell edges.

implicit none
type(fempsgrid),      intent(in) :: grid
integer,              intent(in)  :: igrid, nf, ne
real(kind=kind_real), intent(in)  :: f(nf)
real(kind=kind_real), intent(out) :: df(ne)

integer :: ie1, if1, if2

do ie1 = 1, ne
  if1 = grid%fnxte(ie1,1,igrid)
  if2 = grid%fnxte(ie1,2,igrid)
  df(ie1) = f(if2) - f(if1)
enddo

end subroutine Ddual1

! --------------------------------------------------------------------------------------------------

subroutine buildjlump(grid,oprs)

! Build a lumped version of the j matrix

implicit none
type(fempsgrid), intent(in)    :: grid
type(fempsoprs), intent(inout) :: oprs

! Two choices for lumped J
! ijlump = 1     Diagonal of J
! ijlump = 2     Exact answer when applied to constant input vector
! ijlump = 3     Exact answer when applied to the discrete representation
!                  of a constant scalar (input vector proportional to varea)
integer, parameter :: ijlump = 3
integer :: igrid, iv1, nv
real(kind=kind_real), allocatable :: p1(:), jp1(:)
character(len=2056) :: errormessage

do igrid = 1, grid%ngrids

  if (ijlump == 1) then

    do iv1 = 1, grid%nvert(igrid)
      oprs%jlump(iv1,igrid) = oprs%jstar(iv1,1,igrid)
    enddo

  elseif (ijlump == 2) then

    nv = grid%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = 1.0_kind_real
    call HodgeJ(oprs,p1,jp1,igrid,nv)
    oprs%jlump(1:nv,igrid) = jp1
    deallocate(p1,jp1)

  elseif (ijlump == 3) then

    nv = grid%nvert(igrid)
    allocate(p1(nv), jp1(nv))
    p1 = oprs%varea(1:nv,igrid)
    call HodgeJ(oprs,p1,jp1,igrid,nv)
    oprs%jlump(1:nv,igrid) = jp1/oprs%varea(1:nv,igrid)
    deallocate(p1,jp1)

  else

    write(errormessage,*) 'option ijlump = ',ijlump,' not available in subroutine buildjlump'
    call message(errormessage,fatal)

  endif

enddo


end subroutine buildjlump

! --------------------------------------------------------------------------------------------------

subroutine buildmlump(grid,oprs)

! Build a lumped version of the M matrix

implicit none
type(fempsgrid), intent(in)    :: grid
type(fempsoprs), intent(inout) :: oprs

! Three choices for lumped M
! imlump = 1     Diagonal of M
! imlump = 2     Row sum of absolute values of M
! imlump = 3     Best fit to solid body rotation normal to each edge
integer, parameter :: imlump = 3
integer :: igrid, ie1, ie2, ix, iv0, ne, nv
real(kind=kind_real) :: temp, slon, slat, clon, clat, a1, a2, a3, num, den
real(kind=kind_real), allocatable :: psi1(:), psi2(:), psi3(:),    &
                                     v1(:), v2(:), v3(:), vv(:),   &
                                     mv1(:), mv2(:), mv3(:)
character(len=2056) :: errormessage

do igrid = 1, grid%ngrids

  if (imlump == 1) then

    do ie1 = 1, grid%nedge(igrid)
      oprs%mlump(ie1,igrid) = oprs%mmass(ie1,1,igrid)
    enddo

  elseif (imlump == 2) then

    do ie1 = 1, grid%nedge(igrid)
      temp = 0.0_kind_real
      do ix = 1, oprs%nmsten(ie1,igrid)
        ie2 = oprs%msten(ie1,ix,igrid)
        temp = temp + abs(oprs%mmass(ie1,ix,igrid))
      enddo
      oprs%mlump(ie1,igrid) = temp
    enddo

  elseif (imlump == 3) then

    nv = grid%nvert(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),    &
             v1(ne), v2(ne), v3(ne), vv(ne),  &
             mv1(ne), mv2(ne), mv3(ne))
    ! Set up three solid body rotation fields
    do iv0 = 1, nv
      slon = sin(grid%vlong(iv0,igrid))
      clon = cos(grid%vlong(iv0,igrid))
      slat = sin(grid%vlat(iv0,igrid))
      clat = cos(grid%vlat(iv0,igrid))
      psi1(iv0) = clat*clon
      psi2(iv0) = clat*slon
      psi3(iv0) = slat
    enddo
    call Dprimal1(grid,psi1,v1,igrid,nv,ne)
    call Dprimal1(grid,psi2,v2,igrid,nv,ne)
    call Dprimal1(grid,psi3,v3,igrid,nv,ne)
    call massM(oprs,v1,mv1,igrid,ne)
    call massM(oprs,v2,mv2,igrid,ne)
    call massM(oprs,v3,mv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      den = sqrt(a1*a1 + a2*a2 + a3*a3)
      a1 = a1/den
      a2 = a2/den
      a3 = a3/den
      ! Demand that lumped matrix agrees with full M for this v
      num = a1*mv1(ie1) + a2*mv2(ie1) + a3*mv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      oprs%mlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,vv,mv1,mv2,mv3)

  else

    write(errormessage,*) 'option imlump = ',imlump,' not available in subroutine buildmlump'
    call message(errormessage,fatal)

  endif

enddo

end subroutine buildmlump

! --------------------------------------------------------------------------------------------------

subroutine buildhlump(grid,oprs)

! Build a lumped version of the H matrix

implicit none
type(fempsgrid), intent(in)    :: grid
type(fempsoprs), intent(inout) :: oprs

! Three choices for lumped H
! ihlump = 1     Diagonal of H
! ihlump = 2     Row sum of absolute values of H
! ihlump = 3     Best fit to global divergent flow on dual grid
! ihlump = 4     Best fit to solid body rotation normal to each edge
integer, parameter :: ihlump = 3
integer :: igrid, ie1, ie2, ix, iv0, if0, ne, nv, nf
real(kind=kind_real) :: temp, slon, slat, clon, clat, a1, a2, a3, num, den
real(kind=kind_real), allocatable :: psi1(:), psi2(:), psi3(:),  &
                                     v1(:), v2(:), v3(:),        &
                                     hv1(:), hv2(:), hv3(:)
character(len=2056) :: errormessage

do igrid = 1, grid%ngrids

  if (ihlump == 1) then

    do ie1 = 1, grid%nedge(igrid)
      oprs%hlump(ie1,igrid) = oprs%hstar(ie1,1,igrid)
    enddo

  elseif (ihlump == 2) then

    do ie1 = 1, grid%nedge(igrid)
      temp = 0.0_kind_real
      do ix = 1, oprs%nhsten(ie1,igrid)
        ie2 = oprs%hsten(ie1,ix,igrid)
        temp = temp + abs(oprs%hstar(ie1,ix,igrid))
      enddo
      oprs%hlump(ie1,igrid) = temp
    enddo

  elseif (ihlump == 3) then

    nf = grid%nface(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nf), psi2(nf), psi3(nf),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three global divergent fields
    do if0 = 1, nf
      slon = sin(grid%flong(if0,igrid))
      clon = cos(grid%flong(if0,igrid))
      slat = sin(grid%flat(if0,igrid))
      clat = cos(grid%flat(if0,igrid))
      psi1(if0) = slat
      psi2(if0) = clat*slon
      psi3(if0) = clat*clon
    enddo
    call Ddual1(grid,psi1,v1,igrid,nf,ne)
    call Ddual1(grid,psi2,v2,igrid,nf,ne)
    call Ddual1(grid,psi3,v3,igrid,nf,ne)
    call HodgeH(oprs,v1,hv1,igrid,ne)
    call HodgeH(oprs,v2,hv2,igrid,ne)
    call HodgeH(oprs,v3,hv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      oprs%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  elseif (ihlump == 4) then

    nv = grid%nvert(igrid)
    ne = grid%nedge(igrid)
    allocate(psi1(nv), psi2(nv), psi3(nv),  &
             v1(ne), v2(ne), v3(ne),        &
             hv1(ne), hv2(ne), hv3(ne))
    ! Set up three solid body rotation fields
    do iv0 = 1, nv
      slon = sin(grid%vlong(iv0,igrid))
      clon = cos(grid%vlong(iv0,igrid))
      slat = sin(grid%vlat(iv0,igrid))
      clat = cos(grid%vlat(iv0,igrid))
      psi1(iv0) = slat
      psi2(iv0) = clat*slon
      psi3(iv0) = clat*clon
    enddo
    call Dprimal1(grid,psi1,v1,igrid,nv,ne)
    call Dprimal1(grid,psi2,v2,igrid,nv,ne)
    call Dprimal1(grid,psi3,v3,igrid,nv,ne)
    call HodgeH(oprs,v1,hv1,igrid,ne)
    call HodgeH(oprs,v2,hv2,igrid,ne)
    call HodgeH(oprs,v3,hv3,igrid,ne)
    ! Now loop over edges
    do ie1 = 1, ne
      ! Velocity field that maximizes v(ie1) is
      ! v = a1*v1 + a2*v2 + a3*v3
      ! with
      a1 = v1(ie1)
      a2 = v2(ie1)
      a3 = v3(ie1)
      ! Demand that lumped matrix agrees with full H for this v
      num = a1*hv1(ie1) + a2*hv2(ie1) + a3*hv3(ie1)
      den = a1*v1(ie1)  + a2*v2(ie1)  + a3*v3(ie1)
      oprs%hlump(ie1,igrid) = num/den
    enddo
    deallocate(psi1,psi2,psi3,v1,v2,v3,hv1,hv2,hv3)

  else

    write(errormessage,*) 'option ihlump = ',ihlump,' not available in subroutine buildhlump'
    call message(errormessage,fatal)

  endif

enddo

end subroutine buildhlump

! --------------------------------------------------------------------------------------------------

subroutine buildxminv(grid,oprs)

! Determine stencil and coefficients for a local approximation
! to the inverse of M.
!
! Two options are coded. The first is simply the inverse of
! the diagonal mass lumped approximation to M. The second is
! based on a single underrelaxed Jacobi iteration to the
! inverse of M.

implicit none
type(fempsgrid), intent(in)    :: grid
type(fempsoprs), intent(inout) :: oprs

logical :: llump
integer :: igrid, ie0, ix, ie1
real(kind=kind_real) :: temp, diag
real(kind=kind_real) :: relax

! Underrelaxation coefficient depends on grid
! 0.9 and 1.4 when using mlump in Jacobi
! Check for consistency with massMinv
if (grid%nefmx == 4) then
  ! use sparse approximate inverse on cubed sphere
  llump = .false.
  relax = 0.9_kind_real
else
  ! Use diagonal approximate inverse on hexagonal-icosahedral grid
  llump = .false.
  relax = 1.4_kind_real
endif

if (llump) then

  ! Diagonal operator: stencil size is 1
  oprs%nxmisx = 1
  allocate(oprs%xminvsten(grid%nedgex,oprs%nxmisx,grid%ngrids))
  allocate(oprs%xminv(grid%nedgex,oprs%nxmisx,grid%ngrids))
  do igrid = 1, grid%ngrids
    do ie0 = 1, grid%nedge(igrid)
      ! stencil for edge ie0 is ie0 itself, and coeff is
      ! inverse of the diagonal term of the lumped matrix
      oprs%nxminvsten(ie0,igrid) = 1
      oprs%xminvsten(ie0,1,igrid) = ie0
      oprs%xminv(ie0,1,igrid) = 1.0_kind_real/oprs%mlump(ie0,igrid)
    enddo
  enddo

else

  ! Stencil is the same as for M itself
  oprs%nxmisx = oprs%nmsmx
  allocate(oprs%xminvsten(grid%nedgex,oprs%nxmisx,grid%ngrids))
  allocate(oprs%xminv(grid%nedgex,oprs%nxmisx,grid%ngrids))
  do igrid = 1, grid%ngrids
    do ie0 = 1, grid%nedge(igrid)
      ! stencil for edge ie0 is the same as the stencil for m
      oprs%nxminvsten(ie0,igrid) = oprs%nmsten(ie0,igrid)
      do ix = 1, oprs%nmsten(ie0,igrid)
        ie1 = oprs%msten(ie0,ix,igrid)
        oprs%xminvsten(ie0,ix,igrid) = ie1
        if (ie1 == ie0) then
          diag = 1.0_kind_real + relax
        else
          diag = 0.0_kind_real
        endif
        temp = oprs%mmass(ie0,ix,igrid)/oprs%mlump(ie1,igrid)
        oprs%xminv(ie0,ix,igrid) = (diag - relax*temp)/oprs%mlump(ie0,igrid)
      enddo
    enddo
  enddo

endif

end subroutine buildxminv

! --------------------------------------------------------------------------------------------------

subroutine massL(oprs,f1,f2,igrid,nf)

! Apply the mass matrix L to field f1 to obtain the result f2

implicit none
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, nf
real(kind=kind_real), intent(in)  :: f1(nf)
real(kind=kind_real), intent(out) :: f2(nf)

integer :: if1, if2, ix
real(kind=kind_real) :: temp

do if1 = 1, nf
  temp = 0.0_kind_real
  do ix = 1, oprs%nlsten(if1,igrid)
    if2 = oprs%lsten(if1,ix,igrid)
    temp = temp + f1(if2)*oprs%lmass(if1,ix,igrid)
  enddo
  f2(if1) = temp
enddo

end subroutine massL

! --------------------------------------------------------------------------------------------------

subroutine massM(oprs,f1,f2,igrid,ne)

! Apply the mass matrix M to field f1 to obtain field f2

implicit none
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, ne
real(kind=kind_real), intent(in)  :: f1(ne)
real(kind=kind_real), intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real(kind=kind_real) :: temp

do ie1 = 1, ne
  temp = 0.0_kind_real
  do ix = 1, oprs%nmsten(ie1,igrid)
    ie2 = oprs%msten(ie1,ix,igrid)
    temp = temp + f1(ie2)*oprs%mmass(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine massM

! --------------------------------------------------------------------------------------------------

subroutine approxMinv(oprs,f1,f2,igrid,ne)

! Apply an approximate inverse of the mass matrix M
! to field f1 to obtain field f2

implicit none
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, ne
real(kind=kind_real), intent(in)  :: f1(ne)
real(kind=kind_real), intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real(kind=kind_real) :: temp

do ie1 = 1, ne
  temp = 0.0_kind_real
  do ix = 1, oprs%nxminvsten(ie1,igrid)
    ie2 = oprs%xminvsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*oprs%xminv(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine approxMinv

! --------------------------------------------------------------------------------------------------

subroutine HodgeJ(oprs,f1,f2,igrid,nv)

! Apply the Hodge star J operator that converts dual face
! integrals f1 to vertex values f2 on grid igrid.

implicit none
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, nv
real(kind=kind_real), intent(in)  :: f1(nv)
real(kind=kind_real), intent(out) :: f2(nv)

integer :: iv1, iv2, ix
real(kind=kind_real) :: temp

do iv1 = 1, nv
  temp = 0.0_kind_real
  do ix = 1, oprs%njsten(iv1,igrid)
    iv2 = oprs%jsten(iv1,ix,igrid)
    temp = temp + f1(iv2)*oprs%jstar(iv1,ix,igrid)
  enddo
  f2(iv1) = temp
enddo

end subroutine HodgeJ

! --------------------------------------------------------------------------------------------------

subroutine HodgeH(oprs,f1,f2,igrid,ne)

! Apply the Hodge star H operator that converts dual edge
! integrals (circulations) f1 to primal edge integrals (fluxes) f2
! on grid igrid.

implicit none
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, ne
real(kind=kind_real), intent(in)  :: f1(ne)
real(kind=kind_real), intent(out) :: f2(ne)

integer :: ie1, ie2, ix
real(kind=kind_real) :: temp

do ie1 = 1, ne
  temp = 0.0_kind_real
  do ix = 1, oprs%nhsten(ie1,igrid)
    ie2 = oprs%hsten(ie1,ix,igrid)
    temp = temp + f1(ie2)*oprs%hstar(ie1,ix,igrid)
  enddo
  f2(ie1) = temp
enddo

end subroutine HodgeH

! --------------------------------------------------------------------------------------------------

subroutine massMinv(grid,oprs,f1,f2,igrid,ne,niter)

! Apply the inverse of the mass matrix M to the field f1
! to obtain the result f2

! If niter >= 0 then f2 is assumed to be an initial estimate
! for the solution, and niter further iterations are taken.
! If niter < 0 then f2 is initialized and -niter
! iterations are taken.
! If M is diagonal then there is no need to iterate.

implicit none

type(fempsgrid),      intent(in)    :: grid
type(fempsoprs),      intent(in)    :: oprs
integer,              intent(in)    :: niter
integer,              intent(in)    :: igrid, ne
real(kind=kind_real), intent(in)    :: f1(ne)
real(kind=kind_real), intent(inout) :: f2(ne)

integer :: ie1, iter, miter
real(kind=kind_real) :: temp(ne)
real(kind=kind_real) :: relax

! Underrelaxation coefficient depends on grid
IF (grid%nefmx == 4) THEN
  relax = 0.9_kind_real
ELSE
  relax = 1.4_kind_real
ENDIF

miter = ABS(niter)

IF (niter < 0 .OR. oprs%nmsmx == 1) THEN
  ! First guess based on lumped M
  do ie1 = 1, ne
    f2(ie1) = f1(ie1)/oprs%mlump(ie1,igrid)
  enddo
ENDIF

IF (oprs%nmsmx > 1) THEN
  ! M is not diagonal, so use Jacobi iteration to invert
  do iter = 1, miter
    call massM(oprs,f2,temp,igrid,ne)
    temp = relax*(f1 - temp)
    do ie1 = 1, ne
      f2(ie1) = f2(ie1) + temp(ie1)/oprs%mlump(ie1,igrid)
    enddo
  enddo
ENDIF

end subroutine massMinv

! --------------------------------------------------------------------------------------------------

subroutine restrict(grid,oprs,f1,nf1,f2,nf2,igrid)

! To perform the restriction operation needed for a multigrid solver.
! Restrict field f1 from grid igrid + 1 to grid igrid and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms).

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: nf1, nf2, igrid
real(kind=kind_real), intent(in)  :: f1(nf1)
real(kind=kind_real), intent(out) :: f2(nf2)

integer :: if1, if2, ix
real(kind=kind_real) :: wgt

! Safety check
if (nf2 .ne. grid%nface(igrid)) then
  call message('Wrong size array in subroutine restrict',fatal)
endif

do if2 = 1, nf2
  f2(if2) = 0.0_kind_real
  do ix = 1, oprs%ninj(if2,igrid)
    if1 = oprs%injsten(if2,ix,igrid)
    wgt = oprs%injwgt(if2,ix,igrid)
    f2(if2) = f2(if2) + wgt*f1(if1)
  enddo
enddo

end subroutine restrict

! --------------------------------------------------------------------------------------------------

subroutine prolong(grid,oprs,f2,nf2,f1,nf1,igrid)

! To perform the prolongation operation needed for a multigrid solver.
! Prolong field f2 from grid igrid to grid igrid + 1 and put the
! result in field f2. f1 and f2 are assumed to be area integrals
! (discrete 2-forms); so f2 must be converted to point values
! (zero-form) first and f1 must be converted back to an area-integral
! (2-form) at the end. The prolongation operator is the adjoint of
! the restriction operator, so uses the same stencil and weights.

! *** Currently this operator uses a loop over source entities,
! implying a need for an MPI reduce. The stencil and weights could
! be stored in transpose form to allow a loop over target entities
! and hence avoid an MPI reduce ***

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: nf1, nf2, igrid
real(kind=kind_real), intent(in)  :: f2(nf2)
real(kind=kind_real), intent(out) :: f1(nf1)

integer :: if1, if2, ix, igridp
real(kind=kind_real) :: wgt, f2if2, temp1(nf1), temp2(nf2)

! Safety check
if (nf2 .ne. grid%nface(igrid)) then
  call message('Wrong size array in subroutine prolong',fatal)
endif

igridp = igrid + 1
temp2(1:nf2) = f2(1:nf2)/grid%farea(1:nf2,igrid)
temp1 = 0.0_kind_real
do if2 = 1, nf2
  f2if2 = temp2(if2)
  do ix = 1, oprs%ninj(if2,igrid)
    if1 = oprs%injsten(if2,ix,igrid)
    wgt = oprs%injwgt(if2,ix,igrid)
    temp1(if1) = temp1(if1) + wgt*f2if2
  enddo
enddo
f1(1:nf1) = temp1(1:nf1)*grid%farea(1:nf1,igridp)

end subroutine prolong

! --------------------------------------------------------------------------------------------------

subroutine laplace(grid,oprs,igrid,f,hf)

! To apply the Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(in)  :: f(grid%nface(igrid))
real(kind=kind_real), intent(out) :: hf(grid%nface(igrid))

real(kind=kind_real) :: temp1(grid%nface(igrid)), temp2(grid%nedge(igrid)), temp3(grid%nedge(igrid))
integer :: niter

call massL(oprs,f,temp1,igrid,grid%nface(igrid))
call Ddual1(grid,temp1,temp2,igrid,grid%nface(igrid),grid%nedge(igrid))
niter = -20
call massMinv(grid,oprs,temp2,temp3,igrid,grid%nedge(igrid),niter)
call Dprimal2(grid,oprs,temp3,hf,igrid,grid%nedge(igrid),grid%nface(igrid))

end subroutine laplace

! --------------------------------------------------------------------------------------------------

subroutine inverselaplace(grid,oprs,igrid,hf,f,removemean_in,field_conv)

! To apply the inverse Laplacian operator to the input field hf,
! on grid igrid, the result appearing in the output field f.
! Note f and hf are area integrals (2-forms).

implicit none
type(fempsgrid),                intent(in)     :: grid
type(fempsoprs),                intent(in)     :: oprs
integer,                        intent(in)     :: igrid
real(kind=kind_real),           intent(in)     :: hf(grid%nface(igrid))
real(kind=kind_real),           intent(out)    :: f(grid%nface(igrid))
logical, optional,              intent(in)     :: removemean_in
real(kind=kind_real), optional, intent(inout)  :: field_conv(grid%niter)

integer :: nf, ne, nv, ipass
real(kind=kind_real) :: beta = 1.0_kind_real, fbar
real(kind=kind_real), allocatable, dimension(:) :: ff1, ff2, ff3, ff4, temp1, temp2
logical :: removemean
real(kind=kind_real), allocatable, dimension(:) :: hf_check

removemean = .false.
if (present(removemean_in)) removemean = removemean_in

nf = grid%nface(grid%ngrids)
ne = grid%nedge(grid%ngrids)
nv = grid%nvert(grid%ngrids)

allocate(ff1  (nf))
allocate(ff2  (nf))
allocate(ff3  (nf))
allocate(ff4  (nf))
allocate(temp1(ne))
allocate(temp2(ne))

f = 0.0_kind_real
temp2 = 0.0_kind_real

if (grid%check_convergence) allocate(hf_check(nf))

! Iterate several passes
do ipass = 1, grid%niter

  ! Print
  if (grid%rank == 0) print*, "   FEMPS inverse Laplacian iteration number for root processor: ", ipass

  ! Compute residual based on latest estimate
  call massL(oprs,f,ff2,grid%ngrids,nf)

  call Ddual1(grid,ff2,temp1,grid%ngrids,nf,ne)

  ! Improve the estimate temp2 that we obtained in the previous pass
  call massMinv(grid,oprs,temp1,temp2,grid%ngrids,ne,4)

  call Dprimal2(grid,oprs,temp2,ff3,grid%ngrids,ne,nf)

  ! Residual
  ff4 = hf - ff3

  ! Now solve the approximate Poisson equation
  ! D2 xminv D1bar L psi' = residual
  call mgsolve(grid,oprs,ff3,ff4,grid%ngrids)

  ! And increment best estimate
  ! *** We could think about adding beta*ff3 here and tuning beta for optimal convergence ***
  f = f + beta*ff3

  if (removemean) then
    ! Remove global mean (to ensure unambiguous result)
    fbar = SUM(f)/SUM(grid%farea(:,grid%ngrids))
    f = f - fbar*grid%farea(:,grid%ngrids)
  endif

  ! Optionally track the convergence
  if (grid%check_convergence) then
    if (.not. present(field_conv)) call message("Checking conv provide output argument", fatal)
    call laplace(grid,oprs,igrid,f,hf_check)
    field_conv(ipass) = sqrt(sum((hf_check/grid%farea(:,grid%ngrids)-hf/grid%farea(:,grid%ngrids))**2)/nf)
  endif

enddo

deallocate(ff1,ff2,ff3,ff4,temp1,temp2)

end subroutine inverselaplace

! --------------------------------------------------------------------------------------------------

subroutine xlaplace(grid,oprs,igrid,f,hf)

! To apply the APPROXIMATE Laplacian operator to the input field f,
! on grid igrid, the result appearing in the output field hf.
! Note f and hf are area integrals (2-forms).

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid
real(kind=kind_real), intent(in)  :: f(grid%nface(igrid))
real(kind=kind_real), intent(out) :: hf(grid%nface(igrid))

real(kind=kind_real) :: temp1(grid%nface(igrid)), temp2(grid%nedge(igrid)), temp3(grid%nedge(igrid))

call massL(oprs,f,temp1,igrid,grid%nface(igrid))
call Ddual1(grid,temp1,temp2,igrid,grid%nface(igrid),grid%nedge(igrid))
call approxMinv(oprs,temp2,temp3,igrid,grid%nedge(igrid))
call Dprimal2(grid,oprs,temp3,hf,igrid,grid%nedge(igrid),grid%nface(igrid))

end subroutine xlaplace

! --------------------------------------------------------------------------------------------------

subroutine residual(grid,oprs,f,rhs,res,igrid,nf)

! Compute the residual res in the approximate Poisson equation on grid igrid
! when f is the input field and rhs is the right hand side. Note that
! f, rhs and res are area integrals (2-forms).

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: igrid, nf
real(kind=kind_real), intent(in)  :: f(nf), rhs(nf)
real(kind=kind_real), intent(out) :: res(nf)

call xlaplace(grid,oprs,igrid,f,res)
res = rhs - res

end subroutine residual

! --------------------------------------------------------------------------------------------------

subroutine relax(grid,oprs,f,rhs,igrid,nf,niter)

! To carry out niter Jacobi relaxation iterations for the multigrid
! solver on grid igrid

implicit none
type(fempsgrid),      intent(in)    :: grid
type(fempsoprs),      intent(in)    :: oprs
integer,              intent(in)    :: igrid, nf, niter
real(kind=kind_real), intent(in)    :: rhs(nf)
real(kind=kind_real), intent(inout) :: f(nf)

real(kind=kind_real), allocatable :: res(:), finc(:)
real(kind=kind_real) :: u
integer :: iter

allocate(res(nf), finc(nf))

u = oprs%underrel(igrid)
do iter = 1, niter
  call residual(grid,oprs,f,rhs,res,igrid,nf)
  finc = res/oprs%lapdiag(1:nf,igrid)
  f = f + u*finc
enddo

deallocate(res, finc)

end subroutine relax

! --------------------------------------------------------------------------------------------------

subroutine mgsolve(grid,oprs,phi,rr,ng)

! Multigrid solver for elliptic equation
!
! Dprimal2 xminv Ddual1 L phi = RHS
!
! using a single V-cycle multigrid algorithm.
! Coefficients are contained in module laplace.

implicit none
type(fempsgrid),      intent(in)  :: grid
type(fempsoprs),      intent(in)  :: oprs
integer,              intent(in)  :: ng
real(kind=kind_real), intent(in)  :: rr(grid%nfacex)
real(kind=kind_real), intent(out) :: phi(grid%nfacex)

! Numbers of iterations on coarsest grid and other grids
integer, parameter :: niterc = 10, niter = 2
integer :: nf1, nf2, ne1, ne2, jgrid, jgridp
real(kind=kind_real), allocatable :: ff(:,:), rf(:,:)
real(kind=kind_real) :: temp1(grid%nfacex)

! Allocate space on all grids
allocate(ff(grid%nfacex,grid%ngrids),rf(grid%nfacex,grid%ngrids))

! Initialize solution to zero
phi = 0.0_kind_real

! Initialize rhs on finest grid
rf(:,grid%ngrids) = rr(:)

! Initialize correction to solution on all grids to zero
ff = 0.0_kind_real

! Descending part of V-cycle
do jgrid = grid%ngrids-1, grid%ngrids-ng+1, -1

  jgridp = jgrid + 1
  nf1 = grid%nface(jgridp)
  ne1 = grid%nedge(jgridp)
  nf2 = grid%nface(jgrid)
  ne2 = grid%nedge(jgrid)

  ! Relax on grid jgridp
  call relax(grid,oprs,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,niter)

  ! Calculate residual on jgridp
  call residual(grid,oprs,ff(1,jgridp),rf(1,jgridp),temp1,jgridp,nf1)

  ! Restrict residual to jgrid
  call restrict(grid,oprs,temp1,nf1,rf(1,jgrid),nf2,jgrid)

  ! Set correction first guess to zero on grid jgrid
  ff(1:nf2,jgrid) = 0.0_kind_real

enddo

! Iterate to convergence on coarsest grid
jgrid = grid%ngrids-ng+1
nf1 = grid%nface(jgrid)
ne1 = grid%nedge(jgrid)
ff(1:nf1,jgrid) = 0.0_kind_real
call relax(grid,oprs,ff(1,jgrid),rf(1,jgrid),jgrid,nf1,niterc)

! Ascending part of V-cycle
do jgrid = grid%ngrids-ng+1, grid%ngrids-1

  jgridp = jgrid + 1
  nf1 = grid%nface(jgridp)
  ne1 = grid%nedge(jgridp)
  nf2 = grid%nface(jgrid)
  ne2 = grid%nedge(jgrid)

  ! Prolong correction to jgridp
  call prolong(grid,oprs,ff(1,jgrid),nf2,temp1,nf1,jgrid)

  ! Add correction to solution on jgridp
  ff(1:nf1,jgridp) = ff(1:nf1,jgridp) + temp1(1:nf1)

  ! Relax on grid jgridp
  call relax(grid,oprs,ff(1,jgridp),rf(1,jgridp),jgridp,nf1,niter)

enddo

! Add correction to phi
phi = phi + ff(:,grid%ngrids)

!  ! For diagnostics
!  nf1 = grid%nface(grid%ngrids)
!  ne1 = grid%nedge(grid%ngrids)
!  call residual(grid,oprs,phi,rr,temp1,grid%ngrids,nf1)

deallocate(ff,rf)

end subroutine mgsolve

! --------------------------------------------------------------------------------------------------

end module femps_solve_mod
