! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_utils_mod

use femps_kinds_mod

implicit none

public

interface ll2xyz
   module procedure ll2xyz_sca
   module procedure ll2xyz_vec
end interface

interface xyz2ll
   module procedure xyz2ll_sca
   module procedure xyz2ll_vec
end interface

integer, parameter :: fatal = 1, trace = 2

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine ll2xyz_sca(long,lat,x,y,z)

! To convert longitude and latitude to cartesian coordinates
! on the unit sphere

implicit none
real(kind=kind_real), intent(in)  :: long,lat
real(kind=kind_real), intent(out) :: x, y, z

real(kind=kind_real) :: cln,sln,clt,slt

sln=sin(long)
cln=cos(long)
slt=sin(lat)
clt=cos(lat)

x=cln*clt
y=sln*clt
z=slt

end subroutine ll2xyz_sca

! --------------------------------------------------------------------------------------------------

subroutine ll2xyz_vec(long,lat,x)

implicit none
real(kind=kind_real), intent(in)  :: long,lat
real(kind=kind_real), intent(out) :: x(3)

call ll2xyz_sca(long,lat,x(1),x(2),x(3))

end subroutine ll2xyz_vec

! --------------------------------------------------------------------------------------------------

subroutine xyz2ll_sca(x,y,z,long,lat)

implicit none
real(kind=kind_real), intent(in)  :: x,y,z
real(kind=kind_real), intent(out) :: long,lat

real(kind=kind_real) :: pi,tln,tlt,r

! To convert cartesian coordinates to longitude and latitude
! ----------------------------------------------------------

pi=4.0_kind_real*atan(1.0_kind_real)

if (x.eq.0.0_kind_real) then
  if (y.ge.0.0_kind_real) then
    long=0.5_kind_real*pi
  else
    long=1.5_kind_real*pi
  endif
else
  tln=y/x
  long=atan(tln)
  if (x.lt.0.0_kind_real) then
    long=long+pi
  endif
  if (long.lt.0.0_kind_real) then
    long=long+2.0_kind_real*pi
  endif
endif

r=sqrt(x*x+y*y)
if (r.eq.0.0_kind_real) then
  if (z.gt.0.0_kind_real) then
    lat=0.5_kind_real*pi
  else
    lat=-0.5_kind_real*pi
  endif
else
  tlt=z/r
  lat=atan(tlt)
endif

end subroutine xyz2ll_sca

! --------------------------------------------------------------------------------------------------

subroutine xyz2ll_vec(x,long,lat)

implicit none
real(kind=kind_real), intent(in)  :: x(3)
real(kind=kind_real), intent(out) :: long,lat

call xyz2ll_sca(x(1),x(2),x(3),long,lat)

end subroutine xyz2ll_vec

! --------------------------------------------------------------------------------------------------

subroutine starea2(x0,y0,z0,x1,y1,z1,x2,y2,z2,area)

implicit none
real(kind=kind_real), intent(in)  :: x0, x1, x2
real(kind=kind_real), intent(in)  :: y0, y1, y2
real(kind=kind_real), intent(in)  :: z0, z1, z2
real(kind=kind_real), intent(out) :: area

real(kind=kind_real) :: d0,d1,d2,s,t0,t1,t2,t3

! Calculate the area of the spherical triangle whose corners
! have Cartesian coordinates (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
! The formula below is more robust to roundoff error than the
! better known sum of angle - PI formula

! Distances between pairs of points
call spdist(x0,y0,z0,x1,y1,z1,d2)
call spdist(x1,y1,z1,x2,y2,z2,d0)
call spdist(x2,y2,z2,x0,y0,z0,d1)

! Half perimeter
s=0.5_kind_real*(d0+d1+d2)

! Tangents
t0 = tan(0.5_kind_real*(s-d0))
t1 = tan(0.5_kind_real*(s-d1))
t2 = tan(0.5_kind_real*(s-d2))
t3 = tan(0.5_kind_real*s)

! Area
area = 4.0_kind_real*atan(sqrt(t0*t1*t2*t3))

end subroutine starea2

! --------------------------------------------------------------------------------------------------

subroutine spdist(x1,y1,z1,x2,y2,z2,s)

implicit none
real(kind=kind_real), intent(in)  :: x1, x2
real(kind=kind_real), intent(in)  :: y1, y2
real(kind=kind_real), intent(in)  :: z1, z2
real(kind=kind_real), intent(out) :: s

real(kind=kind_real) :: dx, dy, dz, ad

! Calculate the spherical distance S between two points with Cartesian
! coordinates (X1,Y1,Z1), (X2,Y2,Z2) on the unit sphere

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1
ad = sqrt(dx*dx + dy*dy + dz*dz)
s = 2.0_kind_real*asin(0.5_kind_real*ad)

end subroutine spdist

! --------------------------------------------------------------------------------------------------

subroutine addtab(dim1,dim2,tab,index,tabentry)

implicit none
integer, intent(in)    :: dim1,dim2
integer, intent(inout) :: tab(dim1,dim2)
integer, intent(in)    :: index
integer, intent(in)    :: tabentry

integer :: i
character(len=2056) :: errormessage

! Add an entry to a table
! -----------------------

i=0
100 continue
i=i+1
if (i.gt.dim2) then
  call message('**********')
  call message('table full')
  call message('**********')
  write(errormessage,*) index,tabentry,dim1,dim2
  call message(errormessage)
  write(errormessage,*) tab(index,:)
  call message(errormessage,fatal)
endif
if (tab(index,i).ne.0) goto 100
tab(index,i)=tabentry

end subroutine addtab

! --------------------------------------------------------------------------------------------------

subroutine findinlist(i,list,nl,ix)

implicit none
integer, intent(in)    :: i
integer, intent(in)    :: nl
integer, intent(inout) :: list(nl)
integer, intent(out)   :: ix

character(len=255) :: errormessage

! Find the index ix to the entry i in the input list.
! if i is not already in the list then it is added into
! the first empty (i.e. zero) position.

ix = 1
do while ((list(ix) .ne. i) .and. (list(ix) .ne. 0))
  if (ix == nl) then
    call message('search past end of list in subroutine findinlist.')
    write(errormessage,*) 'i = ',i
    call message(errormessage)
    write(errormessage,*) 'list = ',list
    call message(errormessage,fatal)
  endif
  ix = ix + 1
enddo
list(ix) = i

end subroutine findinlist

! --------------------------------------------------------------------------------------------------

subroutine triangle(x1,x2,x3,l1sq,l2sq,l3sq,area)

implicit none
real(kind=kind_real), intent(in)  :: x1(3), x2(3), x3(3)
real(kind=kind_real), intent(out) :: l1sq, l2sq, l3sq, area

real(kind=kind_real) :: a, b, c

! Compute the squares of the lengths of the sides and the area of
! a triangle defined by three points x1, x2 and x3 in 3d euclidean
! space. Side 1 is opposite vertex 1, etc.

! squares of side by pythagoras
l1sq = (x2(1) - x3(1))**2 + (x2(2) - x3(2))**2 + (x2(3) - x3(3))**2
l2sq = (x3(1) - x1(1))**2 + (x3(2) - x1(2))**2 + (x3(3) - x1(3))**2
l3sq = (x1(1) - x2(1))**2 + (x1(2) - x2(2))**2 + (x1(3) - x2(3))**2

! area by heron's formula
a = sqrt(l1sq)
b = sqrt(l2sq)
c = sqrt(l3sq)
area = 0.25_kind_real*sqrt((a + b + c)*(b + c - a)*(c + a - b)*(a + b - c))

end subroutine triangle

! --------------------------------------------------------------------------------------------------

subroutine nccheck(status,iam)

use netcdf

implicit none
integer,                    intent (in) :: status
character(len=*), optional, intent (in) :: iam

character(len=1024) :: error_descr

if(status /= nf90_noerr) then

  error_descr = "NetCDF error, aborting ... "

  if (present(iam)) then
    error_descr = trim(error_descr)//", "//trim(iam)
  endif

  error_descr = trim(error_descr)//". Error code: "//trim(nf90_strerror(status))

  call message('NCCheck'//trim(error_descr),fatal)

end if

end subroutine nccheck

! ------------------------------------------------------------------------------

subroutine message(note,code)

implicit none
character(len=*),  intent(in) :: note
integer, optional, intent(in) :: code

character(len=2056) :: noteprint

noteprint = note

if (present(code)) then
  if (code == fatal) then
    noteprint = 'FATAL: '//trim(note)
  elseif (code == trace) then
    !noteprint = 'TRACE: '//trim(note)
  endif
endif

!print*, trim(noteprint)

if (present(code)) then
  if (code == fatal) then
    print*, trim(noteprint)
    stop 1
  endif
endif

end subroutine message

! ------------------------------------------------------------------------------

end module femps_utils_mod
