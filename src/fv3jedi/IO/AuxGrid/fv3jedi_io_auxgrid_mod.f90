! (C) Copyright 2022 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_auxgrid_mod

! libs
use netcdf
use mpi

! atlas
use atlas_module,               only: atlas_field, atlas_real, atlas_functionspace, atlas_functionspace_pointcloud

! fckit
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm

! oops
use datetime_mod
use string_utils,               only: swap_name_member
use interpolatorbump_mod,       only: bump_interpolator

! fv3jedi
use fv3jedi_field_mod,          only: fv3jedi_field
use fv3jedi_geom_mod,           only: fv3jedi_geom
use fv3jedi_kinds_mod,          only: kind_int, kind_real
use fv3jedi_io_utils_mod,       only: add_iteration
use fv3jedi_netcdf_utils_mod,   only: nccheck

implicit none
private
public fv3jedi_io_auxgrid

type fv3jedi_io_auxgrid
 type(fckit_mpi_comm) :: comm
 type(bump_interpolator) :: bumpinterp
 integer :: llcomm
 integer :: layout(2)
 logical :: thispe = .false.
 integer :: nx, ny, nxg, nyg
 integer :: npes
 integer :: isc, iec, jsc, jec, npx, npy, npz
 integer :: isd, ied, jsd, jed
 real(kind=kind_real), allocatable :: lons(:)
 real(kind=kind_real), allocatable :: lats(:)
 character(len=24) :: gridtype
 character(len=1024) :: filename
 integer, allocatable :: istart2(:), icount2(:)
 integer, allocatable :: istart3(:), icount3(:)
 integer, allocatable :: istarte(:), icounte(:)
 contains
  procedure :: create
  procedure :: delete
  procedure :: write
  final :: dummy_final
end type fv3jedi_io_auxgrid

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

! Arguments
class(fv3jedi_io_auxgrid),  intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

! Locals
character(len=:), allocatable :: str

! Local copy of geometry things
self%isc  = geom%isc
self%iec  = geom%iec
self%jsc  = geom%jsc
self%jec  = geom%jec
self%isd  = geom%isd
self%ied  = geom%ied
self%jsd  = geom%jsd
self%jed  = geom%jed
self%npx  = geom%npx
self%npy  = geom%npy
self%npz  = geom%npz
self%comm = geom%f_comm

! Set output grid type
self%gridtype = "latlon"
if ( conf%has("gridtype") ) then
    call conf%get_or_die("gridtype", str)
    self%gridtype=str
endif

! Create auxgrid
call create_auxgrid(self, geom)

! Naming convention for the file
! For ensemble methods switch out member template
self%filename = 'Data/fv3jedi.auxgrid.'
call conf%get_or_die("filename",str)
call swap_name_member(conf, str)
call add_iteration(conf, str)
self%filename = str
deallocate(str)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

! Arguments
class(fv3jedi_io_auxgrid), intent(inout) :: self

call delete_auxgrid(self)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine write(self, vdate, fields)

! Arguments
class(fv3jedi_io_auxgrid),  intent(inout) :: self
type(datetime),            intent(in)    :: vdate
type(fv3jedi_field),       intent(in)    :: fields(:)

! Write metadata
! --------------
call write_auxgrid_metadata(self, vdate)

! Write fields
! ------------
call write_auxgrid_fields(self, fields)

end subroutine write

! --------------------------------------------------------------------------------------------------

subroutine create_auxgrid(self, geom)

use fckit_log_module, only: fckit_log

! Arguments
type(fv3jedi_io_auxgrid), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom

! Locals
integer :: color, ix, iy, inode
real(kind=kind_real) :: dx, dy, rad2deg
real(kind=kind_real),allocatable :: slat(:), wlat(:)
real(kind=kind_real),pointer :: lonlat_ptr(:,:)
integer :: i, ierr, jmax
character(len=48) :: debug_msg
type(atlas_field) :: lonlat_field
type(atlas_functionspace) :: afunctionspace_out

! Create output grid and interpolation object for going from cube to output grid
! -------------------------------------------------------------------------------
!Maximum of 12 IO processors
if (self%comm%size() >= 12) then
  self%layout(1) = 12
  self%layout(2) = 1
elseif (self%comm%size() >= 6) then
  self%layout(1) = 6
  self%layout(2) = 1
else
  call abor1_ftn("create_auxgrid error: fewer than 6 npes not anticipated")
endif
self%npes = self%layout(1) * self%layout(2)

!Since the lat lon grid is only for IO purposes it is only
!defined on a subset of PEs - those that will do the writing.
!This is generally more efficient than having many PEs trying
!to write to the same file.

self%nxg = 4*(self%npx - 1)
if ( trim(self%gridtype) == 'latlon' ) then
  self%nyg = 2*(self%npy - 1) + 1
else
  self%nyg = 2*(self%npy - 1)
endif
self%nx = 0
self%ny = 0

color = MPI_UNDEFINED

if (self%comm%rank() <= self%npes-1) then

  !Split communicator
  color = int(self%comm%communicator() / 2)

  self%thispe = .true.

  self%nx = self%nxg / self%layout(1)
  self%ny = self%nyg / self%layout(2)

  allocate(self%lons(self%nx))
  allocate(self%lats(self%ny))

  write(debug_msg,*)'create ouptut grid for ',trim(self%gridtype)
  call fckit_log%debug(debug_msg)

  ! Populate lons and lats arrays for requested auxgrid type
  ! Each processor has subset of lons and all lats.

  dx = 360.0_kind_real / real(self%nxg,kind_real)

  self%lons(1) = dx * self%nx * self%comm%rank()
  do i = 2,self%nx
     self%lons(i) = self%lons(i-1) + dx
  enddo

  if ( trim(self%gridtype) == 'latlon' ) then
     dy = 180.0_kind_real / real(self%nyg-1,kind_real)

     self%lats(1) = -90.0_kind_real
     do i = 2,self%ny
        self%lats(i) = (i-1)*dy + self%lats(1)
     enddo
     self%lats(self%ny) = 90.0_kind_real

  elseif ( trim(self%gridtype) == 'gaussian' ) then
     allocate(slat(self%ny))
     allocate(wlat(self%ny))

     jmax=self%ny
     rad2deg = 180.0_kind_real / (4.0_kind_real * atan(1.0_kind_real))
     call splat(4, jmax, slat, wlat)
     do i = 1,self%ny
        self%lats(i) = -asin(slat(i)) * rad2deg
     enddo

     deallocate(slat)
     deallocate(wlat)
  else
     call abor1_ftn("create_auxgrid error: invalid output gridtype")
  endif

else

  self%nx = 0
  self%ny = 0
  allocate(self%lons(self%nx))
  allocate(self%lats(self%ny))

endif

! Output function space
lonlat_field = atlas_field(name='lonlat',kind=atlas_real(kind_real),shape=(/2,self%nx*self%ny/))
call lonlat_field%data(lonlat_ptr)
inode = 0
do iy=1,self%ny
   do ix=1,self%nx
      inode = inode+1
      lonlat_ptr(1,inode) = self%lons(ix)
      lonlat_ptr(2,inode) = self%lats(iy)
   end do
end do
afunctionspace_out = atlas_functionspace_pointcloud(lonlat_field)

! Initialize bump interpolator
call self%bumpinterp%init(self%comm, geom%afunctionspace, afunctionspace_out, self%npz)

!IO communicator
!call MPI_Comm_split(self%comm%communicator(), color, self%comm%rank(), self%llcomm, ierr)

! NC arrays
allocate(self%istart3(4),self%icount3(4))
allocate(self%istarte(4),self%icounte(4))
allocate(self%istart2(3),self%icount2(3))
self%istart3(1) = self%nx * self%comm%rank() + 1;  self%icount3(1) = self%nx
self%istart3(2) = 1;                                     self%icount3(2) = self%ny
self%istart3(3) = 1;                                     self%icount3(3) = self%npz
self%istart3(4) = 1;                                     self%icount3(4) = 1
self%istart2(1) = self%istart3(1);                     self%icount2(1) = self%icount3(1)
self%istart2(2) = self%istart3(2);                     self%icount2(2) = self%icount3(2)
self%istart2(3) = self%istart3(4);                     self%icount2(3) = self%icount3(4)

self%istarte = self%istarte
self%icounte = self%icount3
self%istarte(3) = self%npz + 1

end subroutine create_auxgrid

! --------------------------------------------------------------------------------------------------

subroutine delete_auxgrid(self)

!Arguments
type(fv3jedi_io_auxgrid), intent(inout) :: self  !< Auxgrid Geometry

if (self%thispe) then
  deallocate(self%lons)
  deallocate(self%lats)

  deallocate ( self%istart2, self%icount2 )
  deallocate ( self%istart3, self%icount3 )
  deallocate ( self%istarte, self%icounte )
endif

call self%bumpinterp%delete()

end subroutine delete_auxgrid

! --------------------------------------------------------------------------------------------------

subroutine write_auxgrid_metadata(self, vdate)

!Arguments
type(fv3jedi_io_auxgrid), intent(inout) :: self
type(datetime),          intent(in)    :: vdate

integer :: date(6)
integer :: idate, isecs
character(len=64)   :: datefile

integer :: ncid, varid(2)
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid

! Current date
call datetime_to_ifs(vdate, idate, isecs)
date(1) = idate/10000
date(2) = idate/100 - date(1)*100
date(3) = idate - (date(1)*10000 + date(2)*100)
date(4) = isecs/3600
date(5) = (isecs - date(4)*3600)/60
date(6) = isecs - (date(4)*3600 + date(5)*60)

! Append filename with the curent datetime from fv3jedi
write(datefile,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') date(1),date(2),date(3),"_",date(4),date(5),date(6)
self%filename = trim(self%filename)//trim(datefile)//trim("z.nc4")

call nccheck( nf90_create( self%filename, ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                           comm = self%comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

! Create dimensions
call nccheck ( nf90_def_dim(ncid, "lon" , self%nxg , x_dimid), "nf90_def_dim lon" )
call nccheck ( nf90_def_dim(ncid, "lat" , self%nyg , y_dimid), "nf90_def_dim lat" )
call nccheck ( nf90_def_dim(ncid, "lev" , self%npz   , z_dimid), "nf90_def_dim lev" )
call nccheck ( nf90_def_dim(ncid, "edge", self%npz+1 , e_dimid), "nf90_def_dim edge")
call nccheck ( nf90_def_dim(ncid, "time", 1          , t_dimid), "nf90_def_dim time")

! Define fields to be written
call nccheck( nf90_def_var(ncid, "lon", NF90_DOUBLE, x_dimid, varid(1)), "nf90_def_var lon" )
call nccheck( nf90_put_att(ncid, varid(1), "long_name", "longitude") )
call nccheck( nf90_put_att(ncid, varid(1), "units", "degrees_east") )

call nccheck( nf90_def_var(ncid, "lat", NF90_DOUBLE, y_dimid, varid(2)), "nf90_def_var lat" )
call nccheck( nf90_put_att(ncid, varid(2), "long_name", "latitude") )
call nccheck( nf90_put_att(ncid, varid(2), "units", "degrees_north") )

call nccheck( nf90_enddef(ncid), "nf90_enddef" )

if (self%thispe) then

  ! Write fields
  call nccheck( nf90_put_var( ncid, varid(1), self%lons, &
                              start = self%istart2(1:1), count = self%icount2(1:1) ), "nf90_put_var lons" )
  call nccheck( nf90_put_var( ncid, varid(2), self%lats, &
                              start = self%istart2(2:2), count = self%icount2(2:2) ), "nf90_put_var lats" )

endif

! Close file
call nccheck( nf90_close(ncid), "nf90_close" )

!Let Auxgrid PEs catch up
call self%comm%barrier()

end subroutine write_auxgrid_metadata

! --------------------------------------------------------------------------------------------------

subroutine write_auxgrid_fields(self, fields)

!Arguments
type(fv3jedi_io_auxgrid), target, intent(inout)   :: self
type(fv3jedi_field),             intent(in)      :: fields(:)

integer :: var, ji, jj, jk, llngrid, ii, i, j, k, n
real(kind=kind_real), allocatable :: llfield(:,:,:)

integer :: ncid, varid
integer :: x_dimid, y_dimid, z_dimid, e_dimid, t_dimid
integer, target  :: dimids3(4), dimids2(3), dimidse(4)
integer, pointer :: dimids(:), istart(:), icount(:)
real(kind_real),allocatable :: array_without_halo(:,:,:)

! Loop over fields
! ----------------
do var = 1, size(fields)

  ! Only certain fields can be written for now
  ! ------------------------------------------
  if ( trim(fields(var)%space) == 'magnitude' .and. &
       trim(fields(var)%horizontal_stagger_location) == 'center' .and. &
       .not. fields(var)%kind == 'integer' ) then

    ! Interpolate the field to the output grid grid
    ! ---------------------------------------------
    allocate(llfield(1:self%nx,1:self%ny,1:self%npz))
    llfield = 0.0_kind_real
    llngrid = self%nx*self%ny

    ! Copy field in array without halo points
    allocate(array_without_halo(self%isc:self%iec,self%jsc:self%jec,1:fields(var)%npz))
    array_without_halo = fields(var)%array(self%isc:self%iec,self%jsc:self%jec,1:fields(var)%npz)

    ! Interpolate
    call self%bumpinterp%apply(array_without_halo, llfield(:,:,1:fields(var)%npz))

    !Open file
    call nccheck( nf90_open( self%filename, ior(NF90_WRITE, NF90_MPIIO), ncid, &
                            comm = self%comm%communicator(), info = MPI_INFO_NULL), "nf90_open" )

    !Dimension ID
    call nccheck( nf90_inq_dimid(ncid, "lon" , x_dimid), "nf90_inq_dimid lon" )
    call nccheck( nf90_inq_dimid(ncid, "lat" , y_dimid), "nf90_inq_dimid lat" )
    call nccheck( nf90_inq_dimid(ncid, "lev" , z_dimid), "nf90_inq_dimid lev" )
    call nccheck( nf90_inq_dimid(ncid, "edge", e_dimid), "nf90_inq_dimid edge" )
    call nccheck( nf90_inq_dimid(ncid, "time", t_dimid), "nf90_inq_dimid time" )

    dimids3 =  (/ x_dimid, y_dimid, z_dimid, t_dimid /)
    dimidse =  (/ x_dimid, y_dimid, e_dimid, t_dimid /)
    dimids2 =  (/ x_dimid, y_dimid, t_dimid /)

    ! Pointer to dimensions based on number of levels
    if (associated(dimids)) nullify(dimids)
    if (associated(istart)) nullify(istart)
    if (associated(icount)) nullify(icount)
    if (fields(var)%npz == self%npz) then
      dimids => dimids3
      istart => self%istart3
      icount => self%icount3
    elseif (fields(var)%npz == 1) then
      dimids => dimids2
      istart => self%istart2
      icount => self%icount2
    elseif (fields(var)%npz == self%npz+1) then
      dimids => dimidse
      istart => self%istarte
      icount => self%icounte
    endif

    if (associated(dimids)) then

      ! Write field to the file
      call nccheck( nf90_def_var( ncid, trim(fields(var)%io_name), NF90_DOUBLE, dimids, varid), &
                    "nf90_def_var"//trim(fields(var)%io_name) )
      call nccheck( nf90_put_att(ncid, varid, "long_name", trim(fields(var)%long_name) ), "nf90_put_att" )
      call nccheck( nf90_put_att(ncid, varid, "units"    , trim(fields(var)%units)     ), "nf90_put_att" )
      call nccheck( nf90_enddef(ncid), "nf90_enddef" )

      if (self%thispe) then
       call nccheck( nf90_put_var( ncid, varid, llfield, start = istart, count = icount), &
                     "nf90_put_var"//trim(fields(var)%io_name) )
      endif

    endif

    ! Close file
    call nccheck( nf90_close(ncid), "nf90_close" )

    !Let Auxgrid PEs catch up
    call self%comm%barrier()

    deallocate(llfield)
    deallocate(array_without_halo)
  endif

enddo

end subroutine write_auxgrid_fields

! --------------------------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_io_auxgrid), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

subroutine splat(idrt,jmax,slat,wlat)
! AUTHOR:
!   IREDELL (NOAA/NWS/EMC) date 96-02-20
!
! Computes cosines of colatitude and gaussian weights
! for one of the following specific global sets of latitudes.
! - Gaussian latitudes (IDRT=4)
! - Equally-spaced latitudes including poles (IDRT=0)
! - Equally-spaced latitudes excluding poles (IDRT=256)
! The gaussian latitudes are located at the zeroes of the
! legendre polynomial of the given order.  these latitudes
! are efficient for reversible transforms from spectral space.
! (About twice as many equally-spaced latitudes are needed.)
! The weights for the equally-spaced latitudes are based on
! Ellsaesser (JAM,1966). (No weight is given the pole point.)
! Note that when analyzing grid to spectral in latitude pairs,
! if an equator point exists, its weight should be halved.
! This version invokes the ibm essl matrix solver.
!
! PROGRAM HISTORY LOG:
! - 1996-02-20  IREDELL
! - 1997-10-20  IREDELL  ADJUST PRECISION
! - 1998-06-11  IREDELL  GENERALIZE PRECISION USING FORTRAN 90 INTRINSIC
! - 1998-12-03  IREDELL  GENERALIZE PRECISION FURTHER
! - 1998-12-03  IREDELL  USES AIX ESSL BLAS CALLS
! - 2009-12-27  DSTARK   updated to switch between ESSL calls on an AIX
!                        platform, and Numerical Recipies calls elsewise.
! - 2010-12-30  SLOVACEK update alignment so preprocessor does not cause
!                        compilation failure
! - 2012-09-01  E.Mirvis & M.Iredell merging & debugging linux errors 
!                        of _d and _8 using generic LU factorization.   
! - 2012-11-05  E.Mirvis generic FFTPACK and LU lapack were removed
!
! ARGUMENTS:
!   IDRT       - INTEGER GRID IDENTIFIER
!                (IDRT=4 FOR GAUSSIAN GRID,
!                 IDRT=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
!                 IDRT=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES)
!   JMAX       - INTEGER NUMBER OF LATITUDES.
!   SLAT (out) - REAL (JMAX) SINES OF LATITUDE.
!   WLAT (out) - REAL (JMAX) GAUSSIAN WEIGHTS.
!
  integer:: idrt,jmax
  REAL(kind=kind_real) SLAT(JMAX),WLAT(JMAX)
  integer:: jh,j,n
  real(kind=kind_real):: r,c
  INTEGER,PARAMETER:: KD=SELECTED_REAL_KIND(15,45)
  REAL(KIND=kind_real):: PK(JMAX/2),PKM1(JMAX/2),PKM2(JMAX/2)
  REAL(KIND=kind_real):: SLATD(JMAX/2),SP,SPMAX,EPS=10.*EPSILON(SP)
  integer,PARAMETER:: JZ=50
  REAL(kind=kind_real) BZ(JZ)
  DATA BZ        / 2.4048255577,  5.5200781103, &
       8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679, &
       21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684, &
       33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132, &
       46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550, &
       58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299, &
       71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
       84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, &
       96.6052679510, 99.7468198587, 102.888374254, 106.029930916, &
       109.171489649, 112.313050280, 115.454612653, 118.596176630, &
       121.737742088, 124.879308913, 128.020877005, 131.162446275, &
       134.304016638, 137.445588020, 140.587160352, 143.728733573, &
       146.870307625, 150.011882457, 153.153458019, 156.295034268 /
  REAL(kind=kind_real):: DLT,D1=1.
  REAL(kind=kind_real) AWORK((JMAX+1)/2,((JMAX+1)/2)),BWORK(((JMAX+1)/2))
  INTEGER:: JHE,JHO,J0=0
  INTEGER IPVT((JMAX+1)/2)
  real(kind=kind_real),PARAMETER :: PI=3.14159265358979

  C=(1.-(2./PI)**2)*0.25

!  GAUSSIAN LATITUDES
IF(IDRT.EQ.4) THEN
   JH=JMAX/2
   JHE=(JMAX+1)/2
   R=1./SQRT((JMAX+0.5)**2+C)
   DO J=1,MIN(JH,JZ)
      SLATD(J)=COS(BZ(J)*R)
   ENDDO
   DO J=JZ+1,JH
      SLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
   ENDDO
   SPMAX=1.
   DO WHILE(SPMAX.GT.EPS)
      SPMAX=0.
      DO J=1,JH
         PKM1(J)=1.
         PK(J)=SLATD(J)
      ENDDO
      DO N=2,JMAX
         DO J=1,JH
            PKM2(J)=PKM1(J)
            PKM1(J)=PK(J)
            PK(J)=((2*N-1)*SLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
         ENDDO
      ENDDO
      DO J=1,JH
         SP=PK(J)*(1.-SLATD(J)**2)/(JMAX*(PKM1(J)-SLATD(J)*PK(J)))
         SLATD(J)=SLATD(J)-SP
         SPMAX=MAX(SPMAX,ABS(SP))
      ENDDO
   ENDDO
!!DIR$ IVDEP
   DO J=1,JH
      SLAT(J)=SLATD(J)
      WLAT(J)=(2.*(1.-SLATD(J)**2))/(JMAX*PKM1(J))**2
      SLAT(JMAX+1-J)=-SLAT(J)
      WLAT(JMAX+1-J)=WLAT(J)
   ENDDO
   IF(JHE.GT.JH) THEN
      SLAT(JHE)=0.
      WLAT(JHE)=2./JMAX**2
      DO N=2,JMAX,2
         WLAT(JHE)=WLAT(JHE)*N**2/(N-1)**2
      ENDDO
   ENDIF
endif
end subroutine splat

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_auxgrid_mod
