! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_interpolation_mod

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_max, fckit_mpi_min

use atlas_module
use fv3jedi_bump_mod,      only: bump_init, bump_apply
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_field_mod,     only: fv3jedi_field, pointer_field, has_field
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_constants_mod, only: rad2deg
use wind_vt_mod,           only: d2a, a2d

use type_bump, only: bump_type
use unstructured_interpolation_mod, only: unstrc_interp

implicit none
private
public :: field2field_interp

! --------------------------------------------------------------------------------------------------

type field2field_interp

  character(len=32) :: interp_type
  type(bump_type) :: bump
  type(unstrc_interp) :: unsinterp
  logical :: need_bump, need_bary
  integer :: nnearest
  contains
    procedure :: create
    procedure :: delete
    procedure :: apply

end type field2field_interp

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, interp_type_in, integer_interp, geom_in, geom_ou)

class(field2field_interp), intent(inout) :: self           ! field2field_interp
character(len=*),          intent(in)    :: interp_type_in ! bump or barycent
logical,                   intent(in)    :: integer_interp ! Need to interpolate integers
type(fv3jedi_geom),        intent(in)    :: geom_in        !Geometry of input grid
type(fv3jedi_geom),        intent(in)    :: geom_ou        !Geometry of output grid

! Locals
integer :: ngrid_ou
character(len=32) :: us_interp_type

self%need_bump = .false.
self%need_bary = .false.

self%interp_type = trim(interp_type_in)
us_interp_type = ''

if (trim(self%interp_type) == 'bump') then
  self%need_bump = .true.
elseif (trim(self%interp_type) == 'barycent') then
  self%need_bary = .true.
  us_interp_type = trim(self%interp_type)
else
  call abor1_ftn("In fv3jedi_interpolation_mod.create: interp_type should be bump or barycent")
endif

if (integer_interp) then
  self%need_bary = .true.
  if (us_interp_type == '') us_interp_type = 'barycent'
endif

! Initialize bump object
! ----------------------
if (self%need_bump) call bump_init(geom_in, geom_ou%ngrid, rad2deg*geom_ou%lat_us, rad2deg*geom_ou%lon_us, self%bump)

! Initialize unstructured interpolation object
! --------------------------------------------
self%nnearest = 4
if (self%need_bary) then
  call self%unsinterp%create( geom_in%f_comm, self%nnearest, trim(us_interp_type), &
                              geom_in%ngrid, rad2deg*geom_in%lat_us, rad2deg*geom_in%lon_us, &
                              geom_ou%ngrid, rad2deg*geom_ou%lat_us, rad2deg*geom_ou%lon_us )
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(field2field_interp), intent(inout) :: self           ! field2field_interp

if (self%need_bump) call self%bump%dealloc()
if (self%need_bary) call self%unsinterp%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine apply(self, nf, geom_in, fields_in, geom_ou, fields_ou)

class(field2field_interp), intent(inout) :: self           !field2field_interp
integer,                   intent(in)    :: nf             !Number of fields
type(fv3jedi_geom),        intent(inout) :: geom_in        !Geometry of input grid
type(fv3jedi_field),       intent(in)    :: fields_in(nf)  !Input fields
type(fv3jedi_geom),        intent(inout) :: geom_ou        !Geometry of output grid
type(fv3jedi_field),       intent(inout) :: fields_ou(nf)  !Output fields

!Locals
integer :: var, i, j, k, n
real(kind=kind_real), allocatable :: field_in(:), field_ou(:), field_ou_2d(:,:)

logical :: do_d2a
real(kind=kind_real), allocatable :: ua_in(:,:,:)
real(kind=kind_real), allocatable :: va_in(:,:,:)
real(kind=kind_real), allocatable :: ud_ou(:,:,:)
real(kind=kind_real), allocatable :: vd_ou(:,:,:)
type(fv3jedi_field), pointer :: ud_in
type(fv3jedi_field), pointer :: vd_in
type(fv3jedi_field), pointer :: ua_ou
type(fv3jedi_field), pointer :: va_ou

! Special case of D-grid winds
! ----------------------------
do_d2a = .false.
if (has_field(fields_in,'ud')) then

  do_d2a = .true.

  call pointer_field(fields_in,'ud',ud_in)
  call pointer_field(fields_in,'vd',vd_in)

  allocate(ua_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,geom_in%npz))
  allocate(va_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,geom_in%npz))

  call d2a(geom_in, ud_in%array, vd_in%array, ua_in, va_in)

  ! Overwrite field with A-Grid for doing interpolation
  ud_in%array(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,:) = ua_in
  vd_in%array(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,:) = va_in

  deallocate(ua_in, va_in)

endif


! Interpolate all fields
! ---------------------
do var = 1,nf

  if (.not. fields_in(var)%integerfield .and. trim(self%interp_type) == 'bump') then

    ! Allocation
    allocate(field_ou_2d(geom_ou%ngrid,1:fields_ou(var)%npz))

    ! Interpolate
    call bump_apply(fields_ou(var)%npz, geom_in, fields_in(var)%array, geom_ou%ngrid, field_ou_2d, self%bump)

    ! Back to structured
    n = 0
    do j = geom_ou%jsc,geom_ou%jec
      do i = geom_ou%isc,geom_ou%iec
          n = n + 1
          fields_ou(var)%array(i,j,1:fields_ou(var)%npz) = field_ou_2d(n,1:fields_ou(var)%npz)
      enddo
    enddo

    ! Release memory
    deallocate(field_ou_2d)

  else

    ! Allocation
    allocate(field_in(geom_in%ngrid))
    allocate(field_ou(geom_ou%ngrid))

    do k = 1, fields_ou(var)%npz

      ! To unstructured
      n = 0
      do j = geom_in%jsc,geom_in%jec
        do i = geom_in%isc,geom_in%iec
            n = n + 1
            field_in(n) = fields_in(var)%array(i,j,k)
        enddo
      enddo

      ! Interpolate
      if (.not. fields_in(var)%integerfield .and. trim(self%interp_type) == 'barycent') then

        call self%unsinterp%apply(field_in, field_ou)

      elseif (fields_in(var)%integerfield) then

        call iinterp_apply(self, geom_in%f_comm, geom_in%ngrid, field_in, geom_ou%ngrid, field_ou)

      endif

      ! Back to structured
      n = 0
      do j = geom_ou%jsc,geom_ou%jec
        do i = geom_ou%isc,geom_ou%iec
            n = n + 1
            fields_ou(var)%array(i,j,k) = field_ou(n)
        enddo
      enddo

    enddo

    ! Release memory
    deallocate(field_in)
    deallocate(field_ou)

  endif

enddo


! Back to D-Grid
! --------------
if (do_d2a) then

  call pointer_field(fields_ou,'ud',ua_ou)
  call pointer_field(fields_ou,'vd',va_ou)

  allocate(ud_ou(geom_ou%isc:geom_ou%iec  ,geom_ou%jsc:geom_ou%jec+1,geom_ou%npz))
  allocate(vd_ou(geom_ou%isc:geom_ou%iec+1,geom_ou%jsc:geom_ou%jec  ,geom_ou%npz))

  call d2a(geom_ou, ua_ou%array(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,:), &
                    va_ou%array(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,:), &
                    ud_ou, vd_ou)

  ! Overwrite field with A-Grid for doing interpolation
  ua_ou%array = ud_ou
  va_ou%array = vd_ou

  deallocate(ud_ou, vd_ou)

endif

end subroutine apply

! --------------------------------------------------------------------------------------------------

subroutine iinterp_apply(self, f_comm, ngrid_in, field_in, ngrid_ou, field_ou)

type(field2field_interp), intent(in)    :: self
type(fckit_mpi_comm),     intent(in)    :: f_comm
integer,                  intent(in)    :: ngrid_in           !Input grid dimensions
real(kind=kind_real),     intent(in)    :: field_in(ngrid_in) !Integer field in
integer,                  intent(in)    :: ngrid_ou           !Output grid dimensions
real(kind=kind_real),     intent(inout) :: field_ou(ngrid_ou) !Integer field out

! Locals
integer :: maxtypel, mintypel, maxtype, mintype
integer :: i, j, k, n, index
real(kind=kind_real), allocatable :: interp_w(:,:)
real(kind=kind_real), allocatable :: field_neighbours(:,:)
real(kind=kind_real), allocatable :: field_types(:)

! Get nearest neighbours
! ----------------------
allocate(field_neighbours(ngrid_ou,self%nnearest))
call self%unsinterp%apply(field_in, field_ou, field_neighbours)


! Global min and max integers in field
! ------------------------------------
maxtypel = int(maxval(field_in))
mintypel = int(minval(field_in))
call f_comm%allreduce(maxtypel,maxtype,fckit_mpi_max())
call f_comm%allreduce(mintypel,mintype,fckit_mpi_min())


! Put weights into field type array and pick max for interpolated value
! ---------------------------------------------------------------------
allocate(field_types(mintype:maxtype))

field_ou = 0.0_kind_real
do i = 1,ngrid_ou
  field_types = 0.0
  do n = 1, self%nnearest
    index = int(field_neighbours(i,n))
    field_types(index) = field_types(index) + self%unsinterp%interp_w(i,n)
  enddo
  field_ou(i) = real(maxloc(field_types,1)+(mintype-1),kind_real)
enddo

end subroutine iinterp_apply

! --------------------------------------------------------------------------------------------------

end module fv3jedi_interpolation_mod
