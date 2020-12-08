! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_interpolation_mod

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm, fckit_mpi_max, fckit_mpi_min

! oops
use unstructured_interpolation_mod, only: unstrc_interp

! saber
use interpolatorbump_mod,         only: bump_interpolator

! fv3jedi
use fv3jedi_kinds_mod,              only: kind_real
use fv3jedi_field_mod,              only: fv3jedi_field, get_field, has_field
use fv3jedi_geom_mod,               only: fv3jedi_geom
use fv3jedi_constants_mod,          only: rad2deg
use wind_vt_mod,                    only: d2a, a2d

implicit none
private
public :: field2field_interp, unsinterp_integer_apply, unsinterp_nearest_apply

! --------------------------------------------------------------------------------------------------

type field2field_interp

  character(len=32) :: interp_type
  type(bump_interpolator) :: bumpinterp
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

! Initialize bump interpolator
! ----------------------------
if (self%need_bump) then
  call self%bumpinterp%init(geom_in%f_comm, afunctionspace_in=geom_in%afunctionspace, afunctionspace_out=geom_ou%afunctionspace, &
     & nl=geom_in%npz)
endif

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

if (self%need_bump) call self%bumpinterp%delete()
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
real(kind=kind_real), allocatable :: field_in(:), field_ou(:)

logical :: do_d2a
real(kind=kind_real), allocatable :: ua(:,:,:)
real(kind=kind_real), allocatable :: va(:,:,:)
real(kind=kind_real), allocatable :: ud(:,:,:)
real(kind=kind_real), allocatable :: vd(:,:,:)
type(fv3jedi_field), pointer :: u_in
type(fv3jedi_field), pointer :: v_in
type(fv3jedi_field), pointer :: u_ou
type(fv3jedi_field), pointer :: v_ou

! Special case of D-grid winds
! ----------------------------
do_d2a = .false.
if (has_field(fields_in,'ud')) then

  do_d2a = .true.

  call get_field(fields_in, 'ud', u_in)
  call get_field(fields_in, 'vd', v_in)

  allocate(ua(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,geom_in%npz))
  allocate(va(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,geom_in%npz))

  call d2a(geom_in, u_in%array, v_in%array, ua, va)

  ! Reallocate without staggering
  deallocate(u_in%array)
  deallocate(v_in%array)
  allocate(u_in%array(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,1:geom_in%npz))
  allocate(v_in%array(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,1:geom_in%npz))

  ! Reallocate output
  call get_field(fields_ou, 'ud', u_ou)
  call get_field(fields_ou, 'vd', v_ou)
  deallocate(u_ou%array)
  deallocate(v_ou%array)
  allocate(u_ou%array(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,1:geom_ou%npz))
  allocate(v_ou%array(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,1:geom_ou%npz))

  ! Overwrite field with A-Grid for doing interpolation
  u_in%array = ua
  v_in%array = va

  ! Change field type
  u_in%space = 'magnitude'
  v_in%space = 'magnitude'

  deallocate(ua, va)

endif


! Interpolate all fields
! ---------------------
do var = 1,nf
  if (trim(self%interp_type) == 'bump') then

    ! Use BUMP
    call self%bumpinterp%apply(fields_in(var)%array(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,1:fields_in(var)%npz), &
       & fields_ou(var)%array(:,:,1:fields_ou(var)%npz), &
       & nn=(fields_in(var)%integerfield .or. trim(fields_in(var)%space) == 'direction'))

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
      if (.not. fields_in(var)%integerfield .and. trim(fields_in(var)%space) == 'magnitude' .and. &
          trim(self%interp_type) == 'barycent') then

        call self%unsinterp%apply(field_in, field_ou)

      elseif (fields_in(var)%integerfield) then

        call unsinterp_integer_apply(self%unsinterp, field_in, field_ou)

      elseif (trim(fields_in(var)%space) == 'direction') then

        call unsinterp_nearest_apply(self%unsinterp, field_in, field_ou)

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

  call get_field(fields_ou, 'ud', u_ou)
  call get_field(fields_ou, 'vd', v_ou)

  allocate(ud(geom_ou%isc:geom_ou%iec  ,geom_ou%jsc:geom_ou%jec+1,geom_ou%npz))
  allocate(vd(geom_ou%isc:geom_ou%iec+1,geom_ou%jsc:geom_ou%jec  ,geom_ou%npz))

  call a2d(geom_ou, u_ou%array, v_ou%array, ud, vd)

  ! Reallocate with staggering
  deallocate(u_ou%array)
  deallocate(v_ou%array)
  allocate(u_ou%array(geom_ou%isc:geom_ou%iec  ,geom_ou%jsc:geom_ou%jec+1,1:geom_ou%npz))
  allocate(v_ou%array(geom_ou%isc:geom_ou%iec+1,geom_ou%jsc:geom_ou%jec  ,1:geom_ou%npz))

  ! Put new D grid back into arrays
  u_ou%array = ud
  v_ou%array = vd

  u_ou%space = 'vector'
  v_ou%space = 'vector'

  deallocate(ud, vd)

endif

end subroutine apply

! --------------------------------------------------------------------------------------------------

subroutine unsinterp_integer_apply(unsinterp, field_in, field_ou)

type(unstrc_interp),      intent(in)    :: unsinterp
real(kind=kind_real),     intent(in)    :: field_in(:) !Integer field in
real(kind=kind_real),     intent(inout) :: field_ou(:) !Integer field out

! Locals
integer :: maxtypel, mintypel, maxtype, mintype, ngrid_ou
integer :: i, j, k, n, index
real(kind=kind_real), allocatable :: interp_w(:,:)
real(kind=kind_real), allocatable :: field_ou_tmp(:)
real(kind=kind_real), allocatable :: field_neighbours(:,:)
real(kind=kind_real), allocatable :: field_types(:)

! Inteprolation of integer fields

! Size of output
! --------------
ngrid_ou = size(field_ou)

! Get nearest neighbours
! ----------------------
allocate(field_neighbours(unsinterp%nn,ngrid_ou))
allocate(field_ou_tmp(ngrid_ou))
call unsinterp%apply(field_in, field_ou_tmp, field_neighbours)

! Global min and max integers in field
! ------------------------------------
maxtypel = int(maxval(field_in))
mintypel = int(minval(field_in))
call unsinterp%comm%allreduce(maxtypel,maxtype,fckit_mpi_max())
call unsinterp%comm%allreduce(mintypel,mintype,fckit_mpi_min())


! Put weights into field type array and pick max for interpolated value
! ---------------------------------------------------------------------
allocate(field_types(mintype:maxtype))

field_ou = 0.0_kind_real
do i = 1,ngrid_ou
  field_types = 0.0
  do n = 1, unsinterp%nn
    index = int(field_neighbours(n,i))
    field_types(index) = field_types(index) + unsinterp%interp_w(n,i)
  enddo
  field_ou(i) = real(maxloc(field_types,1)+(mintype-1),kind_real)
enddo

end subroutine unsinterp_integer_apply

! --------------------------------------------------------------------------------------------------

subroutine unsinterp_nearest_apply(unsinterp, field_in, field_ou)

type(unstrc_interp),      intent(in)    :: unsinterp
real(kind=kind_real),     intent(in)    :: field_in(:) !Integer field in
real(kind=kind_real),     intent(inout) :: field_ou(:) !Integer field out

integer :: n, ngrid_ou
real(kind=kind_real), allocatable :: field_ou_tmp(:)
real(kind=kind_real), allocatable :: field_neighbours(:,:)

! Inteprolation using the nearest neighbour

! Size of output
! --------------
ngrid_ou = size(field_ou)

! Get nearest neighbours
! ----------------------
allocate(field_neighbours(unsinterp%nn,ngrid_ou))
allocate(field_ou_tmp(ngrid_ou))
call unsinterp%apply(field_in, field_ou_tmp, field_neighbours)

! Find nearest neighbour
! ----------------------
do n = 1, ngrid_ou
  field_ou(n) = field_neighbours(minloc(unsinterp%interp_w(:,n),1),n)
enddo


end subroutine unsinterp_nearest_apply

! --------------------------------------------------------------------------------------------------

end module fv3jedi_interpolation_mod
