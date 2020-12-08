! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_getvalues_mod

! fckit
use fckit_mpi_module,               only: fckit_mpi_comm

! oops
use datetime_mod,                   only: datetime
use unstructured_interpolation_mod, only: unstrc_interp

! saber
use interpolatorbump_mod,         only: bump_interpolator

! ufo
use ufo_locs_mod,                   only: ufo_locs, ufo_locs_time_mask
use ufo_geovals_mod,                only: ufo_geovals

! fv3jedi uses
use fv3jedi_constants_mod,          only: rad2deg
use fv3jedi_field_mod,              only: fv3jedi_field, get_field, field_clen, &
                                          long_name_to_fv3jedi_name
use fv3jedi_geom_mod,               only: fv3jedi_geom
use fv3jedi_interpolation_mod,      only: unsinterp_integer_apply, unsinterp_nearest_apply
use fv3jedi_kinds_mod,              only: kind_real

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: fv3jedi_getvalues, fv3jedi_getvalues_base

type, abstract :: fv3jedi_getvalues_base
  integer                 :: isc, iec, jsc, jec, npz, ngrid
  character(len=2048)     :: interp_method
  integer                 :: nnear = 4
  type(bump_interpolator) :: bumpinterp
  type(unstrc_interp)     :: unsinterp
  type(fckit_mpi_comm)    :: comm
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: fill_geovals
    generic, public :: set_trajectory => fill_geovals
    generic, public :: fill_geovals_tl => fill_geovals
end type fv3jedi_getvalues_base

type, extends(fv3jedi_getvalues_base) :: fv3jedi_getvalues
end type fv3jedi_getvalues

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, locs)

class(fv3jedi_getvalues_base), intent(inout) :: self
type(fv3jedi_geom),            intent(in)    :: geom
type(ufo_locs),                intent(in)    :: locs

! Geometry
! --------
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz
self%ngrid = geom%ngrid
self%comm = geom%f_comm

! Initialize bump interpolator
! ----------------------------
self%interp_method = trim(geom%interp_method)
if (trim(self%interp_method) == 'bump') then
  call self%bumpinterp%init(geom%f_comm, afunctionspace_in=geom%afunctionspace, lon_out=locs%lon, lat_out=locs%lat, nl=geom%npz)
endif

! Always create unstructured interpolation as it is used for special case fields, e.g. integers
call self%unsinterp%create( geom%f_comm, self%nnear, 'barycent', &
                            self%ngrid, rad2deg*geom%lat_us, rad2deg*geom%lon_us, &
                            locs%nlocs, locs%lat, locs%lon )

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_getvalues_base), intent(inout) :: self

if (trim(self%interp_method) == 'bump') call self%bumpinterp%delete()

call self%unsinterp%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine fill_geovals(self, geom, fields, t1, t2, locs, geovals)

class(fv3jedi_getvalues_base), intent(inout) :: self
type(fv3jedi_geom),            intent(in)    :: geom
type(fv3jedi_field),           intent(in)    :: fields(:)
type(datetime),                intent(in)    :: t1
type(datetime),                intent(in)    :: t2
type(ufo_locs),                intent(in)    :: locs
type(ufo_geovals),             intent(inout) :: geovals

integer :: gv, n, ji, jj, jlev
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name
logical, allocatable :: time_mask(:)
real(kind=kind_real), allocatable :: field_us(:)
real(kind=kind_real), allocatable :: geovals_all(:,:), geovals_tmp(:)

! Get mask for locations in this time window
! ------------------------------------------
call ufo_locs_time_mask(locs, t1, t2, time_mask)


! Allocate geovals
! ----------------
if (.not. geovals%linit) then
  do gv = 1, geovals%nvar
    geovals%geovals(gv)%nval = fields(gv)%npz
    allocate(geovals%geovals(gv)%vals(geovals%geovals(gv)%nval, geovals%geovals(gv)%nlocs))
    geovals%geovals(gv)%vals = 0.0_kind_real
  enddo
endif
geovals%linit = .true.


! Loop over GeoVaLs
! -----------------
allocate(field_us(self%ngrid))
allocate(geovals_all(locs%nlocs, self%npz+1))
allocate(geovals_tmp(locs%nlocs))

do gv = 1, geovals%nvar

  ! Get GeoVaLs field
  ! -----------------
  call long_name_to_fv3jedi_name(fields, trim(geovals%variables(gv)), fv3jedi_name)
  call get_field(fields, fv3jedi_name, field)

  ! Interpolation
  ! -------------
  geovals_all = 0.0_kind_real

  ! Can optionally interpolate real valued magnitude fields with bump
  ! -----------------------------------------------------------------
  if ( trim(self%interp_method) == 'bump' .and. &
       .not.field%integerfield .and. trim(field%space)=='magnitude' ) then

    ! Interpolate
    call self%bumpinterp%apply(field%array(field%isc:field%iec,field%jsc:field%jec,1:field%npz),geovals_all(:,1:field%npz))

  else ! Otherwise use unstructured interpolation

    do jlev = 1, field%npz
      n = 0
      do jj = field%jsc, field%jec
        do ji = field%isc, field%iec
          n = n + 1
          field_us(n) = field%array(ji, jj, jlev)
        enddo
      enddo

      ! Conditions for integer and directional fields
      ! ---------------------------------------------
      if (.not. field%integerfield .and. trim(field%space)=='magnitude') then
        call self%unsinterp%apply(field_us, geovals_tmp)
      elseif (field%integerfield) then
        call unsinterp_integer_apply(self%unsinterp, field_us, geovals_tmp)
      elseif (trim(field%space)=='direction') then
        call unsinterp_nearest_apply(self%unsinterp, field_us, geovals_tmp)
      else
        call abor1_ftn("fv3jedi_getvalues_mod.fill_geovals: interpolation for this kind of "// &
                       "field is not supported. FieldName: "// trim(field%fv3jedi_name))
      endif
      geovals_all(1:locs%nlocs, jlev) = geovals_tmp(1:locs%nlocs)
    enddo

  endif

  ! Fill GeoVaLs relevant to this window
  ! ------------------------------------
  do n = 1,locs%nlocs
    if (time_mask(n)) geovals%geovals(gv)%vals(1:field%npz, n) = geovals_all(n, 1:field%npz)
  enddo

enddo

deallocate(field_us)
deallocate(geovals_all)
deallocate(geovals_tmp)
deallocate(time_mask)

end subroutine fill_geovals

! --------------------------------------------------------------------------------------------------

end module fv3jedi_getvalues_mod
