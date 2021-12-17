! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_getvalues_mod

use iso_c_binding

! fckit
use fckit_configuration_module,     only: fckit_configuration
use fckit_mpi_module,               only: fckit_mpi_comm

! oops
use datetime_mod,                   only: datetime
use unstructured_interpolation_mod, only: unstrc_interp

! ufo
use ufo_locations_mod
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

subroutine create(self, geom, locs, conf)

class(fv3jedi_getvalues_base), intent(inout) :: self
type(fv3jedi_geom),            intent(in)    :: geom
type(ufo_locations),           intent(in)    :: locs
type(fckit_configuration),     intent(in)    :: conf

real(8), allocatable, dimension(:) :: lons, lats
integer :: nlocs
character(len=:), allocatable :: str

! Geometry
! --------
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz
self%ngrid = geom%ngrid
self%comm = geom%f_comm

nlocs = locs%nlocs()
allocate(lons(nlocs), lats(nlocs))
call locs%get_lons(lons)
call locs%get_lats(lats)

! Get interpolation method
! ------------------------
self%interp_method = geom%interp_method
if (conf%has("interpolation method")) then
    call conf%get_or_die("interpolation method",str)
    self%interp_method = str
    deallocate(str)
endif

! Initialize bump interpolator
! ----------------------------
if (trim(self%interp_method) == 'bump') then
  call abor1_ftn("fv3jedi_getvalues_mod: bump interpolation not supported")
endif

! Always create unstructured interpolation as it is used for special case fields, e.g. integers
call self%unsinterp%create( geom%f_comm, self%nnear, 'barycent', &
                            self%ngrid, rad2deg*geom%lat_us, rad2deg*geom%lon_us, &
                            nlocs, lats, lons )

deallocate(lons, lats)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_getvalues_base), intent(inout) :: self

call self%unsinterp%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine fill_geovals(self, geom, fields, t1, t2, locs, geovals)
implicit none
class(fv3jedi_getvalues_base), intent(inout) :: self
type(fv3jedi_geom),            intent(in)    :: geom
type(fv3jedi_field),           intent(in)    :: fields(:)
type(datetime),                intent(in)    :: t1
type(datetime),                intent(in)    :: t2
type(ufo_locations),           intent(in)    :: locs
type(ufo_geovals),             intent(inout) :: geovals

integer :: gv, n, ji, jj, jlev, nlocs
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name
logical(c_bool), allocatable :: time_mask(:)
real(kind=kind_real), allocatable :: field_us(:)
real(kind=kind_real), allocatable :: geovals_all(:,:)

! Get mask for locations in this time window
! ------------------------------------------
nlocs = locs%nlocs()
allocate(time_mask(nlocs))
call locs%get_timemask(t1, t2, time_mask)

! Loop over GeoVaLs
! -----------------
allocate(field_us(self%ngrid))
allocate(geovals_all(nlocs, self%npz+1))

do gv = 1, geovals%nvar

  ! Get GeoVaLs field
  ! -----------------
  call long_name_to_fv3jedi_name(fields, trim(geovals%variables(gv)), fv3jedi_name)
  call get_field(fields, fv3jedi_name, field)

  ! Interpolation
  ! -------------

  ! Can optionally interpolate real valued magnitude fields with bump
  ! -----------------------------------------------------------------
  if ( trim(self%interp_method) == 'bump' .and. trim(field%interp_type) == "default") then
    call abor1_ftn("fv3jedi_getvalues_mod: bump interpolation not supported")
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
      if ( trim(field%interp_type) == "default" ) then
        call self%unsinterp%apply(field_us, geovals_all(:, jlev))
      elseif ( trim(field%interp_type) == "integer" ) then
        call unsinterp_integer_apply(self%unsinterp, field_us, geovals_all(:, jlev))
      elseif ( trim(field%interp_type) == "nearest" ) then
        call unsinterp_nearest_apply(self%unsinterp, field_us, geovals_all(:, jlev))
      else
        call abor1_ftn("fv3jedi_getvalues_mod.fill_geovals: interpolation for this kind of "// &
                       "field is not supported. FieldName: "// trim(field%fv3jedi_name)// &
                       "interp_type: "//trim(field%interp_type))
      endif
    enddo

  endif

  ! Fill GeoVaLs relevant to this window
  ! ------------------------------------
  do n = 1,nlocs
    if (time_mask(n)) geovals%geovals(gv)%vals(1:field%npz, n) = geovals_all(n, 1:field%npz)
  enddo

enddo

deallocate(field_us)
deallocate(geovals_all)
deallocate(time_mask)
geovals%linit = .true.

end subroutine fill_geovals

! --------------------------------------------------------------------------------------------------

end module fv3jedi_getvalues_mod
