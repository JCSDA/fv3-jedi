! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_lineargetvalues_mod

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
use fv3jedi_field_mod,              only: fv3jedi_field, get_field, field_clen, &
                                          long_name_to_fv3jedi_name
use fv3jedi_geom_mod,               only: fv3jedi_geom
use fv3jedi_getvalues_mod,          only: fv3jedi_getvalues_base
use fv3jedi_increment_mod,          only: fv3jedi_increment
use fv3jedi_kinds_mod,              only: kind_real

implicit none

private
public :: fv3jedi_lineargetvalues

type, extends(fv3jedi_getvalues_base) :: fv3jedi_lineargetvalues
  contains
    procedure, public :: fill_geovals_ad
end type fv3jedi_lineargetvalues

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fill_geovals_ad(self, geom, fields, t1, t2, locs, geovals)

class(fv3jedi_lineargetvalues), intent(inout) :: self
type(fv3jedi_geom),             intent(in)    :: geom
type(fv3jedi_field),            intent(inout) :: fields(:)
type(datetime),                 intent(in)    :: t1
type(datetime),                 intent(in)    :: t2
type(ufo_locs),                 intent(in)    :: locs
type(ufo_geovals),              intent(in)    :: geovals

integer :: gv, n, ji, jj, jlev
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name
logical, allocatable :: time_mask(:)
real(kind=kind_real), allocatable :: field_us(:)
real(kind=kind_real), allocatable :: geovals_all(:,:), geovals_tmp(:)

! Get mask for locations in this time window
! ------------------------------------------
call ufo_locs_time_mask(locs, t1, t2, time_mask)


! Allocate intermediate variables
! -------------------------------
allocate(field_us(self%ngrid))
allocate(geovals_all(locs%nlocs, self%npz+1))
allocate(geovals_tmp(locs%nlocs))
field_us = 0.0_kind_real
geovals_all = 0.0_kind_real
geovals_tmp = 0.0_kind_real


! Loop over GeoVaLs
! -----------------
do gv = 1, geovals%nvar

  ! Get GeoVaLs field
  ! -----------------
  call long_name_to_fv3jedi_name(fields, trim(geovals%variables(gv)), fv3jedi_name)
  call get_field(fields, fv3jedi_name, field)

  ! Adjoint of fill GeoVaLs relevant to this window
  ! -----------------------------------------------
  do n = 1,locs%nlocs
    if (time_mask(n)) geovals_all(n, 1:field%npz) = geovals%geovals(gv)%vals(1:field%npz, n)
  enddo

  ! Adjoint of interpolation
  ! ------------------------
  if ( trim(self%interp_method) == 'bump' .and. &
       .not.field%integerfield .and. trim(field%space)=='magnitude' ) then

    ! Interpolate
    call self%bumpinterp%apply_ad(geovals_all(:,1:field%npz),field%array(field%isc:field%iec,field%jsc:field%jec,1:field%npz))

  else ! Otherwise use unstructured interpolation

    do jlev = 1, field%npz

      geovals_tmp(1:locs%nlocs) = geovals_all(1:locs%nlocs, jlev)

      ! Conditions for integer and directional fields
      ! ---------------------------------------------
      if (.not. field%integerfield .and. trim(field%space)=='magnitude') then
        call self%unsinterp%apply_ad(field_us, geovals_tmp)
      else
        call abor1_ftn("fv3jedi_getvalues_mod.fill_geovals: interpolation for this kind of "// &
                       "field is not supported. FieldName: "// trim(field%fv3jedi_name))
      endif

      n = 0
      do jj = field%jsc, field%jec
        do ji = field%isc, field%iec
          n = n + 1
          field%array(ji, jj, jlev) = field_us(n)
        enddo
      enddo

    enddo

  endif

enddo

deallocate(time_mask)
deallocate(field_us)
deallocate(geovals_all)
deallocate(geovals_tmp)

end subroutine fill_geovals_ad

! ------------------------------------------------------------------------------

end module fv3jedi_lineargetvalues_mod
