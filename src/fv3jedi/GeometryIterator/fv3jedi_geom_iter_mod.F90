!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_geom_iter_mod

  use iso_c_binding
  use kinds
  use fv3jedi_geom_mod, only: fv3jedi_geom
  use fv3jedi_constants_mod, only: rad2deg

  implicit none

  private
  public :: fv3jedi_geom_iter
  public :: fv3jedi_geom_iter_registry
  public :: fv3jedi_geom_iter_setup, fv3jedi_geom_iter_clone, fv3jedi_geom_iter_equals
  public :: fv3jedi_geom_iter_current, fv3jedi_geom_iter_next
  public :: fv3jedi_geom_iter_orography

  type :: fv3jedi_geom_iter
    type(fv3jedi_geom), pointer :: geom => null() !< Geometry
    integer :: iindex = 1  !< index e.g. lat(iindex,jindex)
    integer :: jindex = 1  !<
  end type fv3jedi_geom_iter

#define LISTED_TYPE fv3jedi_geom_iter

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: fv3jedi_geom_iter_registry

contains

  ! ------------------------------------------------------------------------------
  ! Public
  ! ------------------------------------------------------------------------------

  !> Linked list implementation
#include "oops/util/linkedList_c.f"

  ! ------------------------------------------------------------------------------
  !> Setup for the geometry iterator
  subroutine fv3jedi_geom_iter_setup(self, geom, iindex, jindex)

    ! Passed variables
    type(fv3jedi_geom_iter),     intent(inout) :: self !< Geometry iterator
    type(fv3jedi_geom), pointer, intent(   in) :: geom !< Geometry
    integer,                  intent(   in) :: iindex, jindex  !< Index

    ! Associate geometry
    self%geom => geom

    ! Define iindex/jindex for local tile
    self%iindex = iindex
    self%jindex = jindex

  end subroutine fv3jedi_geom_iter_setup

  ! ------------------------------------------------------------------------------
  !> Clone for the geometry iterator
  subroutine fv3jedi_geom_iter_clone(self, other)

    ! Passed variables
    type(fv3jedi_geom_iter), intent(inout) :: self  !< Geometry iterator
    type(fv3jedi_geom_iter), intent(   in) :: other !< Other geometry iterator

    ! Associate geometry
    self%geom => other%geom

    ! Copy iindex/jindex
    self%iindex = other%iindex
    self%jindex = other%jindex

  end subroutine fv3jedi_geom_iter_clone

  ! ------------------------------------------------------------------------------
  !> Check for the geometry iterator equality
  subroutine fv3jedi_geom_iter_equals(self, other, equals)

    ! Passed variables
    type(fv3jedi_geom_iter), intent( in) :: self   !< Geometry iterator
    type(fv3jedi_geom_iter), intent( in) :: other  !< Other geometry iterator
    integer,            intent(out) :: equals !< Equality flag

    ! Initialization
    equals = 0

    ! Check equality
    if (associated(self%geom, other%geom) .and. (self%iindex==other%iindex) .and. (self%jindex==other%jindex)) equals = 1

  end subroutine fv3jedi_geom_iter_equals

  ! ------------------------------------------------------------------------------
  !> Get geometry iterator current lat/lon
  subroutine fv3jedi_geom_iter_current(self, lon, lat)

    ! Passed variables
    type(fv3jedi_geom_iter), intent( in) :: self !< Geometry iterator
    real(kind_real),    intent(out) :: lat  !< Latitude
    real(kind_real),    intent(out) :: lon  !< Longitude

    ! Check iindex/jindex
    if (self%iindex == -1 .AND. self%jindex == -1) then
      ! special case of {-1,-1} means end of the grid
      lat = self%geom%grid_lat(self%geom%iec,self%geom%jec)
      lon = self%geom%grid_lon(self%geom%iec,self%geom%jec)
    elseif (self%iindex < self%geom%isc .OR. self%iindex > self%geom%iec .OR. &
            self%jindex < self%geom%jsc .OR. self%jindex > self%geom%jec) then
      ! outside of the grid
      call abor1_ftn('fv3jedi_geom_iter_current: iterator out of bounds')
    else
      ! inside of the grid
      lat = self%geom%grid_lat(self%iindex,self%jindex)
      lon = self%geom%grid_lon(self%iindex,self%jindex)
    endif

    !convert to degrees from radians
    lat = rad2deg*lat
    lon = rad2deg*lon

  end subroutine fv3jedi_geom_iter_current

  ! ------------------------------------------------------------------------------
  !> Get geometry iterator current lat/lon
  subroutine fv3jedi_geom_iter_orography(self, oro)

    ! Passed variables
    type(fv3jedi_geom_iter), intent( in) :: self !< Geometry iterator
    real(kind_real),    intent(out) :: oro  !< Orography

    ! Check iindex/jindex
    if (self%iindex == -1 .AND. self%jindex == -1) then
      ! special case of {-1,-1} means end of the grid
      oro = self%geom%orography(self%geom%iec,self%geom%jec,1)
    elseif (self%iindex < self%geom%isc .OR. self%iindex > self%geom%iec .OR. &
            self%jindex < self%geom%jsc .OR. self%jindex > self%geom%jec) then
      ! outside of the grid
      call abor1_ftn('fv3jedi_geom_iter_orography: iterator out of bounds')
    else
      ! inside of the grid
      oro = self%geom%orography(self%iindex,self%jindex,1)
    endif

  end subroutine fv3jedi_geom_iter_orography

  ! ------------------------------------------------------------------------------
  !> Update geometry iterator to next point
  subroutine fv3jedi_geom_iter_next(self)

    ! Passed variables
    type(fv3jedi_geom_iter), intent(inout) :: self !< Geometry iterator
    integer :: iindex, jindex

    iindex = self%iindex
    jindex = self%jindex

    ! increment by 1
    if (iindex.lt.self%geom%iec) then
      iindex = iindex + 1
    elseif (iindex.eq.self%geom%iec) then
      iindex = self%geom%isc
      jindex = jindex + 1
    end if

    if (jindex > self%geom%jec) then
        iindex=-1
        jindex=-1
    end if

    self%iindex = iindex
    self%jindex = jindex

  end subroutine fv3jedi_geom_iter_next
  ! ------------------------------------------------------------------------------

end module fv3jedi_geom_iter_mod
