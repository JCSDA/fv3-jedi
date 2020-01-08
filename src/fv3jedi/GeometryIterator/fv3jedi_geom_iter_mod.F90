!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_geom_iter_mod

  use iso_c_binding
  use kinds
  use fv3jedi_geom_mod, only: fv3jedi_geom

  implicit none

  private
  public :: fv3jedi_geom_iter
  public :: fv3jedi_geom_iter_registry
  public :: fv3jedi_geom_iter_setup, fv3jedi_geom_iter_clone, fv3jedi_geom_iter_equals
  public :: fv3jedi_geom_iter_current, fv3jedi_geom_iter_next

  type :: fv3jedi_geom_iter
    type(fv3jedi_geom), pointer :: geom => null() !< Geometry
    integer :: iind = 1  !< index e.g. lat(iind,jind)
    integer :: jind = 1  !< 
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
  subroutine fv3jedi_geom_iter_setup(self, geom, iind, jind)

    ! Passed variables
    type(fv3jedi_geom_iter),     intent(inout) :: self !< Geometry iterator
    type(fv3jedi_geom), pointer, intent(   in) :: geom !< Geometry
    integer,                  intent(   in) :: iind, jind  !< Index

    ! Associate geometry
    self%geom => geom

    ! Define iind/jind for local tile
    self%iind = iind
    self%jind = jind

  end subroutine fv3jedi_geom_iter_setup

  ! ------------------------------------------------------------------------------
  !> Clone for the geometry iterator
  subroutine fv3jedi_geom_iter_clone(self, other)

    ! Passed variables
    type(fv3jedi_geom_iter), intent(inout) :: self  !< Geometry iterator
    type(fv3jedi_geom_iter), intent(   in) :: other !< Other geometry iterator

    ! Associate geometry
    self%geom => other%geom

    ! Copy iind/jind
    self%iind = other%iind
    self%jind = other%jind

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
    if (associated(self%geom, other%geom) .and. (self%iind==other%iind) .and. (self%jind==other%jind)) equals = 1

  end subroutine fv3jedi_geom_iter_equals

  ! ------------------------------------------------------------------------------
  !> Get geometry iterator current lat/lon
  subroutine fv3jedi_geom_iter_current(self, lon, lat)

    ! Passed variables
    type(fv3jedi_geom_iter), intent( in) :: self !< Geometry iterator
    real(kind_real),    intent(out) :: lat  !< Latitude
    real(kind_real),    intent(out) :: lon  !< Longitude

    ! Check iind/jind
    if (self%iind == -1 .AND. self%jind == -1) then
      ! special case of {-1,-1} means end of the grid
      lat = self%geom%grid_lat(self%geom%iec,self%geom%jec)
      lon = self%geom%grid_lon(self%geom%iec,self%geom%jec) 
    elseif (self%iind < self%geom%isc .OR. self%iind > self%geom%iec .OR. &
            self%jind < self%geom%jsc .OR. self%jind > self%geom%jec) then
      ! outside of the grid
      call abor1_ftn('fv3jedi_geom_iter_current: iterator out of bounds')
    else
      ! inside of the grid
      lat = self%geom%grid_lat(self%iind,self%jind)
      lon = self%geom%grid_lon(self%iind,self%jind)
    endif

  end subroutine fv3jedi_geom_iter_current

  ! ------------------------------------------------------------------------------
  !> Update geometry iterator to next point
  subroutine fv3jedi_geom_iter_next(self)

    ! Passed variables
    type(fv3jedi_geom_iter), intent(inout) :: self !< Geometry iterator
    integer :: iind, jind

    iind = self%iind
    jind = self%jind

    ! increment by 1
    if (iind.lt.self%geom%iec) then 
      iind = iind + 1
    elseif (iind.eq.self%geom%iec) then
      iind = self%geom%isc
      jind = jind + 1
    end if

    if (jind > self%geom%jec) then
        iind=-1
        jind=-1
    end if

    self%iind = iind
    self%jind = jind

  end subroutine fv3jedi_geom_iter_next
  ! ------------------------------------------------------------------------------

end module fv3jedi_geom_iter_mod
