!
! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_geom_iter_interface

  use iso_c_binding
  use kinds
  use fv3jedi_geom_iter_mod
  use fv3jedi_geom_interface_mod,  only : fv3jedi_geom_registry
  use fv3jedi_geom_mod, only: fv3jedi_geom

  implicit none

  private

contains

  ! ------------------------------------------------------------------------------
  !> Setup geometry iterator
  subroutine fv3jedi_geom_iter_setup_c(c_key_self, c_key_geom, &
                                       c_iindex, c_jindex, c_kindex) & 
                                  bind(c, name='fv3jedi_geom_iter_setup_f90')

    ! Passed variables
    integer(c_int), intent(inout) :: c_key_self !< Geometry iterator
    integer(c_int), intent(   in) :: c_key_geom !< Geometry
    integer(c_int), intent(   in) :: c_iindex   !< Index
    integer(c_int), intent(   in) :: c_jindex   !< Index
    integer(c_int), intent(   in) :: c_kindex   !< Index

    ! Local variables
    type(fv3jedi_geom_iter),     pointer :: self
    type(fv3jedi_geom),          pointer :: geom

    ! Interface
    call fv3jedi_geom_iter_registry%init()
    call fv3jedi_geom_iter_registry%add(c_key_self)
    call fv3jedi_geom_iter_registry%get(c_key_self, self)
    call fv3jedi_geom_registry%get(c_key_geom, geom)

    ! Call Fortran
    call fv3jedi_geom_iter_setup(self, geom, c_iindex, c_jindex, c_kindex)

  end subroutine fv3jedi_geom_iter_setup_c

  ! ------------------------------------------------------------------------------
  !> Clone geometry iterator
  subroutine fv3jedi_geom_iter_clone_c(c_key_self, c_key_other) bind(c, name='fv3jedi_geom_iter_clone_f90')

    ! Passed variables
    integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
    integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator

    ! Local variables
    type(fv3jedi_geom_iter), pointer :: self, other

    ! Interface
    call fv3jedi_geom_iter_registry%get(c_key_other, other)
    call fv3jedi_geom_iter_registry%init()
    call fv3jedi_geom_iter_registry%add(c_key_self)
    call fv3jedi_geom_iter_registry%get(c_key_self, self)

    ! Call Fortran
    call fv3jedi_geom_iter_clone(self, other)

  end subroutine fv3jedi_geom_iter_clone_c

  ! ------------------------------------------------------------------------------
  !> Delete geometry iterator
  subroutine fv3jedi_geom_iter_delete_c(c_key_self) bind(c, name='fv3jedi_geom_iter_delete_f90')

      ! Passed variables
      integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

      ! Clear interface
      call fv3jedi_geom_iter_registry%remove(c_key_self)

  end subroutine fv3jedi_geom_iter_delete_c

  ! ------------------------------------------------------------------------------
  !> Check geometry iterator equality
  subroutine fv3jedi_geom_iter_equals_c(c_key_self, c_key_other, c_equals) bind(c, name='fv3jedi_geom_iter_equals_f90')

    ! Passed variables
    integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
    integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator
    integer(c_int), intent(inout) :: c_equals    !< Equality flag

    ! Local variables
    type(fv3jedi_geom_iter),pointer :: self,other

    ! Interface
    call fv3jedi_geom_iter_registry%get(c_key_self, self)
    call fv3jedi_geom_iter_registry%get(c_key_other, other)

    ! Call Fortran
    call fv3jedi_geom_iter_equals(self, other, c_equals)

  end subroutine fv3jedi_geom_iter_equals_c

  ! ------------------------------------------------------------------------------
  !> Get geometry iterator current lon/lat/vCoord
  subroutine fv3jedi_geom_iter_current_c(c_key_self, c_lon, c_lat, c_vCoord) &
                                    bind(c, name='fv3jedi_geom_iter_current_f90')

    ! Passed variables
    integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
    real(c_double), intent(inout) :: c_lat      !< Latitude
    real(c_double), intent(inout) :: c_lon      !< Longitude
    real(c_double), intent(inout) :: c_vCoord   !< Vertical Coordinator

    ! Local variables
    type(fv3jedi_geom_iter), pointer :: self

    ! Interface
    call fv3jedi_geom_iter_registry%get(c_key_self, self)

    ! Call Fortran
    call fv3jedi_geom_iter_current(self, c_lon, c_lat, c_vCoord)

  end subroutine fv3jedi_geom_iter_current_c

  ! ------------------------------------------------------------------------------
  !> Get geometry iterator current orography
  subroutine fv3jedi_geom_iter_orography_c(c_key_self, c_oro) bind(c, name='fv3jedi_geom_iter_orography_f90')

    ! Passed variables
    integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
    real(c_double), intent(inout) :: c_oro      !< Orography

    ! Local variables
    type(fv3jedi_geom_iter), pointer :: self

    ! Interface
    call fv3jedi_geom_iter_registry%get(c_key_self, self)

    ! Call Fortran
    call fv3jedi_geom_iter_orography(self, c_oro)

  end subroutine fv3jedi_geom_iter_orography_c

  ! ------------------------------------------------------------------------------
  !> Update geometry iterator to next point
  subroutine fv3jedi_geom_iter_next_c(c_key_self) bind(c, name='fv3jedi_geom_iter_next_f90')

    ! Passed variables
    integer(c_int), intent(in) :: c_key_self !< Geometry iterator

    ! Local variables
    type(fv3jedi_geom_iter), pointer :: self

    ! Interface
    call fv3jedi_geom_iter_registry%get(c_key_self, self)

    ! Call Fortran
    call fv3jedi_geom_iter_next(self)

  end subroutine fv3jedi_geom_iter_next_c

end module fv3jedi_geom_iter_interface
