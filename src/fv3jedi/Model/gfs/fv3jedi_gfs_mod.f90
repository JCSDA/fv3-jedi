! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_gfs_mod

use iso_c_binding
use config_mod
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment

use fv3gfs_cap_mod, only: FV3SS => SetServices
use esmf

implicit none
private

public :: gfs_model
public :: gfs_create
public :: gfs_delete
public :: gfs_initialize
public :: gfs_step
public :: gfs_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: gfs_model
  integer :: gfs_goes_here
end type gfs_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine gfs_create(self, geom, c_conf)

implicit none
type(c_ptr),        intent(in)    :: c_conf
type(gfs_model),   intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom

end subroutine gfs_create

! ------------------------------------------------------------------------------

subroutine gfs_delete(self)

implicit none
type(gfs_model), intent(inout) :: self

end subroutine gfs_delete

! ------------------------------------------------------------------------------

subroutine gfs_initialize(self, state)

implicit none
type(gfs_model), target :: self
type(fv3jedi_state)      :: state

end subroutine gfs_initialize

! ------------------------------------------------------------------------------

subroutine gfs_step(self, state, vdate)

implicit none
type(gfs_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate !< Valid datetime after step

call state_to_gfs( state, self )

!call advance gfs

call gfs_to_state( self, state )

end subroutine gfs_step

! ------------------------------------------------------------------------------

subroutine gfs_finalize(self, state)

implicit none
type(gfs_model), target :: self
type(fv3jedi_state)      :: state

end subroutine gfs_finalize

! ------------------------------------------------------------------------------

subroutine state_to_gfs( state, self )

implicit none
type(fv3jedi_state), intent(in)    :: state
type(gfs_model),    intent(inout) :: self

end subroutine state_to_gfs

! ------------------------------------------------------------------------------

subroutine gfs_to_state( self, state )

implicit none
type(gfs_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state

end subroutine gfs_to_state

! ------------------------------------------------------------------------------

end module fv3jedi_gfs_mod
