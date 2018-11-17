! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_geos_mod

use iso_c_binding
use config_mod
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 

implicit none
private

public :: geos_model
public :: geos_create
public :: geos_delete
public :: geos_initialize
public :: geos_step
public :: geos_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: geos_model
  integer :: geos_goes_here
end type geos_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine geos_create(self, geom, c_conf)

implicit none
type(c_ptr),        intent(in)    :: c_conf
type(geos_model),   intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom

end subroutine geos_create

! ------------------------------------------------------------------------------

subroutine geos_delete(self)

implicit none
type(geos_model), intent(inout) :: self

end subroutine geos_delete

! ------------------------------------------------------------------------------

subroutine geos_initialize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

end subroutine geos_initialize

! ------------------------------------------------------------------------------

subroutine geos_step(self, state, vdate)

implicit none
type(geos_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate !< Valid datetime after step

call state_to_geos( state, self )

!call advance geos

call geos_to_state( self, state )

end subroutine geos_step

! ------------------------------------------------------------------------------

subroutine geos_finalize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

end subroutine geos_finalize

! ------------------------------------------------------------------------------

subroutine state_to_geos( state, self )

implicit none
type(fv3jedi_state), intent(in)    :: state
type(geos_model),    intent(inout) :: self
 
end subroutine state_to_geos

! ------------------------------------------------------------------------------

subroutine geos_to_state( self, state )

implicit none
type(geos_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state
 
end subroutine geos_to_state

! ------------------------------------------------------------------------------

end module fv3jedi_geos_mod
