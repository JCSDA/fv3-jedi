! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_model_mod

use kinds
use iso_c_binding
use config_mod
use duration_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state

implicit none
private

public :: fv3jedi_model
public :: model_create
public :: model_delete
public :: model_initialize
public :: model_step
public :: model_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: fv3jedi_model
 integer :: gfs_gfs_here
end type fv3jedi_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine model_create(self, geom, c_conf)

implicit none
type(fv3jedi_model), intent(inout) :: self    !Model type
type(fv3jedi_geom)   intent(in)    :: geom    !Geometry
type(c_ptr),         intent(in)    :: c_conf  !User configuration

!Create the GFS model

end subroutine model_create

! ------------------------------------------------------------------------------

subroutine model_delete(self)

implicit none
type(fv3jedi_model), intent(inout) :: self    !Model type

!Delete the GFS model

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_initialize(self, state)

implicit none
type(fv3jedi_model), intent(inout) :: self    !Model type
type(fv3jedi_state), intent(in)    :: state   !JEDI state fields

!Initialize the GFS model

end subroutine model_initialize

! ------------------------------------------------------------------------------

subroutine model_step(self, state, vdate)

implicit none
type(fv3jedi_model), intent(inout) :: self    !Model type
type(fv3jedi_state), intent(inout) :: state   !JEDI state fields
type(datetime),      intent(in)    :: vdate   !Current time

call state_to_gfs(state,self)

!Propagate GFS model by one time step

call gfs_to_state(state,self)

end subroutine model_step

! ------------------------------------------------------------------------------

subroutine model_finalize(self, state)

implicit none
type(fv3jedi_model), intent(inout) :: self    !Model type
type(fv3jedi_state), intent(in)    :: state   !JEDI state fields

!Finalize GFS model

end subroutine model_finalize

! ------------------------------------------------------------------------------

subroutine state_to_gfs(state,self)

implicit none
type(fv3jedi_state), intent(in)    :: state
type(fv3jedi_model), intent(inout) :: self

!Go from JEDI state to GFS fields

end subroutine state_to_gfs

! ------------------------------------------------------------------------------

subroutine gfs_to_state(self,state)

implicit none
type(fv3jedi_model), intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state

!Go from GFS fields to JEDI state

end subroutine gfs_to_state

! ------------------------------------------------------------------------------

end module fv3jedi_model_mod
