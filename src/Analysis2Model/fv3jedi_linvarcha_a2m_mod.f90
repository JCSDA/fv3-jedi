! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_linvarcha_a2m_mod

use iso_c_binding
use config_mod
use fv3jedi_kinds_mod

use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment

implicit none
private

public :: fv3jedi_linvarcha_a2m
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

type :: fv3jedi_linvarcha_a2m
 integer :: dummy
end type fv3jedi_linvarcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, bg, fg, geom, c_conf)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_state),         intent(in)    :: bg
type(fv3jedi_state),         intent(in)    :: fg
type(fv3jedi_geom),          intent(in)    :: geom
type(c_ptr),                 intent(in)    :: c_conf

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xctl,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xctl
type(fv3jedi_increment),     intent(inout) :: xmod

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xmod,xctl)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xctl

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xmod,xctl)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xctl

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xctl,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(in)    :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fv3jedi_increment),     intent(in)    :: xctl
type(fv3jedi_increment),     intent(inout) :: xmod

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_mod
