! (C) Copyright 2018-2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_varcha_a2m_mod

use iso_c_binding
use config_mod
use fv3jedi_kinds_mod

use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state

implicit none
private

public :: fv3jedi_varcha_a2m
public :: create
public :: delete
public :: changevar
public :: changevarinverse

type :: fv3jedi_varcha_a2m
  integer :: dummy
end type fv3jedi_varcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, c_conf)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(c_ptr),              intent(in)    :: c_conf

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self

end subroutine delete

! ------------------------------------------------------------------------------

subroutine changevar(self,geom,xana,xmod)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xmod

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xana,xmod)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xmod

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_a2m_mod
