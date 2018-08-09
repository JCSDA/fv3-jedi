! 
! (C) Copyright 2017-2018  UCAR. 
!  
! This software is licensed under the terms of the Apache Licence Version 2.0 
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! 

!> Fortran module to handle variables for the FV3JEDI model

module fv3jedi_vars_mod

use iso_c_binding
use config_mod

implicit none
private
public :: fv3jedi_vars, fv3jedi_vars_create

! ------------------------------------------------------------------------------
!> Fortran derived type to represent fv3jedi model variables

type :: fv3jedi_vars
  integer :: nv
  character(len=100), allocatable :: fldnames(:) !< Variable identifiers
end type fv3jedi_vars

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine fv3jedi_vars_create(c_vars,self)
implicit none

type(c_ptr),        intent(in)    :: c_vars
type(fv3jedi_vars), intent(inout) :: self

character(len=1023) :: varlist
integer :: n

!Get long comma seperated list of variables
varlist = config_get_string(c_vars,len(varlist),"variables")

!Count number of commas
self%nv = 1 + count(transfer(varlist, 'a', len(varlist)) == ",")

!Allocate array to hold variable names
allocate(self%fldnames(self%nv))

!Place list into the 
read(varlist,*) self%fldnames

end subroutine fv3jedi_vars_create

! ------------------------------------------------------------------------------

end module fv3jedi_vars_mod
