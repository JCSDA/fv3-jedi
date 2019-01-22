! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_field_mod

use fv3jedi_kinds_mod, only: kind_real

implicit none

private
public :: fv3jedi_field

!Field type
type :: fv3jedi_field
 logical :: lalloc = .false.
 character(len=32) :: short_name = "null"   !Short name (to match file name)
 character(len=64) :: fv3jedi_name = "null" !Common name
 character(len=64) :: long_name = "null"    !More descriptive name
 character(len=32) :: units = "null"        !Units for the field
 integer :: staggerloc   !Middle, corners, east, south, etc
 integer :: isc, iec, jsc, jec, npz
 real(kind=kind_real), allocatable :: field(:,:,:)
 contains
  procedure :: allocate_field
  procedure :: equals
  generic :: assignment(=) => equals
  procedure :: deallocate_field
endtype fv3jedi_field

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine allocate_field(self,isc,iec,jsc,jec,npz,short_name,long_name,&
                          fv3jedi_name,units,staggerloc)

implicit none
class(fv3jedi_field), intent(inout) :: self
integer, intent(in) :: isc,iec,jsc,jec,npz
character(len=*), optional, intent(in) :: short_name
character(len=*), optional, intent(in) :: long_name
character(len=*), optional, intent(in) :: fv3jedi_name
character(len=*), optional, intent(in) :: units
integer, optional, intent(in) :: staggerloc 

self%isc = isc
self%iec = iec
self%jsc = jsc
self%jec = jec
self%npz = npz

if(.not.self%lalloc) &
allocate(self%field(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
self%lalloc = .true.

if (present(short_name)) self%short_name = trim(short_name)
if (present(long_name)) self%long_name = trim(long_name)
if (present(fv3jedi_name)) self%fv3jedi_name = trim(fv3jedi_name)
if (present(units)) self%units = trim(units)
if (present(staggerloc)) self%staggerloc = staggerloc

end subroutine allocate_field

! ------------------------------------------------------------------------------

subroutine deallocate_field(self)

implicit none
class(fv3jedi_field), intent(inout) :: self

if(self%lalloc) deallocate(self%field)
self%lalloc = .false.

end subroutine deallocate_field

! ------------------------------------------------------------------------------

subroutine equals(self,rhs)

implicit none
class(fv3jedi_field), intent(inout) :: self
type (fv3jedi_field), intent(in)    :: rhs

call self%allocate_field( rhs%isc,rhs%iec,rhs%jsc,rhs%jec,rhs%npz, &
                          short_name=rhs%short_name, &
                          long_name=rhs%long_name, &
                          fv3jedi_name=rhs%fv3jedi_name, &
                          units=rhs%units, &
                          staggerloc=rhs%staggerloc)

self%field = rhs%field

end subroutine equals

! ------------------------------------------------------------------------------

end module fv3jedi_field_mod
