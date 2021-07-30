! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fields_metadata_mod

use iso_c_binding

use fckit_C_interop_module, only: c_ptr_to_string

use string_f_c_mod, only: f_c_string

implicit none

private
public fields_metadata, field_metadata

integer, parameter :: clen = 2048 ! If this changes, change below and in FieldsMetadata.interface.cc

type fields_metadata
 private
 type(c_ptr) :: ptr = c_null_ptr
 contains
   procedure, public :: get_field
end type fields_metadata

type field_metadata
 character(len=clen) :: field_io_name
 character(len=clen) :: field_name
 character(len=clen) :: array_kind
 integer :: levels
 character(len=clen) :: long_name
 character(len=clen) :: space
 character(len=clen) :: stagger_loc
 logical :: tracer
 character(len=clen) :: units
 character(len=clen) :: interp_type
 character(len=clen) :: io_file
end type field_metadata

interface fields_metadata
  module procedure create
end interface

! --------------------------------------------------------------------------------------------------

interface
  subroutine c_fields_metadata_get_field(ptr, field_io_name, field_name, array_kind, levels, &
                                         long_name, space, stagger_loc, tracer, units, interp_type, &
                                         io_file) bind(c, name='fields_metadata_get_field_f')
    use iso_c_binding
    integer, parameter :: clen = 2048
    type(c_ptr), value :: ptr
    character(len=1, kind=c_char), intent(in) :: field_io_name(clen)
    character(len=1, kind=c_char), intent(inout) :: field_name(clen)
    character(len=1, kind=c_char), intent(inout) :: array_kind(clen)
    integer(kind=c_int),           intent(inout) :: levels
    character(len=1, kind=c_char), intent(inout) :: long_name(clen)
    character(len=1, kind=c_char), intent(inout) :: space(clen)
    character(len=1, kind=c_char), intent(inout) :: stagger_loc(clen)
    logical(c_bool),               intent(inout) :: tracer
    character(len=1, kind=c_char), intent(inout) :: units(clen)
    character(len=1, kind=c_char), intent(inout) :: interp_type(clen)
    character(len=1, kind=c_char), intent(inout) :: io_file(clen)
  end subroutine c_fields_metadata_get_field
end interface

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

function create(c_ptr_this) result(this)

type(c_ptr), value    :: c_ptr_this
type(fields_metadata) :: this

this%ptr = c_ptr_this

end function create

! --------------------------------------------------------------------------------------------------

function get_field(self, field_io_name_in) result(field)

class(fields_metadata), intent(in)    :: self
character(len=*),       intent(in)    :: field_io_name_in
type(field_metadata) :: field

character(len=1, kind=c_char), allocatable :: field_io_name(:)

! String pointers
character(len=1, kind=c_char), allocatable :: field_name(:)
character(len=1, kind=c_char), allocatable :: array_kind(:)
character(len=1, kind=c_char), allocatable :: long_name(:)
character(len=1, kind=c_char), allocatable :: space(:)
character(len=1, kind=c_char), allocatable :: stagger_loc(:)
character(len=1, kind=c_char), allocatable :: units(:)
character(len=1, kind=c_char), allocatable :: interp_type(:)
character(len=1, kind=c_char), allocatable :: io_file(:)

character(len=clen, kind=c_char) :: field_name_

integer(c_int)  :: levels
logical(c_bool) :: tracer

integer :: n, iolen

! field_io_name that can be passed to c++
iolen = len(trim(field_io_name_in))
allocate(field_io_name(clen))
field_io_name = ''
do n = 1,iolen
  field_io_name(n) = field_io_name_in(n:n)
enddo
field_io_name(iolen+1) = c_null_char

! Allocate outputs
allocate(field_name(clen))
allocate(array_kind(clen))
allocate(long_name(clen))
allocate(space(clen))
allocate(stagger_loc(clen))
allocate(units(clen))
allocate(interp_type(clen))
allocate(io_file(clen))

field_name = c_null_char
array_kind = c_null_char
long_name = c_null_char
space = c_null_char
stagger_loc = c_null_char
units = c_null_char
interp_type = c_null_char
io_file = c_null_char

! Get information from C++ object
call c_fields_metadata_get_field(self%ptr, field_io_name, field_name, array_kind, levels, &
                                 long_name, space, stagger_loc, tracer, units, interp_type, io_file )

! FieldNameIO
field%field_io_name = ''
do n = 1,clen
  if (field_io_name(n) == c_null_char) exit
  field%field_io_name(n:n) = field_io_name(n)
enddo

! FieldName
field%field_name = ''
do n = 1,clen
  if (field_name(n) == c_null_char) exit
  field%field_name(n:n) = field_name(n)
enddo

! ArrayKind
field%array_kind = ''
do n = 1,clen
  if (array_kind(n) == c_null_char) exit
  field%array_kind(n:n) = array_kind(n)
enddo

! Levels
field%levels = levels

! LongName
field%long_name = ''
do n = 1,clen
  if (long_name(n) == c_null_char) exit
  field%long_name(n:n) = long_name(n)
enddo

! Space
field%space = ''
do n = 1,clen
  if (space(n) == c_null_char) exit
  field%space(n:n) = space(n)
enddo

! StaggerLoc
field%stagger_loc = ''
do n = 1,clen
  if (stagger_loc(n) == c_null_char) exit
  field%stagger_loc(n:n) = stagger_loc(n)
enddo

! Tracer
field%tracer = tracer

! Units
field%units = ''
do n = 1,clen
  if (units(n) == c_null_char) exit
  field%units(n:n) = units(n)
enddo

! Interpolation Type
field%interp_type = ''
do n = 1,clen
  if (interp_type(n) == c_null_char) exit
  field%interp_type(n:n) = interp_type(n)
enddo

! IO file name
field%io_file = ''
do n = 1,clen
  if (io_file(n) == c_null_char) exit
  field%io_file(n:n) = io_file(n)
enddo

end function get_field

! --------------------------------------------------------------------------------------------------

end module fields_metadata_mod
