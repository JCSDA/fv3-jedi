! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fields_metadata_mod

use iso_c_binding

use fckit_C_interop_module, only: c_ptr_to_string

use string_f_c_mod, only: f_c_string, c_f_string

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
 character(len=clen) :: long_name
 character(len=clen) :: short_name
 character(len=clen) :: units
 character(len=clen) :: kind
 logical :: tracer
 character(len=clen) :: horizontal_stagger_location
 integer :: levels
 character(len=clen) :: space
 character(len=clen) :: io_name
 character(len=clen) :: io_file
 character(len=clen) :: interpolation_type
 character(len=clen) :: interpolation_source_point_mask
end type field_metadata

interface fields_metadata
  module procedure create
end interface

! --------------------------------------------------------------------------------------------------

interface
  subroutine c_fields_metadata_get_field(ptr, longshortio_name, long_name, short_name, &
                                         units, kindd, tracer, horizontal_stagger_location, &
                                         levels, space, io_name, io_file, interpolation_type, &
                                         interpolation_source_point_mask) &
                                         bind(c, name='fields_metadata_get_field_f')
    use iso_c_binding
    integer, parameter :: clen = 2048
    type(c_ptr), value :: ptr
    character(len=1, kind=c_char), intent(in) :: longshortio_name(clen)
    character(len=1, kind=c_char) :: long_name(clen)
    character(len=1, kind=c_char) :: short_name(clen)
    character(len=1, kind=c_char) :: units(clen)
    character(len=1, kind=c_char) :: kindd(clen)
    logical(c_bool) :: tracer
    character(len=1, kind=c_char) :: horizontal_stagger_location(clen)
    integer(kind=c_int) :: levels
    character(len=1, kind=c_char) :: space(clen)
    character(len=1, kind=c_char) :: io_name(clen)
    character(len=1, kind=c_char) :: io_file(clen)
    character(len=1, kind=c_char) :: interpolation_type(clen)
    character(len=1, kind=c_char) :: interpolation_source_point_mask(clen)
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

function get_field(self, longshortio_name_in) result(field)

class(fields_metadata), intent(in) :: self
character(len=*),       intent(in) :: longshortio_name_in
type(field_metadata) :: field

character(len=1, kind=c_char), allocatable :: longshortio_name(:)

! Returned from c++
character(len=1, kind=c_char), allocatable :: long_name(:)
character(len=1, kind=c_char), allocatable :: short_name(:)
character(len=1, kind=c_char), allocatable :: units(:)
character(len=1, kind=c_char), allocatable :: kindd(:)
logical(c_bool) :: tracer
character(len=1, kind=c_char), allocatable :: horizontal_stagger_location(:)
integer(kind=c_int) :: levels
character(len=1, kind=c_char), allocatable :: space(:)
character(len=1, kind=c_char), allocatable :: io_name(:)
character(len=1, kind=c_char), allocatable :: io_file(:)
character(len=1, kind=c_char), allocatable :: interpolation_type(:)
character(len=1, kind=c_char), allocatable :: interpolation_source_point_mask(:)

integer :: n, longshortiolen

! longshortio_name that can be passed to c++
longshortiolen = len(trim(longshortio_name_in))
allocate(longshortio_name(clen))
longshortio_name = ''
do n = 1,longshortiolen
  longshortio_name(n) = longshortio_name_in(n:n)
enddo
longshortio_name(longshortiolen+1) = c_null_char

! Allocate string outputs
allocate(long_name(clen))
allocate(short_name(clen))
allocate(units(clen))
allocate(kindd(clen))
allocate(horizontal_stagger_location(clen))
allocate(space(clen))
allocate(io_name(clen))
allocate(io_file(clen))
allocate(interpolation_type(clen))
allocate(interpolation_source_point_mask(clen))

long_name = c_null_char
short_name = c_null_char
units = c_null_char
kindd = c_null_char
horizontal_stagger_location = c_null_char
space = c_null_char
io_name = c_null_char
io_file = c_null_char
interpolation_type = c_null_char
interpolation_source_point_mask = c_null_char

! Get information from C++ object
call c_fields_metadata_get_field(self%ptr, longshortio_name, long_name, short_name, units, kindd, &
                                 tracer, horizontal_stagger_location, levels, space, io_name, &
                                 io_file, interpolation_type, interpolation_source_point_mask)

! Copy non string
field%tracer = tracer
field%levels = levels

! Copy string
call c_f_string(long_name, field%long_name)
call c_f_string(short_name, field%short_name)
call c_f_string(units, field%units)
call c_f_string(kindd, field%kind)
call c_f_string(horizontal_stagger_location, field%horizontal_stagger_location)
call c_f_string(space, field%space)
call c_f_string(io_name, field%io_name)
call c_f_string(io_file, field%io_file)
call c_f_string(interpolation_type, field%interpolation_type)
call c_f_string(interpolation_source_point_mask, field%interpolation_source_point_mask)

end function get_field

! --------------------------------------------------------------------------------------------------

end module fields_metadata_mod
