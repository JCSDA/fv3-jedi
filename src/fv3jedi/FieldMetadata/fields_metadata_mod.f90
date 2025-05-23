! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fields_metadata_mod

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_char, c_size_t, &
       & c_int, c_int32_t, c_int64_t, c_float, c_double, c_bool, c_null_char

use fckit_C_interop_module, only: c_ptr_to_string

use string_f_c_mod, only: f_c_string, c_f_string

implicit none

private
public fields_metadata, field_metadata

integer, public, parameter :: clen = 2048 ! If this changes, change below and in
                                          ! FieldsMetadata.interface.cc

type fields_metadata
 private
 type(c_ptr) :: ptr = c_null_ptr
 contains
   procedure, public :: get_field_metadata
end type fields_metadata

type field_metadata
 character(len=clen) :: long_name
 character(len=clen) :: units
 character(len=clen) :: kind
 logical :: tracer
 integer :: levels
 character(len=clen) :: space
end type field_metadata

interface fields_metadata
  module procedure create
end interface

! --------------------------------------------------------------------------------------------------

interface
  subroutine c_get_field_metadata(ptr, long_name, units, kindd, tracer, levels, space) &
                                  bind(c, name='get_field_metadata_f')
    use iso_c_binding
    integer, parameter :: clen = 2048
    type(c_ptr), value :: ptr
    character(len=1, kind=c_char), intent(in) :: long_name(clen)
    character(len=1, kind=c_char) :: units(clen)
    character(len=1, kind=c_char) :: kindd(clen)
    logical(c_bool) :: tracer
    integer(kind=c_int) :: levels
    character(len=1, kind=c_char) :: space(clen)
  end subroutine c_get_field_metadata
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

function get_field_metadata(self, long_name_in) result(fmd)

class(fields_metadata), intent(in) :: self
character(len=*),       intent(in) :: long_name_in
type(field_metadata) :: fmd

! Long name that is sent to c++
character(len=1, kind=c_char), allocatable :: long_name(:)

! Returned from c++
character(len=1, kind=c_char), allocatable :: units(:)
character(len=1, kind=c_char), allocatable :: kindd(:)
logical(c_bool) :: tracer
integer(kind=c_int) :: levels
character(len=1, kind=c_char), allocatable :: space(:)

integer :: n, longlen

! long_name that can be passed to c++
longlen = len(trim(long_name_in))
allocate(long_name(clen))
long_name = ''
do n = 1,longlen
  long_name(n) = long_name_in(n:n)
enddo
long_name(longlen+1) = c_null_char

! Allocate string outputs
allocate(units(clen))
allocate(kindd(clen))
allocate(space(clen))

units = c_null_char
kindd = c_null_char
space = c_null_char

! Get information from C++ object
call c_get_field_metadata(self%ptr, long_name, units, kindd, tracer, levels, space)

! Copy non string
fmd%tracer = tracer
fmd%levels = levels

! Copy string
fmd%long_name = long_name_in
call c_f_string(units, fmd%units)
call c_f_string(kindd, fmd%kind)
call c_f_string(space, fmd%space)

end function get_field_metadata

! --------------------------------------------------------------------------------------------------

end module fields_metadata_mod
