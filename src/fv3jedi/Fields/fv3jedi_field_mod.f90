! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_field_mod

! fckit
use fckit_mpi_module

! oops
use datetime_mod
use oops_variables_mod

! fv3jedi
use fv3jedi_kinds_mod, only: kind_real

implicit none
private
public :: fv3jedi_field, has_field, get_field, put_field, checksame, copy_subset, &
          long_name_to_fv3jedi_name

! These are the same methods as used in fv3jedi_fields but with argument being a list of individual
! fields instead of the fv3jedi_fields class
interface get_field
  module procedure get_field_return_type_pointer
  module procedure get_field_return_array_pointer
  module procedure get_field_return_array_allocatable
end interface

integer, parameter, public :: field_clen = 2048

! --------------------------------------------------------------------------------------------------

!Field type (individual field)
type :: fv3jedi_field
 logical :: lalloc = .false.
 character(len=field_clen) :: short_name = "null"   !Short name (to match file name)
 character(len=field_clen) :: fv3jedi_name = "null" !Common name
 character(len=field_clen) :: long_name = "null"    !More descriptive name
 character(len=field_clen) :: units = "null"        !Units for the field
 character(len=field_clen) :: io_file = "null"      !Which restart to read/write if not the default
 logical                   :: tracer = .false.      !Whether field is classified as tracer (pos. def.)
 character(len=field_clen) :: space                 !One of vector, magnitude, direction
 character(len=field_clen) :: staggerloc            !One of center, eastwest, northsouth, corner
 integer :: isc, iec, jsc, jec, npz
 real(kind=kind_real), allocatable :: array(:,:,:)
 logical :: integerfield = .false.
endtype fv3jedi_field

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

logical function has_field(fields, field_name, field_index)

type(fv3jedi_field), intent(in)  :: fields(:)
character(len=*),    intent(in)  :: field_name
integer, optional,   intent(out) :: field_index

integer :: var

has_field = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(field_name)) then
    has_field = .true.
    if (present(field_index)) field_index = var
    exit
  endif
enddo

end function has_field

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_type_pointer(fields, field_name, field)

type(fv3jedi_field), target,  intent(in)  :: fields(:)
character(len=*),             intent(in)  :: field_name
type(fv3jedi_field), pointer, intent(out) :: field

integer :: var
logical :: found

if(associated(field)) nullify(field)

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(field_name)) then
    field => fields(var)
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.get_field_return_type_pointer: field "&
                                //trim(field_name)//" not found in fields")

end subroutine get_field_return_type_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_pointer(fields, field_name, field)

type(fv3jedi_field), target,   intent(in)    :: fields(:)
character(len=*),              intent(in)    :: field_name
real(kind=kind_real), pointer, intent(out)   :: field(:,:,:)

integer :: var
logical :: found

if(associated(field)) nullify(field)

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(field_name)) then
    field => fields(var)%array
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.get_field_return_array_pointer: field "&
                                //trim(field_name)//" not found in fields")

end subroutine get_field_return_array_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_allocatable(fields, field_name, field)

type(fv3jedi_field),               intent(in)  :: fields(:)
character(len=*),                  intent(in)  :: field_name
real(kind=kind_real), allocatable, intent(out) :: field(:,:,:)

integer :: var
logical :: found, boundsmatch

found = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(field_name)) then

    if (.not. allocated(field)) then
      ! If not allocated allocate
    else
      ! If allocated check bounds
      boundsmatch = lbound(field,1) == fields(var)%isc .and. ubound(field,1) == fields(var)%iec .and. &
                    lbound(field,2) == fields(var)%jsc .and. ubound(field,2) == fields(var)%jec .and. &
                    lbound(field,3) == 1               .and. ubound(field,3) == fields(var)%npz
      if (.not.boundsmatch) call abor1_ftn("get_field_return_array_allocatable: field "//&
                                           trim(field_name)//" bounds mismatch")
    endif

    ! Copy the field
    field = fields(var)%array

    ! Set found flag
    found = .true.

    exit
  endif
enddo

if (.not.found) call abor1_ftn("get_field_return_array_allocatable: field "//trim(field_name)//&
                               " not found in fields")

end subroutine get_field_return_array_allocatable


! --------------------------------------------------------------------------------------------------

subroutine put_field(fields, field_name, field)

type(fv3jedi_field),               intent(inout) :: fields(:)
character(len=*),                  intent(in)    :: field_name
real(kind=kind_real), allocatable, intent(in)    :: field(:,:,:)

integer :: var
logical :: found, boundsmatch

if (.not. allocated(field)) call abor1_ftn("put_field: field "//trim(field_name)//" not allocated")

found = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(field_name)) then

    ! Check for matching bounds
    boundsmatch = lbound(field,1) == fields(var)%isc .and. ubound(field,1) == fields(var)%iec .and. &
                  lbound(field,2) == fields(var)%jsc .and. ubound(field,2) == fields(var)%jec .and. &
                  lbound(field,3) == 1               .and. ubound(field,3) == fields(var)%npz
    if (.not.boundsmatch) call abor1_ftn("put_field: field "//trim(field_name)//" bounds mismatch")

    ! Copy the field
    fields(var)%array = field

    ! Set found flag
    found = .true.

    exit
  endif
enddo

if (.not.found) call abor1_ftn("put_field: field "//trim(field_name)//" not found in fields")

end subroutine put_field

! --------------------------------------------------------------------------------------------------

subroutine checksame(fields1, fields2, calling_method)

implicit none
type(fv3jedi_field), intent(in) :: fields1(:)
type(fv3jedi_field), intent(in) :: fields2(:)
character(len=*),    intent(in) :: calling_method

integer :: var

if (size(fields1) .ne. size(fields2)) then
  call abor1_ftn(trim(calling_method)//"(checksame): Different number of fields")
endif

do var = 1,size(fields1)
  if (fields1(var)%fv3jedi_name .ne. fields2(var)%fv3jedi_name) then
      call abor1_ftn(trim(calling_method)//"(checksame): field "//trim(fields1(var)%fv3jedi_name)//&
                     " not in the equivalent position in the right hand side")
  endif
enddo

end subroutine checksame

! ------------------------------------------------------------------------------

subroutine copy_subset(field_in, field_ou, not_copied)

implicit none
type(fv3jedi_field),                             intent(in)    :: field_in(:)
type(fv3jedi_field),                             intent(inout) :: field_ou(:)
character(len=field_clen), allocatable, optional, intent(out)   :: not_copied(:)

integer :: var
character(len=field_clen) :: not_copied_(size(field_ou))
integer :: num_not_copied

! Loop over fields and copy if existing in both
num_not_copied = 0
do var = 1, size(field_ou)
  if (has_field(field_in, field_ou(var)%fv3jedi_name )) then
    call get_field(field_in, field_ou(var)%fv3jedi_name, field_ou(var)%array)
  else
    num_not_copied = num_not_copied + 1
    not_copied_(num_not_copied) = field_ou(var)%fv3jedi_name
  endif
enddo

! Send back list of variables not retrivable from field_in
if (present(not_copied) .and. num_not_copied > 0) then
  allocate(not_copied(num_not_copied))
  not_copied(1:num_not_copied) = not_copied_(1:num_not_copied)
endif

end subroutine copy_subset

! --------------------------------------------------------------------------------------------------

subroutine long_name_to_fv3jedi_name(fields, long_name, fv3jedi_name)

type(fv3jedi_field), intent(in)  :: fields(:)
character(len=*),    intent(in)  :: long_name
character(len=*),    intent(out) :: fv3jedi_name

integer :: n

do n = 1, size(fields)
  if (trim(long_name) == trim(fields(n)%long_name)) then
    fv3jedi_name = trim(fields(n)%fv3jedi_name)
    return
  endif
enddo

! Try with increment_of_ prepended to long_name
do n = 1, size(fields)
  if ("increment_of_"//trim(long_name) == trim(fields(n)%long_name)) then
    fv3jedi_name = trim(fields(n)%fv3jedi_name)
    return
  endif
enddo

call abor1_ftn("fv3jedi_field_mod.long_name_to_fv3jedi_name long_name "//trim(long_name)//&
               " not found in fields.")

end subroutine long_name_to_fv3jedi_name

! --------------------------------------------------------------------------------------------------

end module fv3jedi_field_mod
