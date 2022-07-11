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
use fv3jedi_kinds_mod,   only: kind_real
use fields_metadata_mod, only: field_metadata

implicit none
private
public :: fv3jedi_field, hasfield, get_field, put_field, checksame, copy_subset, create_field

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
 character(len=field_clen) :: long_name                       ! Field long name
 character(len=field_clen) :: short_name                      ! Field short name
 character(len=field_clen) :: units                           ! Field units
 character(len=field_clen) :: kind                            ! Data kind, real, integer etc (always allocate real data)
 logical                   :: tracer                          ! Whether field is tracer or not
 character(len=field_clen) :: horizontal_stagger_location     ! Stagger location in horizontal
 character(len=field_clen) :: space                           ! Vector, magnitude, direction
 character(len=field_clen) :: io_name                         ! Name used for IO
 character(len=field_clen) :: io_file                         ! File used for IO
 character(len=field_clen) :: interpolation_type              ! Type of interpolation to use
 integer :: isc, iec, jsc, jec, npz
 real(kind=kind_real), allocatable :: array(:,:,:)
 type(fckit_mpi_comm) :: comm                       ! Communicator
endtype fv3jedi_field

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create_field(self, fmd, comm)

type(fv3jedi_field),  intent(inout) :: self
type(field_metadata), intent(in)    :: fmd
type(fckit_mpi_comm), intent(in)    :: comm

! Check that the names in the field meta data are not longer than expected
! ------------------------------------------------------------------------
if (len(trim(fmd%long_name)) > field_clen) &
  call abor1_ftn("fv3jedi_field.create: " //trim(fmd%long_name)// " too long")
if (len(trim(fmd%short_name)) > field_clen) &
  call abor1_ftn("fv3jedi_field.create: " //trim(fmd%short_name)// " too long")
if (len(trim(fmd%io_name)) > field_clen) &
  call abor1_ftn("fv3jedi_field.create: " //trim(fmd%io_name)// " too long")

! Copy metadata
! -------------
self%long_name = fmd%long_name
self%short_name = fmd%short_name
self%units = fmd%units
self%kind = fmd%kind
self%tracer = fmd%tracer
self%horizontal_stagger_location = fmd%horizontal_stagger_location
self%npz = fmd%levels
self%space = fmd%space
self%io_name = fmd%io_name
self%io_file = fmd%io_file
self%interpolation_type = fmd%interpolation_type

! Allocate the field array data
! -----------------------------
if(.not.self%lalloc) then

  if (trim(self%horizontal_stagger_location) == 'center') then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  elseif (trim(self%horizontal_stagger_location) == 'northsouth') then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
  elseif (trim(self%horizontal_stagger_location) == 'eastwest') then
    allocate(self%array(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
  endif

  ! Initialize to zero and set allocated
  self%array = 0.0_kind_real
  self%lalloc = .true.

endif

! Ensure the interpolation type is consistent with other metadata
! ---------------------------------------------------------------
if (self%interpolation_type == 'default') then
  if ( trim(self%kind)=='integer' ) then
    self%interpolation_type  = 'integer'
  elseif ( trim(self%space) == 'direction' ) then
    self%interpolation_type  = 'nearest'
  endif
endif

! Communicator
! ------------
self%comm = comm

end subroutine create_field

! --------------------------------------------------------------------------------------------------

logical function hasfield(fields, field_name, field_index)

type(fv3jedi_field), intent(in)  :: fields(:)
character(len=*),    intent(in)  :: field_name
integer, optional,   intent(out) :: field_index

integer :: var

hasfield = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%short_name) == trim(field_name) .or. &
       trim(fields(var)%long_name) == trim(field_name) ) then
    hasfield = .true.
    if (present(field_index)) field_index = var
    exit
  endif
enddo

end function hasfield

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
  if ( trim(fields(var)%short_name) == trim(field_name) .or. &
       trim(fields(var)%long_name) == trim(field_name) ) then
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
  if ( trim(fields(var)%short_name) == trim(field_name) .or. &
       trim(fields(var)%long_name) == trim(field_name) ) then
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
  if ( trim(fields(var)%short_name) == trim(field_name) .or. &
       trim(fields(var)%long_name) == trim(field_name) ) then

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
  if ( trim(fields(var)%short_name) == trim(field_name) .or. &
       trim(fields(var)%long_name) == trim(field_name) ) then

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

subroutine print_fields_debug(fields1, fields2)

type(fv3jedi_field), intent(in) :: fields1(:)
type(fv3jedi_field), intent(in) :: fields2(:)

integer :: var

! Only root prints the fields
if (fields1(1)%comm%rank() == 0) then

  print*, 'Number of fields: ', size(fields1), "and", size(fields2)

  ! Print list of fields in fields1
  print*, "List of fields in fields1:"
  do var = 1,size(fields1)
    print*, trim(fields1(var)%short_name)
  enddo

  ! Print list of fields in fields2
  print*, "List of fields in fields2:"
  do var = 1,size(fields2)
    print*, 'fields2:', var, trim(fields2(var)%short_name)
  enddo

endif

end subroutine print_fields_debug

! --------------------------------------------------------------------------------------------------

subroutine checksame(fields1, fields2, calling_method)

implicit none
type(fv3jedi_field), intent(in) :: fields1(:)
type(fv3jedi_field), intent(in) :: fields2(:)
character(len=*),    intent(in) :: calling_method

integer :: var

if (size(fields1) .ne. size(fields2)) then
  if (fields1(1)%comm%rank() == 0) print*, 'fv3jedi.fields checksame different number of fields'
  call print_fields_debug(fields1, fields2)
  call abor1_ftn(trim(calling_method)//"(checksame): Different number of fields")
endif

do var = 1,size(fields1)
  if (fields1(var)%short_name .ne. fields2(var)%short_name) then
    if (fields1(1)%comm%rank() == 0) print*, 'fv3jedi.fields checksame positional differences'
    call print_fields_debug(fields1, fields2)
    call abor1_ftn(trim(calling_method)//"(checksame): field "//trim(fields1(var)%short_name)//&
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
  if (hasfield(field_in, field_ou(var)%short_name )) then
    call get_field(field_in, field_ou(var)%short_name, field_ou(var)%array)
  else
    num_not_copied = num_not_copied + 1
    not_copied_(num_not_copied) = field_ou(var)%short_name
  endif
enddo

! Send back list of variables not retrivable from field_in
if (present(not_copied) .and. num_not_copied > 0) then
  allocate(not_copied(num_not_copied))
  not_copied(1:num_not_copied) = not_copied_(1:num_not_copied)
endif

end subroutine copy_subset

! --------------------------------------------------------------------------------------------------

end module fv3jedi_field_mod
