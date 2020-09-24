! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_field_mod

use fckit_mpi_module
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

implicit none

private
public :: fv3jedi_field, &
          has_field, &
          pointer_field, &
          pointer_field_array, &
          copy_field_array, &
          allocate_copy_field_array, &
          fields_rms, &
          field_getminmaxrms, &
          checksame, &
          copy_subset, &
          long_name_to_fv3jedi_name, &
          zero_fields

integer, parameter, public :: field_clen = 2048

!Field type
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
 contains
  procedure :: allocate_field
  procedure :: equals
  generic :: assignment(=) => equals
  procedure :: deallocate_field
endtype fv3jedi_field

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine allocate_field(self,isc,iec,jsc,jec,npz,short_name,long_name,&
                          fv3jedi_name,units,io_file,space,staggerloc,tracer,integerfield)

implicit none
class(fv3jedi_field), target,  intent(inout) :: self
integer,                       intent(in)    :: isc,iec,jsc,jec,npz
character(len=*),              intent(in)    :: short_name
character(len=*),              intent(in)    :: long_name
character(len=*),              intent(in)    :: fv3jedi_name
character(len=*),              intent(in)    :: units
character(len=*),              intent(in)    :: io_file
character(len=*),              intent(in)    :: space
character(len=*),              intent(in)    :: staggerloc
logical,                       intent(in)    :: tracer
logical,                       intent(in)    :: integerfield

self%isc = isc
self%iec = iec
self%jsc = jsc
self%jec = jec
self%npz = npz

if (len(trim(fv3jedi_name)) > field_clen) &
  call abor1_ftn("fv3jedi_field_mod.allocate_field: fv3jedi_name " //trim(fv3jedi_name)// &
                 " should not be longer than 2048 characters")

if(.not.self%lalloc) then

  if (trim(staggerloc) == 'center') then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  elseif (trim(staggerloc) == 'northsouth') then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
  elseif (trim(staggerloc) == 'eastwest') then
    allocate(self%array(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
  endif

endif

self%lalloc = .true.

self%short_name   = trim(short_name)
self%long_name    = trim(long_name)
self%fv3jedi_name = trim(fv3jedi_name)
self%units        = trim(units)
self%io_file      = trim(io_file)
self%space        = space
self%staggerloc   = staggerloc
self%tracer       = tracer
self%integerfield = integerfield

end subroutine allocate_field

! --------------------------------------------------------------------------------------------------

subroutine deallocate_field(self)

implicit none
class(fv3jedi_field), intent(inout) :: self

if(self%lalloc) deallocate(self%array)
self%lalloc = .false.

end subroutine deallocate_field

! --------------------------------------------------------------------------------------------------

subroutine equals(self,rhs)

implicit none
class(fv3jedi_field), intent(inout) :: self
type (fv3jedi_field), intent(in)    :: rhs

call self%allocate_field( rhs%isc,rhs%iec,rhs%jsc,rhs%jec,rhs%npz, &
                          short_name   = rhs%short_name, &
                          long_name    = rhs%long_name, &
                          fv3jedi_name = rhs%fv3jedi_name, &
                          units        = rhs%units, &
                          io_file      = rhs%io_file, &
                          space        = rhs%space, &
                          staggerloc   = rhs%staggerloc, &
                          tracer       = rhs%tracer, &
                          integerfield = rhs%integerfield)

self%array = rhs%array

end subroutine equals

! --------------------------------------------------------------------------------------------------

logical function has_field(fields, fv3jedi_name, index)

type(fv3jedi_field), target,  intent(in)  :: fields(:)
character(len=*),             intent(in)  :: fv3jedi_name
integer, optional,            intent(out) :: index

integer :: var

has_field = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    has_field = .true.
    if (present(index)) index = var
    exit
  endif
enddo

end function has_field

! --------------------------------------------------------------------------------------------------

subroutine allocate_copy_field_array(fields, fv3jedi_name, field_array)

type(fv3jedi_field),               intent(in)  :: fields(:)
character(len=*),                  intent(in)  :: fv3jedi_name
real(kind=kind_real), allocatable, intent(out) :: field_array(:,:,:)

integer :: var
logical :: found

if(allocated(field_array)) deallocate(field_array)

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    allocate(field_array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,fields(var)%npz))
    field_array = fields(var)%array
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.allocate_copy_field_array: field "&
                                //trim(fv3jedi_name)//" not found in fields")

end subroutine allocate_copy_field_array

! --------------------------------------------------------------------------------------------------

subroutine copy_field_array(fields, fv3jedi_name, field_array)

type(fv3jedi_field),  intent(in)  :: fields(:)
character(len=*),     intent(in)  :: fv3jedi_name
real(kind=kind_real), intent(out) :: field_array(:,:,:)

integer :: var
logical :: found

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    field_array = fields(var)%array
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.copy_field_array: field "&
                                //trim(fv3jedi_name)//" not found in fields")

end subroutine copy_field_array

! --------------------------------------------------------------------------------------------------

subroutine pointer_field(fields, fv3jedi_name, field_pointer)

type(fv3jedi_field), target,  intent(in)  :: fields(:)
character(len=*),             intent(in)  :: fv3jedi_name
type(fv3jedi_field), pointer, intent(out) :: field_pointer

integer :: var
logical :: found

if(associated(field_pointer)) nullify(field_pointer)

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    field_pointer => fields(var)
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.pointer_field: field "&
                                //trim(fv3jedi_name)//" not found in fields")

end subroutine pointer_field

! --------------------------------------------------------------------------------------------------

subroutine pointer_field_array(fields, fv3jedi_name, array_pointer)

type(fv3jedi_field), target,   intent(in)    :: fields(:)
character(len=*),              intent(in)    :: fv3jedi_name
real(kind=kind_real), pointer, intent(out)   :: array_pointer(:,:,:)

integer :: var
logical :: found

if(associated(array_pointer)) nullify(array_pointer)

found = .false.
do var = 1,size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    array_pointer => fields(var)%array
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field.pointer_field_array: field "&
                                //trim(fv3jedi_name)//" not found in fields")

end subroutine pointer_field_array

! --------------------------------------------------------------------------------------------------

subroutine fields_rms(nf,fields,rms, f_comm)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
real(kind=kind_real), intent(inout) :: rms
type(fckit_mpi_comm), intent(in)    :: f_comm

integer :: i, j, k, ii, iisum, var
real(kind=kind_real) :: zz

zz = 0.0_kind_real
ii = 0

do var = 1,nf

  do k = 1,fields(var)%npz
    do j = fields(var)%jsc,fields(var)%jec
      do i = fields(var)%isc,fields(var)%iec
        zz = zz + fields(var)%array(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo

enddo

!Get global values
call f_comm%allreduce(zz,rms,fckit_mpi_sum())
call f_comm%allreduce(ii,iisum,fckit_mpi_sum())
rms = sqrt(rms/real(iisum,kind_real))

end subroutine fields_rms

! --------------------------------------------------------------------------------------------------

subroutine field_getminmaxrms(fields, field_num, field_name, minmaxrms, f_comm)

implicit none
type(fv3jedi_field),  intent(in)    :: fields(:)
integer,              intent(in)    :: field_num
character(len=*),     intent(inout) :: field_name
real(kind=kind_real), intent(inout) :: minmaxrms(3)
type(fckit_mpi_comm), intent(in)    :: f_comm

integer :: isc, iec, jsc, jec, npz
real(kind=kind_real) :: tmp(3), gs3, gs3g

isc = fields(field_num)%isc
iec = fields(field_num)%iec
jsc = fields(field_num)%jsc
jec = fields(field_num)%jec
npz = fields(field_num)%npz

field_name = fields(field_num)%short_name

! Compute global sum over the field
gs3 = real((iec-isc+1)*(jec-jsc+1)*npz, kind_real)
call f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

! Min/Max/SumSquares
tmp(1) = minval(fields(field_num)%array(isc:iec,jsc:jec,1:npz))
tmp(2) = maxval(fields(field_num)%array(isc:iec,jsc:jec,1:npz))
tmp(3) =    sum(fields(field_num)%array(isc:iec,jsc:jec,1:npz)**2)

! Get global min/max/sum
call f_comm%allreduce(tmp(1), minmaxrms(1), fckit_mpi_min())
call f_comm%allreduce(tmp(2), minmaxrms(2), fckit_mpi_max())
call f_comm%allreduce(tmp(3), minmaxrms(3), fckit_mpi_sum())

! SumSquares to rms
minmaxrms(3) = sqrt(minmaxrms(3)/gs3g)

end subroutine field_getminmaxrms

! --------------------------------------------------------------------------------------------------

subroutine checksame(self,other,method)

implicit none
type(fv3jedi_field), intent(in) :: self(:)
type(fv3jedi_field), intent(in) :: other(:)
character(len=*),    intent(in) :: method

integer :: var

if (size(self) .ne. size(other)) then
  call abor1_ftn(trim(method)//"(checksame): Different number of fields")
endif

do var = 1,size(self)
  if (self(var)%fv3jedi_name .ne. other(var)%fv3jedi_name) then
      call abor1_ftn(trim(method)//"(checksame): field "//trim(self(var)%fv3jedi_name)//&
                     " not in the equivalent position in the right hand side")
  endif
enddo

end subroutine checksame

! ------------------------------------------------------------------------------

subroutine copy_subset(field_in,field_ou,not_copied)

implicit none
type(fv3jedi_field),                             intent(in)    :: field_in(:)
type(fv3jedi_field),                             intent(inout) :: field_ou(:)
character(len=field_clen), allocatable, optional, intent(out)   :: not_copied(:)

integer :: var
character(len=field_clen) :: not_copied_(10000)
integer :: num_not_copied

! Loop over fields and copy if existing in both
num_not_copied = 0
do var = 1, size(field_ou)
  if (has_field(field_in, field_ou(var)%fv3jedi_name )) then
    call copy_field_array(field_in, field_ou(var)%fv3jedi_name, field_ou(var)%array)
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

subroutine zero_fields(fields)

type(fv3jedi_field), intent(inout) :: fields(:)

integer :: var

do var = 1, size(fields)
  fields(var)%array = 0.0_kind_real
enddo

end subroutine zero_fields

! --------------------------------------------------------------------------------------------------

end module fv3jedi_field_mod
