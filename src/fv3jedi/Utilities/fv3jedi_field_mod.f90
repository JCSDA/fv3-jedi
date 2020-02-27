! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_field_mod

use fckit_mpi_module
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use mpp_domains_mod,   only: east, north, center

implicit none

private
public :: fv3jedi_field, &
          has_field, &
          pointer_field, &
          pointer_field_array, &
          copy_field_array, &
          allocate_copy_field_array, &
          fields_rms, &
          fields_gpnorm, &
          fields_print, &
          checksame, &
          flip_array_vertical, &
          copy_subset

!Field type
type :: fv3jedi_field
 logical :: lalloc = .false.
 character(len=32) :: short_name = "null"   !Short name (to match file name)
 character(len=10) :: fv3jedi_name = "null" !Common name
 character(len=64) :: long_name = "null"    !More descriptive name
 character(len=32) :: units = "null"        !Units for the field
 logical :: tracer = .false.
 integer :: staggerloc   !Middle, corners, east, south, etc
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
                          fv3jedi_name,units,staggerloc,tracer,integerfield)

implicit none
class(fv3jedi_field), target,  intent(inout) :: self
integer,                       intent(in)    :: isc,iec,jsc,jec,npz
character(len=*),              intent(in)    :: short_name
character(len=*),              intent(in)    :: long_name
character(len=*),              intent(in)    :: fv3jedi_name
character(len=*),              intent(in)    :: units
integer,                       intent(in)    :: staggerloc
logical, optional,             intent(in)    :: tracer
logical, optional,             intent(in)    :: integerfield

self%isc = isc
self%iec = iec
self%jsc = jsc
self%jec = jec
self%npz = npz

if (len(fv3jedi_name) > 10) &
call abor1_ftn("fv3jedi_field_mod.allocate_field: fv3jedi_name should not be longer than ten characters")

if(.not.self%lalloc) then

  if (staggerloc == center) then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  elseif (staggerloc == north) then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
  elseif (staggerloc == east) then
    allocate(self%array(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
  endif

endif

self%lalloc = .true.

self%short_name   = trim(short_name)
self%long_name    = trim(long_name)
self%fv3jedi_name = trim(fv3jedi_name)
self%units        = trim(units)
self%staggerloc   = staggerloc

if (present(tracer)) then
  self%tracer = tracer
else
  self%tracer = .false.
endif

! Fields that are e.g. types
if (present(integerfield)) then
  self%integerfield = integerfield
else
  self%integerfield = .false.
endif

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
                          short_name=rhs%short_name, &
                          long_name=rhs%long_name, &
                          fv3jedi_name=rhs%fv3jedi_name, &
                          units=rhs%units, &
                          staggerloc=rhs%staggerloc, &
                          tracer = rhs%tracer)

self%array = rhs%array

end subroutine equals

! --------------------------------------------------------------------------------------------------

logical function has_field(fields, fv3jedi_name)

type(fv3jedi_field), target,  intent(in)  :: fields(:)
character(len=*),             intent(in)  :: fv3jedi_name

integer :: var

has_field = .false.
do var = 1, size(fields)
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    has_field = .true.
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

subroutine fields_gpnorm(nf, fields, pstat, f_comm)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
real(kind=kind_real), intent(inout) :: pstat(3, nf)
type(fckit_mpi_comm), intent(in)    :: f_comm

integer :: var
real(kind=kind_real) :: tmp(3),  gs3, gs3g

do var = 1,nf

  gs3 = real((fields(var)%iec-fields(var)%isc+1)*(fields(var)%jec-fields(var)%jsc+1)*fields(var)%npz, kind_real)
  call f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

  tmp(1) = minval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
  tmp(2) = maxval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
  tmp(3) =    sum(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz)**2)

  call f_comm%allreduce(tmp(1),pstat(1,var),fckit_mpi_min())
  call f_comm%allreduce(tmp(2),pstat(2,var),fckit_mpi_max())
  call f_comm%allreduce(tmp(3),pstat(3,var),fckit_mpi_sum())
  pstat(3,var) = sqrt(pstat(3,var)/gs3g)

enddo

end subroutine fields_gpnorm

! --------------------------------------------------------------------------------------------------

subroutine fields_print(nf, fields, name, f_comm)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
character(len=*),     intent(in)    :: name
type(fckit_mpi_comm), intent(in)    :: f_comm

integer :: var
real(kind=kind_real) :: tmp(3), pstat(3), gs3, gs3g
character(len=34) :: printname

integer :: ngrid, sngrid

ngrid = (fields(1)%iec-fields(1)%isc+1)*(fields(1)%iec-fields(1)%isc+1)
call f_comm%allreduce(ngrid,sngrid,fckit_mpi_sum())
sngrid = nint(sqrt(real(sngrid,kind_real)/6.0_kind_real))

printname = "|     "//trim(name)//" print"

if (f_comm%rank() == 0) then
  write(*,"(A70)") "----------------------------------------------------------------------"
  write(*,"(A34)") printname
  write(*,"(A70)") "----------------------------------------------------------------------"
  write(*,"(A70)") " "
  write(*,"(A27,I5)") "    Cube sphere face size: ", sngrid
  write(*,"(A70)") " "
endif

do var = 1,nf

  gs3 = real((fields(var)%iec-fields(var)%isc+1)*(fields(var)%jec-fields(var)%jsc+1)*fields(var)%npz, kind_real)
  call f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

  tmp(1) = minval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
  tmp(2) = maxval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
  tmp(3) =    sum(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz)**2)

  call f_comm%allreduce(tmp(1),pstat(1),fckit_mpi_min())
  call f_comm%allreduce(tmp(2),pstat(2),fckit_mpi_max())
  call f_comm%allreduce(tmp(3),pstat(3),fckit_mpi_sum())
  pstat(3) = sqrt(pstat(3)/gs3g)

  if (f_comm%rank() == 0) write(*,"(A10,A6,ES14.7,A6,ES14.7,A6,ES14.7)") &
                                   trim(fields(var)%short_name),&
                                   "| Min=",real(pstat(1),4),&
                                   ", Max=",real(pstat(2),4),&
                                   ", RMS=",real(pstat(3),4)

enddo

if (f_comm%rank() == 0) &
  write(*,"(A70)") "---------------------------------------------------------------------"

end subroutine fields_print

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

! --------------------------------------------------------------------------------------------------

subroutine flip_array_vertical(nf,fields)

implicit none
integer,             intent(in)    :: nf
type(fv3jedi_field), intent(inout) :: fields(nf)

integer :: n, lev_in, lev_out
real(kind=kind_real), allocatable :: array_tmp(:,:,:)

do n = 1, nf

  if (fields(n)%npz > 1) then

    if (fields(n)%staggerloc == center) then
      allocate(array_tmp(fields(n)%isc:fields(n)%iec,fields(n)%jsc:fields(n)%jec,1:fields(n)%npz))
    elseif (fields(n)%staggerloc == north) then
      allocate(array_tmp(fields(n)%isc:fields(n)%iec,fields(n)%jsc:fields(n)%jec+1,1:fields(n)%npz))
    elseif (fields(n)%staggerloc == east) then
      allocate(array_tmp(fields(n)%isc:fields(n)%iec+1,fields(n)%jsc:fields(n)%jec,1:fields(n)%npz))
    endif

    do lev_in = 1,fields(n)%npz

      lev_out = fields(n)%npz-lev_in+1
      array_tmp(:,:,lev_out) = fields(n)%array(:,:,lev_in)

    enddo

    fields(n)%array = array_tmp

    deallocate(array_tmp)

  endif

enddo

end subroutine flip_array_vertical

! ------------------------------------------------------------------------------

subroutine copy_subset(rhs,lhs,not_copied)

implicit none
type(fv3jedi_field),                      intent(in)    :: rhs(:)
type(fv3jedi_field),                      intent(inout) :: lhs(:)
character(len=10), allocatable, optional, intent(out)   :: not_copied(:)

integer :: var
character(len=10) :: not_copied_(10000)
integer :: num_not_copied

! Loop over fields and copy if existing in both
num_not_copied = 0
do var = 1, size(lhs)
  if (has_field(rhs, lhs(var)%fv3jedi_name )) then
    call copy_field_array(rhs, lhs(var)%fv3jedi_name, lhs(var)%array)
  else
    num_not_copied = num_not_copied + 1
    not_copied_(num_not_copied) = lhs(var)%fv3jedi_name
  endif
enddo

! Send back list of variables not retrivable from rhs
if (present(not_copied) .and. num_not_copied > 0) then
  allocate(not_copied(num_not_copied))
  not_copied(1:num_not_copied) = not_copied_(1:num_not_copied)
endif

end subroutine copy_subset

! --------------------------------------------------------------------------------------------------

end module fv3jedi_field_mod
