! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_field_mod

use fckit_mpi_module
use fv3jedi_kinds_mod, only: kind_real
use mpp_domains_mod,   only: east, north, center

implicit none

private
public :: fv3jedi_field, get_field, &
          fields_rms, fields_gpnorm, fields_print, &
          checksame, flip_array_vertical

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
 contains
  procedure :: allocate_field
  procedure :: array_pointer
  procedure :: equals
  generic :: assignment(=) => equals
  procedure :: deallocate_field
endtype fv3jedi_field

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine allocate_field(self,isc,iec,jsc,jec,npz,short_name,long_name,&
                          fv3jedi_name,units,staggerloc,tracer,arraypointer)

implicit none
class(fv3jedi_field), target,  intent(inout) :: self
integer,                       intent(in)    :: isc,iec,jsc,jec,npz
character(len=*),              intent(in)    :: short_name
character(len=*),              intent(in)    :: long_name
character(len=*),              intent(in)    :: fv3jedi_name
character(len=*),              intent(in)    :: units
integer,                       intent(in)    :: staggerloc
real(kind=kind_real), pointer, intent(inout) :: arraypointer(:,:,:)
logical, optional,             intent(in)    :: tracer

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
    allocate(arraypointer(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
  elseif (staggerloc == north) then
    allocate(self%array(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
    allocate(arraypointer(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
  elseif (staggerloc == east) then
    allocate(self%array(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
    allocate(arraypointer(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
  endif

endif

arraypointer => self%array

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

end subroutine allocate_field

! ------------------------------------------------------------------------------

subroutine array_pointer(self,arraypointer)

implicit none
class(fv3jedi_field), target,  intent(in)    :: self
real(kind=kind_real), pointer, intent(inout) :: arraypointer(:,:,:)

if (self%staggerloc == center) then
  allocate(arraypointer(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
elseif (self%staggerloc == north) then
  allocate(arraypointer(self%isc:self%iec,self%jsc:self%jec+1,1:self%npz))
elseif (self%staggerloc == east) then
  allocate(arraypointer(self%isc:self%iec+1,self%jsc:self%jec,1:self%npz))
endif

arraypointer => self%array

end subroutine array_pointer

! ------------------------------------------------------------------------------

subroutine deallocate_field(self)

implicit none
class(fv3jedi_field), intent(inout) :: self

if(self%lalloc) deallocate(self%array)
self%lalloc = .false.

end subroutine deallocate_field

! ------------------------------------------------------------------------------

subroutine equals(self,rhs)

implicit none
class(fv3jedi_field), intent(inout) :: self
type (fv3jedi_field), intent(in)    :: rhs

real(kind=kind_real), pointer :: dummy_pointer(:,:,:) => null()

call self%allocate_field( rhs%isc,rhs%iec,rhs%jsc,rhs%jec,rhs%npz, &
                          short_name=rhs%short_name, &
                          long_name=rhs%long_name, &
                          fv3jedi_name=rhs%fv3jedi_name, &
                          units=rhs%units, &
                          staggerloc=rhs%staggerloc, &
                          tracer = rhs%tracer, &
                          arraypointer = dummy_pointer)

self%array = rhs%array

nullify(dummy_pointer)

end subroutine equals

! ------------------------------------------------------------------------------

subroutine get_field(nf,fields,fv3jedi_name,field_pointer)

integer,                      intent(in)  :: nf
type(fv3jedi_field), target,  intent(in)  :: fields(nf)
character(len=10),            intent(in)  :: fv3jedi_name
type(fv3jedi_field), pointer, intent(out) :: field_pointer

integer :: var
logical :: found

found = .false.
do var = 1,nf
  if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
    field_pointer => fields(var)
    found = .true.
    exit
  endif
enddo

if (.not.found) call abor1_ftn("fv3jedi_field get_field. Field "&
                                //fv3jedi_name//" not found in field array")

end subroutine get_field

! ------------------------------------------------------------------------------

subroutine fields_rms(nf,fields,rms)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
real(kind=kind_real), intent(inout) :: rms

integer :: i, j, k, ii, iisum, var
real(kind=kind_real) :: zz
type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

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

! ------------------------------------------------------------------------------

subroutine fields_gpnorm(nf, fields, pstat)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
real(kind=kind_real), intent(inout) :: pstat(3, nf)

integer :: var
real(kind=kind_real) :: tmp(3),  gs3, gs3g
type(fckit_mpi_comm) :: f_comm

f_comm = fckit_mpi_comm()

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

! ------------------------------------------------------------------------------

subroutine fields_print(nf, fields, name)

implicit none
integer,              intent(in)    :: nf
type(fv3jedi_field),  intent(in)    :: fields(nf)
character(len=*),     intent(in)    :: name

integer :: var
real(kind=kind_real) :: tmp(3), pstat(3), gs3, gs3g
type(fckit_mpi_comm) :: f_comm
character(len=34) :: printname

integer :: ngrid, sngrid

f_comm = fckit_mpi_comm()

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

! ------------------------------------------------------------------------------

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

end module fv3jedi_field_mod
