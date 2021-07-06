! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fields_mod

! fckit
use fckit_mpi_module
use fckit_configuration_module

! oops
use datetime_mod
use oops_variables_mod
use string_utils

! fv3jedi
use fv3jedi_field_mod,         only: fv3jedi_field, field_clen, checksame, get_field, put_field, hasfield
use fv3jedi_geom_mod,          only: fv3jedi_geom
use fv3jedi_interpolation_mod, only: field2field_interp
use fv3jedi_io_gfs_mod,        only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod,       only: fv3jedi_io_geos
use fv3jedi_io_latlon_mod,     only: fv3jedi_llgeom
use fv3jedi_kinds_mod,         only: kind_real
use fields_metadata_mod,       only: field_metadata

implicit none
private
public :: fv3jedi_fields

! --------------------------------------------------------------------------------------------------

! Fields type (base for State and Increment)
type :: fv3jedi_fields

  integer :: isc, iec, jsc, jec, npx, npy, npz, nf             ! Geometry convenience
  integer :: calendar_type, date_init(6)                       ! FMS style datetime
  type(fckit_mpi_comm) :: f_comm                               ! Communicator
  type(fv3jedi_field), allocatable :: fields(:)                ! Array of fields

  contains

    ! Methods needed by both state and increment classes
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: copy
    procedure, public :: zero
    procedure, public :: norm
    procedure, public :: change_resol
    procedure, public :: minmaxrms
    procedure, public :: read
    procedure, public :: write
    procedure, public :: accumul
    procedure, public :: serialize
    procedure, public :: deserialize

    ! Public array/field accessor functions
    procedure, public :: has_field => has_field_
    generic,   public :: get_field => get_field_return_type_pointer, &
                                      get_field_return_array_pointer, &
                                      get_field_return_array_allocatable
    procedure, public :: put_field => put_field_

    ! Private array/field accessor functions
    procedure, private :: get_field_return_type_pointer
    procedure, private :: get_field_return_array_pointer
    procedure, private :: get_field_return_array_allocatable

endtype fv3jedi_fields

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, vars, increment)

  class(fv3jedi_fields), intent(inout) :: self
  type(fv3jedi_geom),    intent(in)    :: geom
  type(oops_variables),  intent(in)    :: vars
  logical, optional,     intent(in)    :: increment

  integer :: var, fc
  type(field_metadata) :: fmd
  logical :: is_increment, field_fail

  ! Allocate fields structure
  ! -------------------------
  self%nf = vars%nvars()
  allocate(self%fields(self%nf))

  ! Check if this is an increment rather than state
  ! -----------------------------------------------
  if (.not. present(increment)) then
    is_increment = .false.
  else
    is_increment = increment
  endif

  ! Loop through and allocate actual fields
  ! ---------------------------------------
  fc = 0
  do var = 1, vars%nvars()

    fmd = geom%fields%get_field(trim(vars%variable(var)))

    ! Uptick counter
    fc=fc+1;

    self%fields(fc)%isc = geom%isc
    self%fields(fc)%iec = geom%iec
    self%fields(fc)%jsc = geom%jsc
    self%fields(fc)%jec = geom%jec
    self%fields(fc)%npz = fmd%levels

    field_fail = len(trim(fmd%field_name)) > field_clen
    if (field_fail) call abor1_ftn("fv3jedi_fields.create: " //trim(fmd%field_name)// " too long")

    if(.not.self%fields(fc)%lalloc) then

      if (trim(fmd%stagger_loc) == 'center') then
        allocate(self%fields(fc)%array(geom%isc:geom%iec,geom%jsc:geom%jec,1:fmd%levels))
      elseif (trim(fmd%stagger_loc) == 'northsouth') then
        allocate(self%fields(fc)%array(geom%isc:geom%iec,geom%jsc:geom%jec+1,1:fmd%levels))
      elseif (trim(fmd%stagger_loc) == 'eastwest') then
        allocate(self%fields(fc)%array(geom%isc:geom%iec+1,geom%jsc:geom%jec,1:fmd%levels))
      endif

    endif

    self%fields(fc)%lalloc = .true.

    self%fields(fc)%short_name   = trim(fmd%field_io_name)
    if (is_increment) then
      self%fields(fc)%long_name    = "increment_of_"//trim(fmd%long_name)
    else
      self%fields(fc)%long_name    = trim(fmd%long_name)
    endif
    self%fields(fc)%fv3jedi_name = trim(fmd%field_name)
    self%fields(fc)%units        = trim(fmd%units)
    self%fields(fc)%io_file      = trim(fmd%io_file)
    self%fields(fc)%space        = trim(fmd%space)
    self%fields(fc)%staggerloc   = trim(fmd%stagger_loc)
    self%fields(fc)%tracer       = fmd%tracer
    self%fields(fc)%integerfield = trim(fmd%array_kind)=='integer'

  enddo

  ! Check field count
  if (fc .ne. self%nf) call abor1_ftn("fv3jedi_fields_mod create: fc does not equal self%nf")

  ! Initialize all arrays to zero
  call self%zero()

  ! Copy some geometry for convenience
  self%isc    = geom%isc
  self%iec    = geom%iec
  self%jsc    = geom%jsc
  self%jec    = geom%jec
  self%npx    = geom%npx
  self%npy    = geom%npy
  self%npz    = geom%npz

  ! Pointer to fv3jedi communicator
  self%f_comm = geom%f_comm

  ! Check winds
  field_fail = self%has_field('ua') .and. .not.self%has_field('va')
  if (field_fail) call abor1_ftn("fv3jedi_fields.create: found A-Grid u but not v")
  field_fail = .not.self%has_field('ua') .and. self%has_field('va')
  if (field_fail) call abor1_ftn("fv3jedi_fields.create: found A-Grid v but not u")
  field_fail = self%has_field('ud') .and. .not.self%has_field('vd')
  if (field_fail) call abor1_ftn("fv3jedi_fields.create: found D-Grid u but not v")
  field_fail = .not.self%has_field('ud') .and. self%has_field('vd')
  if (field_fail) call abor1_ftn("fv3jedi_fields.create: found D-Grid v but not u")

  !Check User's choice of ozone variables.
  field_fail = self%has_field('o3mr') .and. self%has_field('o3ppmv')
  if (field_fail) call abor1_ftn("fv3jedi_fields.create: o3mr and o3ppmv created")

endsubroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_fields), intent(inout) :: self

integer :: var

do var = 1, self%nf
  if(self%fields(var)%lalloc) deallocate(self%fields(var)%array)
  self%fields(var)%lalloc = .false.
enddo
deallocate(self%fields)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine copy(self, other)

class(fv3jedi_fields), intent(inout) :: self
class(fv3jedi_fields),  intent(in)    :: other

integer :: var

call checksame(self%fields, other%fields, "fv3jedi_fields_mod.copy")

do var = 1, self%nf
  self%fields(var)%array = other%fields(var)%array
enddo

self%calendar_type = other%calendar_type
self%date_init = other%date_init

end subroutine copy

! --------------------------------------------------------------------------------------------------

subroutine zero(self)

class(fv3jedi_fields), intent(inout) :: self

integer :: var

do var = 1, self%nf
  self%fields(var)%array = 0.0_kind_real
enddo

endsubroutine zero

! --------------------------------------------------------------------------------------------------

subroutine norm(self, normout)

class(fv3jedi_fields), intent(inout) :: self
real(kind=kind_real),  intent(out)   :: normout

integer :: i, j, k, ii, iisum, var
real(kind=kind_real) :: zz

zz = 0.0_kind_real
ii = 0

do var = 1, self%nf

  do k = 1, self%fields(var)%npz
    do j = self%fields(var)%jsc, self%fields(var)%jec
      do i = self%fields(var)%isc, self%fields(var)%iec
        zz = zz + self%fields(var)%array(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo

enddo

!Get global values
call self%f_comm%allreduce(zz,normout,fckit_mpi_sum())
call self%f_comm%allreduce(ii,iisum,fckit_mpi_sum())
normout = sqrt(normout/real(iisum,kind_real))

endsubroutine norm

! --------------------------------------------------------------------------------------------------

subroutine change_resol(self, geom, other, geom_other)

implicit none
class(fv3jedi_fields), intent(inout) :: self
type(fv3jedi_geom),    intent(inout) :: geom
class(fv3jedi_fields), intent(in)    :: other
type(fv3jedi_geom),    intent(inout) :: geom_other

! Interpolation
integer :: var
type(field2field_interp) :: interp
logical :: integer_interp = .false.

call checksame(self%fields, other%fields, "fv3jedi_fields_mod.change_resol")

if ((other%iec-other%isc+1)-(self%iec-self%isc+1) == 0) then

  call self%copy(other)

else

  ! Check if integer interp needed
  do var = 1, self%nf
    if (other%fields(var)%integerfield) integer_interp = .true.
  enddo

  call interp%create(geom%interp_method, integer_interp, geom_other, geom)
  call interp%apply(self%nf, geom_other, other%fields, geom, self%fields)
  call interp%delete()

  self%calendar_type = other%calendar_type
  self%date_init = other%date_init

endif

end subroutine change_resol

! --------------------------------------------------------------------------------------------------

subroutine minmaxrms(self, field_num, field_name, minmaxrmsout)

class(fv3jedi_fields), intent(inout) :: self
integer,               intent(in)    :: field_num
character(len=*),      intent(inout) :: field_name
real(kind=kind_real),  intent(out)   :: minmaxrmsout(3)

integer :: isc, iec, jsc, jec, npz
real(kind=kind_real) :: tmp(3), gs3, gs3g

isc = self%fields(field_num)%isc
iec = self%fields(field_num)%iec
jsc = self%fields(field_num)%jsc
jec = self%fields(field_num)%jec
npz = self%fields(field_num)%npz

field_name = self%fields(field_num)%short_name

! Compute global sum over the field
gs3 = real((iec-isc+1)*(jec-jsc+1)*npz, kind_real)
call self%f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

! Min/Max/SumSquares
tmp(1) = minval(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz))
tmp(2) = maxval(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz))
tmp(3) =    sum(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz)**2)

! Get global min/max/sum
call self%f_comm%allreduce(tmp(1), minmaxrmsout(1), fckit_mpi_min())
call self%f_comm%allreduce(tmp(2), minmaxrmsout(2), fckit_mpi_max())
call self%f_comm%allreduce(tmp(3), minmaxrmsout(3), fckit_mpi_sum())

! SumSquares to rms
minmaxrmsout(3) = sqrt(minmaxrmsout(3)/gs3g)

endsubroutine minmaxrms

! --------------------------------------------------------------------------------------------------

subroutine read(self, geom, conf, vdate)

class(fv3jedi_fields),     intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf
type(datetime),            intent(inout) :: vdate

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos

character(len=10) :: filetype
character(len=:), allocatable :: str


! IO type
call conf%get_or_die("filetype",str)
filetype = str
deallocate(str)


if (trim(filetype) == 'gfs') then

  call gfs%setup_conf(conf)
  call gfs%read_meta(vdate, self%calendar_type, self%date_init)
  call gfs%read_fields(self%fields, geom%domain, geom%npz)

elseif (trim(filetype) == 'geos') then

  call geos%setup_conf(geom, conf)
  call geos%read_meta(geom, vdate, self%calendar_type, self%date_init, self%fields)
  call geos%read_fields(geom, self%fields)
  call geos%delete()

else

  call abor1_ftn("fv3jedi_fields_mod.read: restart type not supported")

endif

endsubroutine read

! --------------------------------------------------------------------------------------------------

subroutine write(self, geom, conf, vdate)

class(fv3jedi_fields),     intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf
type(datetime),            intent(inout) :: vdate

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos
type(fv3jedi_llgeom)  :: latlon

character(len=10) :: filetype
character(len=:), allocatable :: str

! IO type
call conf%get_or_die("filetype",str)
filetype = str
deallocate(str)

if (trim(filetype) == 'gfs') then

  call gfs%setup_conf(conf)
  call gfs%setup_date(vdate)
  call gfs%write(geom%domain, self%fields, vdate, self%calendar_type, self%date_init)

elseif (trim(filetype) == 'geos') then

  call geos%setup_conf(geom, conf)
  call geos%setup_date(vdate)
  call geos%write(geom, self%fields, vdate)
  call geos%delete()

elseif (trim(filetype) == 'latlon') then

  call latlon%setup_conf(geom)
  call latlon%setup_date(vdate)
  call latlon%write(geom, self%fields, conf, vdate)

else

    call abor1_ftn("fv3jedi_fields.write: restart type not supported")

endif

endsubroutine write

! --------------------------------------------------------------------------------------------------

subroutine accumul(self, zz, rhs)

implicit none
class(fv3jedi_fields), intent(inout) :: self
real(kind=kind_real),  intent(in)    :: zz
class(fv3jedi_fields), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_fields.accumul")

do var = 1, self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine accumul

! --------------------------------------------------------------------------------------------------

subroutine serialize(self,vsize,vect_inc)

implicit none

! Passed variables
class(fv3jedi_fields), intent(in)  :: self
integer,               intent(in)  :: vsize
real(kind_real),       intent(out) :: vect_inc(vsize)

! Local variables
integer :: ind, var, i, j, k

! Initialize
ind = 0

! Copy
do var = 1, self%nf
  do k = 1,self%fields(var)%npz
    do j = self%fields(var)%jsc,self%fields(var)%jec
      do i = self%fields(var)%isc,self%fields(var)%iec
        ind = ind + 1
        vect_inc(ind) = self%fields(var)%array(i, j, k)
      enddo
    enddo
  enddo
enddo

end subroutine serialize

! --------------------------------------------------------------------------------------------------

subroutine deserialize(self, vsize, vect_inc, index)

implicit none

! Passed variables
class(fv3jedi_fields), intent(inout) :: self
integer,               intent(in) :: vsize
real(kind_real),       intent(in) :: vect_inc(vsize)
integer,               intent(inout) :: index

! Local variables
integer :: ind, var, i, j, k

! Copy
do var = 1, self%nf
  do k = 1,self%fields(var)%npz
    do j = self%fields(var)%jsc,self%fields(var)%jec
      do i = self%fields(var)%isc,self%fields(var)%iec
        self%fields(var)%array(i, j, k) = vect_inc(index + 1)
        index = index + 1
      enddo
    enddo
  enddo
enddo

end subroutine deserialize

! --------------------------------------------------------------------------------------------------

logical function has_field_(self, field_name, field_index)

class(fv3jedi_fields), intent(in)  :: self
character(len=*),      intent(in)  :: field_name
integer, optional,     intent(out) :: field_index

has_field_ = hasfield(self%fields, field_name, field_index)

end function has_field_

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_type_pointer(self, field_name, field)

class(fv3jedi_fields), target, intent(in)    :: self
character(len=*),              intent(in)    :: field_name
type(fv3jedi_field), pointer,  intent(inout) :: field

call get_field(self%fields, field_name, field)

endsubroutine get_field_return_type_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_pointer(self, field_name, field)

class(fv3jedi_fields), target, intent(in)    :: self
character(len=*),              intent(in)    :: field_name
real(kind=kind_real), pointer, intent(inout) :: field(:,:,:)

call get_field(self%fields, field_name, field)

endsubroutine get_field_return_array_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_allocatable(self, field_name, field)

class(fv3jedi_fields), target,     intent(in)    :: self
character(len=*),                  intent(in)    :: field_name
real(kind=kind_real), allocatable, intent(inout) :: field(:,:,:)

call get_field(self%fields, field_name, field)

endsubroutine get_field_return_array_allocatable

! --------------------------------------------------------------------------------------------------

subroutine put_field_(self, field_name, field)

class(fv3jedi_fields), target,     intent(inout) :: self
character(len=*),                  intent(in)    :: field_name
real(kind=kind_real), allocatable, intent(in)    :: field(:,:,:)

call put_field(self%fields, field_name, field)

endsubroutine put_field_

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fields_mod
