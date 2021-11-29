! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fields_mod

! atlas
use atlas_module, only: atlas_field, atlas_fieldset, atlas_real

! fckit
use fckit_mpi_module
use fckit_configuration_module

! oops
use datetime_mod
use missing_values_mod
use oops_variables_mod
use string_utils

! fv3jedi
use fv3jedi_field_mod,         only: fv3jedi_field, field_clen, checksame, get_field, put_field, &
                                     hasfield, short_name_to_fv3jedi_name
use fv3jedi_geom_mod,          only: fv3jedi_geom
use fv3jedi_interpolation_mod, only: field2field_interp
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
    procedure, public :: accumul
    procedure, public :: serialize
    procedure, public :: deserialize
    procedure, public :: set_atlas
    procedure, public :: to_atlas
    procedure, public :: from_atlas

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

subroutine create(self, geom, vars)

  class(fv3jedi_fields), intent(inout) :: self
  type(fv3jedi_geom),    intent(in)    :: geom
  type(oops_variables),  intent(in)    :: vars

  integer :: var, fc
  type(field_metadata) :: fmd
  logical :: field_fail

  ! Allocate fields structure
  ! -------------------------
  self%nf = vars%nvars()
  allocate(self%fields(self%nf))


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
    self%fields(fc)%long_name    = trim(fmd%long_name)
    self%fields(fc)%fv3jedi_name = trim(fmd%field_name)
    self%fields(fc)%units        = trim(fmd%units)
    self%fields(fc)%io_file      = trim(fmd%io_file)
    self%fields(fc)%space        = trim(fmd%space)
    self%fields(fc)%staggerloc   = trim(fmd%stagger_loc)
    self%fields(fc)%tracer       = fmd%tracer
    self%fields(fc)%integerfield = trim(fmd%array_kind)=='integer'
    ! hard-coding options for direction and integer fields, otherwise use user spec.
    if ( self%fields(fc)%integerfield ) then
      self%fields(fc)%interp_type  = 'integer'
    elseif (  trim(self%fields(fc)%space) == 'direction' ) then
      self%fields(fc)%interp_type  = 'nearest'
    else
      self%fields(fc)%interp_type  = trim(fmd%interp_type)
    endif

  enddo

  ! Check field count
  if (fc .ne. self%nf) call abor1_ftn("fv3jedi_fields_mod create: fc does not equal self%nf")

  ! Initialize all arrays to zero
  call self%zero()

  ! Initialize datetime to zero
  self%calendar_type = 0
  self%date_init(:) = 0

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

  !Dummy FMS style datetime
  self%calendar_type = 0
  self%date_init = 0

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
        if (self%fields(var)%array(i,j,k)/=missing_value(0.0_kind_real)) then
           zz = zz + self%fields(var)%array(i,j,k)**2
           ii = ii + 1
        endif
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
    if (trim(other%fields(var)%interp_type) == "integer") integer_interp = .true.
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
tmp(1) = minval(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz), &
 & mask=self%fields(field_num)%array(isc:iec,jsc:jec,1:npz)/=missing_value(0.0_kind_real))
tmp(2) = maxval(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz), &
 & mask=self%fields(field_num)%array(isc:iec,jsc:jec,1:npz)/=missing_value(0.0_kind_real))
tmp(3) =    sum(self%fields(field_num)%array(isc:iec,jsc:jec,1:npz)**2, &
 & mask=self%fields(field_num)%array(isc:iec,jsc:jec,1:npz)/=missing_value(0.0_kind_real))

! Get global min/max/sum
call self%f_comm%allreduce(tmp(1), minmaxrmsout(1), fckit_mpi_min())
call self%f_comm%allreduce(tmp(2), minmaxrmsout(2), fckit_mpi_max())
call self%f_comm%allreduce(tmp(3), minmaxrmsout(3), fckit_mpi_sum())

! SumSquares to rms
minmaxrmsout(3) = sqrt(minmaxrmsout(3)/gs3g)

endsubroutine minmaxrms

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

subroutine set_atlas(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(in)    :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(inout) :: afieldset

integer :: jvar, npz
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name


do jvar = 1,vars%nvars()

  ! Get fv3-jedi field
  call short_name_to_fv3jedi_name(self%fields, trim(vars%variable(jvar)), fv3jedi_name)
  call self%get_field(fv3jedi_name, field)

  ! Create the field in the fieldset if it doesn't already exists
  if (.not.afieldset%has_field(field%long_name)) then

    ! Variable dimension
    npz = field%npz
    if (npz==1) npz = 0

    ! Create Atlas field
    afield = geom%afunctionspace%create_field( name=field%long_name, kind=atlas_real(kind_real), &
                                               levels=npz )

    ! Add to fieldsets
    call afieldset%add(afield)

    ! Release pointer
    call afield%final()

  endif

enddo

end subroutine set_atlas

! --------------------------------------------------------------------------------------------------

subroutine to_atlas(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(in)    :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(inout) :: afieldset

integer :: jvar, npz, jl
real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name

do jvar = 1,vars%nvars()

  ! Get fv3-jedi field
  call short_name_to_fv3jedi_name(self%fields, trim(vars%variable(jvar)), fv3jedi_name)
  call self%get_field(fv3jedi_name, field)

  ! Variable dimension
  npz = field%npz
  if (npz==1) npz = 0

  ! Get/create Atlas field
  if (afieldset%has_field(field%long_name)) then

    ! Get Atlas field
    afield = afieldset%field(field%long_name)

  else

    ! Create field
    afield = geom%afunctionspace%create_field( name=field%long_name, kind=atlas_real(kind_real), &
                                               levels=npz )

    ! Add to fieldsets
    call afieldset%add(afield)

  endif

  if (npz==0) then

    ! Get pointer to Atlas data
    call afield%data(real_ptr_1)

    ! Copy the data
    real_ptr_1 = reshape(field%array(geom%isc:geom%iec,geom%jsc:geom%jec,1),(/geom%ngrid/))

  else

    ! Get pointer to Atlas data
    call afield%data(real_ptr_2)

    ! Copy the data
    do jl=1, npz
      real_ptr_2(jl,:) = reshape(field%array(geom%isc:geom%iec,geom%jsc:geom%jec,jl),(/geom%ngrid/))
    enddo

  endif

  ! Release pointer
  call afield%final()

enddo

end subroutine to_atlas

! --------------------------------------------------------------------------------------------------

subroutine from_atlas(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(inout) :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(in)    :: afieldset

integer :: jvar, npz, jl
real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
character(len=field_clen) :: fv3jedi_name


do jvar = 1,vars%nvars()

  ! Get fv3-jedi field
  call short_name_to_fv3jedi_name(self%fields, trim(vars%variable(jvar)), fv3jedi_name)
  call self%get_field(fv3jedi_name, field)

  ! Variable dimension
  npz = field%npz
  if (npz==1) npz = 0

  ! Get Atlas field
  afield = afieldset%field(field%long_name)

  if (npz==0) then

    ! Get pointer to Atlas field data
    call afield%data(real_ptr_1)

    ! Copy to fv3-jedi field
    field%array(geom%isc:geom%iec,geom%jsc:geom%jec,1) = reshape(real_ptr_1, &
                                                       & (/geom%iec-geom%isc+1,geom%jec-geom%jsc+1/))
  else

    ! Get pointer to Atlas field data
    call afield%data(real_ptr_2)

    ! Copy to fv3-jedi field
    do jl=1, npz
      field%array(geom%isc:geom%iec,geom%jsc:geom%jec,jl) = reshape(real_ptr_2(jl,:), &
                                                          & (/geom%iec-geom%isc+1,geom%jec-geom%jsc+1/))
    enddo

  endif

  ! Release pointer
  call afield%final()

enddo

end subroutine from_atlas

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
