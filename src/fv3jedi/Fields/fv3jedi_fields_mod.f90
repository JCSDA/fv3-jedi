! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fields_mod

! atlas
use atlas_module, only: atlas_field, atlas_fieldset, atlas_real, atlas_metadata

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
                                     hasfield, create_field
use fv3jedi_geom_mod,          only: fv3jedi_geom
use fv3jedi_interpolation_mod, only: field2field_interp
use fv3jedi_kinds_mod,         only: kind_real
use fields_metadata_mod,       only: field_metadata

use mpp_domains_mod, only: mpp_update_domains, mpp_update_domains_ad

implicit none
private
public :: fv3jedi_fields

! --------------------------------------------------------------------------------------------------

! Fields type (base for State and Increment)
type :: fv3jedi_fields

  integer :: isc, iec, jsc, jec, npx, npy, npz, nf             ! Geometry convenience
  type(fckit_mpi_comm) :: f_comm                               ! Communicator
  type(fv3jedi_field), allocatable :: fields(:)                ! Array of fields
  type(datetime) :: time

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
    procedure, public :: to_fieldset
    procedure, public :: to_fieldset_ad
    procedure, public :: from_fieldset
    procedure, public :: update_fields

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
  logical :: field_fail

  ! Allocate fields structure
  ! -------------------------
  self%nf = vars%nvars()
  allocate(self%fields(self%nf))

  ! Loop through and allocate actual fields
  ! ---------------------------------------
  fc = 0
  do var = 1, vars%nvars()

    ! Uptick counter
    fc=fc+1;

    ! Copy geometry
    self%fields(fc)%isc = geom%isc
    self%fields(fc)%iec = geom%iec
    self%fields(fc)%jsc = geom%jsc
    self%fields(fc)%jec = geom%jec

    ! Set this fields meta data
    call create_field(self%fields(fc), geom%fields%get_field(trim(vars%variable(var))), geom%f_comm)

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

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_fields), intent(inout) :: self

integer :: var

do var = 1, size(self%fields)
  if (self%fields(var)%lalloc) deallocate(self%fields(var)%array)
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
    if (trim(other%fields(var)%interpolation_type) == "integer") integer_interp = .true.
  enddo

  call interp%create(geom%interp_method, integer_interp, geom_other, geom)
  call interp%apply(self%nf, geom_other, other%fields, geom, self%fields)
  call interp%delete()

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

field_name = self%fields(field_num)%long_name

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

! Fills a FieldSet with model data, including halos
! - exchanges halo data between fv3 grid tiles
! - stores fv3-jedi field (local + halo) data in atlas FieldSet as 1d fields
subroutine to_fieldset(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(in)    :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(inout) :: afieldset

integer :: jvar, jl
real(kind=kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
type(atlas_metadata) :: meta

real(kind=kind_real), allocatable :: tmp_fv3_data(:,:,:)
integer :: max_npz

! To allocate tmp_fv3_data just once, we need the max number of levels
max_npz = -1
do jvar = 1,vars%nvars()
  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)
  max_npz = max(max_npz, field%npz)
end do

allocate(tmp_fv3_data(geom%isd:geom%ied, geom%jsd:geom%jed, max_npz))

do jvar = 1,vars%nvars()

  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)

  ! Fill tmp buffer with local data
  ! Note: if field%npz < max_npz, then some of the memory allocated is unused
  tmp_fv3_data = 0.0_kind_real
  tmp_fv3_data(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz) = &
    field%array(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz)

  ! Perform halo exchange to fill halos -- this uses MPI to exchange halo data
  call mpp_update_domains(tmp_fv3_data, geom%domain)

  ! Get/create atlas field
  if (afieldset%has_field(field%long_name)) then
    afield = afieldset%field(field%long_name)
  else
    afield = geom%afunctionspace_incl_halo%create_field( name=field%long_name, &
                                                         kind=atlas_real(kind_real), &
                                                         levels=field%npz )
    call afieldset%add(afield)
  endif

  ! Copy fv3-jedi field data into atlas field
  ! Note the jedi interp grid excludes fv3's corner halos at cubed-sphere grid corners, which
  ! results in irregularly-shaped data. We pack this irregular data into a 1d atlas field array.
  call afield%data(real_ptr)
  do jl=1, field%npz
    call geom%trim_fv3_grid_to_jedi_interp_grid(tmp_fv3_data(:,:,jl), real_ptr(jl,:))
  enddo

  meta = afield%metadata()
  call meta%set('interp_type', trim(field%interpolation_type))
  ! Set atlas::Field metadata for interp mask IFF the user specifically requested a non-default mask
  if (trim(field%interpolation_source_point_mask) .ne. 'default') then
    call meta%set('interp_source_point_mask', trim(field%interpolation_source_point_mask))
  end if

  ! Release pointer
  call afield%final()

enddo

deallocate(tmp_fv3_data)

end subroutine to_fieldset

! --------------------------------------------------------------------------------------------------

! (Adjoint of) Fills a FieldSet with model data, including halos
subroutine to_fieldset_ad(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(inout) :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(in)    :: afieldset

integer :: jvar, jl
real(kind=kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
character(len=1024) :: errmsg

real(kind=kind_real), allocatable :: tmp_fv3_data(:,:,:)
integer :: max_npz
integer :: ngrid

ngrid = geom%ngrid_including_halo()

! To allocate tmp_fv3_data just once, we need the max number of levels
max_npz = -1
do jvar = 1,vars%nvars()
  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)
  max_npz = max(max_npz, field%npz)
end do

allocate(tmp_fv3_data(geom%isd:geom%ied, geom%jsd:geom%jed, max_npz))

do jvar = 1,vars%nvars()

  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)

  ! Get atlas field
  afield = afieldset%field(field%long_name)

  ! Sanity check that afield has expected size, i.e., includes halo points
  if (ngrid /= afield%size()/field%npz) then
    write (errmsg,*) "Bad afield size in to_fieldset_ad: expect ngrid = ", ngrid, &
                     ", but got afield ngrid = ", afield%size()/field%npz
    call abor1_ftn(trim(errmsg))
  end if

  ! (Adjoint of) Copy fv3-jedi field data into atlas field
  ! Note: this includes reshaping of the 1d afield data (local + halo) into 2d fv3-jedi field
  tmp_fv3_data = 0.0
  call afield%data(real_ptr)
  do jl=1,field%npz
    call geom%trim_fv3_grid_to_jedi_interp_grid_ad(tmp_fv3_data(:,:,jl), real_ptr(jl,:))
  enddo

  ! (Adjoint of) Perform halo exchange
  call mpp_update_domains_ad(tmp_fv3_data, geom%domain)

  ! (Adjoint of) Fill tmp buffer with local data
  field%array(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz) = &
    field%array(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz) &
    + tmp_fv3_data(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz)

  ! Release pointer
  call afield%final()

enddo

deallocate(tmp_fv3_data)

end subroutine to_fieldset_ad

! --------------------------------------------------------------------------------------------------

! Fills model data from an atlas FieldSet
! - copies local portion of packed (local + halo) data into fv3-jedi field
! - Note 1: assumes the local data occupies the first contiguous portion of the packed data
! - Note 2: does NOT copy the halo data; halo data values are out-of-sync after this call
subroutine from_fieldset(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(inout) :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(in)    :: afieldset

integer :: jvar, jl
real(kind=kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: afield
type(fv3jedi_field), pointer :: field
character(len=1024) :: errmsg

integer :: ngrid

ngrid = geom%ngrid_including_halo()

do jvar = 1,vars%nvars()

  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)

  ! Get Atlas field
  afield = afieldset%field(field%long_name)

  ! Sanity check that afield has expected size, i.e., includes halo points
  if (ngrid /= afield%size()/field%npz) then
    write (errmsg,*) "Bad afield size in from_fieldset: expect ngrid = ", ngrid, &
                     ", but got afield ngrid = ", afield%size()/field%npz
    call abor1_ftn(trim(errmsg))
  end if

  ! Reshape the first portion of afield, the local data, into 2d fv3-jedi field
  call afield%data(real_ptr)
  do jl=1,field%npz
    field%array(geom%isc:geom%iec, geom%jsc:geom%jec, jl) = &
      reshape(real_ptr(jl, 1:geom%ngrid), (/geom%iec-geom%isc+1, geom%jec-geom%jsc+1/))
  enddo

  ! Release pointer
  call afield%final()

enddo

end subroutine from_fieldset

! --------------------------------------------------------------------------------------------------

subroutine update_fields(self, geom, new_vars)

implicit none

! Passed variables
class(fv3jedi_fields), intent(inout) :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: new_vars

type(fv3jedi_field), allocatable :: fields_tmp(:)
integer :: f, findex
type(field_metadata) :: fmd

! Allocate temporary array to hold fields
allocate(fields_tmp(new_vars%nvars()))

! Loop over and move fields or allocate new fields
do f = 1, new_vars%nvars()

  fields_tmp(f)%isc = geom%isc
  fields_tmp(f)%iec = geom%iec
  fields_tmp(f)%jsc = geom%jsc
  fields_tmp(f)%jec = geom%jec

  fmd = geom%fields%get_field(trim(new_vars%variable(f)))

  if (self%has_field(trim(fmd%short_name), findex)) then

    ! If already allocated then move to temporary
    call move_alloc(self%fields(findex)%array, fields_tmp(f)%array)
    self%fields(findex)%lalloc = .false.
    fields_tmp(f)%lalloc = .true.

  endif

  ! Create field in temporary (allocate will be skipped if moved)
  call create_field(fields_tmp(f), fmd, geom%f_comm)

enddo

! Move the temporary array back to self
call self%delete()
call move_alloc(fields_tmp, self%fields)

! Update number of fields
self%nf = size(self%fields)

end subroutine update_fields

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
