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
use fv3jedi_kinds_mod,         only: kind_real
use fields_metadata_mod,       only: field_metadata, clen
use wind_vt_mod,               only: a_to_d

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
  integer :: ntracers

  contains

    ! Methods needed by both state and increment classes
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: copy
    procedure, public :: zero
    procedure, public :: norm
    procedure, public :: minmaxrms
    procedure, public :: accumul
    procedure, public :: serialize
    procedure, public :: deserialize
    procedure, public :: to_fieldset
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

  integer :: var, fc, nfmd
  type(field_metadata) :: fmd
  logical :: field_fail

  ! Allocate fields structure
  ! -------------------------
  self%nf = vars%nvars()
  allocate(self%fields(self%nf))

  ! Loop through and allocate actual fields
  ! ---------------------------------------
  fc = 0
  self%ntracers = 0
  do var = 1, vars%nvars()

    ! Uptick counter
    fc=fc+1;

    ! Copy geometry
    self%fields(fc)%isc = geom%isc
    self%fields(fc)%iec = geom%iec
    self%fields(fc)%jsc = geom%jsc
    self%fields(fc)%jec = geom%jec

    ! Set this fields meta data
    call create_field(self%fields(fc), geom%fmd%get_field_metadata(trim(vars%variable(var))), &
                      geom%f_comm)

    ! count number of tracers using tracer flag
    if (self%fields(fc)%tracer) then
      self%ntracers = self%ntracers + 1
    end if
  enddo

  ! Check field count
  if (fc .ne. self%nf) call abor1_ftn("fv3jedi_fields_mod.create: fc does not equal self%nf")

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

self%ntracers = other%ntracers

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

! Fills a FieldSet with model data (follows "option 2" from guidelines: halos are set to 0)
! - fills fv3-jedi owned data in owned portion of atlas::FieldSet as 1d array
! - fills 0 in halo portion of atlas::FieldSet
subroutine to_fieldset(self, geom, vars, afieldset)

implicit none
class(fv3jedi_fields), intent(in)    :: self
type(fv3jedi_geom),    intent(in)    :: geom
type(oops_variables),  intent(in)    :: vars
type(atlas_fieldset),  intent(inout) :: afieldset

integer :: jvar, jl
type(fv3jedi_field), pointer :: field
real(kind=kind_real), pointer :: ptr(:,:)
type(atlas_field) :: afield
type(atlas_metadata) :: meta
character(len=:), allocatable :: conf_str

! Variables needed to handle boundary-condition grid points in a regional grid
integer :: max_npz, ngrid, dummy_ntris, dummy_nquads
real(kind=kind_real), allocatable :: tmp(:,:,:)

! To allocate tmp just once, we need the max number of levels
max_npz = -1
do jvar = 1, vars%nvars()
  call self%get_field(trim(vars%variable(jvar)), field)
  max_npz = max(max_npz, field%npz)
end do

! Set up regional-grid variables
if (geom%ntiles == 1) then
  call geom%get_num_nodes_and_elements(ngrid, dummy_ntris, dummy_nquads)
  allocate(tmp(geom%isd:geom%ied, geom%jsd:geom%jed, max_npz))
end if

do jvar = 1, vars%nvars()

  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)

  ! Get/create atlas field
  if (afieldset%has_field(field%long_name)) then
    afield = afieldset%field(field%long_name)
  else
    afield = geom%afunctionspace%create_field( name=field%long_name, &
                                               kind=atlas_real(kind_real), &
                                               levels=field%npz )
    call afieldset%add(afield)
  endif

  ! Copy fv3-jedi's owned field data into atlas field
  ! - For a global model, this is easy: copy the owned data.
  ! - For a regional model, extra steps are needed because atlas and fv3 have different definitions
  !   of the halo. The boundary condition points that surround a regional grid are halos in fv3,
  !   but they are "owned" boundary condition points in atlas. Here we're trying to fill the owned
  !   data for atlas, so the regional cases has to copy the fv3 owned data AND some fv3 halo data
  !   from the boundary condition regions into the atlas array.
  call afield%data(ptr)
  if (geom%ntiles == 6) then
    do jl=1, field%npz
      ptr(jl, 1:geom%ngrid) = reshape(field%array(:, :, jl), (/geom%ngrid/))
    end do
    ptr(:, geom%ngrid+1:) = 0.0_kind_real  ! Assign 0 at the atlas-generated halo points
  else if (geom%ntiles == 1) then
    tmp = 0.0_kind_real
    ! Fill the owned data at the center of the grid patch
    tmp(geom%isc:geom%iec, geom%jsc:geom%jec, 1:field%npz) = field%array(:, :, :)
    ! Manually extend data out into BC regions if BCs are detected. We do this by copying the
    ! outermost row of points from the owned region into the BC region, filling the innermost row
    ! of the BC region (1 point deep is enough for atlas).
    ! TODO(): Long term, would be better to grab the actual BC's when running the model.
    ! First, do the sides,
    if (self%isc == 1) then  ! left edge
      tmp(geom%isc-1, geom%jsc:geom%jec, 1:field%npz) = &
        tmp(geom%isc, geom%jsc:geom%jec, 1:field%npz)
    end if
    if (self%iec == self%npx-1) then  ! right edge
      tmp(geom%iec+1, geom%jsc:geom%jec, 1:field%npz) = &
        tmp(geom%iec, geom%jsc:geom%jec, 1:field%npz)
    end if
    if (self%jsc == 1) then  ! bottom edge
      tmp(geom%isc:geom%iec, geom%jsc-1, 1:field%npz) = &
        tmp(geom%isc:geom%iec, geom%jsc, 1:field%npz)
    end if
    if (self%jec == self%npy-1) then  ! top edge
      tmp(geom%isc:geom%iec, geom%jec+1, 1:field%npz) = &
        tmp(geom%isc:geom%iec, geom%jec, 1:field%npz)
    end if
    ! Then, propagate corners diagonally,
    if (self%isc == 1 .and. self%jsc == 1) then  ! bottom-left corner
      tmp(geom%isc-1, geom%jsc-1, 1:field%npz) = tmp(geom%isc, geom%jsc, 1:field%npz)
    end if
    if (self%iec == self%npx-1 .and. self%jsc == 1) then  ! bottom-right corner
      tmp(geom%iec+1, geom%jsc-1, 1:field%npz) = tmp(geom%iec, geom%jsc, 1:field%npz)
    end if
    if (self%isc == 1 .and. self%jec == self%npy-1) then  ! top-left corner
      tmp(geom%isc-1, geom%jec+1, 1:field%npz) = tmp(geom%isc, geom%jec, 1:field%npz)
    end if
    if (self%iec == self%npx-1 .and. self%jec == self%npy-1) then  ! top-right corner
      tmp(geom%iec+1, geom%jec+1, 1:field%npz) = tmp(geom%iec, geom%jec, 1:field%npz)
    end if

    ! First, assign the entire atlas Field to 0. Then, use the subroutine fv3_nodes_to_atlas_nodes
    ! to fill the owned and BC points into the atlas array in the correct order
    ptr = 0.0_kind_real
    do jl=1, field%npz
      call geom%fv3_nodes_to_atlas_nodes(tmp(:, :, jl), ptr(jl, 1:ngrid))
    end do
  end if

  ! Mark halos as being out-of-date
  call afield%set_dirty(.true.)

  ! Set the interpolation strategy based on field properties
  meta = afield%metadata()
  call meta%set('interp_type', 'default')
  if ( trim(field%kind)=='integer' ) then
    call meta%set('interp_type', 'integer')
  elseif ( trim(field%space) == 'direction' ) then
    call meta%set('interp_type', 'nearest')
  endif

  ! Further set the interpolation strategy based on user specification
  if (geom%field_interp_methods%has(trim(field%long_name))) then
    call geom%field_interp_methods%get_or_die(trim(field%long_name), conf_str)
    call meta%set('interp_type', conf_str)
  endif

  ! Set interpolation masks based on user specification
  if (geom%field_masks%has(trim(field%long_name))) then
    call geom%field_masks%get_or_die(trim(field%long_name), conf_str)
    call meta%set('mask', trim(conf_str))
  endif

  ! Set the nearest 3d level for 2d variables
  call meta%set('nearest 3d level', 'bottom')

  ! Release pointer
  call afield%final()
enddo

! Clean up
if (allocated(tmp)) deallocate(tmp)
call afield%final()
call meta%final()

end subroutine to_fieldset

! --------------------------------------------------------------------------------------------------

! Fills model data from an atlas FieldSet (follows "option 2" from guidelines: halos are set to 0)
! - copies owned portion of atlas array into fv3-jedi field
! - sets fv3-jedi field halos to 0 (they are out-of-sync after this call)
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

do jvar = 1,vars%nvars()

  ! Get field
  call self%get_field(trim(vars%variable(jvar)), field)

  ! Set to zero (could be more efficient by doing this only in halos...)
  field%array = 0.0_kind_real

  ! Get Atlas field
  afield = afieldset%field(field%long_name)

  ! Reshape the first portion of afield, the local data, into 2d fv3-jedi field
  call afield%data(real_ptr)
  do jl=1,field%npz
    field%array(geom%isc:geom%iec, geom%jsc:geom%jec, jl) = &
      reshape(real_ptr(jl, 1:geom%ngrid), (/geom%iec-geom%isc+1, geom%jec-geom%jsc+1/))
  enddo
enddo
call afield%final()

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

  fmd = geom%fmd%get_field_metadata(trim(new_vars%variable(f)))

  if (self%has_field(trim(fmd%long_name), findex)) then

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

subroutine get_field_return_type_pointer(self, field_name, field, ok_to_fail_in)

class(fv3jedi_fields), target, intent(in)    :: self
character(len=*),              intent(in)    :: field_name
type(fv3jedi_field), pointer,  intent(inout) :: field
logical, optional,             intent(in)    :: ok_to_fail_in

integer :: field_index
logical :: ok_to_fail, have_field

! Set whether the get can fail
ok_to_fail = .false.
if (present(ok_to_fail_in)) ok_to_fail = ok_to_fail_in

! Check field exists
have_field = hasfield(self%fields, field_name)

! Fail conditions
if (.not.have_field) then
  if (.not.ok_to_fail) then
    call abor1_ftn("fv3jedi_fields_mod.get_field: field " // trim(field_name) &
                   // " not found")
  else
    return
  endif
endif

! Get the actual field
call get_field(self%fields, field_name, field)

endsubroutine get_field_return_type_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_pointer(self, field_name, field, ok_to_fail_in)

class(fv3jedi_fields), target, intent(in)    :: self
character(len=*),              intent(in)    :: field_name
real(kind=kind_real), pointer, intent(inout) :: field(:,:,:)
logical, optional,             intent(in)    :: ok_to_fail_in

integer :: field_index
logical :: ok_to_fail, have_field

! Set whether the get can fail
ok_to_fail = .false.
if (present(ok_to_fail_in)) ok_to_fail = ok_to_fail_in

! Check field exists
have_field = hasfield(self%fields, field_name)

! Fail conditions
if (.not.have_field) then
  if (.not.ok_to_fail) then
    call abor1_ftn("fv3jedi_fields_mod.get_field: field " // trim(field_name) &
                   // " not found")
  else
    return
  endif
endif

! Get the actual field
call get_field(self%fields, field_name, field)

endsubroutine get_field_return_array_pointer

! --------------------------------------------------------------------------------------------------

subroutine get_field_return_array_allocatable(self, field_name, field, ok_to_fail_in)

class(fv3jedi_fields), target,     intent(in)    :: self
character(len=*),                  intent(in)    :: field_name
real(kind=kind_real), allocatable, intent(inout) :: field(:,:,:)
logical, optional,                 intent(in)    :: ok_to_fail_in

integer :: field_index
logical :: ok_to_fail, have_field

! Set whether the get can fail
ok_to_fail = .false.
if (present(ok_to_fail_in)) ok_to_fail = ok_to_fail_in

! Check field exists
have_field = hasfield(self%fields, field_name)

! Fail conditions
if (.not.have_field) then
  if (.not.ok_to_fail) then
    call abor1_ftn("fv3jedi_fields_mod.get_field: field " // trim(field_name) &
                   // " not found")
  else
    return
  endif
endif

! Get the actual field
call get_field(self%fields, field_name, field)

endsubroutine get_field_return_array_allocatable

! --------------------------------------------------------------------------------------------------

subroutine put_field_(self, field_name, field)

class(fv3jedi_fields), target,     intent(inout) :: self
character(len=*),                  intent(in)    :: field_name
real(kind=kind_real), allocatable, intent(in)    :: field(:,:,:)

integer :: field_index

call put_field(self%fields, field_name, field)

endsubroutine put_field_

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fields_mod
