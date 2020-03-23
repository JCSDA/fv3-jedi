! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_increment_mod

use atlas_module
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use oops_variables_mod, only: oops_variables

use random_mod
use fckit_mpi_module

use fields_metadata_mod, only: field_metadata

use fv3jedi_field_mod
use fv3jedi_constants_mod,       only: constoz, cp, alhl, rgas
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_iter_mod,       only: fv3jedi_geom_iter
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_interpolation_mod,   only: field2field_interp
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_state_utils_mod,     only: fv3jedi_state
use fv3jedi_getvalues_mod,       only: getvalues_tl, getvalues_ad

use wind_vt_mod, only: d2a

use mpp_domains_mod, only: mpp_global_sum, bitwise_efp_sum, center, east, north, center

implicit none
private
public :: fv3jedi_increment, create, delete, zeros, random, set_atlas, to_atlas, from_atlas, copy, &
          self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          dot_prod, diff_incr, &
          read_file, write_file, gpnorm, rms, &
          change_resol, getvalues_tl, getvalues_ad, &
          dirac, jnormgrad, &
          fv3jedi_increment_serialize, fv3jedi_increment_deserialize, &
          fv3jedi_getpoint, fv3jedi_setpoint, &
          increment_print

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(oops_variables),    intent(in)    :: vars

integer :: var, fc
type(field_metadata) :: fmd

! Allocate fields structure
! -------------------------
self%nf = vars%nvars()
allocate(self%fields(self%nf))

! Loop through and allocate main state fields
! -------------------------------------------
fc = 0
do var = 1, vars%nvars()

  fmd = geom%fields%get_field(trim(vars%variable(var)))

  fc=fc+1;
  call self%fields(fc)%allocate_field(geom%isc, geom%iec, geom%jsc, geom%jec, &
                                      fmd%levels, &
                                      short_name   = trim(fmd%field_io_name), &
                                      long_name    = "increment_of_"//trim(fmd%long_name), &
                                      fv3jedi_name = trim(fmd%field_name), &
                                      units        = fmd%units, &
                                      space        = trim(fmd%space), &
                                      staggerloc   = trim(fmd%stagger_loc), &
                                      tracer       = fmd%tracer, &
                                      integerfield = trim(fmd%array_kind)=='integer')

enddo

if (fc .ne. self%nf) &
call abor1_ftn("fv3jedi_increment_mod.create: fc does not equal self%nf")

self%hydrostatic = .true.
if (has_field(self%fields, 'delz') .and. has_field(self%fields, 'w')) self%hydrostatic = .false.

! Initialize all domain arrays to zero
call zeros(self)

! For convenience
self%isc    = geom%isc
self%iec    = geom%iec
self%jsc    = geom%jsc
self%jec    = geom%jec
self%isd    = geom%isd
self%ied    = geom%ied
self%jsd    = geom%jsd
self%jed    = geom%jed
self%npx    = geom%npx
self%npy    = geom%npy
self%npz    = geom%npz
self%ntile  = geom%ntile
self%ntiles = geom%ntiles

! Pointer to fv3jedi communicator
self%f_comm = geom%f_comm

! Check winds
if (has_field(self%fields, 'ua') .and. .not.has_field(self%fields, 'va')) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid u but not v")
if (.not.has_field(self%fields, 'ua') .and. has_field(self%fields, 'va')) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid v but not u")
if (has_field(self%fields, 'ud') .and. .not.has_field(self%fields, 'vd')) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid u but not v")
if (.not.has_field(self%fields, 'ud') .and. has_field(self%fields, 'vd')) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid v but not u")

self%have_agrid = .false.
self%have_dgrid = .false.
if (has_field(self%fields, 'ua')) self%have_agrid = .true.
if (has_field(self%fields, 'ud')) self%have_dgrid = .true.

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

!Deallocate fields
do var = 1, self%nf
  call self%fields(var)%deallocate_field()
enddo
deallocate(self%fields)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 0.0_kind_real
enddo

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine ones(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 1.0_kind_real
enddo

end subroutine ones

! ------------------------------------------------------------------------------

subroutine random(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var
integer, parameter :: rseed = 7

do var = 1,self%nf
  call normal_distribution(self%fields(var)%array, 0.0_kind_real, 1.0_kind_real, rseed)
enddo

end subroutine random

! ------------------------------------------------------------------------------

subroutine set_atlas(self, geom, vars, vdate, afieldset)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(oops_variables),    intent(in)    :: vars
type(datetime),          intent(in)    :: vdate
type(atlas_fieldset),    intent(inout) :: afieldset

integer :: jvar, jf, jl
logical :: var_found
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

do jvar = 1,vars%nvars()
  var_found = .false.
  do jf = 1,self%nf
    if (trim(vars%variable(jvar))==trim(self%fields(jf)%short_name)) then
      ! Get or create field
      fieldname = trim(vars%variable(jvar))//'_'//sdate
      if (afieldset%has_field(trim(fieldname))) then
        ! Get field
        afield = afieldset%field(trim(fieldname))
      else
        ! Create field
        afield = geom%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%fields(jvar)%npz)

        ! Add field
        call afieldset%add(afield)
      endif

      ! Release pointer
      call afield%final()

      ! Set flag
      var_found = .true.
      exit
    end if
  end do
  if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
end do

end subroutine set_atlas

! ------------------------------------------------------------------------------

subroutine to_atlas(self, geom, vars, vdate, afieldset)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(oops_variables),    intent(in)    :: vars
type(datetime),          intent(in)    :: vdate
type(atlas_fieldset),    intent(inout) :: afieldset

integer :: jvar, jf, jl
real(kind=kind_real), pointer :: real_ptr_2(:,:)
logical :: var_found
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

do jvar = 1,vars%nvars()
  var_found = .false.
  do jf = 1,self%nf
    if (trim(vars%variable(jvar))==trim(self%fields(jf)%short_name)) then
      ! Get or create field
      fieldname = trim(vars%variable(jvar))//'_'//sdate
      if (afieldset%has_field(trim(fieldname))) then
        ! Get field
        afield = afieldset%field(trim(fieldname))
      else
        ! Create field
        afield = geom%afunctionspace%create_field(name=trim(fieldname),kind=atlas_real(kind_real),levels=self%fields(jvar)%npz)

        ! Add field
        call afieldset%add(afield)
      endif

      ! Copy data
      call afield%data(real_ptr_2)
      do jl=1,self%fields(jf)%npz
        real_ptr_2(jl,:) = pack(self%fields(jf)%array(geom%isc:geom%iec,geom%jsc:geom%jec,jl),.true.)
      enddo

      ! Release pointer
      call afield%final()

      ! Set flag
      var_found = .true.
      exit
    end if
  end do
  if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
end do

end subroutine to_atlas

! ------------------------------------------------------------------------------

subroutine from_atlas(self, geom, vars, vdate, afieldset)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(oops_variables),    intent(in)    :: vars
type(datetime),          intent(in)    :: vdate
type(atlas_fieldset),    intent(in)    :: afieldset

integer :: jvar, jf, jl
real(kind=kind_real), pointer :: real_ptr_2(:,:)
logical :: umask(geom%isc:geom%iec,geom%jsc:geom%jec),var_found
character(len=20) :: sdate
character(len=1024) :: fieldname
type(atlas_field) :: afield

! Set date
call datetime_to_string(vdate,sdate)

! Initialization
umask = .true.

do jvar = 1,vars%nvars()
  var_found = .false.
  do jf = 1,self%nf
    if (trim(vars%variable(jvar))==trim(self%fields(jf)%short_name)) then
      ! Get field
      fieldname = trim(vars%variable(jvar))//'_'//sdate
      afield = afieldset%field(trim(fieldname))

      ! Copy data
      call afield%data(real_ptr_2)
      do jl=1,self%fields(jf)%npz
        self%fields(jf)%array(geom%isc:geom%iec,geom%jsc:geom%jec,jl) = unpack(real_ptr_2(jl,:), &
      & umask,self%fields(jf)%array(geom%isc:geom%iec,geom%jsc:geom%jec,jl))
      enddo

      ! Release pointer
      call afield%final()

      ! Set flag
      var_found = .true.
      exit
    end if
  end do
  if (.not.var_found) call abor1_ftn('variable '//trim(vars%variable(jvar))//' not found in increment')
end do

end subroutine from_atlas

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

do var = 1, self%nf
  self%fields(var) = rhs%fields(var)
enddo

self%calendar_type  = rhs%calendar_type
self%date_init      = rhs%date_init

end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_add")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + rhs%fields(var)%array
enddo

end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_schur")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array * rhs%fields(var)%array
enddo

end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.self_sub")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array - rhs%fields(var)%array
enddo

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz

integer :: var

do var = 1,self%nf
  self%fields(var)%array = zz * self%fields(var)%array
enddo

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy_inc(self,zz,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.axpy_inc")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine axpy_inc

! ------------------------------------------------------------------------------

subroutine axpy_state(self,zz,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz
type(fv3jedi_state),     intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.axpy_state")

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine axpy_state

! ------------------------------------------------------------------------------

subroutine dot_prod(self,other,zprod)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(fv3jedi_increment), intent(in)    :: other
real(kind=kind_real),    intent(inout) :: zprod

real(kind=kind_real) :: zp
integer :: i,j,k
integer :: var

call checksame(self%fields,other%fields,"fv3jedi_increment_mod.dot_prod")

zp=0.0_kind_real
do var = 1,self%nf
  do k = 1,self%fields(var)%npz
    do j = self%jsc,self%jec
      do i = self%isc,self%iec
        zp = zp + self%fields(var)%array(i,j,k) * other%fields(var)%array(i,j,k)
      enddo
    enddo
  enddo
enddo

!Get global dot product
call self%f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
!if (self%f_comm%rank() == 0) print*, "Dot product test result: ", zprod

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine diff_incr(self,x1,x2,geom)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_state),     intent(in)    :: x1
type(fv3jedi_state),     intent(in)    :: x2
type(fv3jedi_geom),      intent(inout) :: geom

integer :: var
type(fv3jedi_field), pointer :: x1p, x2p
real(kind=kind_real), allocatable :: x1_ua(:,:,:), x1_va(:,:,:)
real(kind=kind_real), allocatable :: x2_ua(:,:,:), x2_va(:,:,:)
real(kind=kind_real), pointer :: x1_ud(:,:,:), x1_vd(:,:,:)
real(kind=kind_real), pointer :: x2_ud(:,:,:), x2_vd(:,:,:)
real(kind=kind_real), pointer :: x1_delp(:,:,:)
real(kind=kind_real), pointer :: x2_delp(:,:,:)
real(kind=kind_real), allocatable :: x1_ps(:,:,:), x2_ps(:,:,:)

! Make sure two states have same fields and in same position
call checksame(x1%fields,x2%fields,"fv3jedi_state_mod.diff_incr x1 vs x2")

!D-Grid to A-Grid (if needed)
if (self%have_agrid) then
  allocate(x1_ua(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  allocate(x1_va(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  allocate(x2_ua(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  allocate(x2_va(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  if (x1%have_agrid) then
    call copy_field_array(x1%fields, 'ua', x1_ua)
    call copy_field_array(x1%fields, 'va', x1_va)
    call copy_field_array(x2%fields, 'ua', x2_ua)
    call copy_field_array(x2%fields, 'va', x2_va)
  elseif (x1%have_dgrid) then
    call pointer_field_array(x1%fields, 'ud', x1_ud)
    call pointer_field_array(x1%fields, 'vd', x1_vd)
    call pointer_field_array(x2%fields, 'ud', x2_ud)
    call pointer_field_array(x2%fields, 'vd', x2_vd)
    call d2a(geom, x1_ud, x1_vd, x1_ua, x1_va)
    call d2a(geom, x2_ud, x2_vd, x2_ua, x2_va)
  else
    call abor1_ftn("fv3jedi_state_mod.diff_incr: no way to determine A grid winds")
  endif
endif

!delp to ps
if (has_field(self%fields, 'ps')) then

  allocate(x1_ps(x1%isc:x1%iec,x1%jsc:x1%jec,1))
  allocate(x2_ps(x2%isc:x2%iec,x2%jsc:x2%jec,1))

  if (has_field(x1%fields, 'delp')) then
    call pointer_field_array(x1%fields, 'delp', x1_delp)
    x1_ps(:,:,1) = sum(x1_delp,3)
  elseif (has_field(x1%fields, 'ps')) then
    call copy_field_array(x1%fields, 'ps', x1_ps)
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x1")
  endif

  if (has_field(x2%fields, 'delp')) then
    call pointer_field_array(x2%fields, 'delp', x2_delp)
    x2_ps(:,:,1) = sum(x2_delp,3)
  elseif (has_field(x2%fields, 'ps')) then
    call copy_field_array(x2%fields, 'ps', x2_ps)
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x2")
  endif

endif

do var = 1,self%nf

  !A-Grid winds can be a special case
  if (self%fields(var)%fv3jedi_name == 'ua') then

    self%fields(var)%array = x1_ua - x2_ua

  elseif (self%fields(var)%fv3jedi_name == 'va') then

    self%fields(var)%array = x1_va - x2_va

  !Ps can be a special case
  elseif (self%fields(var)%fv3jedi_name == 'ps') then

    self%fields(var)%array = x1_ps - x2_ps

  else

    !Get pointer to states
    call pointer_field(x1%fields, self%fields(var)%fv3jedi_name, x1p)
    call pointer_field(x2%fields, self%fields(var)%fv3jedi_name, x2p)

    !inc = state - state
    self%fields(var)%array = x1p%array - x2p%array

    !Nullify pointers
    nullify(x1p,x2p)

  endif

enddo

if (allocated(x1_ua)) deallocate(x1_ua)
if (allocated(x1_va)) deallocate(x1_va)
if (allocated(x2_ua)) deallocate(x2_ua)
if (allocated(x2_va)) deallocate(x2_va)
if (allocated(x1_ps)) deallocate(x1_ps)
if (allocated(x2_ps)) deallocate(x2_ps)

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(self,geom,rhs,geom_rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: rhs
type(fv3jedi_geom),      intent(inout) :: geom_rhs

! Interpolation
integer :: var
type(field2field_interp) :: interp
logical :: integer_interp = .false.

call checksame(self%fields,rhs%fields,"fv3jedi_increment_mod.change_resol")

if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1)==0) then

  call copy(self, rhs)

else

  ! Check if integer interp needed
  do var = 1, self%nf
    if (rhs%fields(var)%integerfield) integer_interp = .true.
  enddo

  call interp%create(geom%interp_method, integer_interp, geom_rhs, geom)
  call interp%apply(self%nf, geom_rhs, rhs%fields, geom, self%fields)
  call interp%delete()

  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init

endif

end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(geom, self, c_conf, vdate)

  implicit none
  type(fv3jedi_geom),      intent(inout) :: geom     !< Geometry
  type(fv3jedi_increment), intent(inout) :: self     !< Increment
  type(c_ptr),             intent(in)    :: c_conf   !< Configuration
  type(datetime),          intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype
  integer :: psinfile
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filetype",str)
  filetype = str
  deallocate(str)

  psinfile = 0
  if (f_conf%has("psinfile")) then
    call f_conf%get_or_die("psinfile",psinfile)
  endif

  if (trim(filetype) == 'gfs') then

    call gfs%setup(f_conf,psinfile = psinfile)
    call gfs%read_meta(geom, vdate, self%calendar_type, self%date_init)
    call gfs%read_fields(geom, self%fields)

  elseif (trim(filetype) == 'geos') then

    call geos%setup(geom, self%fields, vdate, 'read', f_conf)
    call geos%read_meta(geom, vdate, self%calendar_type, self%date_init)
    call geos%read_fields(geom, self%fields)
    call geos%delete()

  else

     call abor1_ftn("fv3jedi_increment_mod.read: restart type not supported")

  endif

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, self, c_conf, vdate)

  implicit none

  type(fv3jedi_geom),      intent(inout) :: geom     !< Geometry
  type(fv3jedi_increment), intent(in)    :: self     !< Increment
  type(c_ptr),             intent(in)    :: c_conf   !< Configuration
  type(datetime),          intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filetype",str)
  filetype = str
  deallocate(str)

  if (trim(filetype) == 'gfs') then

    call gfs%setup(f_conf)
    call gfs%write_all(geom, self%fields, vdate, self%calendar_type, self%date_init)

  elseif (trim(filetype) == 'geos') then

    call geos%setup(geom, self%fields, vdate, 'write', f_conf)
    call geos%write_all(geom, self%fields, vdate)
    call geos%delete()

  else

     call abor1_ftn("fv3jedi_increment_mod.write: restart type not supported")

  endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine increment_print(self)

implicit none
type(fv3jedi_increment), intent(in) :: self

call fields_print(self%nf, self%fields, "Increment", self%f_comm)

end subroutine increment_print

! ------------------------------------------------------------------------------

subroutine gpnorm(self, nf, pstat)

implicit none
type(fv3jedi_increment), intent(in)    :: self
integer,                 intent(in)    :: nf
real(kind=kind_real),    intent(inout) :: pstat(3, nf)

if (nf .ne. self%nf) then
  call abor1_ftn("fv3jedi_increment.gpnorm: nf passed in does not match expeted nf")
endif

call fields_gpnorm(nf, self%fields, pstat, self%f_comm)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_increment), intent(in)  :: self
real(kind=kind_real),    intent(out) :: prms

call fields_rms(self%nf, self%fields, prms, self%f_comm)

end subroutine rms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf, geom)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(c_ptr),             intent(in)    :: c_conf
type(fv3jedi_geom),      intent(in)    :: geom

integer :: ndir,idir,var

integer, allocatable :: ixdir(:),iydir(:),ildir(:),itdir(:)
character(len=32), allocatable :: ifdir(:)

logical :: found
type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str
character(kind=c_char,len=255), allocatable :: ifdir_array(:)
integer(c_size_t),parameter :: csize = 255

f_conf = fckit_configuration(c_conf)

! Get Diracs positions
call f_conf%get_or_die("ndir",ndir)

allocate(ixdir(ndir))
allocate(iydir(ndir))
allocate(ildir(ndir))
allocate(itdir(ndir))

if ((f_conf%get_size("ixdir")/=ndir) .or. &
    (f_conf%get_size("iydir")/=ndir) .or. &
    (f_conf%get_size("ildir")/=ndir) .or. &
    (f_conf%get_size("itdir")/=ndir) .or. &
    (f_conf%get_size("ifdir")/=ndir)) &
  call abor1_ftn("fv3jedi_increment_mod.diracL=: dimension inconsistency")

call f_conf%get_or_die("ixdir",ixdir)
call f_conf%get_or_die("iydir",iydir)
call f_conf%get_or_die("ildir",ildir)
call f_conf%get_or_die("itdir",itdir)

call f_conf%get_or_die("ifdir",csize,ifdir_array)
ifdir = ifdir_array
deallocate(ifdir_array)

! Setup Diracs
call zeros(self)

! only u, v, T, ps and tracers allowed
do idir=1,ndir

  found = .false.

  ! is specified grid point, tile number on this processor
  if (geom%ntile == itdir(idir) .and. &
    ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
    iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then

    ! If so, perturb desired increment and level
    do var = 1,self%nf
      if (trim(self%fields(var)%fv3jedi_name) == trim(ifdir(idir))) then
        found = .true.
        self%fields(var)%array(ixdir(idir),iydir(idir),ildir) = 1.0
      endif
    enddo

    if (.not.found) call abor1_ftn("fv3jedi_increment_mod.dirac: dirac not found")

  endif

enddo

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine jnormgrad(self,geom,ref,c_conf)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_state),     intent(in)    :: ref !To linearize around if nl
type(c_ptr),             intent(in)    :: c_conf

integer :: i,j,k
integer :: isc,iec,jsc,jec,npz
real(kind=kind_real), allocatable :: cellweight(:,:,:), ref_ps(:,:)

real(kind=kind_real) :: global_area

real(kind=kind_real) :: Ufac
real(kind=kind_real) :: Tfac, Tref
real(kind=kind_real) :: qfac, qeps
real(kind=kind_real) :: pfac, pref

real(kind=kind_real), pointer :: ua(:,:,:)
real(kind=kind_real), pointer :: va(:,:,:)
real(kind=kind_real), pointer :: t (:,:,:)
real(kind=kind_real), pointer :: q (:,:,:)
real(kind=kind_real), pointer :: ps(:,:,:)

real(kind=kind_real), pointer :: ref_ua  (:,:,:)
real(kind=kind_real), pointer :: ref_va  (:,:,:)
real(kind=kind_real), pointer :: ref_t   (:,:,:)
real(kind=kind_real), pointer :: ref_q   (:,:,:)
real(kind=kind_real), pointer :: ref_delp(:,:,:)

type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str

!Code to compute a vector norm for an increment, e.g. the energy norm for FSOI

isc = self%isc
iec = self%iec
jsc = self%jsc
jec = self%jec
npz = self%npz

! Constants
! ---------
f_conf = fckit_configuration(c_conf)

call f_conf%get_or_die("Tref",tref)
call f_conf%get_or_die("qepsilon",qeps)
call f_conf%get_or_die("pref",pref)

Ufac = 0.5_kind_real
Tfac = 0.5_kind_real*cp/tref
qfac = 0.5_kind_real*qeps*alhl*alhl/(cp*tref)
pfac = 0.5_kind_real*Rgas*tref/pref**2

! Compute grid weighting based on volume
! --------------------------------------

global_area = mpp_global_sum(geom%domain, geom%area, flags=bitwise_efp_sum)

call pointer_field_array(ref%fields, 'ua'  , ref_ua)
call pointer_field_array(ref%fields, 'va'  , ref_va)
call pointer_field_array(ref%fields, 't'   , ref_t)
call pointer_field_array(ref%fields, 'q'   , ref_q)
call pointer_field_array(ref%fields, 'delp', ref_delp)

allocate(ref_ps(isc:iec,jsc:jec))
ref_ps = sum(ref_delp,3)

allocate(cellweight(isc:iec,jsc:jec,1:npz))
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      cellweight(i,j,k) = (ref_delp(i,j,k)/ref_ps(i,j)) * geom%area(i,j)/global_area
    enddo
  enddo
enddo

!ua
call pointer_field_array(self%fields, 'ua', ua)
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      ua(i,j,k) = Ufac * 2.0_kind_real * ref_ua(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!va
call pointer_field_array(self%fields, 'va', va)
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      va(i,j,k) = Ufac * 2.0_kind_real * ref_va(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!T
call pointer_field_array(self%fields, 't', t)
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      t(i,j,k) = Tfac * 2.0_kind_real * ref_t(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!q
call pointer_field_array(self%fields, 'q', q)
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      q(i,j,k) = qfac * 2.0_kind_real * ref_q(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!ps
if (has_field(self%fields, 'ps')) then
  call pointer_field_array(self%fields, 'ps', ps)
  do j = jsc,jec
    do i = isc,iec
      ps(i,j,1) = pfac * 2.0_kind_real * ref_ps (i,j) * cellweight(i,j,npz) &
                                          / (ref_delp(i,j,npz)/ref_ps(i,j))
    enddo
  enddo
else
  call abor1_ftn("fv3jedi_increment_mod.jnormgrad: Ps must be in the increment")
endif

deallocate(cellweight)
deallocate(ref_ps)

end subroutine jnormgrad

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_serialize(self,vsize,vect_inc)

implicit none

! Passed variables
type(fv3jedi_increment),intent(in) :: self       !< Increment
integer,intent(in) :: vsize                      !< Size
real(kind_real),intent(out) :: vect_inc(vsize)   !< Vector

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
      end do
    end do
  end do
end do



end subroutine fv3jedi_increment_serialize

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_deserialize(self,vsize,vect_inc,index)

implicit none

! Passed variables
type(fv3jedi_increment),intent(inout) :: self         !< Increment
integer,intent(in) :: vsize                           !< Size
real(kind_real),intent(in) :: vect_inc(vsize)         !< Vector
integer,intent(inout) :: index                        !< Index

! Local variables
integer :: ind, var, i, j, k

! Copy
do var = 1, self%nf
  do k = 1,self%fields(var)%npz
    do j = self%fields(var)%jsc,self%fields(var)%jec
      do i = self%fields(var)%isc,self%fields(var)%iec
        self%fields(var)%array(i, j, k) = vect_inc(index + 1)
        index = index + 1
      end do
    end do
  end do
end do


end subroutine fv3jedi_increment_deserialize

! ------------------------------------------------------------------------------
subroutine fv3jedi_getpoint(self, geoiter, values)

implicit none

type(fv3jedi_increment),           intent(   in) :: self
type(fv3jedi_geom_iter),           intent(   in) :: geoiter
real(kind=kind_real),              intent(inout) :: values(:)

integer :: var, nz, ii

ii = 0
do var = 1,self%nf
  nz = self%fields(var)%npz
  values(ii+1:ii+nz) = self%fields(var)%array(geoiter%iind, geoiter%jind,:)
  ii = ii + nz
enddo

end subroutine fv3jedi_getpoint

! ------------------------------------------------------------------------------
subroutine fv3jedi_setpoint(self, geoiter, values)

implicit none

! Passed variables
type(fv3jedi_increment),        intent(inout) :: self
type(fv3jedi_geom_iter),        intent(   in) :: geoiter
real(kind=kind_real),           intent(   in) :: values(:)

integer :: var, nz, ii

ii = 0
do var = 1,self%nf
  nz = self%fields(var)%npz
  self%fields(var)%array(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+nz)
  ii = ii + nz
enddo

end subroutine fv3jedi_setpoint

! ------------------------------------------------------------------------------
end module fv3jedi_increment_mod
