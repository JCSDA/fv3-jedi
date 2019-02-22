! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_increment_mod

use iso_c_binding
use config_mod
use datetime_mod

use random_mod
use fckit_mpi_module
use unstructured_grid_mod

use fv3jedi_field_mod,           only: fv3jedi_field, get_field, fields_rms, fields_print, fields_gpnorm
use fv3jedi_constants_mod,       only: rad2deg, constoz, cp, alhl, rgas
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs 
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_state_utils_mod,     only: fv3jedi_state
use fv3jedi_vars_mod,            only: fv3jedi_vars
use fv3jedi_getvalues_mod,       only: getvalues_tl, getvalues_ad

use mpp_domains_mod, only: mpp_global_sum, bitwise_efp_sum, center, east, north, center

implicit none
private
public :: fv3jedi_increment, create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, rms, &
          change_resol, getvalues_tl, getvalues_ad, &
          ug_coord, increment_to_ug, increment_from_ug, dirac, jnormgrad, &
          increment_print

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_vars),      intent(in)    :: vars

integer :: var, vcount

! Total fields
! ------------
self%nf = vars%nv

! Allocate fields structure
! -------------------------
allocate(self%fields(self%nf))

! Loop through and allocate main increment fields
! -----------------------------------------------
vcount = 0
do var = 1, vars%nv
  select case (trim(vars%fldnames(var)))
    case("u","ud")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_d_grid_u_wind', &
           fv3jedi_name = 'ud', units = 'm s-1', staggerloc = north, arraypointer = self%ud)
    case("v","vd")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_d_grid_v_wind', &
           fv3jedi_name = 'vd', units = 'm s-1', staggerloc = east, arraypointer = self%vd)
    case("ua")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_eastward_wind', &
           fv3jedi_name = 'ua', units = 'm s-1', staggerloc = center, arraypointer = self%ua)
    case("va")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_northward_wind', &
           fv3jedi_name = 'va', units = 'm s-1', staggerloc = center, arraypointer = self%va)
    case("t","T")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_air_temperature', &
           fv3jedi_name = 't', units = 'K', staggerloc = center, arraypointer = self%t)
    case("ps")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
           short_name = vars%fldnames(var), long_name = 'increment_of_surface_pressure', &
           fv3jedi_name = 'ps', units = 'Pa', staggerloc = center, arraypointer = self%ps)
    case("q","sphum")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_specific_humidity', &
           fv3jedi_name = 'q', units = 'kg kg-1', staggerloc = center, arraypointer = self%q)
    case("qi","ice_wat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_cloud_liquid_ice', &
           fv3jedi_name = 'qi', units = 'kg kg-1', staggerloc = center, arraypointer = self%qi)
    case("ql","liq_wat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_cloud_liquid_ice_water', &
           fv3jedi_name = 'ql', units = 'kg kg-1', staggerloc = center, arraypointer = self%ql)
    case("o3","o3mr")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_ozone_mass_mixing_ratio', &
           fv3jedi_name = 'o3', units = 'kg kg-1', staggerloc = center, arraypointer = self%o3)
    case("psi")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_stream_function', &
           fv3jedi_name = 'psi', units = 'm+2 s', staggerloc = center, arraypointer = self%psi)
    case("chi")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_velocity_potential', &
           fv3jedi_name = 'chi', units = 'm+2 s', staggerloc = center, arraypointer = self%chi)
    case("tv")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_virtual_temperature', &
           fv3jedi_name = 'tv', units = 'K', staggerloc = center, arraypointer = self%tv)
    case("rh")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_relative_humidity', &
           fv3jedi_name = 'rh', units = '1', staggerloc = center, arraypointer = self%rh)
    case("delp","DELP")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_pressure_thickness', &
           fv3jedi_name = 'delp', units = 'Pa', staggerloc = center, arraypointer = self%delp)
    case("w","W")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_vertical_wind', &
           fv3jedi_name = 'w', units = 'm s-1', staggerloc = center, arraypointer = self%w)
    case("delz","DZ")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_layer_thickness', &
           fv3jedi_name = 'delz', units = 'm', staggerloc = center, arraypointer = self%delz)
      !Fields not in increment
    case default 
      call abor1_ftn("fv3jedi_increment_mod.create: unknown variable "//trim(vars%fldnames(var)))

  end select

enddo

if (vcount .ne. self%nf) &
call abor1_ftn("fv3jedi_increment_mod.create: vcount does not equal self%nf")

self%hydrostatic = .true.
if (associated(self%w) .and. associated(self%delz)) self%hydrostatic = .false.

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
self%f_comm = fckit_mpi_comm()

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

!Nullify pointers
if (associated(self%ud  )) nullify(self%ud  )
if (associated(self%vd  )) nullify(self%vd  )
if (associated(self%ua  )) nullify(self%ua  )
if (associated(self%va  )) nullify(self%va  )
if (associated(self%t   )) nullify(self%t   )
if (associated(self%ps  )) nullify(self%ps  )
if (associated(self%delp)) nullify(self%delp)
if (associated(self%w   )) nullify(self%w   )
if (associated(self%delz)) nullify(self%delz)
if (associated(self%q   )) nullify(self%q   )
if (associated(self%qi  )) nullify(self%qi  )
if (associated(self%ql  )) nullify(self%ql  )
if (associated(self%o3  )) nullify(self%o3  )
if (associated(self%psi )) nullify(self%psi )
if (associated(self%chi )) nullify(self%chi )
if (associated(self%tv  )) nullify(self%tv  )
if (associated(self%rh  )) nullify(self%rh  )

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

subroutine copy(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

logical :: found
integer :: self_var, rhs_var

self%isc            = rhs%isc
self%iec            = rhs%iec
self%jsc            = rhs%jsc
self%jec            = rhs%jec
self%isd            = rhs%isd
self%ied            = rhs%ied
self%jsd            = rhs%jsd
self%jed            = rhs%jed
self%npx            = rhs%npx
self%npy            = rhs%npy
self%npz            = rhs%npz
self%ntiles         = rhs%ntiles
self%ntile          = rhs%ntile
self%hydrostatic    = rhs%hydrostatic
self%tladphystrj    = rhs%tladphystrj
self%calendar_type  = rhs%calendar_type
self%date_init      = rhs%date_init   

!Copy the individual fields
if (.not.allocated(self%fields)) then

  !Direct copy of one increment to another
  self%nf = rhs%nf
  allocate(self%fields(self%nf))

  do self_var = 1, self%nf

    !Field copy
    self%fields(self_var) = rhs%fields(self_var)

  enddo
  
else

  !Increment copy, potentialy with differnt fields
  do self_var = 1, self%nf
    found = .false.
    do rhs_var = 1, rhs%nf
      if (trim(self%fields(self_var)%fv3jedi_name) == trim(rhs%fields(rhs_var)%fv3jedi_name)) then
        self%fields(self_var) = rhs%fields(rhs_var)
        found = .true.
        exit
      endif
    enddo
    if (.not.found) call abor1_ftn("fv3jedi_increment_mod.copy: Field "//&
                    trim(self%fields(self_var)%fv3jedi_name)//" not found in increment being copied from." )

  enddo

endif

!Set pointers
do self_var = 1, self%nf
  select case (trim(self%fields(self_var)%fv3jedi_name))
    case("u","ud")
      call self%fields(self_var)%array_pointer(self%ud)
    case("v","vd")
      call self%fields(self_var)%array_pointer(self%vd)
    case("ua")
      call self%fields(self_var)%array_pointer(self%ua)
    case("va")
      call self%fields(self_var)%array_pointer(self%va)
    case("t","T")
      call self%fields(self_var)%array_pointer(self%t)
    case("delp","DELP")
      call self%fields(self_var)%array_pointer(self%delp)
    case("ps")
      call self%fields(self_var)%array_pointer(self%ps)
    case("w","W")
      call self%fields(self_var)%array_pointer(self%w)
    case("delz","DZ")
      call self%fields(self_var)%array_pointer(self%delz)
    case("q","sphum")
      call self%fields(self_var)%array_pointer(self%q)
    case("qi","ice_wat")
      call self%fields(self_var)%array_pointer(self%qi)
    case("ql","liq_wat")
      call self%fields(self_var)%array_pointer(self%ql)
    case("o3","o3mr")
      call self%fields(self_var)%array_pointer(self%o3)
    case("psi")
      call self%fields(self_var)%array_pointer(self%psi)
    case("chi")
      call self%fields(self_var)%array_pointer(self%chi)
    case("tv")
      call self%fields(self_var)%array_pointer(self%tv)
    case("rh")
      call self%fields(self_var)%array_pointer(self%rh)
    case default 
      call abor1_ftn("fv3jedi_increment_mod.copy: unknown variable "//trim(self%fields(self_var)%fv3jedi_name))
  end select
enddo

end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self,rhs)

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

call checksame(self,rhs)

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

call checksame(self,rhs)

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

call checksame(self,rhs)

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
type(fv3jedi_field), pointer :: field_pointer

do var = 1,self%nf

  !Get pointer to equivalent state
  call get_field(rhs%nf,rhs%fields,self%fields(var)%fv3jedi_name,field_pointer)

  !inc = inc + zz state
  self%fields(var)%array = self%fields(var)%array + zz * field_pointer%array

  !Clear pointer
  nullify(field_pointer)

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

call checksame(self,other)

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
if (self%f_comm%rank() == 0) print*, "Dot product test result: ", zprod

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: check, var

check = (rhs%iec-rhs%isc+1) - (self%iec-self%isc+1)

if (check==0) then
  call checksame(self,rhs)
  do var = 1,self%nf
    self%fields(var)%array = self%fields(var)%array + rhs%fields(var)%array
  enddo
else
   call abor1_ftn("fv3jedi_increment_mod.add_incr not implemented for low res increment yet")
endif

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(self,x1,x2)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_state),     intent(in)    :: x1
type(fv3jedi_state),     intent(in)    :: x2

real(kind=kind_real), allocatable :: x1_ps(:,:), x2_ps(:,:)
integer :: var, check
type(fv3jedi_field), pointer :: x1p, x2p

check = (x1%iec-x1%isc+1) - (x2%iec-x2%isc+1)

call zeros(self)
if (check==0) then

  do var = 1,self%nf
  
    !Get pointer to states
    call get_field(x1%nf,x1%fields,self%fields(var)%fv3jedi_name,x1p)
    call get_field(x2%nf,x2%fields,self%fields(var)%fv3jedi_name,x2p)

    !inc = state - state
    self%fields(var)%array = x1p%array - x2p%array

    !Nullify pointers
    nullify(x1p,x2p)

  enddo

else

   call abor1_ftn("fv3jedi_increment_mod.diff_incr: not implemented for low res increment yet")

endif

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: check

check = (rhs%iec-rhs%isc+1) - (self%iec-self%isc+1)

if (check==0) then
   call copy(self, rhs)
else
   call abor1_ftn("fv3jedi_increment_mod.change_resol: not implmeneted yet")
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
  character(len=255) :: filename

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%read_meta(geom, vdate, self%calendar_type, self%date_init)
    call gfs%read_fields(geom, self%fields)
  elseif (trim(filetype) == 'geos') then
    filename = config_get_string(c_conf,len(filename),"filename")
    call geos%create(geom, 'read', filename)
    call geos%read_time(vdate)
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

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%write_all(geom, self%fields, vdate, self%calendar_type, self%date_init)
  elseif (trim(filetype) == 'geos') then
    call geos%create(geom, 'write')
    call geos%write_all(geom, self%fields, c_conf, vdate)
    call geos%delete()
  else
     call abor1_ftn("fv3jedi_increment_mod.write: restart type not supported")
  endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine increment_print(self)

implicit none
type(fv3jedi_increment), intent(in) :: self

call fields_print(self%nf, self%fields, "Increment")

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

call fields_gpnorm(nf, self%fields, pstat)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_increment), intent(in)  :: self
real(kind=kind_real),    intent(out) :: prms

call fields_rms(self%nf,self%fields,prms)

end subroutine rms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf, geom)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(c_ptr),             intent(in)    :: c_conf
type(fv3jedi_geom),      intent(in)    :: geom

integer :: ndir,idir,ildir,itiledir,var
integer,allocatable :: ixdir(:),iydir(:)
character(len=3) :: idirchar
character(len=32) :: ifdir
logical :: found

! Get Diracs positions
ndir = config_get_int(c_conf,"ndir")
allocate(ixdir(ndir))
allocate(iydir(ndir))

do idir=1,ndir
   write(idirchar,'(i3)') idir
   ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
   iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
end do
ildir = config_get_int(c_conf,"ildir")
ifdir = config_get_string(c_conf,len(ifdir),"ifdir")
itiledir = config_get_int(c_conf,"itiledir")

! Check
if (ndir<1) call abor1_ftn("fv3jedi_increment_mod.dirac: non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%npx)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid ixdir")
endif
if (any(iydir<1).or.any(iydir>geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid iydir")
endif
if ((ildir<1).or.(ildir>self%npz)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac invalid ildir")
endif
if ((itiledir<1).or.(itiledir>6)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid itiledir")
endif

! Setup Diracs
call zeros(self)

! only u, v, T, ps and tracers allowed
do idir=1,ndir

   ! is specified grid point, tile number on this processor
   if (geom%ntile == itiledir .and. &
       ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
       iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then
       ! If so, perturb desired increment and level
       found = .false.
       do var = 1,self%nf
         if (trim(self%fields(var)%fv3jedi_name) == trim(ifdir)) then
           found = .true.
           self%fields(var)%array(ixdir(idir),iydir(idir),ildir) = 1.0
         endif
       enddo
       if (.not.found) call abor1_ftn("fv3jedi_increment_mod.dirac:: field "&
                                       //trim(ifdir)//" not found")
   endif
end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine ug_size(self, ug)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug

! Set number of grids
if (ug%colocated==1) then
   ! Colocatd
   ug%ngrid = 1
else
   ! Not colocated
   ug%ngrid = 1
   call abor1_ftn("fv3jedi_increment_mod.ug_size: Uncolocated grids not coded yet")
end if

! Allocate grid instances
if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

! Set local number of points
ug%grid(1)%nmga = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1) 

! Set number of levels
ug%grid(1)%nl0 = self%npz

! Set number of variables
ug%grid(1)%nv = self%nf 

! Set number of timeslots
ug%grid(1)%nts = 1

end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine ug_coord(self, ug, colocated, geom)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug
integer,                 intent(in)    :: colocated
type(fv3jedi_geom),      intent(in)    :: geom

integer :: imga,jx,jy,jl
real(kind=kind_real) :: sigmaup,sigmadn

! Copy colocated
ug%colocated = colocated

! Define size
call ug_size(self, ug)

! Allocate unstructured grid coordinates
call allocate_unstructured_grid_coord(ug)

if (ug%colocated==1) then
  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      ug%grid(1)%lon(imga) = rad2deg*geom%grid_lon(jx,jy)
      ug%grid(1)%lat(imga) = rad2deg*geom%grid_lat(jx,jy)
      ug%grid(1)%area(imga) = geom%area(jx,jy)
      do jl=1,self%npz
        sigmaup = geom%ak(jl+1)/101300.0+geom%bk(jl+1) ! si are now sigmas
        sigmadn = geom%ak(jl  )/101300.0+geom%bk(jl  )
        ug%grid(1)%vunit(imga,jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
        ug%grid(1)%lmask(imga,jl) = .true.
      enddo
    enddo
  enddo 
endif

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine increment_to_ug(self, ug, colocated)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug
integer,                 intent(in)    :: colocated

integer :: var,imga,jx,jy,jl

! Copy colocated
ug%colocated = colocated

! Define size
call ug_size(self, ug)

! Allocate unstructured grid increment
call allocate_unstructured_grid_field(ug)

! Copy increment

ug%grid(1)%fld = 0.0_kind_real

if (ug%colocated==1) then

  do var = 1,self%nf
    imga = 0
    do jy=self%jsc,self%jec
      do jx=self%isc,self%iec
        imga = imga+1
        do jl=1,self%fields(var)%npz
          ug%grid(1)%fld(imga,jl,var,1) = self%fields(var)%array(jx,jy,jl)
        enddo
      enddo
    enddo
  enddo

endif

end subroutine increment_to_ug

! -----------------------------------------------------------------------------

subroutine increment_from_ug(self, ug)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(unstructured_grid), intent(in)    :: ug

integer :: imga,jx,jy,jl,var

! Copy increment

if (ug%colocated==1) then

  do var = 1,self%nf
    imga = 0
    do jy=self%jsc,self%jec
      do jx=self%isc,self%iec
        imga = imga+1
        do jl=1,self%fields(var)%npz
          self%fields(var)%array(jx,jy,jl) = ug%grid(1)%fld(imga,jl,var,1)
        enddo
      enddo
    enddo
  enddo

endif

end subroutine increment_from_ug

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

!Code to compute a vector norm for an increment, e.g. the energy norm for FSOI

isc = self%isc
iec = self%iec
jsc = self%jsc
jec = self%jec
npz = self%npz

! Constants
! ---------
tref = config_get_real(c_conf,"Tref")
qeps = config_get_real(c_conf,"qepsilon")
pref = config_get_real(c_conf,"pref")

Ufac = 0.5_kind_real
Tfac = 0.5_kind_real*cp/tref
qfac = 0.5_kind_real*qeps*alhl*alhl/(cp*tref)
pfac = 0.5_kind_real*Rgas*tref/pref**2

! Compute grid weighting based on volume
! --------------------------------------

global_area = mpp_global_sum(geom%domain, geom%area, flags=bitwise_efp_sum)

allocate(ref_ps(isc:iec,jsc:jec))
ref_ps = sum(ref%delp,3)

allocate(cellweight(isc:iec,jsc:jec,1:npz))
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      cellweight(i,j,k) = (ref%delp(i,j,k)/ref_ps(i,j)) * geom%area(i,j)/global_area
    enddo
  enddo
enddo

!ua
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%ua(i,j,k) = Ufac * 2.0_kind_real * ref%ua(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!va
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%va(i,j,k) = Ufac * 2.0_kind_real * ref%va(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!T
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%t(i,j,k) = Tfac * 2.0_kind_real * ref%t(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!q
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%q(i,j,k) = qfac * 2.0_kind_real * ref%q(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!ps
if (associated(self%ps)) then
  do j = jsc,jec
    do i = isc,iec
      self%ps(i,j,1) = pfac * 2.0_kind_real * ref_ps (i,j) * cellweight(i,j,npz) &
                                          / (ref%delp(i,j,npz)/ref_ps(i,j))
    enddo
  enddo
else
  call abor1_ftn("fv3jedi_increment_mod.jnormgrad: Ps must be in the increment")
endif

deallocate(cellweight)
deallocate(ref_ps)

end subroutine jnormgrad

! ------------------------------------------------------------------------------

subroutine checksame(self,other)

implicit none
type(fv3jedi_increment), intent(in) :: self
type(fv3jedi_increment), intent(in) :: other

integer :: var

if (self%nf .ne. other%nf) then
  call abor1_ftn("fv3jedi_increment_mod.checksame: Different number of fields")
endif

do var = 1,self%nf
  if (self%fields(var)%fv3jedi_name .ne. other%fields(var)%fv3jedi_name) then
      call abor1_ftn("fv3jedi_increment_mod.checksame: field "//trim(self%fields(var)%fv3jedi_name)//&
                     " not in the equivalent position in the right hand side")
  endif
enddo

end subroutine checksame

! ------------------------------------------------------------------------------

end module fv3jedi_increment_mod
