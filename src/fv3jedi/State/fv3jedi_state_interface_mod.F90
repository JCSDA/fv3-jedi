! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_state_interface_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module,     only: fckit_configuration

! atlas
use atlas_module, only: atlas_fieldset

! oops
use datetime_mod
use duration_mod
use oops_variables_mod

! fv3jedi
use fv3jedi_field_mod,               only: field_clen
use fv3jedi_kinds_mod,               only: kind_real
use fv3jedi_geom_mod,                only: fv3jedi_geom
use fv3jedi_geom_interface_mod,      only: fv3jedi_geom_registry
use fv3jedi_increment_mod,           only: fv3jedi_increment, fv3jedi_increment_registry
use fv3jedi_state_mod,               only: fv3jedi_state

private
public :: fv3jedi_state_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_state

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_state_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_create_c(c_key_self, c_key_geom, c_vars, c_time) &
           bind(c,name='fv3jedi_state_create_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_key_geom !< Geometry
type(c_ptr), value, intent(in) :: c_vars     !< List of all variables
type(c_ptr), value, intent(in) :: c_time     !< Datetime

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables)         :: vars

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%init()
call fv3jedi_state_registry%add(c_key_self)
call fv3jedi_state_registry%get(c_key_self,self)

vars = oops_variables(c_vars)

! Create Fortran pointer to datetime
call c_f_datetime(c_time, self%time)

call self%create(geom, vars)

end subroutine fv3jedi_state_create_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_delete_c(c_key_self) bind(c,name='fv3jedi_state_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)

call self%delete()

call fv3jedi_state_registry%remove(c_key_self)

end subroutine fv3jedi_state_delete_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_zero_c(c_key_self) bind(c,name='fv3jedi_state_zero_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)
call self%zero()

end subroutine fv3jedi_state_zero_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_state_copy_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_state), pointer :: self
type(fv3jedi_state), pointer :: other
call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_rhs,other)

call self%copy(other)

end subroutine fv3jedi_state_copy_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='fv3jedi_state_axpy_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_state), pointer :: self
type(fv3jedi_state), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%accumul(zz,rhs)

end subroutine fv3jedi_state_axpy_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_add_increment_c(c_key_self,c_key_rhs,c_key_geom) &
           bind(c,name='fv3jedi_state_add_increment_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
integer(c_int), intent(in) :: c_key_geom

type(fv3jedi_state), pointer :: self
type(fv3jedi_increment), pointer :: rhs
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)
call fv3jedi_geom_registry%get(c_key_geom,geom)

call self%add_increment(rhs%fields,geom)

end subroutine fv3jedi_state_add_increment_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_change_resol_c(c_key_state,c_key_geom,c_key_rhs,c_key_geom_rhs) &
           bind(c,name='fv3jedi_state_change_resol_f90')

implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_rhs
integer(c_int), intent(in) :: c_key_geom_rhs

type(fv3jedi_state), pointer :: self, other
type(fv3jedi_geom),  pointer :: geom, geom_other

call fv3jedi_state_registry%get(c_key_state,self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_rhs, other)
call fv3jedi_geom_registry%get(c_key_geom_rhs, geom_other)

call self%change_resol(geom, other, geom_other)

end subroutine fv3jedi_state_change_resol_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_analytic_init_c(c_key_state, c_key_geom, c_conf) &
           bind(c,name='fv3jedi_state_analytic_init_f90')

implicit none
integer(c_int),     intent(in)    :: c_key_state
integer(c_int),     intent(in)    :: c_key_geom
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_state), pointer :: f_self
type(fv3jedi_geom), pointer :: f_geom
type(fckit_configuration) :: f_conf

! Linked list
! -----------
call fv3jedi_geom_registry%get(c_key_geom, f_geom)
call fv3jedi_state_registry%get(c_key_state, f_self)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)

! Call implementation
! -------------------
call f_self%analytic_IC(f_geom, f_conf)

end subroutine fv3jedi_state_analytic_init_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_norm_c(c_key_state, prms) bind(c,name='fv3jedi_state_norm_f90')

implicit none
integer(c_int), intent(in) :: c_key_state
real(c_double), intent(inout) :: prms

type(fv3jedi_state), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_state_registry%get(c_key_state,self)

call self%norm(zz)

prms = zz

end subroutine fv3jedi_state_norm_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_getnfieldsncube_c(c_key_self, c_number_fields, c_cube_size) &
           bind(c,name='fv3jedi_state_getnfieldsncube_f90')

implicit none
integer(c_int), intent(in)  :: c_key_self
integer(c_int), intent(out) :: c_number_fields
integer(c_int), intent(out) :: c_cube_size

type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)

c_number_fields = self%nf
c_cube_size = self%npx-1

end subroutine fv3jedi_state_getnfieldsncube_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_getminmaxrms_c(c_key_self, c_f_num, c_f_name_len, c_f_name, c_minmaxrms ) &
           bind(c,name='fv3jedi_state_getminmaxrms_f90')

implicit none
integer(c_int),               intent(in)    :: c_key_self
integer(c_int),               intent(in)    :: c_f_num
integer(c_int),               intent(in)    :: c_f_name_len
character(len=1,kind=c_char), intent(inout) :: c_f_name(c_f_name_len + 1)
real(c_double),               intent(inout) :: c_minmaxrms(3)

type(fv3jedi_state), pointer :: self
character(len=field_clen) :: field_name
integer :: n, trunc_name_len

call fv3jedi_state_registry%get(c_key_self,self)

call self%minmaxrms(c_f_num, field_name, c_minmaxrms)

! logic from oops f_c_string, but without allocation of c string array
trunc_name_len = min(len_trim(field_name), c_f_name_len)
do n = 1,trunc_name_len
  c_f_name(n) = field_name(n:n)
enddo

! if field_name is shorter than C char array, pad with spaces before adding null terminator
do n = trunc_name_len+1,c_f_name_len
  c_f_name(n) = ' '
enddo
c_f_name(c_f_name_len+1) = c_null_char

end subroutine fv3jedi_state_getminmaxrms_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_to_fieldset_c(c_key_self, c_key_geom, c_vars, c_afieldset) &
 & bind (c,name='fv3jedi_state_to_fieldset_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call fv3jedi_state_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%to_fieldset(geom, vars, afieldset)

end subroutine fv3jedi_state_to_fieldset_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_from_fieldset_c(c_key_self, c_key_geom, c_vars, c_afieldset) &
 & bind (c,name='fv3jedi_state_from_fieldset_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call fv3jedi_state_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%from_fieldset(geom, vars, afieldset)

end subroutine fv3jedi_state_from_fieldset_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_synchronize_interface_fields_c(c_key_self, c_key_geom) &
 & bind (c,name='fv3jedi_state_synchronize_interface_fields_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_state_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call self%synchronize_interface_fields(geom)

end subroutine fv3jedi_state_synchronize_interface_fields_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_set_interface_fields_outofdate_c(c_key_self, c_outofdate) &
 & bind (c,name='fv3jedi_state_set_interface_fields_outofdate_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
logical(c_bool), intent(in) :: c_outofdate

type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self, self)

self%interface_fields_are_out_of_date = c_outofdate

end subroutine fv3jedi_state_set_interface_fields_outofdate_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_sersize_c(c_key_self,inc_size) bind(c,name='fv3jedi_state_sersize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: inc_size
type(fv3jedi_state), pointer :: self
integer var, i, j, k

call fv3jedi_state_registry%get(c_key_self, self)

inc_size = 0
do var = 1, self%nf
  inc_size = inc_size + (self%fields(var)%iec-self%fields(var)%isc+1)*&
                        (self%fields(var)%jec-self%fields(var)%jsc+1)*&
                         self%fields(var)%npz
end do

end subroutine fv3jedi_state_sersize_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_serialize_c(c_key_self,c_vsize,c_vect_inc) &
           bind(c,name='fv3jedi_state_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self           !< State
integer(c_int),intent(in) :: c_vsize              !< Size
real(c_double),intent(out) :: c_vect_inc(c_vsize) !< Vector

type(fv3jedi_state),pointer :: self

call fv3jedi_state_registry%get(c_key_self, self)
! Call Fortran
call self%serialize(c_vsize,c_vect_inc)

end subroutine fv3jedi_state_serialize_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_deserialize_c(c_key_self,c_vsize,c_vect_inc,c_index) &
           bind(c,name='fv3jedi_state_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< State
integer(c_int),intent(in) :: c_vsize             !< Size
real(c_double),intent(in) :: c_vect_inc(c_vsize) !< Vector
integer(c_int), intent(inout):: c_index          !< Index

type(fv3jedi_state),pointer :: self

call fv3jedi_state_registry%get(c_key_self, self)

! Call Fortran
call self%deserialize(c_vsize,c_vect_inc,c_index)


end subroutine fv3jedi_state_deserialize_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_update_fields_c(c_key_state, c_key_geom, c_vars) &
  bind(c,name='fv3jedi_state_update_fields_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_state
integer(c_int), intent(in)     :: c_key_geom !< Geometry
type(c_ptr), value, intent(in) :: c_vars     !< List of variables

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables)         :: vars

call fv3jedi_state_registry%get(c_key_state,self)
call fv3jedi_geom_registry%get(c_key_geom, geom)

vars = oops_variables(c_vars)
call self%update_fields(geom, vars)

end subroutine fv3jedi_state_update_fields_c

! --------------------------------------------------------------------------------------------------

end module fv3jedi_state_interface_mod
