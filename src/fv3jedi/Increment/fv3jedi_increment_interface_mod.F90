! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_increment_interface_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! atlas
use atlas_module, only: atlas_fieldset

! oops
use datetime_mod
use duration_mod
use oops_variables_mod

! fv3jedi
use fv3jedi_field_mod,           only: field_clen
use fv3jedi_geom_iter_mod,       only: fv3jedi_geom_iter, fv3jedi_geom_iter_registry
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry
use fv3jedi_increment_mod,       only: fv3jedi_increment, fv3jedi_increment_registry
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_state_mod,           only: fv3jedi_state

implicit none
private
public :: fv3jedi_increment_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_create_c(c_key_self, c_key_geom, c_vars, c_time) &
           bind(c,name='fv3jedi_increment_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom !< Geometry
type(c_ptr), value, intent(in)    :: c_vars     !< List of variables
type(c_ptr), value, intent(in)    :: c_time     !< Datetime

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(fckit_configuration)    :: f_conf

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%init()
call fv3jedi_increment_registry%add(c_key_self)
call fv3jedi_increment_registry%get(c_key_self, self)

vars = oops_variables(c_vars)

! Create Fortran pointer to datetime
call c_f_datetime(c_time, self%time)

call self%create(geom, vars)

end subroutine fv3jedi_increment_create_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_delete_c(c_key_self) bind(c,name='fv3jedi_increment_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)

call self%delete()

call fv3jedi_increment_registry%remove(c_key_self)

end subroutine fv3jedi_increment_delete_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_zero_c(c_key_self) bind(c,name='fv3jedi_increment_zero_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)
call self%zero()

end subroutine fv3jedi_increment_zero_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_ones_c(c_key_self) bind(c,name='fv3jedi_increment_ones_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)
call self%ones()

end subroutine fv3jedi_increment_ones_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_dirac_c(c_key_self, c_conf, c_key_geom) &
           bind(c,name='fv3jedi_increment_dirac_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_conf
integer(c_int),     intent(in) :: c_key_geom

type(fv3jedi_increment), pointer :: f_self
type(fv3jedi_geom),      pointer :: f_geom
type(fckit_configuration)        :: f_conf

! Linked list
! -----------
call fv3jedi_increment_registry%get(c_key_self, f_self)
call fv3jedi_geom_registry%get(c_key_geom, f_geom)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)

! Call implementation
! -------------------
call f_self%dirac(f_conf, f_geom)

end subroutine fv3jedi_increment_dirac_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_random_c(c_key_self) bind(c,name='fv3jedi_increment_random_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)
call self%random()

end subroutine fv3jedi_increment_random_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_update_fields_c(c_key_self, c_key_geom, c_vars) &
  bind(c,name='fv3jedi_increment_update_fields_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self
integer(c_int), intent(in)     :: c_key_geom !< Geometry
type(c_ptr), value, intent(in) :: c_vars     !< List of variables

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),      pointer :: geom
type(oops_variables)             :: vars

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)

vars = oops_variables(c_vars)
call self%update_fields(geom, vars)

end subroutine fv3jedi_increment_update_fields_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_set_atlas_c(c_key_self, c_key_geom, c_vars, c_afieldset) &
 & bind (c,name='fv3jedi_increment_set_atlas_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%set_atlas(geom, vars, afieldset)

end subroutine fv3jedi_increment_set_atlas_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_to_atlas_c(c_key_self, c_key_geom, c_vars, c_afieldset) &
 & bind (c,name='fv3jedi_increment_to_atlas_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%to_atlas(geom, vars, afieldset)

end subroutine fv3jedi_increment_to_atlas_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_from_atlas_c(c_key_self, c_key_geom, c_vars, c_afieldset) &
 & bind (c,name='fv3jedi_increment_from_atlas_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_variables) :: vars
type(atlas_fieldset) :: afieldset

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%from_atlas(geom, vars, afieldset)

end subroutine fv3jedi_increment_from_atlas_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_copy_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: other
call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_rhs,other)

call self%copy(other)

end subroutine fv3jedi_increment_copy_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_add_c(c_key_self,c_key_rhs) &
           bind(c,name='fv3jedi_increment_self_add_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self%self_add(rhs)

end subroutine fv3jedi_increment_self_add_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_schur_c(c_key_self,c_key_rhs) &
           bind(c,name='fv3jedi_increment_self_schur_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self%self_schur(rhs)

end subroutine fv3jedi_increment_self_schur_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_sub_c(c_key_self,c_key_rhs) &
           bind(c,name='fv3jedi_increment_self_sub_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self%self_sub(rhs)

end subroutine fv3jedi_increment_self_sub_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_mul_c(c_key_self,c_zz) &
           bind(c,name='fv3jedi_increment_self_mul_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(fv3jedi_increment), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self, self)
zz = c_zz

call self%self_mul(zz)

end subroutine fv3jedi_increment_self_mul_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_axpy_inc_c(c_key_self,c_zz,c_key_rhs) &
           bind(c,name='fv3jedi_increment_axpy_inc_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%accumul(zz,rhs)

end subroutine fv3jedi_increment_axpy_inc_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_axpy_state_c(c_key_self,c_zz,c_key_rhs) &
           bind(c,name='fv3jedi_increment_axpy_state_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_state), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_rhs,rhs)
zz = c_zz

call self%accumul(zz,rhs)
!call self%axpy_state(zz,rhs%fields)

end subroutine fv3jedi_increment_axpy_state_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_dot_prod_c(c_key_inc1,c_key_inc2,c_prod) &
           bind(c,name='fv3jedi_increment_dot_prod_f90')

implicit none
integer(c_int), intent(in)    :: c_key_inc1, c_key_inc2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(fv3jedi_increment), pointer :: self, inc2

call fv3jedi_increment_registry%get(c_key_inc1,self)
call fv3jedi_increment_registry%get(c_key_inc2,inc2)

call self%dot_prod(inc2,zz)

c_prod = zz

end subroutine fv3jedi_increment_dot_prod_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2,c_key_geom) &
           bind(c,name='fv3jedi_increment_diff_incr_f90')

implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
integer(c_int), intent(in) :: c_key_geom

type(fv3jedi_increment), pointer :: self
type(fv3jedi_state), pointer :: x1
type(fv3jedi_state), pointer :: x2
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_increment_registry%get(c_key_lhs,self)
call fv3jedi_state_registry%get(c_key_x1,x1)
call fv3jedi_state_registry%get(c_key_x2,x2)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call self%diff_incr(x1%fields,x2%fields,geom)

end subroutine fv3jedi_increment_diff_incr_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_change_resol_c(c_key_inc,c_key_geom,c_key_rhs,c_key_geom_rhs) &
           bind(c,name='fv3jedi_increment_change_resol_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_rhs
integer(c_int), intent(in) :: c_key_geom_rhs

type(fv3jedi_increment), pointer :: self, other
type(fv3jedi_geom),  pointer :: geom, geom_other

call fv3jedi_increment_registry%get(c_key_inc,self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_rhs,other)
call fv3jedi_geom_registry%get(c_key_geom_rhs, geom_other)

call self%change_resol(geom, other, geom_other)

end subroutine fv3jedi_increment_change_resol_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_norm_c(c_key_inc, prms) bind(c,name='fv3jedi_increment_norm_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
real(c_double), intent(inout) :: prms

type(fv3jedi_increment), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_inc,self)

call self%norm(zz)

prms = zz

end subroutine fv3jedi_increment_norm_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_sizes_c(c_key_self,inc_size) bind(c,name='fv3jedi_increment_sizes_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: inc_size
type(fv3jedi_increment), pointer :: self
integer var, i, j, k

call fv3jedi_increment_registry%get(c_key_self, self)

inc_size = 0
do var = 1, self%nf
  inc_size = inc_size + (self%fields(var)%iec-self%fields(var)%isc+1)*&
                        (self%fields(var)%jec-self%fields(var)%jsc+1)*&
                         self%fields(var)%npz
end do

end subroutine fv3jedi_increment_sizes_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_serialize_c(c_key_self,c_vsize,c_vect_inc) &
           bind(c,name='fv3jedi_increment_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self           !< Increment
integer(c_int),intent(in) :: c_vsize              !< Size
real(c_double),intent(out) :: c_vect_inc(c_vsize) !< Vector

type(fv3jedi_increment),pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)
! Call Fortran
call self%serialize(c_vsize,c_vect_inc)

end subroutine fv3jedi_increment_serialize_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_deserialize_c(c_key_self,c_vsize,c_vect_inc,c_index) &
           bind(c,name='fv3jedi_increment_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Increment
integer(c_int),intent(in) :: c_vsize             !< Size
real(c_double),intent(in) :: c_vect_inc(c_vsize) !< Vector
integer(c_int), intent(inout):: c_index          !< Index

type(fv3jedi_increment),pointer :: self

call fv3jedi_increment_registry%get(c_key_self, self)

! Call Fortran
call self%deserialize(c_vsize,c_vect_inc,c_index)


end subroutine fv3jedi_increment_deserialize_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_getpoint_c(c_key_self, c_key_iter, values, values_len) &
           bind(c,name='fv3jedi_increment_getpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self           !< Increment
integer(c_int), intent(in) :: c_key_iter
integer(c_int), intent(in) :: values_len
real(c_double), intent(inout) :: values(values_len)

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom_iter), pointer :: iter

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_iter_registry%get(c_key_iter,iter)

call self%getpoint(iter, values)

end subroutine fv3jedi_increment_getpoint_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_setpoint_c(c_key_self, c_key_iter, values, values_len) &
           bind(c,name='fv3jedi_increment_setpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self           !< Increment
integer(c_int), intent(in)   :: c_key_iter
integer(c_int), intent(in)   :: values_len
real(c_double), intent(in)   :: values(values_len)

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom_iter), pointer :: iter

call fv3jedi_increment_registry%get(c_key_self, self)
call fv3jedi_geom_iter_registry%get(c_key_iter,iter)

call self%setpoint(iter, values)

end subroutine fv3jedi_increment_setpoint_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_getnfieldsncube_c(c_key_self, c_number_fields, c_cube_size) &
           bind(c,name='fv3jedi_increment_getnfieldsncube_f90')

implicit none
integer(c_int), intent(in)  :: c_key_self
integer(c_int), intent(out) :: c_number_fields
integer(c_int), intent(out) :: c_cube_size

type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

c_number_fields = self%nf
c_cube_size = self%npx-1

end subroutine fv3jedi_increment_getnfieldsncube_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_increment_getminmaxrms_c(c_key_self, c_f_num, c_f_name_len, c_f_name, &
                                            c_minmaxrms ) &
           bind(c,name='fv3jedi_increment_getminmaxrms_f90')

implicit none
integer(c_int),               intent(in)    :: c_key_self
integer(c_int),               intent(in)    :: c_f_num
integer(c_int),               intent(in)    :: c_f_name_len
character(len=1,kind=c_char), intent(inout) :: c_f_name(c_f_name_len)
real(c_double),               intent(inout) :: c_minmaxrms(3)

type(fv3jedi_increment), pointer :: self
character(len=field_clen) :: field_name
integer :: n

call fv3jedi_increment_registry%get(c_key_self,self)

call self%minmaxrms(c_f_num, field_name, c_minmaxrms)

do n = 1,c_f_name_len
  c_f_name(n) = field_name(n:n)
enddo

end subroutine fv3jedi_increment_getminmaxrms_c

! --------------------------------------------------------------------------------------------------

end module fv3jedi_increment_interface_mod
