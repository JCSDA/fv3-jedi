! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_constants_mod

! Bindings to c++ code
use iso_c_binding

! fv3jedi uses
use fv3jedi_kinds_mod, only: kind_real

implicit none
private
public constant

! --------------------------------------------------------------------------------------------------

! Interface to get the constant
interface
  subroutine get_constant_c(constant_name, constant_value) bind(c, name='getConstantF')
    use iso_c_binding
    integer, parameter :: clen = 2048
    character(len=1, kind=c_char), intent(in) :: constant_name(clen)
    real(kind=c_double) :: constant_value
  end subroutine get_constant_c
end interface

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

function constant(constant_name) result(constant_value)

    ! Arguments
    character(len=*), intent(in) :: constant_name
    real(kind=kind_real) :: constant_value

    ! Locals
    real(kind=c_double) :: constant_value_c
    real(kind=kind_real) :: constant_value_tmp

    ! Get constant from cpp code
    call get_constant_c([constant_name, C_NULL_CHAR], constant_value_c)
    constant_value = real(constant_value_c, kind=kind_real)

end function

! --------------------------------------------------------------------------------------------------

end module fv3jedi_constants_mod
