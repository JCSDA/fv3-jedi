! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_setup_mod

use iso_c_binding

use string_f_c_mod
use fckit_mpi_module, only: fckit_mpi_comm

use mpp_mod, only: mpp_init, mpp_clock_id, mpp_clock_begin
use fms_mod, only: fms_init!, fms_end

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3jedi_setup(c_conf, c_comm) bind(c,name='fv3jedi_setup_f')

implicit none

type(c_ptr),        intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_comm

type(fckit_mpi_comm) :: f_comm
integer :: initClock

! Fortran configuration
! ---------------------
f_comm = fckit_mpi_comm(c_comm)

call fms_init(localcomm=f_comm%communicator())
call mpp_init()

initClock = mpp_clock_id( 'Initialization' )
call mpp_clock_begin (initClock) !nesting problem

call fms_init

end subroutine fv3jedi_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_finalize() bind(c,name='fv3jedi_finalize_f')

implicit none

!call fms_end

end subroutine fv3jedi_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_setup_mod
