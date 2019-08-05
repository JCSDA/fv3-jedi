! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_setup_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration

use mpp_mod,         only: mpp_init, mpp_exit
use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
use fms_io_mod,      only: fms_io_init, fms_io_exit

use mpp_domains_mod, only: mpp_domains_set_stack_size
use fckit_mpi_module,only: fckit_mpi_comm
use fms_mod,         only: fms_init

implicit none
private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3jedi_setup(c_conf) bind(c,name='fv3jedi_setup_f')

implicit none

type(c_ptr), intent(in) :: c_conf
integer :: stackmax = 4000000
type(fckit_mpi_comm) :: f_comm
type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm()

call mpp_init(localcomm=f_comm%communicator())
call mpp_domains_init
call fms_io_init
call fms_init()

if (f_conf%has("stackmax")) then
  call f_conf%get_or_die("stackmax",stackmax)
endif

call mpp_domains_set_stack_size(stackmax)

end subroutine fv3jedi_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_finalize() bind(c,name='fv3jedi_finalize_f')

implicit none

call fms_io_exit
call mpp_domains_exit
!call mpp_exit

end subroutine fv3jedi_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_setup_mod
