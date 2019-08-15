! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_covariance_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom

implicit none

!> Fortran derived type to hold configuration data for the background/model covariance
type :: fv3jedi_covar
  integer :: nothing_yet
end type fv3jedi_covar

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

!> Setup for the model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine fv3jedi_covar_setup(self, geom, c_conf)

implicit none
type(fv3jedi_covar), intent(inout) :: self    !< Covariance structure
type(c_ptr), intent(in)            :: c_conf  !< Configuration
type(fv3jedi_geom), intent(in)     :: geom    !< Geometry

type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)

end subroutine fv3jedi_covar_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_covar_delete(self)

implicit none
type(fv3jedi_covar), intent(inout) :: self  !< Covariance structure

end subroutine fv3jedi_covar_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)), where C is 3d covariance matrix

subroutine fv3jedi_covar_sqrt_inv_mult(self, xctl, xincr)

implicit none
type(fv3jedi_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(fv3jedi_increment), intent(in)    :: xincr

end subroutine fv3jedi_covar_sqrt_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)) - Adjoint

subroutine fv3jedi_covar_sqrt_inv_mult_ad(self, xctl, xincr)

implicit none
type(fv3jedi_covar), intent(in)    :: self
type(fv3jedi_increment), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine fv3jedi_covar_sqrt_inv_mult_ad

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C), where C is a 3d covariance matrix

subroutine fv3jedi_covar_sqrt_mult(self, xincr, xctl)

implicit none
type(fv3jedi_covar), intent(in)    :: self
type(fv3jedi_increment), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine fv3jedi_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine fv3jedi_covar_sqrt_mult_ad(self, xincr, xctl)

implicit none
type(fv3jedi_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(fv3jedi_increment), intent(in)    :: xincr

end subroutine fv3jedi_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module fv3jedi_covariance_mod