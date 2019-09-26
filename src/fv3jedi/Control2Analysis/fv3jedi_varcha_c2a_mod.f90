! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_varcha_c2a_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod

use fckit_log_module, only : fckit_log
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm

use fv3jedi_kinds_mod,   only: kind_real
use fv3jedi_geom_mod,    only: fv3jedi_geom
use fv3jedi_state_mod,   only: fv3jedi_state
use fv3jedi_io_gfs_mod,  only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos

use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use wind_vt_mod

use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod
use femps_fv3_mod

implicit none
private

public :: fv3jedi_varcha_c2a
public :: create
public :: delete
public :: changevar
public :: changevarinverse

type :: fv3jedi_varcha_c2a
  type(fempsgrid) :: grid
  type(fempsoprs) :: oprs
  type(fckit_mpi_comm) :: f_comm
end type fv3jedi_varcha_c2a

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
type(fv3jedi_varcha_c2a),  intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf

integer :: ngrids, niter
character(len=:), allocatable :: str
character(len=2055) :: path2fv3gridfiles

! Grid and operators for the femps Poisson solver
! -----------------------------------------------


! Pointer to fv3jedi communicator
self%f_comm = fckit_mpi_comm()

! Solver only works on one PE
if (self%f_comm%rank() == 0) then

  ! Configuration
  call conf%get_or_die("femps_iterations",niter)
  call conf%get_or_die("femps_ngrids",ngrids)
  call conf%get_or_die("femps_path2fv3gridfiles",str); path2fv3gridfiles = str

  call self%grid%setup('cs',ngrids=ngrids,cube=geom%npx-1,niter=niter)
  call fv3grid_to_ugrid(self%grid,path2fv3gridfiles)

  ! Build the connectivity and extra geom
  call self%grid%build_cs(1,1)

  ! Perform all the setup
  call preliminary(self%grid,self%oprs)

endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self

if (self%f_comm%rank() == 0) then
  call self%oprs%delete()
  call self%grid%delete()
endif

end subroutine delete

! ------------------------------------------------------------------------------

subroutine changevar(self,geom,xctl,xana)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xctl
type(fv3jedi_state),      intent(inout) :: xana

! psi and chi to A-grid u and v
! -----------------------------
call psichi_to_udvd(geom,xctl%psi,xctl%chi,xana%ud,xana%vd)

! For testing Poisson solver recover vorticity and divergence
! -----------------------------------------------------------
if (associated(xana%vort) .and. associated(xana%divg)) then
  call psichi_to_vortdivg(geom,self%grid,self%oprs,xctl%psi,xctl%chi,xana%vort,xana%divg)
endif

! Copied variables
! ----------------
xana%t    = xctl%t
xana%delp = xctl%delp
xana%q    = xctl%q
xana%qi   = xctl%qi
xana%ql   = xctl%ql
xana%o3   = xctl%o3

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xana,xctl)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xctl

real(kind=kind_real) :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

! A-grid winds to stream func and velocity potential
! --------------------------------------------------
if (associated(xctl%vort) .and. associated(xctl%divg)) then
  call udvd_to_psichi(geom,self%grid,self%oprs,xana%ud,xana%vd,xctl%psi,xctl%chi,xctl%vort,xctl%divg)
else
  call udvd_to_psichi(geom,self%grid,self%oprs,xana%ud,xana%vd,xctl%psi,xctl%chi)
endif

! Temperature to virtual temperature
! ----------------------------------
call t_to_tv(geom,xana%t,xana%q,xctl%tv)

! Specific humidity to relative humidity
! --------------------------------------
call get_qsat(geom,xana%delp,xana%t,xana%q,qsat)
call q_to_rh(geom,qsat,xana%q,xctl%rh)

! Surface pressure
! ----------------
xctl%ps(:,:,1) = sum(xana%delp,3)

! Copied variables
! ----------------
xctl%t    = xana%t
xctl%delp = xana%delp
xctl%q    = xana%q
xctl%qi   = xana%qi
xctl%ql   = xana%ql
xctl%o3   = xana%o3

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_c2a_mod
