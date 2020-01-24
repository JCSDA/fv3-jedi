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
  integer :: lprocs
  integer, allocatable :: lev_start(:), lev_final(:)
end type fv3jedi_varcha_c2a

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
type(fv3jedi_varcha_c2a),  intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf

integer :: ngrids, niter, lprocs, lstart
logical :: check_convergence
character(len=:), allocatable :: str
character(len=2055) :: path2fv3gridfiles

integer :: n, levs_per_proc

! Grid and operators for the femps Poisson solver
! -----------------------------------------------

! Configuration
call conf%get_or_die("femps_iterations",niter)
call conf%get_or_die("femps_ngrids",ngrids)
call conf%get_or_die("femps_path2fv3gridfiles",str); path2fv3gridfiles = str
if( .not. conf%get('femps_levelprocs',lprocs) ) then
  lprocs = -1
endif
if( .not. conf%get('femps_checkconvergence',check_convergence) ) then
  check_convergence = .false.
endif

! Processors that will do the work
! --------------------------------
lprocs = min(lprocs,geom%f_comm%size())
lprocs = min(lprocs,geom%npz)

if (lprocs == -1) then
  self%lprocs = min(geom%npz,geom%f_comm%size())
else
  self%lprocs = lprocs
endif

if (geom%f_comm%rank() == 0 ) print*, "Running femps with ", self%lprocs ," processors."

allocate(self%lev_start(self%lprocs))
allocate(self%lev_final(self%lprocs))

if (self%lprocs == geom%npz) then
  do n = 1,self%lprocs
    self%lev_start(n) = n
    self%lev_final(n) = n
  enddo
else
  levs_per_proc = floor(real(geom%npz,kind_real)/real(self%lprocs,kind_real))
  lstart = 0
  do n = 1,self%lprocs
    self%lev_start(n) = lstart+1
    self%lev_final(n) = self%lev_start(n) + levs_per_proc - 1
    if (n .le. mod(geom%npz, self%lprocs)) self%lev_final(n) = self%lev_final(n) + 1
    lstart = self%lev_final(n)
  enddo
endif

if (self%lev_final(self%lprocs) .ne. geom%npz) &
  call abor1_ftn("fv3jedi_varcha_c2a_mod.create: last level not equal to number of levels.")

! Processors doing the work need grid and operators
if (geom%f_comm%rank() < self%lprocs ) then

  call self%grid%setup('cs',ngrids=ngrids,cube=geom%npx-1,niter=niter,&
                       comm = geom%f_comm%communicator(), &
                       rank = geom%f_comm%rank(), &
                       csize = geom%f_comm%size(), &
                       check_convergence = check_convergence )

  call fv3grid_to_ugrid(self%grid,path2fv3gridfiles)

  ! Build the connectivity and extra geom
  call self%grid%build_cs(1,1)

  ! Perform all the setup
  call preliminary(self%grid,self%oprs)

  ! Partial delete of operators not needed
  call self%oprs%pdelete()

endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self

call self%oprs%delete()
call self%grid%delete()

deallocate(self%lev_start,self%lev_final)

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
  call psichi_to_vortdivg(geom,self%grid,self%oprs,xctl%psi,xctl%chi,self%lprocs,self%lev_start,&
                          self%lev_final,xana%vort,xana%divg)
endif

! Copied variables
! ----------------
xana%t    = xctl%t
xana%delp = xctl%delp
xana%q    = xctl%q
xana%qi   = xctl%qi
xana%ql   = xctl%ql
xana%o3   = xctl%o3

! Copy calendar infomation
xana%calendar_type = xctl%calendar_type
xana%date_init = xctl%date_init

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xana,xctl)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xctl

real(kind=kind_real) :: qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

! D-grid winds to stream func and velocity potential
! --------------------------------------------------
if (associated(xctl%vort) .and. associated(xctl%divg)) then
  call udvd_to_psichi(geom,self%grid,self%oprs,xana%ud,xana%vd,xctl%psi,xctl%chi,&
                      self%lprocs, self%lev_start, self%lev_final, &
                      xctl%vort,xctl%divg)
else
  call udvd_to_psichi(geom,self%grid,self%oprs,xana%ud,xana%vd,xctl%psi,xctl%chi,&
                      self%lprocs, self%lev_start, self%lev_final)
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

! Copy calendar infomation
xctl%calendar_type = xana%calendar_type
xctl%date_init = xana%date_init

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_c2a_mod
