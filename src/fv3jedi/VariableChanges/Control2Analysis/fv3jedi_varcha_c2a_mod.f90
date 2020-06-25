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

use fv3jedi_field_mod, only: copy_subset, has_field, pointer_field_array

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

real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_chi

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_vd
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_vort
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_divg

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xctl%fields,xana%fields)

! psi and chi to A-grid u and v
! -----------------------------
call pointer_field_array(xctl%fields, 'psi', xctl_psi)
call pointer_field_array(xctl%fields, 'chi', xctl_chi)
call pointer_field_array(xana%fields, 'ud' , xana_ud)
call pointer_field_array(xana%fields, 'vd' , xana_vd)
call psichi_to_udvd(geom, xctl_psi, xctl_chi, xana_ud, xana_vd)

! For testing Poisson solver recover vorticity and divergence
! -----------------------------------------------------------
if (has_field(xana%fields, 'vort') .and. has_field(xana%fields, 'divg')) then
  call pointer_field_array(xana%fields, 'vort', xana_vort)
  call pointer_field_array(xana%fields, 'divg', xana_divg)
  call psichi_to_vortdivg(geom, self%grid, self%oprs, xctl_psi, xctl_chi, self%lprocs, self%lev_start,&
                          self%lev_final, xana_vort, xana_divg)
endif

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

real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_vort
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_divg
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_rh
real(kind=kind_real), pointer, dimension(:,:,:) :: xctl_ps

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_vd
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_q
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_delp

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xana%fields,xctl%fields)

! D-grid winds to stream func and velocity potential
! --------------------------------------------------
if (has_field(xctl%fields, 'vort') .and. has_field(xctl%fields, 'divg')) then
  call pointer_field_array(xana%fields, 'ud' , xana_ud)
  call pointer_field_array(xana%fields, 'vd' , xana_vd)
  call pointer_field_array(xctl%fields, 'psi', xctl_psi)
  call pointer_field_array(xctl%fields, 'chi', xctl_chi)
  call pointer_field_array(xctl%fields, 'vort', xctl_vort)
  call pointer_field_array(xctl%fields, 'divg', xctl_divg)
  call udvd_to_psichi(geom, self%grid, self%oprs, xana_ud, xana_vd, xctl_psi, xctl_chi, &
                      self%lprocs, self%lev_start, self%lev_final, &
                      xctl_vort, xctl_divg)
else
  call pointer_field_array(xana%fields, 'ud' , xana_ud)
  call pointer_field_array(xana%fields, 'vd' , xana_vd)
  call pointer_field_array(xctl%fields, 'psi', xctl_psi)
  call pointer_field_array(xctl%fields, 'chi', xctl_chi)
  call udvd_to_psichi(geom, self%grid, self%oprs, xana_ud, xana_vd, xctl_psi, xctl_chi, &
                      self%lprocs, self%lev_start, self%lev_final)
endif

! Temperature to virtual temperature
! ----------------------------------
call pointer_field_array(xana%fields, 't' , xana_t)
call pointer_field_array(xana%fields, 'sphum' , xana_q)
call pointer_field_array(xctl%fields, 'tv', xctl_tv)
call t_to_tv(geom, xana_t, xana_q, xctl_tv)

! Specific humidity to relative humidity
! --------------------------------------
call pointer_field_array(xana%fields, 'delp', xana_delp)
call pointer_field_array(xana%fields, 't'   , xana_t)
call pointer_field_array(xana%fields, 'sphum'   , xana_q)
call pointer_field_array(xctl%fields, 'rh'  , xctl_rh)
call get_qsat(geom, xana_delp, xana_t, xana_q, qsat)
call q_to_rh(geom, qsat, xana_q, xctl_rh)

! Surface pressure
! ----------------
call pointer_field_array(xana%fields, 'delp', xana_delp)
call pointer_field_array(xctl%fields, 'ps'  , xctl_ps)
xctl_ps(:,:,1) = sum(xana_delp,3)

! Copy calendar infomation
xctl%calendar_type = xana%calendar_type
xctl%date_init = xana%date_init

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_c2a_mod
