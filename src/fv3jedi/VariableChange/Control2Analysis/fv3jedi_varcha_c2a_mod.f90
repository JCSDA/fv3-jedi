! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_varcha_c2a_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use datetime_mod

! femps
use femps_grid_mod
use femps_operators_mod
use femps_testgrid_mod
use femps_solve_mod
use femps_fv3_mod

! fv3jedi
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state
use fv3jedi_field_mod,     only: copy_subset, field_clen

use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use wind_vt_mod

implicit none
private

public :: fv3jedi_varcha_c2a
public :: create
public :: delete
public :: changevar
public :: changevarinverse

type :: fv3jedi_varcha_c2a
  logical :: skip_femps_init
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
character(len=:), allocatable :: str
character(len=2055) :: path2fv3gridfiles

integer :: n, levs_per_proc
logical :: tmp

! Grid and operators for the femps Poisson solver
! -----------------------------------------------
self%skip_femps_init = .false.
if (conf%has('skip femps initialization')) then
   call conf%get_or_die('skip femps initialization',self%skip_femps_init)
end if

if (.not. self%skip_femps_init) then

  ! Configuration
  call conf%get_or_die('femps_iterations',niter)
  call conf%get_or_die('femps_ngrids',ngrids)
  call conf%get_or_die('femps_path2fv3gridfiles',str); path2fv3gridfiles = str
  if( .not. conf%get('femps_levelprocs',lprocs) ) then
    lprocs = -1
  endif
  if( .not. conf%get('femps_checkconvergence',self%grid%check_convergence) ) then
    self%grid%check_convergence = .false.
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

  if (geom%f_comm%rank() == 0 ) print*, 'Running femps with ', self%lprocs, ' processors.'

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
    call abor1_ftn('fv3jedi_varcha_c2a_mod.create: last level not equal to number of levels.')

  ! Processors doing the work need grid and operators
  if (geom%f_comm%rank() < self%lprocs ) then

    if (geom%f_comm%rank() == 0 ) print*, 'Creating FEMPS grid object'
    call self%grid%setup('cs',ngrids=ngrids,cube=geom%npx-1,niter=niter,&
                         comm = geom%f_comm%communicator(), &
                         rank = geom%f_comm%rank(), &
                         csize = geom%f_comm%size() )

    if (geom%f_comm%rank() == 0 ) print*, 'Creating FEMPS grid hierarchy from files'
    call fv3grid_to_ugrid(self%grid,path2fv3gridfiles)

    ! Build the connectivity and extra geom
    if (geom%f_comm%rank() == 0 ) print*, 'Creating FEMPS cubed-sphere connectivity'
    call self%grid%build_cs(1,1)

    ! Perform all the setup
    if (geom%f_comm%rank() == 0 ) print*, 'Creating FEMPS static operators'
    call preliminary(self%grid,self%oprs)

    ! Partial delete of operators not needed
    if (geom%f_comm%rank() == 0 ) print*, 'FEMPS partial deallocate'
    call self%oprs%pdelete()

  endif

endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self

if (.not. self%skip_femps_init) then
  call self%oprs%delete()
  call self%grid%delete()
endif

if (allocated(self%lev_start)) deallocate(self%lev_start)
if (allocated(self%lev_final)) deallocate(self%lev_final)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine changevar(self,geom,xctl,xana)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xctl
type(fv3jedi_state),      intent(inout) :: xana

integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Stream function/velocity potential
logical :: have_vodi, have_uava
real(kind=kind_real), pointer     ::   psi(:,:,:)     ! Stream function
real(kind=kind_real), pointer     ::   chi(:,:,:)     ! Velocity potentail
real(kind=kind_real), allocatable ::    ud(:,:,:)     ! D-grid u wind
real(kind=kind_real), allocatable ::    vd(:,:,:)     ! D-grid v wind
real(kind=kind_real), allocatable ::  vort(:,:,:)     ! Vorticity
real(kind=kind_real), allocatable ::  divg(:,:,:)     ! Divergence
real(kind=kind_real), allocatable ::    ua(:,:,:)     ! D-grid u wind
real(kind=kind_real), allocatable ::    va(:,:,:)     ! D-grid v wind

! Pressure
logical :: have_delp
real(kind=kind_real), pointer     ::  ps  (:,:,:)     ! Pressure thickness
real(kind=kind_real), allocatable ::  delp(:,:,:)     ! Pressure thickness

! Clouds
logical :: have_cld4
real(kind=kind_real), pointer     :: qi   (:,:,:)     ! Cloud liquid ice
real(kind=kind_real), pointer     :: ql   (:,:,:)     ! Cloud liquid water
real(kind=kind_real), pointer     :: qilsf(:,:,:)     ! Fraction ice large scale
real(kind=kind_real), pointer     :: qicnf(:,:,:)     ! Fraction ice convective
real(kind=kind_real), allocatable :: qils (:,:,:)     ! Cloud liquid ice large scale
real(kind=kind_real), allocatable :: qicn (:,:,:)     ! Cloud liquid ice convective
real(kind=kind_real), allocatable :: qlls (:,:,:)     ! Cloud liquid water large scale
real(kind=kind_real), allocatable :: qlcn (:,:,:)     ! Cloud liquid water convective

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xctl%fields, xana%fields, fields_to_do)

! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return

! Wind variables
! --------------
have_uava = .false.
have_vodi = .false.
if (xctl%has_field('air_horizontal_streamfunction') .and. &
    xctl%has_field('air_horizontal_velocity_potential')) then
  call xctl%get_field('air_horizontal_streamfunction', psi)
  call xctl%get_field('air_horizontal_velocity_potential', chi)
  allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
  allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
  call psichi_to_udvd(geom, psi, chi, ud, vd)
  if (.not. self%skip_femps_init) then
    allocate(vort(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
    allocate(divg(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
    call psichi_to_vortdivg(geom, self%grid, self%oprs, psi, chi, vort, divg, self%lprocs, &
                            self%lev_start, self%lev_final)
    have_vodi = .true.
  endif
  allocate(ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call d_to_a(geom, ud, vd, ua, va)
  have_uava = .true.
endif

! Pressure
! --------
have_delp = .false.
if (xctl%has_field('air_pressure_at_surface')) then
  call xctl%get_field('air_pressure_at_surface', ps)
  allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call ps_to_delp(geom, ps, delp)
  have_delp = .true.
endif

! Clouds
! ------
have_cld4 = .false.
if ( xctl%has_field('cloud_liquid_ice') .and. &
     xctl%has_field('cloud_liquid_water') .and. &
     xctl%has_field('fraction_of_large_scale_cloud_that_is_ice') .and. &
     xctl%has_field('fraction_of_convective_cloud_that_is_ice')) then
  call xctl%get_field('cloud_liquid_ice', qi)
  call xctl%get_field('cloud_liquid_water', ql)
  call xctl%get_field('fraction_of_large_scale_cloud_that_is_ice', qilsf)
  call xctl%get_field('fraction_of_convective_cloud_that_is_ice', qicnf)
  allocate(qils(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(qicn(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(qlls(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(qlcn(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call q2_to_q4(geom, qi, ql, qilsf, qicnf, qils, qicn, qlls, qlcn)
  have_cld4 = .true.
endif

! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call xana%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ('eastward_wind')

    if (.not. have_uava) call field_fail(fields_to_do(f))
    field_ptr = ua

  case ('northward_wind')

    if (.not. have_uava) call field_fail(fields_to_do(f))
    field_ptr = va

  case ('air_upward_absolute_vorticity')

    if (.not. have_vodi) call field_fail(fields_to_do(f))
    field_ptr = vort

  case ('air_horizontal_divergence')

    if (.not. have_vodi) call field_fail(fields_to_do(f))
    field_ptr = divg

  case ('mass_fraction_of_large_scale_cloud_ice_water')

    if (.not. have_cld4) call field_fail(fields_to_do(f))
    field_ptr = qils

  case ('mass_fraction_of_convective_cloud_ice_water')

    if (.not. have_cld4) call field_fail(fields_to_do(f))
    field_ptr = qicn

  case ('mass_fraction_of_large_scale_cloud_liquid_water')

    if (.not. have_cld4) call field_fail(fields_to_do(f))
    field_ptr = qlls

  case ('mass_fraction_of_convective_cloud_liquid_water')

    if (.not. have_cld4) call field_fail(fields_to_do(f))
    field_ptr = qlcn

  case ('air_pressure_thickness')

    if (.not. have_delp) call field_fail(fields_to_do(f))
    field_ptr = delp

  end select

enddo

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xana,xctl)

implicit none
type(fv3jedi_varcha_c2a), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xctl

integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Stream function/velocity potential
logical :: have_pcvd
real(kind=kind_real), pointer     ::    ua(:,:,:)     ! D-grid u wind
real(kind=kind_real), pointer     ::    va(:,:,:)     ! D-grid v wind
real(kind=kind_real), allocatable ::    ud(:,:,:)     ! D-grid u wind
real(kind=kind_real), allocatable ::    vd(:,:,:)     ! D-grid v wind
real(kind=kind_real), allocatable ::   psi(:,:,:)     ! Stream function
real(kind=kind_real), allocatable ::   chi(:,:,:)     ! Velocity potentail
real(kind=kind_real), allocatable ::  vort(:,:,:)     ! Vorticity
real(kind=kind_real), allocatable ::  divg(:,:,:)     ! Divergence

! Pressure
logical :: have_ps
real(kind=kind_real), allocatable ::  ps  (:,:,:)     ! Pressure thickness
real(kind=kind_real), pointer     ::  delp(:,:,:)     ! Pressure thickness

! Humidity
logical :: have_rhum
real(kind=kind_real), pointer     ::     t(:,:,:)     ! Temperature
real(kind=kind_real), pointer     ::     q(:,:,:)     ! Specific humidity
real(kind=kind_real), allocatable ::  qsat(:,:,:)     ! Saturation specific humidity
real(kind=kind_real), allocatable ::    rh(:,:,:)     ! 'relative_humidity' humidity

! Clouds
logical :: have_qiql, have_cfrc
real(kind=kind_real), allocatable :: qi   (:,:,:)     ! Cloud liquid ice
real(kind=kind_real), allocatable :: ql   (:,:,:)     ! Cloud liquid water
real(kind=kind_real), pointer     :: qils (:,:,:)     ! Cloud liquid ice large scale
real(kind=kind_real), pointer     :: qicn (:,:,:)     ! Cloud liquid ice convective
real(kind=kind_real), pointer     :: qlls (:,:,:)     ! Cloud liquid water large scale
real(kind=kind_real), pointer     :: qlcn (:,:,:)     ! Cloud liquid water convective
real(kind=kind_real), allocatable :: qilsf(:,:,:)     ! Fraction ice large scale
real(kind=kind_real), allocatable :: qicnf(:,:,:)     ! Fraction ice convective


! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xana%fields, xctl%fields, fields_to_do)

! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return

! Wind variables
! --------------
have_pcvd = .false.
if (xana%has_field('eastward_wind') .and. xana%has_field('northward_wind')) then
  call xana%get_field('eastward_wind', ua)
  call xana%get_field('northward_wind', va)
  allocate(ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz))
  allocate(vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz))
  call a_to_d(geom, ua, va, ud, vd)
  allocate( psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate( chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(vort(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(divg(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call udvd_to_vortdivg(geom, ud, vd, vort, divg)
  call vortdivg_to_psichi(geom, self%grid, self%oprs, vort, divg, psi, chi, &
                          self%lprocs, self%lev_start, self%lev_final)
  have_pcvd = .true.
endif

! Surface pressure
! ----------------
have_ps = .false.
if (xana%has_field('air_pressure_thickness')) then
  call xana%get_field('air_pressure_thickness', delp)
  allocate(ps(geom%isc:geom%iec,geom%jsc:geom%jec,1))
  ps(:,:,1) = sum(delp, dim=3) + geom%ptop
  have_ps = .true.
endif

! Humidity
! --------
have_rhum = .false.
if (xana%has_field('water_vapor_mixing_ratio_wrt_moist_air') &
    .and. xana%has_field('air_temperature') .and. xana%has_field('air_pressure_thickness')) then
  call xana%get_field('air_pressure_thickness', delp)
  call xana%get_field('air_temperature', t)
  allocate(qsat(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(  rh(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call xana%get_field('water_vapor_mixing_ratio_wrt_moist_air', q)
  call get_qsat(geom, delp, t, q, qsat)
  call q_to_rh(geom, qsat, q, rh)
  have_rhum = .true.
endif

! Clouds
! ------
have_qiql = .false.
have_cfrc = .false.
if (xana%has_field('cloud_liquid_ice') .and. xana%has_field('cloud_liquid_water')) then
  call xana%get_field('cloud_liquid_ice', qi)
  call xana%get_field('cloud_liquid_water', ql)
  have_qiql = .true.
elseif (xana%has_field('mass_fraction_of_large_scale_cloud_ice_water') .and. &
        xana%has_field('mass_fraction_of_convective_cloud_ice_water') .and. &
        xana%has_field('mass_fraction_of_large_scale_cloud_liquid_water') .and. &
        xana%has_field('mass_fraction_of_convective_cloud_liquid_water')) then
  call xana%get_field('mass_fraction_of_large_scale_cloud_ice_water', qils)
  call xana%get_field('mass_fraction_of_convective_cloud_ice_water', qicn)
  call xana%get_field('mass_fraction_of_large_scale_cloud_liquid_water', qlls)
  call xana%get_field('mass_fraction_of_convective_cloud_liquid_water', qlcn)
  allocate(qi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(ql(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(qilsf(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  allocate(qicnf(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
  call q4_to_q2(geom, qils, qicn, qlls, qlcn, qi, ql, qilsf, qicnf)
  have_cfrc = .true.
  have_qiql = .true.
endif


! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call xctl%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ('air_horizontal_streamfunction')

    if (.not. have_pcvd) call field_fail(fields_to_do(f))
    field_ptr = psi

  case ('air_horizontal_velocity_potential')

    if (.not. have_pcvd) call field_fail(fields_to_do(f))
    field_ptr = chi

  case ('air_upward_absolute_vorticity')

    if (.not. have_pcvd) call field_fail(fields_to_do(f))
    field_ptr = vort

  case ('air_horizontal_divergence')

    if (.not. have_pcvd) call field_fail(fields_to_do(f))
    field_ptr = divg

  case ('relative_humidity')

    if (.not. have_rhum) call field_fail(fields_to_do(f))
    field_ptr = rh

  case ('cloud_liquid_ice')

    if (.not. have_qiql) call field_fail(fields_to_do(f))
    field_ptr = qi

  case ('cloud_liquid_water')

    if (.not. have_qiql) call field_fail(fields_to_do(f))
    field_ptr = ql

  case ('fraction_of_large_scale_cloud_that_is_ice')

    if (.not. have_cfrc) call field_fail(fields_to_do(f))
    field_ptr = qilsf

  case ('fraction_of_convective_cloud_that_is_ice')

    if (.not. have_cfrc) call field_fail(fields_to_do(f))
    field_ptr = qicnf

  case ('air_pressure_at_surface')

    if (.not. have_ps) call field_fail(fields_to_do(f))
    field_ptr = ps

  end select

enddo

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_c2a_mod
