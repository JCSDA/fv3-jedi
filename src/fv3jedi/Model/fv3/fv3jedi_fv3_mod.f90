! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3_mod

! oops uses
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

!! fms uses
!use constants_mod,      only: cp_air, rdgas, grav, rvgas, kappa, omega, pi
!use field_manager_mod,  only: MODEL_ATMOS
!use fms_mod,            only: set_domain, nullify_domain
!use mpp_domains_mod,    only: mpp_update_domains, mpp_get_boundary, DGRID_NE
!use time_manager_mod,   only: time_type, time_type_to_real, set_time, get_time, set_date, &
!                              operator(-), JULIAN, set_calendar_type
!use tracer_manager_mod, only: get_tracer_index, get_number_tracers, NO_TRACER
!
!! fv3 uses
!use fv_arrays_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
!use fv_control_mod,     only: fv_init, fv_end, ngrids
!use fv_dynamics_mod,    only: fv_dynamics
!use fv_grid_utils_mod,  only: grid_utils_end
!use fv_nesting_mod,     only: twoway_nesting
!use fv_regional_mod,    only: a_step, p_step, current_time_in_seconds, read_new_bc_data

! fv3-jedi uses
use fv3jedi_kinds_mod,  only: kind_real
use fv3jedi_geom_mod,   only: fv3jedi_geom
use fv3jedi_state_mod,  only: fv3jedi_state
use fv3jedi_field_mod,  only: pointer_field_array

implicit none
private
public :: fv3_model

! --------------------------------------------------------------------------------------------------

type :: fv3_model
  !type(fv_atmos_type), allocatable :: Atm(:)      !<FV3 structure
  !integer :: p_split, mytile
  !integer :: npx, npy, npz, ncnst, nq             !<Grid, domains
  !integer :: isc, iec, jsc, jec                   !<Grid, compute region
  !integer :: isd, ied, jsd, jed                   !<Grid, with halo
  !type(time_type) :: Time, Time_step_atmos
  !real(kind=kind_real) :: dt_atmos
  !logical, allocatable :: grids_on_this_pe(:)
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: initialize
    procedure, public :: step
    procedure, public :: finalize
end type fv3_model

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
class(fv3_model), target,       intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

!!Locals
!integer :: i, j, n, dtsecs
!character(len=20) :: ststep
!type(duration) :: dtstep
!real(kind=kind_real) :: dt
!character(len=:), allocatable :: str
!real(kind_real) :: f_coriolis_angle
!
!! Model time step
!! ---------------
!call conf%get_or_die("tstep",str)
!ststep = str
!deallocate(str)
!dtstep = trim(ststep)
!dtsecs = duration_seconds(dtstep)
!self%dt_atmos = real(duration_seconds(dtstep),kind_real)
!self%Time_step_atmos = set_time(dtsecs)
!
!
!! Call initialize
!! ---------------
!self%p_split = 1
!call fv_init( self%Atm, self%dt_atmos, self%grids_on_this_pe, self%p_split )
!
!do n = 1, ngrids
!  if (self%grids_on_this_pe(n)) self%mytile = n
!enddo
!
!a_step = 0
!if(self%Atm(self%mytile)%flagstruct%warm_start) then
!  a_step = nint(current_time_in_seconds/self%dt_atmos)
!endif
!
!! Convenience
!! -----------
!self%npx   = self%Atm(self%mytile)%npx
!self%npy   = self%Atm(self%mytile)%npy
!self%npz   = self%Atm(self%mytile)%npz
!self%ncnst = self%Atm(self%mytile)%ncnst
!
!self%isc = self%Atm(self%mytile)%bd%isc
!self%iec = self%Atm(self%mytile)%bd%iec
!self%jsc = self%Atm(self%mytile)%bd%jsc
!self%jec = self%Atm(self%mytile)%bd%jec
!
!self%isd = self%isc - self%Atm(self%mytile)%bd%ng
!self%ied = self%iec + self%Atm(self%mytile)%bd%ng
!self%jsd = self%jsc - self%Atm(self%mytile)%bd%ng
!self%jed = self%jec + self%Atm(self%mytile)%bd%ng
!
!self%nq = self%Atm(self%mytile)%ncnst-self%Atm(self%mytile)%flagstruct%pnats
!
!! Set ptop, ak, bk
!! ----------------
!self%Atm(self%mytile)%ak = geom%ak
!self%Atm(self%mytile)%bk = geom%bk
!self%Atm(self%mytile)%ptop = geom%ptop
!
!! Assert grid similarity with Geometry TODO does this work for regional?
!if (geom%isc .ne. self%isc .or. geom%iec .ne. self%iec .or. &
!    geom%jsc .ne. self%jsc .or. geom%jec .ne. self%jec) then
!  call abor1_ftn("fv3jedi_fv3_mod: dimension mismatch between Geometry and the model")
!endif
!
!call set_calendar_type(JULIAN)
!
!
!! Things normally done by fv_restart
!! ----------------------------------
!f_coriolis_angle = 0.0_kind_real
!
!!fC and f0
!if (self%Atm(self%mytile)%flagstruct%grid_type == 4) then
!   self%Atm(self%mytile)%gridstruct%fC(:,:) = 2.0_kind_real*omega*sin(self%Atm(self%mytile)%flagstruct%deglat/180.0_kind_real*pi)
!   self%Atm(self%mytile)%gridstruct%f0(:,:) = 2.0_kind_real*omega*sin(self%Atm(self%mytile)%flagstruct%deglat/180.0_kind_real*pi)
!else
!   if (f_coriolis_angle == -999.0_kind_real) then
!      self%Atm(self%mytile)%gridstruct%fC(:,:) = 0.0_kind_real
!      self%Atm(self%mytile)%gridstruct%f0(:,:) = 0.0_kind_real
!   else
!      do j=self%Atm(self%mytile)%bd%jsd,self%Atm(self%mytile)%bd%jed+1
!         do i=self%Atm(self%mytile)%bd%isd,self%Atm(self%mytile)%bd%ied+1
!            self%Atm(self%mytile)%gridstruct%fC(i,j) = 2.0_kind_real*omega*( -COS(self%Atm(self%mytile)%gridstruct%grid(i,j,1))*&
!                                           COS(self%Atm(self%mytile)%gridstruct%grid(i,j,2))*SIN(f_coriolis_angle) + &
!                                           SIN(self%Atm(self%mytile)%gridstruct%grid(i,j,2))*COS(f_coriolis_angle) )
!         enddo
!      enddo
!      do j=self%Atm(self%mytile)%bd%jsd,self%Atm(self%mytile)%bd%jed
!         do i=self%Atm(self%mytile)%bd%isd,self%Atm(self%mytile)%bd%ied
!            self%Atm(self%mytile)%gridstruct%f0(i,j) = 2.0_kind_real*omega*( -COS(self%Atm(self%mytile)%gridstruct%agrid(i,j,1))*&
!                                           COS(self%Atm(self%mytile)%gridstruct%agrid(i,j,2))*SIN(f_coriolis_angle) + &
!                                           SIN(self%Atm(self%mytile)%gridstruct%agrid(i,j,2))*COS(f_coriolis_angle) )
!         enddo
!      enddo
!   endif
!endif
!
!!Pointer to self when not nested
!if (.not. self%Atm(self%mytile)%gridstruct%nested) self%Atm(self%mytile)%parent_grid => self%Atm(self%mytile)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3_model), intent(inout) :: self

integer :: n

!call grid_utils_end()

!do n = 1, size(self%Atm(:))
!   call deallocate_fv_atmos_type(self%Atm(n))
!end do

!deallocate(self%Atm)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state, vdate)

implicit none
class(fv3_model),    intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state
type(datetime),      intent(inout) :: vdate

integer :: field, tindex
character(len=20) :: vdatestrz
character(len=19) :: vdatestr

!! Check that all the tracers in fv3-jedi are in field_table
!! ---------------------------------------------------------
!do field = 1,state%nf
!  if (state%fields(field)%tracer) then
!    ! Get tracer index from the fv3 model side
!    tindex = get_tracer_index(MODEL_ATMOS, trim(state%fields(field)%fv3jedi_name))
!    ! Abort if this tracer was not put into field_table
!    if (tindex == NO_TRACER) call abor1_ftn("fv3_model: tracer "//&
!                                      trim(state%fields(field)%fv3jedi_name)//" not in field_table")
!  endif
!enddo
!
!! Set the time
!! ------------
!call datetime_to_string(vdate, vdatestrz)
!vdatestr(1:19) = vdatestrz(1:19)
!vdatestr(11:11) = ' '
!self%Time = set_date(vdatestr)
!
!! Internal init time
!! ------------------
!self%Atm(self%mytile)%Time_init = self%Time
!
!! Regional
!! --------
!current_time_in_seconds = time_type_to_real( self%Time - self%Atm(self%mytile)%Time_init )
!a_step = 0
!if(self%Atm(self%mytile)%flagstruct%warm_start) then
!  !TODO TimeInit is not correct. Need to pass the proper time from coupler.res
!  a_step = nint(current_time_in_seconds/self%dt_atmos)
!endif
!
!! Set domain
!call set_domain ( self%Atm(self%mytile)%domain )

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, vdate)

implicit none
class(fv3_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(inout) :: vdate

integer :: n, seconds, days, n_split_loc, psc
real(8) :: zvir

!! Mytile
!! ------
!n = self%mytile
!
!!Copy from traj to the fv3 structure
!!-----------------------------------
!call state_to_fv3(state, self%Atm(n))
!
!! nsplits
!! -------
!call get_time (self%Time, seconds,  days)
!if (seconds < nint(3600*self%Atm(n)%flagstruct%fhouri) .and. self%Atm(n)%flagstruct%fac_n_spl > 1.0) then
!  n_split_loc = nint(self%Atm(n)%flagstruct%n_split * self%Atm(n)%flagstruct%fac_n_spl)
!else
!  n_split_loc = self%Atm(n)%flagstruct%n_split
!endif
!
!! Constants
!! ---------
!zvir = rvgas/rdgas - 1.
!
!! Regional
!! --------
!a_step = a_step + 1
!if(self%Atm(n)%flagstruct%regional)then
!  call read_new_bc_data(self%Atm(n), self%Time, self%Time_step_atmos, self%p_split, &
!                        self%isd, self%ied, self%jsd, self%jed )
!endif
!
!!Propagate FV3
!!-------------
!do psc=1,abs(self%p_split)
!
!  print*, 'u', maxval(self%Atm(n)%u), minval(self%Atm(n)%u)
!  print*, 'v', maxval(self%Atm(n)%v), minval(self%Atm(n)%v)
!  print*, 't', maxval(self%Atm(n)%pt), minval(self%Atm(n)%pt)
!  print*, 'q', maxval(self%Atm(n)%q(:,:,:,1)), minval(self%Atm(n)%q(:,:,:,1))
!  print*, 'p', maxval(self%Atm(n)%delp), minval(self%Atm(n)%delp)
!
!  p_step = psc
!  call fv_dynamics( self%npx, self%npy, self%npz, self%nq, self%Atm(n)%ng,                  &
!                    self%dt_atmos/real(abs(self%p_split)),                                  &
!                    self%Atm(n)%flagstruct%consv_te, self%Atm(n)%flagstruct%fill,           &
!                    self%Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,              &
!                    self%Atm(n)%ptop, self%Atm(n)%ks, self%nq,                              &
!                    n_split_loc, self%Atm(n)%flagstruct%q_split,                            &
!                    self%Atm(n)%u,    self%Atm(n)%v,     self%Atm(n)%w,  self%Atm(n)%delz,  &
!                    self%Atm(n)%flagstruct%hydrostatic,                                     &
!                    self%Atm(n)%pt  , self%Atm(n)%delp,  self%Atm(n)%q,  self%Atm(n)%ps,    &
!                    self%Atm(n)%pe,   self%Atm(n)%pk,    self%Atm(n)%peln,                  &
!                    self%Atm(n)%pkz,  self%Atm(n)%phis,  self%Atm(n)%q_con,                 &
!                    self%Atm(n)%omga, self%Atm(n)%ua,    self%Atm(n)%va, self%Atm(n)%uc,    &
!                    self%Atm(n)%vc,   self%Atm(n)%ak,    self%Atm(n)%bk, self%Atm(n)%mfx,   &
!                    self%Atm(n)%mfy , self%Atm(n)%cx,    self%Atm(n)%cy, self%Atm(n)%ze0,   &
!                    self%Atm(n)%flagstruct%hybrid_z,                                        &
!                    self%Atm(n)%gridstruct,  self%Atm(n)%flagstruct,                        &
!                    self%Atm(n)%neststruct,  self%Atm(n)%idiag, self%Atm(n)%bd,             &
!                    self%Atm(n)%parent_grid, self%Atm(n)%domain, self%Atm(n)%diss_est       )
!
!  if (ngrids > 1 .and. (psc < self%p_split .or. self%p_split < 0)) then
!    call twoway_nesting(self%Atm, ngrids, self%grids_on_this_pe, zvir)
!  endif
!
!enddo !p_split
!
!
!!Copy from fv3 back to traj structure
!!------------------------------------
!call fv3_to_state(self%Atm(self%mytile), state)

end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine finalize(self, state)

implicit none
class(fv3_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state

! MPP nulify
! ----------
!call nullify_domain ( )

end subroutine finalize

! --------------------------------------------------------------------------------------------------

!subroutine state_to_fv3( state, Atm )
!
!implicit none
!type(fv3jedi_state), intent(in)    :: state
!type(FV_Atmos_type), intent(inout) :: Atm
!
!integer :: i, j, k, field, tindex
!integer :: isc, iec, jsc, jec, npz
!integer :: isd, ied, jsd, jed
!logical :: update_winds
!real(kind=kind_real), pointer, dimension(:,:,:) :: field_ptr
!
!real(kind_real), allocatable, dimension(:,:) :: ebuffery    !<Halo holder
!real(kind_real), allocatable, dimension(:,:) :: nbufferx    !<Halo holder
!real(kind_real), allocatable, dimension(:,:) :: wbuffery    !<Halo holder
!real(kind_real), allocatable, dimension(:,:) :: sbufferx    !<Halo holder
!
!! Convenience
!! -----------
!isc = Atm%bd%isc
!iec = Atm%bd%iec
!jsc = Atm%bd%jsc
!jec = Atm%bd%jec
!isd = Atm%bd%isd
!ied = Atm%bd%ied
!jsd = Atm%bd%jsd
!jed = Atm%bd%jed
!npz = Atm%npz
!
!update_winds = .false.
!
!! Loop through fields on fv3-jedi side and copy
!!----------------------------------------------
!do field = 1,state%nf
!
!  ! Get pointer to fv3-jedi field
!  call pointer_field_array(state%fields, state%fields(field)%fv3jedi_name, field_ptr)
!
!  ! Tracer names match between systems so easy copy
!  if (state%fields(field)%tracer) then
!
!    ! Get index on fv3 model side and fill
!    tindex = get_tracer_index(MODEL_ATMOS, trim(state%fields(field)%fv3jedi_name))
!    Atm%q(isc:iec,jsc:jec,:,tindex) = field_ptr(isc:iec,jsc:jec,:)
!
!  else
!
!    select case (trim(state%fields(field)%fv3jedi_name))
!
!    case ("ud")
!      Atm%u(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!      update_winds = .true.
!    case ("vd")
!      Atm%v(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!      update_winds = .true.
!    case ("ua")
!      Atm%ua(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!    case ("va")
!      Atm%va(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!    case ("t")
!      Atm%pt(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!    case ("delp")
!      Atm%delp(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!      print*, 'delp', maxval(Atm%delp), minval(Atm%delp), Atm%ptop
!      call compute_fv3_pressures( isc, iec, jsc, jec, isd, ied, jsd, jed, &
!                                  npz, kappa, Atm%ptop, &
!                                  Atm%delp, Atm%pe, Atm%pk, Atm%pkz, Atm%peln )
!    case ("phis")
!      Atm%phis(isc:iec,jsc:jec) = field_ptr(isc:iec,jsc:jec,1)
!    case ("w")
!      Atm%w(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!    case ("delz")
!      Atm%delz(isc:iec,jsc:jec,:) = field_ptr(isc:iec,jsc:jec,:)
!
!    end select
!
!  endif
!
!enddo
!
!
!!Update edges of d-grid winds
!!----------------------------
!if (update_winds) then
!  allocate(wbuffery(jsc:jec,npz))
!  allocate(sbufferx(isc:iec,npz))
!  allocate(ebuffery(jsc:jec,npz))
!  allocate(nbufferx(isc:iec,npz))
!
!  call mpp_get_boundary(Atm%u, Atm%v, Atm%domain, &
!                        wbuffery=wbuffery, ebuffery=ebuffery, &
!                        sbufferx=sbufferx, nbufferx=nbufferx, &
!                        gridtype=DGRID_NE, complete=.true. )
!  do k=1,npz
!     do i=isc,iec
!        Atm%u(i,jec+1,k) = nbufferx(i,k)
!     enddo
!  enddo
!  do k=1,npz
!     do j=jsc,jec
!        Atm%v(iec+1,j,k) = ebuffery(j,k)
!     enddo
!  enddo
!endif
!
!! Fill phi halos
!! --------------
!call mpp_update_domains(Atm%phis, Atm%domain, complete=.true.)
!
!end subroutine state_to_fv3
!
!! --------------------------------------------------------------------------------------------------
!
!subroutine fv3_to_state( Atm, state )
!
!implicit none
!type(FV_Atmos_type), intent(in)    :: Atm
!type(fv3jedi_state), intent(inout) :: state
!
!integer :: field, tindex
!integer :: isc, iec, jsc, jec, npz
!integer :: isd, ied, jsd, jed
!
!real(kind=kind_real), pointer, dimension(:,:,:) :: field_ptr
!
!! Convenience
!! -----------
!isc = Atm%bd%isc
!iec = Atm%bd%iec
!jsc = Atm%bd%jsc
!jec = Atm%bd%jec
!isd = Atm%bd%isd
!ied = Atm%bd%ied
!jsd = Atm%bd%jsd
!jed = Atm%bd%jed
!npz = Atm%npz
!
!! Loop through fields on fv3-jedi side and copy
!!----------------------------------------------
!do field = 1,state%nf
!
!  ! Get pointer to fv3-jedi field
!  call pointer_field_array(state%fields, state%fields(field)%fv3jedi_name, field_ptr)
!
!  ! Tracer names match between systems
!  if (state%fields(field)%tracer) then
!
!    ! Get index on fv3 model side and fill
!    tindex = get_tracer_index(MODEL_ATMOS, trim(state%fields(field)%fv3jedi_name))
!    field_ptr(isc:iec,jsc:jec,:) = real(Atm%q(isc:iec,jsc:jec,:,tindex),kind_real)
!
!  else
!
!    select case (trim(state%fields(field)%fv3jedi_name))
!
!    case ("ud")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%u(isc:iec,jsc:jec,:),kind_real)
!    case ("vd")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%v(isc:iec,jsc:jec,:),kind_real)
!    case ("ua")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%ua(isc:iec,jsc:jec,:),kind_real)
!    case ("va")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%va(isc:iec,jsc:jec,:),kind_real)
!    case ("t")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%pt(isc:iec,jsc:jec,:),kind_real)
!    case ("delp")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%delp(isc:iec,jsc:jec,:),kind_real)
!    case ("phis")
!      field_ptr(isc:iec,jsc:jec,1) = real(Atm%phis(isc:iec,jsc:jec),kind_real)
!    case ("w")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%w(isc:iec,jsc:jec,:),kind_real)
!    case ("delz")
!      field_ptr(isc:iec,jsc:jec,:) = real(Atm%delz(isc:iec,jsc:jec,:),kind_real)
!
!    end select
!
!  endif
!
!enddo
!
!end subroutine fv3_to_state
!
!! --------------------------------------------------------------------------------------------------
!
!subroutine compute_fv3_pressures_r4( is, ie, js, je, isd, ied, jsd, jed, npz, &
!                                     kappa, ptop, delp, pe, pk, pkz, peln )
!
!implicit none
!integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
!real(8), intent(in) :: kappa, ptop
!real(4), intent(in) :: delp(isd:ied, jsd:jed, npz)
!real(4), intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
!real(4), intent(out) :: pk(is:ie, js:je, npz+1)
!real(4), intent(out) :: peln(is:ie, npz+1, js:je)
!real(4), intent(out) :: pkz(is:ie, js:je, npz)
!integer :: i, j, k
!
!pe(:, :, :) = 0.0
!pe(:, 1, :) = real(ptop,4)
!do k=2,npz+1
!  do j=js,je
!    do i=is,ie
!      pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
!    end do
!  end do
!end do
!
!do k=1,npz+1
!  do j=js,je
!    do i=is,ie
!      peln(i, k, j) = log(pe(i, k, j))
!    end do
!  end do
!end do
!
!do k=1,npz+1
!  do j=js,je
!    do i=is,ie
!      pk(i, j, k) = exp(real(kappa,4)*peln(i, k, j))
!    end do
!  end do
!end do
!
!do k=1,npz
!  do j=js,je
!    do i=is,ie
!      pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(real(kappa,4)*(peln(i, k+1, j)-peln(i, k, j)))
!    end do
!  end do
!end do
!
!end subroutine compute_fv3_pressures_r4
!
!! --------------------------------------------------------------------------------------------------
!
!subroutine compute_fv3_pressures_r8( is, ie, js, je, isd, ied, jsd, jed, npz, &
!                                  kappa, ptop, delp, pe, pk, pkz, peln )
!
!implicit none
!integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
!real(8), intent(in) :: kappa, ptop
!real(8), intent(in) :: delp(isd:ied, jsd:jed, npz)
!real(8), intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
!real(8), intent(out) :: pk(is:ie, js:je, npz+1)
!real(8), intent(out) :: peln(is:ie, npz+1, js:je)
!real(8), intent(out) :: pkz(is:ie, js:je, npz)
!integer :: i, j, k
!
!pe(:, :, :) = 0.0
!pe(:, 1, :) = ptop
!do k=2,npz+1
!  do j=js,je
!    do i=is,ie
!      pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
!    end do
!  end do
!end do
!
!do k=1,npz+1
!  do j=js,je
!    do i=is,ie
!      peln(i, k, j) = log(pe(i, k, j))
!    end do
!  end do
!end do
!
!do k=1,npz+1
!  do j=js,je
!    do i=is,ie
!      pk(i, j, k) = exp(kappa*peln(i, k, j))
!    end do
!  end do
!end do
!
!do k=1,npz
!  do j=js,je
!    do i=is,ie
!      pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1, j)-peln(i, k, j)))
!    end do
!  end do
!end do
!
!end subroutine compute_fv3_pressures_r8

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fv3_mod
