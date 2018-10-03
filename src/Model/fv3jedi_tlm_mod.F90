! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_tlm_mod

use iso_c_binding
use config_mod
use duration_mod
use fv3jedi_geom_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 
use fv3jedi_traj_mod
use fv3jedi_constants
use kinds

use fv_arrays_mod,  only: fv_atmos_type
use mpp_mod,        only: mpp_pe, mpp_root_pe 
use mpp_domains_mod, only: mpp_update_domains, mpp_get_boundary, DGRID_NE
use pressure_vt_mod, only: compute_fv3_pressures, compute_fv3_pressures_tlm, compute_fv3_pressures_bwd
use fms_mod, only: set_domain, nullify_domain

use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index

use fv_arrays_nlm_mod, only: fv_atmos_pert_type

implicit none
private

public :: fv3jedi_tlm 
public :: tlm_create
public :: tlm_delete
public :: tlm_initialize_tl
public :: tlm_initialize_ad
public :: tlm_step_tl
public :: tlm_step_ad
public :: tlm_finalize_tl
public :: tlm_finalize_ad

! ------------------------------------------------------------------------------

!> Fortran derived type to hold tlm definition
type:: fv3jedi_tlm
  real(kind=kind_real)                         :: DT                  !<TLM big timestep
  real(kind_real), allocatable, dimension(:,:) :: ebuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: nbufferx            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: wbuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: sbufferx            !<Halo holder
  type(fv_atmos_type), allocatable             :: FV_Atm(:)           !<Main FV3 construct 
  type(fv_atmos_pert_type), allocatable        :: FV_AtmP(:)          !<Main FV3 construct perturbation variables
  logical, allocatable                         :: grids_on_this_pe(:) !<FV3 record
  integer                                      :: p_split = 1         !<FV3 record
  integer                                      :: isc,iec,jsc,jec     !<Convenience pointer to grid
  integer                                      :: isd,ied,jsd,jed     !<Convenience pointer to grid
  integer                                      :: npz                 !<Convenience
  logical                                      :: hydrostatic         !<Convenience
  integer                                      :: cp_dyn_ind          !<Module index for checkpointing
  integer                                      :: update_dgridwind=1  !<Update the fv3 pressures each time step
  integer                                      :: update_pressures=1  !<Update the fv3 pressures each time step
  integer                                      :: init_tlmadm         !<Is the TLM/ADM needed in this instance
  integer                                      :: linmodtest          !<Testing the linear model
  integer                                      :: ti_q                !<Tracer index for specific humidity
  integer                                      :: ti_qi               !<Tracer index for cloud ice water
  integer                                      :: ti_ql               !<Tracer index for cloud liquid water
  integer                                      :: ti_o3               !<Tracer index for ozone
end type fv3jedi_tlm

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine tlm_create(tlm, geom, c_conf)

use fv_control_mod, only: fv_init, pelist_all

use fv_control_nlm_mod, only: fv_init_pert
use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter

implicit none
type(c_ptr), intent(in)     :: c_conf !< pointer to object of class Config
type(fv3jedi_tlm), target :: tlm  ! should I put intent on these?
type(fv3jedi_geom)          :: geom

character(len=20) :: ststep
type(duration) :: dtstep

integer :: i,j
integer :: tmp

!TLM time step
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
tlm%DT = real(duration_seconds(dtstep),kind_real)

!Call to fv_init
call fv_init(tlm%FV_Atm, tlm%DT, tlm%grids_on_this_pe, tlm%p_split)
deallocate(pelist_all)

!Compute grid must be same as geometry
if ( (geom%isc .ne. tlm%FV_Atm(1)%bd%isc) .or. (geom%iec .ne. tlm%FV_Atm(1)%bd%iec) .or. &
     (geom%jsc .ne. tlm%FV_Atm(1)%bd%jsc) .or. (geom%jec .ne. tlm%FV_Atm(1)%bd%jec) .or. &
     (geom%npz .ne. tlm%FV_Atm(1)%npz) ) then
   call abor1_ftn("fv3jedi tlm: compute areas for geometry and tlm do not agree")
endif

!Copy of grid info for convenience
tlm%isc = tlm%FV_Atm(1)%bd%isc
tlm%iec = tlm%FV_Atm(1)%bd%iec
tlm%jsc = tlm%FV_Atm(1)%bd%jsc
tlm%jec = tlm%FV_Atm(1)%bd%jec
tlm%isd = tlm%FV_Atm(1)%bd%isd
tlm%ied = tlm%FV_Atm(1)%bd%ied
tlm%jsd = tlm%FV_Atm(1)%bd%jsd
tlm%jed = tlm%FV_Atm(1)%bd%jed
tlm%npz = tlm%FV_Atm(1)%npz
tlm%hydrostatic = tlm%FV_Atm(1)%flagstruct%hydrostatic

!Halo holders for domain grid
allocate(tlm%wbuffery(tlm%FV_Atm(1)%bd%jsc:tlm%FV_Atm(1)%bd%jec,tlm%FV_Atm(1)%npz))
allocate(tlm%sbufferx(tlm%FV_Atm(1)%bd%isc:tlm%FV_Atm(1)%bd%iec,tlm%FV_Atm(1)%npz))
allocate(tlm%ebuffery(tlm%FV_Atm(1)%bd%jsc:tlm%FV_Atm(1)%bd%jec,tlm%FV_Atm(1)%npz))
allocate(tlm%nbufferx(tlm%FV_Atm(1)%bd%isc:tlm%FV_Atm(1)%bd%iec,tlm%FV_Atm(1)%npz))

!Set ptop, ak, bk to be same as geometry
tlm%FV_Atm(1)%ak = geom%ak
tlm%FV_Atm(1)%bk = geom%bk
tlm%FV_Atm(1)%ptop = geom%ptop

!Tracer indexes
tlm%ti_q  = 1 !get_tracer_index (MODEL_ATMOS, 'sphum')
tlm%ti_ql = 2 !get_tracer_index (MODEL_ATMOS, 'liq_wat')
tlm%ti_qi = 3 !get_tracer_index (MODEL_ATMOS, 'ice_wat')
tlm%ti_o3 = 4 !get_tracer_index (MODEL_ATMOS, 'o3mr')

!Always allocate w, delz, q_con for now
deallocate(tlm%FV_Atm(1)%w)
deallocate(tlm%FV_Atm(1)%delz)
deallocate(tlm%FV_Atm(1)%q_con)
allocate  ( tlm%FV_Atm(1)%w (tlm%FV_Atm(1)%bd%isd:tlm%FV_Atm(1)%bd%ied,tlm%FV_Atm(1)%bd%jsd:tlm%FV_Atm(1)%bd%jed,&
                               tlm%FV_Atm(1)%flagstruct%npz) )
allocate  ( tlm%FV_Atm(1)%delz (tlm%FV_Atm(1)%bd%isd:tlm%FV_Atm(1)%bd%ied,tlm%FV_Atm(1)%bd%jsd:tlm%FV_Atm(1)%bd%jed,&
                               tlm%FV_Atm(1)%flagstruct%npz) )
allocate  ( tlm%FV_Atm(1)%q_con(tlm%FV_Atm(1)%bd%isd:tlm%FV_Atm(1)%bd%ied,tlm%FV_Atm(1)%bd%jsd:tlm%FV_Atm(1)%bd%jed,&
                               tlm%FV_Atm(1)%flagstruct%npz) )
tlm%FV_Atm(1)%w = 0.0
tlm%FV_Atm(1)%delz = 0.0
tlm%FV_Atm(1)%q_con = 0.0

!fC and f0
if (tlm%FV_Atm(1)%flagstruct%grid_type == 4) then
   tlm%FV_Atm(1)%gridstruct%fC(:,:) = 2.*omega*sin(tlm%FV_Atm(1)%flagstruct%deglat/180.*pi)
   tlm%FV_Atm(1)%gridstruct%f0(:,:) = 2.*omega*sin(tlm%FV_Atm(1)%flagstruct%deglat/180.*pi)
else
   if (f_coriolis_angle == -999) then
      tlm%FV_Atm(1)%gridstruct%fC(:,:) = 0.0
      tlm%FV_Atm(1)%gridstruct%f0(:,:) = 0.0
   else
      do j=tlm%FV_Atm(1)%bd%jsd,tlm%FV_Atm(1)%bd%jed+1
         do i=tlm%FV_Atm(1)%bd%isd,tlm%FV_Atm(1)%bd%ied+1
            tlm%FV_Atm(1)%gridstruct%fC(i,j) = 2.*omega*( -COS(tlm%FV_Atm(1)%gridstruct%grid(i,j,1))*&
                                           COS(tlm%FV_Atm(1)%gridstruct%grid(i,j,2))*SIN(f_coriolis_angle) + &
                                           SIN(tlm%FV_Atm(1)%gridstruct%grid(i,j,2))*COS(f_coriolis_angle) )
         enddo
      enddo
      do j=tlm%FV_Atm(1)%bd%jsd,tlm%FV_Atm(1)%bd%jed
         do i=tlm%FV_Atm(1)%bd%isd,tlm%FV_Atm(1)%bd%ied
            tlm%FV_Atm(1)%gridstruct%f0(i,j) = 2.*omega*( -COS(tlm%FV_Atm(1)%gridstruct%agrid(i,j,1))*&
                                           COS(tlm%FV_Atm(1)%gridstruct%agrid(i,j,2))*SIN(f_coriolis_angle) + &
                                           SIN(tlm%FV_Atm(1)%gridstruct%agrid(i,j,2))*COS(f_coriolis_angle) )
         enddo
      enddo
   endif
endif

if (config_element_exists(c_conf,"update_dgridwind")) tlm%update_dgridwind = config_get_int(c_conf,"update_dgridwind")
if (config_element_exists(c_conf,"update_pressures")) tlm%update_pressures = config_get_int(c_conf,"update_pressures")

!Pointer to self when not nested
if (.not. tlm%FV_Atm(1)%gridstruct%nested) tlm%FV_Atm(1)%parent_grid => tlm%FV_Atm(1)

!Harwire some flags
tlm%FV_Atm(1)%flagstruct%reproduce_sum = .false.
tlm%FV_Atm(1)%flagstruct%fill = .false.
tlm%FV_Atm(1)%flagstruct%fv_debug = .false.
tlm%FV_Atm(1)%flagstruct%adiabatic = .false.
tlm%FV_Atm(1)%flagstruct%do_sat_adj = .false.
tlm%FV_Atm(1)%flagstruct%breed_vortex_inline = .false.

tlm%init_tlmadm = 0
if (config_element_exists(c_conf,"init_tlmadm")) &
  tlm%init_tlmadm = config_get_int(c_conf,"init_tlmadm")

tlm%linmodtest = 0
if (config_element_exists(c_conf,"linmodtest")) &
  tlm%linmodtest = config_get_int(c_conf,"linmodtest")

if (tlm%linmodtest == 1) tlm%init_tlmadm = 1

if (tlm%init_tlmadm == 1) then

   !Initialize perturbation variables and read config
   call fv_init_pert(tlm%FV_Atm,tlm%FV_AtmP)
      
   !Set up controls for iterative checkpointing
   !-------------------------------------------
   
   !Global
   cp_iter_controls%cp_i  = 0! config_get_int (c_conf,"CP_i")
   cp_iter_controls%cp_nt = 4!config_get_int (c_conf,"CP_nt")
   cp_iter_controls%cp_gb = -0.1!config_get_real(c_conf,"CP_gb")
   cp_iter_controls%cp_nm = 1!config_get_int (c_conf,"CP_nm")
   call initialize_cp_iter
   
   if (cp_iter_controls%cp_i .ne. 0) then
   
      !Dynamics
      tlm%cp_dyn_ind = config_get_int (c_conf,"CP_dyn_ind")
      cp_iter(tlm%cp_dyn_ind)%my_name(1:3) = 'dyn'
      
      cp_iter(tlm%cp_dyn_ind)%cp_test = .false.
      tmp = config_get_int (c_conf,"CP_dyn_test")
      if (tmp==1) cp_iter(tlm%cp_dyn_ind)%cp_test = .true.
      
      cp_iter(tlm%cp_dyn_ind)%cp_rep = .false.
      tmp = config_get_int (c_conf,"CP_dyn_rep")
      if (tmp==1) cp_iter(tlm%cp_dyn_ind)%cp_test = .true.
      
      !Hardwire these for now
      cp_iter(tlm%cp_dyn_ind)%check_st_control = .false.
      cp_iter(tlm%cp_dyn_ind)%check_st_integer = .false.
      cp_iter(tlm%cp_dyn_ind)%check_st_real_r4 = .false.
      cp_iter(tlm%cp_dyn_ind)%check_st_real_r8 = .false.
      
      cp_iter(tlm%cp_dyn_ind)%test_dim_st_control = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_st_integer = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_st_real_r4 = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_st_real_r8 = 0
      
      cp_iter(tlm%cp_dyn_ind)%test_dim_cp_control = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_cp_integer = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_cp_real_r4 = 0
      cp_iter(tlm%cp_dyn_ind)%test_dim_cp_real_r8 = 0
   
   endif

endif

end subroutine tlm_create

! ------------------------------------------------------------------------------

subroutine tlm_delete(self)

use fv_arrays_mod, only: deallocate_fv_atmos_type

use fv_arrays_nlm_mod, only: deallocate_fv_atmos_pert_type
use tapenade_iter, only: cp_iter_controls, finalize_cp_iter

implicit none
type(fv3jedi_tlm) :: self

deallocate(self%ebuffery)
deallocate(self%wbuffery)
deallocate(self%nbufferx)
deallocate(self%sbufferx)

call deallocate_fv_atmos_type(self%FV_Atm(1))
deallocate(self%FV_Atm)

if (self%init_tlmadm == 1) then
   call deallocate_fv_atmos_pert_type(self%FV_AtmP(1))
   deallocate(self%FV_AtmP)
   
   if (cp_iter_controls%cp_i .ne. 0) call finalize_cp_iter
endif

end subroutine tlm_delete

! ------------------------------------------------------------------------------

subroutine tlm_initialize_ad(self, inc)

implicit none
type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)     :: inc

end subroutine tlm_initialize_ad

! ------------------------------------------------------------------------------

subroutine tlm_initialize_tl(self, inc)

implicit none
type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)     :: inc

end subroutine tlm_initialize_tl

! ------------------------------------------------------------------------------

subroutine tlm_step_ad(geom, self, inc, traj)

use fv_dynamics_adm_mod, only: fv_dynamics_fwd, fv_dynamics_bwd
use mpp_domains_mod, only: mpp_get_boundary_ad
use tapenade_iter,   only: cp_iter_controls, cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray

implicit none

type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)     :: inc
type(fv3jedi_traj)    :: traj
type(fv3jedi_geom)          :: geom

type(fv_atmos_type), pointer :: FV_Atm(:)
type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
integer :: i,j,k


if (mpp_pe() == mpp_root_pe()) print*, 'Step adjoint'

! Convenience pointer to the main FV_Atm structure
! ------------------------------------------------
FV_Atm => self%FV_Atm
FV_AtmP => self%FV_AtmP


! Get up the trajectory for this time step 
! ----------------------------------------
call traj_get( traj, &
               self%isc,self%iec,self%jsc,self%jec, &
               self%isd,self%ied,self%jsd,self%jed, &
               self%npz, self%hydrostatic, &
               self%FV_Atm(1)%u,self%FV_Atm(1)%v,self%FV_Atm(1)%pt,self%FV_Atm(1)%delp,&
               self%FV_Atm(1)%q(:,:,:,self%ti_q ),self%FV_Atm(1)%q(:,:,:,self%ti_qi),&
               self%FV_Atm(1)%q(:,:,:,self%ti_ql),self%FV_Atm(1)%q(:,:,:,self%ti_o3),&
               self%FV_Atm(1)%w,self%FV_Atm(1)%delz,self%FV_Atm(1)%phis )


! Fill phi halos
! --------------
call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)


! Zero parts of the trajectory that are recomputed or outputs
! -----------------------------------------------------------
FV_Atm(1)%pe    = 0.0
FV_Atm(1)%peln  = 0.0
FV_Atm(1)%pk    = 0.0
FV_Atm(1)%pkz   = 0.0
FV_Atm(1)%ua    = 0.0
FV_Atm(1)%va    = 0.0
FV_Atm(1)%uc    = 0.0
FV_Atm(1)%vc    = 0.0
FV_Atm(1)%omga  = 0.0
FV_Atm(1)%mfx   = 0.0
FV_Atm(1)%mfy   = 0.0
FV_Atm(1)%cx    = 0.0
FV_Atm(1)%cy    = 0.0
FV_Atm(1)%ze0   = 0.0
FV_Atm(1)%q_con = 0.0


! Update edges of traj d-grid winds
! ---------------------------------
call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                      wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                      sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                      gridtype=DGRID_NE, complete=.true. )
do k=1,self%npz
   do i=self%isc,self%iec
      FV_Atm(1)%u(i,self%jec+1,k) = self%nbufferx(i,k)
   enddo
enddo
do k=1,self%npz
   do j=self%jsc,self%jec
      FV_Atm(1)%v(self%iec+1,j,k) = self%ebuffery(j,k)
   enddo
enddo


!Compute the other pressure variables needed by FV3
!--------------------------------------------------
if (self%update_pressures == 1) then
   call compute_fv3_pressures( self%isc, self%iec, self%jsc, self%jec, self%isd, self%ied, self%jsd, self%jed, &
                               self%npz, kappa, FV_Atm(1)%ptop, &
                               FV_Atm(1)%delp, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%pkz, FV_Atm(1)%peln )
endif


! MPP set domain
! --------------
call set_domain(FV_Atm(1)%domain)


! Initilize the module level checkpointing
! ----------------------------------------
if (cp_iter_controls%cp_i .ne. 0) then
   call cp_mod_ini(self%cp_dyn_ind)
endif


! Forward sweep of the dynamics with saving of checkpoints for use in backward sweep
! ----------------------------------------------------------------------------------
if (cp_iter_controls%cp_i <= 3) then

   call fv_dynamics_fwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                        self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                        FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                        cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                        FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                        FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w,                                                           &
                        FV_Atm(1)%delz, FV_Atm(1)%flagstruct%hydrostatic,                                                &
                        FV_Atm(1)%pt, FV_Atm(1)%delp,                                                                    &
                        FV_Atm(1)%q, FV_Atm(1)%ps, FV_Atm(1)%pe,                                                         &
                        FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                                                     &
                        FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                                                 &
                        FV_Atm(1)%ua, FV_Atm(1)%va,                                                                      &
                        FV_Atm(1)%uc, FV_Atm(1)%vc,                                                                      &
                        FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                        FV_Atm(1)%mfx, FV_Atm(1)%mfy,                                                                    &
                        FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,                                                       &
                        FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
                        FV_AtmP(1)%flagstruct,                                                                           &
                        FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )

   if (cp_iter_controls%cp_i .ne. 0) then
      !Push end of timestep trajectory to stack
      call PUSHREALARRAY(FV_Atm(1)%u   ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%v   ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%w   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%delz,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%pt  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%delp,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%q   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz*FV_Atm(1)%ncnst)
      call PUSHREALARRAY(FV_Atm(1)%ps  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
      call PUSHREALARRAY(FV_Atm(1)%pe  ,(self%iec-self%isc+3)*(self%jec-self%jsc+3)*(self%npz+1))
      call PUSHREALARRAY(FV_Atm(1)%pk  ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
      call PUSHREALARRAY(FV_Atm(1)%peln,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
      call PUSHREALARRAY(FV_Atm(1)%pkz ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%phis,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
      call PUSHREALARRAY(FV_Atm(1)%omga,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%ua  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%va  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%uc  ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%vc  ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%mfx ,(self%iec-self%isc+2)*(self%jec-self%jsc+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%mfy ,(self%iec-self%isc+1)*(self%jec-self%jsc+2)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%cx  ,(self%iec-self%isc+2)*(self%jed-self%jsd+1)*self%npz)
      call PUSHREALARRAY(FV_Atm(1)%cy  ,(self%ied-self%isd+1)*(self%jec-self%jsc+2)*self%npz)
      !Trick checkpoint schemes into not considering these superfluous checkpoints,
      !about to recover with the pop anyway.
      FV_Atm(1)%u    = 2.0_kind_real*FV_Atm(1)%u
      FV_Atm(1)%v    = 2.0_kind_real*FV_Atm(1)%v
      FV_Atm(1)%w    = 2.0_kind_real*FV_Atm(1)%w
      FV_Atm(1)%delz = 2.0_kind_real*FV_Atm(1)%delz
      FV_Atm(1)%pt   = 2.0_kind_real*FV_Atm(1)%pt
      FV_Atm(1)%delp = 2.0_kind_real*FV_Atm(1)%delp
      FV_Atm(1)%q    = 2.0_kind_real*FV_Atm(1)%q
      FV_Atm(1)%ps   = 2.0_kind_real*FV_Atm(1)%ps
      FV_Atm(1)%pe   = 2.0_kind_real*FV_Atm(1)%pe
      FV_Atm(1)%pk   = 2.0_kind_real*FV_Atm(1)%pk
      FV_Atm(1)%peln = 2.0_kind_real*FV_Atm(1)%peln
      FV_Atm(1)%pkz  = 2.0_kind_real*FV_Atm(1)%pkz
      FV_Atm(1)%phis = 2.0_kind_real*FV_Atm(1)%phis
      FV_Atm(1)%omga = 2.0_kind_real*FV_Atm(1)%omga
      FV_Atm(1)%ua   = 2.0_kind_real*FV_Atm(1)%ua
      FV_Atm(1)%va   = 2.0_kind_real*FV_Atm(1)%va
      FV_Atm(1)%uc   = 2.0_kind_real*FV_Atm(1)%uc
      FV_Atm(1)%vc   = 2.0_kind_real*FV_Atm(1)%vc
      FV_Atm(1)%mfx  = 2.0_kind_real*FV_Atm(1)%mfx
      FV_Atm(1)%mfy  = 2.0_kind_real*FV_Atm(1)%mfy
      FV_Atm(1)%cx   = 2.0_kind_real*FV_Atm(1)%cx
      FV_Atm(1)%cy   = 2.0_kind_real*FV_Atm(1)%cy
   endif

endif


! Checkpoint mid point, reset counters etc
! ----------------------------------------
if (cp_iter_controls%cp_i .ne. 0) then 
   call cp_mod_mid
endif

if (cp_iter_controls%cp_i .ne. 0) then
   !Populate end of timestep trajectory from stack
   call POPREALARRAY(FV_Atm(1)%cy  ,(self%ied-self%isd+1)*(self%jec-self%jsc+2)*self%npz)
   call POPREALARRAY(FV_Atm(1)%cx  ,(self%iec-self%isc+2)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%mfy ,(self%iec-self%isc+1)*(self%jec-self%jsc+2)*self%npz)
   call POPREALARRAY(FV_Atm(1)%mfx ,(self%iec-self%isc+2)*(self%jec-self%jsc+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%vc  ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
   call POPREALARRAY(FV_Atm(1)%uc  ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%va  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%ua  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%omga,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%phis,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
   call POPREALARRAY(FV_Atm(1)%pkz ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%peln,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
   call POPREALARRAY(FV_Atm(1)%pk  ,(self%iec-self%isc+1)*(self%jec-self%jsc+1)*(self%npz+1))
   call POPREALARRAY(FV_Atm(1)%pe  ,(self%iec-self%isc+3)*(self%jec-self%jsc+3)*(self%npz+1))
   call POPREALARRAY(FV_Atm(1)%ps  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1))
   call POPREALARRAY(FV_Atm(1)%q   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz*FV_Atm(1)%ncnst)
   call POPREALARRAY(FV_Atm(1)%delp,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%pt  ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%delz,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%w   ,(self%ied-self%isd+1)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%v   ,(self%ied-self%isd+2)*(self%jed-self%jsd+1)*self%npz)
   call POPREALARRAY(FV_Atm(1)%u   ,(self%ied-self%isd+1)*(self%jed-self%jsd+2)*self%npz)
endif


! Make sure everything is zero
! ----------------------------
call zero_pert_vars(FV_AtmP(1))


! Copy to increment variables
! ---------------------------
call tlm_to_inc_ad(geom,inc,self)


! Backward adjoint sweep of the dynamics
! --------------------------------------
call fv_dynamics_bwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                     self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                     FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                     cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                     FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                     FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
                     FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
                     FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
                     FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
                     FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
                     FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
                     FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
                     FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
                     FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                     FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
                     FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                      &
                     FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
                     FV_AtmP(1)%flagstruct,                                                                           &
                     FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )


!Adjoint of compute the other pressure variables needed by FV3
!-------------------------------------------------------------
if (self%update_pressures == 1) then
   call compute_fv3_pressures_bwd( self%isc, self%iec, self%jsc, self%jec, &
                                   self%isd, self%ied, self%jsd, self%jed, &
                                   self%npz, kappa, FV_Atm(1)%ptop, &
                                   FV_Atm(1)%delp, FV_AtmP(1)%delpp, &
                                   FV_Atm(1)%pe, FV_AtmP(1)%pep, &
                                   FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                   FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, &
                                   FV_Atm(1)%peln, FV_AtmP(1)%pelnp )
endif


!Edge of pert always needs to be filled
!--------------------------------------
self%nbufferx = 0.0_kind_real
do k=1,self%npz
   do i=self%isc,self%iec
      self%nbufferx(i,k) = FV_AtmP(1)%up(i,self%jec+1,k)
   enddo
enddo
self%ebuffery = 0.0_kind_real
do k=1,self%npz
   do j=self%jsc,self%jec
      self%ebuffery(j,k) = FV_AtmP(1)%vp(self%iec+1,j,k)
   enddo
enddo

call mpp_get_boundary_ad( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                          wbuffery=self%wbuffery, ebuffery=self%ebuffery, sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                          gridtype=DGRID_NE, complete=.true. )


! MPP nulify
! ----------
call nullify_domain()


! Copy back to increment 
! ----------------------
call inc_to_tlm_ad(geom,self,inc)


! Make sure everything is zero
! ----------------------------
call zero_pert_vars(FV_AtmP(1))


! Finalize iterative step and get ready for next iteration
! --------------------------------------------------------
if (cp_iter_controls%cp_i .ne. 0) then
   call cp_mod_end
endif

end subroutine tlm_step_ad

! ------------------------------------------------------------------------------

subroutine tlm_step_tl(geom, self, inc, traj)

use fv_dynamics_tlm_mod, only: fv_dynamics_tlm

implicit none
type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)         :: inc
type(fv3jedi_traj)    :: traj
type(fv3jedi_geom)          :: geom

type(fv_atmos_type), pointer :: FV_Atm(:)
type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
integer :: i,j,k


if (mpp_pe() == mpp_root_pe()) print*, 'Step tangent linear model'


!Convenience pointer to the main FV_Atm structure
!------------------------------------------------
FV_Atm => self%FV_Atm
FV_AtmP => self%FV_AtmP


!Get up the trajectory for this time step 
!----------------------------------------
call traj_get( traj, &
               self%isc,self%iec,self%jsc,self%jec, &
               self%isd,self%ied,self%jsd,self%jed, &
               self%npz, self%hydrostatic, &
               self%FV_Atm(1)%u,self%FV_Atm(1)%v,self%FV_Atm(1)%pt,self%FV_Atm(1)%delp,&
               self%FV_Atm(1)%q(:,:,:,self%ti_q ),self%FV_Atm(1)%q(:,:,:,self%ti_qi),&
               self%FV_Atm(1)%q(:,:,:,self%ti_ql),self%FV_Atm(1)%q(:,:,:,self%ti_o3),&
               self%FV_Atm(1)%w,self%FV_Atm(1)%delz,self%FV_Atm(1)%phis )


! Fill phi halos
! --------------
call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)


! Zero parts of the trajectory that are recomputed or outputs
! -----------------------------------------------------------
FV_Atm(1)%pe    = 0.0
FV_Atm(1)%peln  = 0.0
FV_Atm(1)%pk    = 0.0
FV_Atm(1)%pkz   = 0.0
FV_Atm(1)%ua    = 0.0
FV_Atm(1)%va    = 0.0
FV_Atm(1)%uc    = 0.0
FV_Atm(1)%vc    = 0.0
FV_Atm(1)%omga  = 0.0
FV_Atm(1)%mfx   = 0.0
FV_Atm(1)%mfy   = 0.0
FV_Atm(1)%cx    = 0.0
FV_Atm(1)%cy    = 0.0
FV_Atm(1)%ze0   = 0.0
FV_Atm(1)%q_con = 0.0


!Update edges of d-grid winds
!----------------------------
call mpp_get_boundary( FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                       wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                       sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                       gridtype=DGRID_NE, complete=.true. )
do k=1,self%npz
   do i=self%isc,self%iec
      FV_Atm(1)%u(i,self%jec+1,k) = self%nbufferx(i,k)
   enddo
enddo
do k=1,self%npz
   do j=self%jsc,self%jec
      FV_Atm(1)%v(self%iec+1,j,k) = self%ebuffery(j,k)
   enddo
enddo


! Make sure everything is zero
! ----------------------------
call zero_pert_vars(FV_AtmP(1))


!Copy to tlm variables
!---------------------
call inc_to_tlm_tl(geom,inc,self)


!Edge of pert always needs to be filled
!--------------------------------------
call mpp_get_boundary( FV_AtmP(1)%up, FV_AtmP(1)%vp, FV_Atm(1)%domain, &
                       wbuffery=self%wbuffery, ebuffery=self%ebuffery, &
                       sbufferx=self%sbufferx, nbufferx=self%nbufferx, &
                       gridtype=DGRID_NE, complete=.true. )
do k=1,self%npz
   do i=self%isc,self%iec
      FV_AtmP(1)%up(i,self%jec+1,k) = self%nbufferx(i,k)
   enddo
enddo
do k=1,self%npz
   do j=self%jsc,self%jec
      FV_AtmP(1)%vp(self%iec+1,j,k) = self%ebuffery(j,k)
   enddo
enddo


!Compute the other pressure variables needed by FV3
!--------------------------------------------------
if (self%update_pressures == 1) then
   call compute_fv3_pressures_tlm( self%isc, self%iec, self%jsc, self%jec, &
                                   self%isd, self%ied, self%jsd, self%jed, &
                                   self%npz, kappa, FV_Atm(1)%ptop, &
                                   FV_Atm(1)%delp, FV_AtmP(1)%delpp, &
                                   FV_Atm(1)%pe, FV_AtmP(1)%pep, &
                                   FV_Atm(1)%pk, FV_AtmP(1)%pkp, &
                                   FV_Atm(1)%pkz, FV_AtmP(1)%pkzp, &
                                   FV_Atm(1)%peln, FV_AtmP(1)%pelnp )
endif


! MPP set domain
! --------------
call set_domain(FV_Atm(1)%domain)


!Step TLM one time step
!----------------------
call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
                     self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
                     FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
                     cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
                     FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
                     FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
                     FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
                     FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
                     FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
                     FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
                     FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
                     FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
                     FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
                     FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
                     FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
                     FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                       &
                     FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, &
                     FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )


! MPP nulify
! ----------
call nullify_domain()


! Copy back to increment
! ----------------------
call tlm_to_inc_tl(geom,self,inc)


! Make sure everything is zero
! ----------------------------
call zero_pert_vars(FV_AtmP(1))

end subroutine tlm_step_tl

! ------------------------------------------------------------------------------

subroutine tlm_finalize_ad(self, inc)

implicit none
type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)     :: inc

end subroutine tlm_finalize_ad

! ------------------------------------------------------------------------------

subroutine tlm_finalize_tl(self, inc)

implicit none
type(fv3jedi_tlm), target :: self
type(fv3jedi_increment)     :: inc

end subroutine tlm_finalize_tl

! ------------------------------------------------------------------------------

subroutine inc_to_tlm_tl(geom,inc,self)

use wind_vt_mod, only: a2d

implicit none
type(fv3jedi_geom), intent(inout)   :: geom
type(fv3jedi_increment), intent(in) :: inc
type(fv3jedi_tlm), intent(inout)  :: self

integer :: isc,iec,jsc,jec
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

allocate(ud(isc:iec  ,jsc:jec+1,1:geom%npz))
allocate(vd(isc:iec+1,jsc:jec  ,1:geom%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

self%FV_AtmP(1)%up    = 0.0
self%FV_AtmP(1)%vp    = 0.0
self%FV_AtmP(1)%ptp   = 0.0
self%FV_AtmP(1)%delpp = 0.0
self%FV_AtmP(1)%qp    = 0.0
self%FV_AtmP(1)%wp    = 0.0
self%FV_AtmP(1)%delzp = 0.0

call a2d(geom, inc%ua(isc:iec,jsc:jec  ,:), inc%va(isc:iec  ,jsc:jec,:), &
                   ud(isc:iec,jsc:jec+1,:),     vd(isc:iec+1,jsc:jec,:))

self%FV_AtmP(1)%up   (isc:iec,jsc:jec,:)            = ud      (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)            = vd      (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)            = inc%t   (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)            = inc%delp(isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_q ) = inc%q   (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_qi) = inc%qi  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_ql) = inc%ql  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_o3) = inc%o3  (isc:iec,jsc:jec,:)
if (.not. inc%hydrostatic) then
   self%FV_AtmP(1)%delzp(isc:iec,jsc:jec,:) = inc%delz(isc:iec,jsc:jec,:)
   self%FV_AtmP(1)%wp   (isc:iec,jsc:jec,:) = inc%w   (isc:iec,jsc:jec,:)
endif

deallocate(ud,vd)

end subroutine inc_to_tlm_tl

! ------------------------------------------------------------------------------

subroutine tlm_to_inc_tl(geom,self,inc)

use wind_vt_mod, only: d2a

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_tlm),     intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

inc%ua   = 0.0
inc%va   = 0.0
inc%t    = 0.0
inc%delp = 0.0
inc%q    = 0.0
inc%qi   = 0.0
inc%ql   = 0.0
inc%o3   = 0.0
if (.not. inc%hydrostatic) then
   inc%delz = 0.0
   inc%w    = 0.0
endif

call d2a(geom, self%FV_AtmP(1)%up, self%FV_AtmP(1)%vp, self%FV_AtmP(1)%uap, self%FV_AtmP(1)%vap)

inc%ua  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%uap  (isc:iec,jsc:jec,:)
inc%va  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%vap  (isc:iec,jsc:jec,:)
inc%t   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)
inc%delp(isc:iec,jsc:jec,:) = self%FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
inc%q   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_q )
inc%qi  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_qi)
inc%ql  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_ql)
inc%o3  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_o3)
if (.not. inc%hydrostatic) then
  inc%delz(isc:iec,jsc:jec,:) = self%FV_AtmP(1)%delzp(isc:iec,jsc:jec,:)
  inc%w   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%wp   (isc:iec,jsc:jec,:)
endif

end subroutine tlm_to_inc_tl

! ------------------------------------------------------------------------------

subroutine tlm_to_inc_ad(geom,inc,self)

use wind_vt_mod, only: d2a_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_tlm),     intent(inout) :: self

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

self%FV_AtmP(1)%up    = 0.0
self%FV_AtmP(1)%vp    = 0.0
self%FV_AtmP(1)%uap   = 0.0
self%FV_AtmP(1)%vap   = 0.0
self%FV_AtmP(1)%ptp   = 0.0
self%FV_AtmP(1)%delpp = 0.0
self%FV_AtmP(1)%qp    = 0.0
self%FV_AtmP(1)%wp    = 0.0
self%FV_AtmP(1)%delzp = 0.0

self%FV_AtmP(1)%uap   (isc:iec,jsc:jec,:)            = inc%ua  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%vap   (isc:iec,jsc:jec,:)            = inc%va  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)            = inc%t   (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)            = inc%delp(isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_q ) = inc%q   (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_qi) = inc%qi  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_ql) = inc%ql  (isc:iec,jsc:jec,:)
self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_o3) = inc%o3  (isc:iec,jsc:jec,:)
if (.not. inc%hydrostatic) then
   self%FV_AtmP(1)%delzp(isc:iec,jsc:jec,:) = inc%delz(isc:iec,jsc:jec,:)
   self%FV_AtmP(1)%wp   (isc:iec,jsc:jec,:) = inc%w   (isc:iec,jsc:jec,:)
endif

call d2a_ad(geom, self%FV_AtmP(1)%up, self%FV_AtmP(1)%vp, self%FV_AtmP(1)%uap, self%FV_AtmP(1)%vap)

self%FV_AtmP(1)%uap(isc:iec,jsc:jec,:) = 0.0
self%FV_AtmP(1)%vap(isc:iec,jsc:jec,:) = 0.0

end subroutine tlm_to_inc_ad

! ------------------------------------------------------------------------------

subroutine inc_to_tlm_ad(geom,self,inc)

use wind_vt_mod, only: a2d_ad

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_tlm),     intent(inout) :: self

integer :: isc,iec,jsc,jec
real(kind=kind_real), allocatable, dimension(:,:,:) :: ud,vd

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

allocate(ud(isc:iec  ,jsc:jec+1,1:geom%npz))
allocate(vd(isc:iec+1,jsc:jec  ,1:geom%npz))
ud = 0.0_kind_real
vd = 0.0_kind_real

inc%ua   = 0.0
inc%va   = 0.0
inc%t    = 0.0
inc%delp = 0.0
inc%q    = 0.0
inc%qi   = 0.0
inc%ql   = 0.0
inc%o3   = 0.0
if (.not. inc%hydrostatic) then
   inc%delz = 0.0
   inc%w    = 0.0
endif

ud      (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%up   (isc:iec,jsc:jec,:)
vd      (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%vp   (isc:iec,jsc:jec,:)
inc%t   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%ptp  (isc:iec,jsc:jec,:)
inc%delp(isc:iec,jsc:jec,:) = self%FV_AtmP(1)%delpp(isc:iec,jsc:jec,:)
inc%q   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_q )
inc%qi  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_qi)
inc%ql  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_ql)
inc%o3  (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%qp   (isc:iec,jsc:jec,:,self%ti_o3)
if (.not. inc%hydrostatic) then
  inc%delz(isc:iec,jsc:jec,:) = self%FV_AtmP(1)%delzp(isc:iec,jsc:jec,:)
  inc%w   (isc:iec,jsc:jec,:) = self%FV_AtmP(1)%wp   (isc:iec,jsc:jec,:)
endif

!Convert A to D
call a2d_ad(geom, inc%ua(isc:iec,jsc:jec  ,:), inc%va(isc:iec  ,jsc:jec,:), &
                      ud(isc:iec,jsc:jec+1,:),     vd(isc:iec+1,jsc:jec,:))

deallocate(ud,vd)

end subroutine inc_to_tlm_ad

! ------------------------------------------------------------------------------

subroutine get_phi_from_state(state, self)

implicit none
type(fv3jedi_state) :: state
type(fv3jedi_tlm) :: self

!To zero halo
self%FV_Atm(1)%phis = 0.0

!Get compute domain from state
self%FV_Atm(1)%phis(self%isc:self%iec,self%jsc:self%jec) = state%phis(self%isc:self%iec,self%jsc:self%jec)

!Fill halos
call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)

end subroutine get_phi_from_state

! ------------------------------------------------------------------------------

subroutine zero_pert_vars(FV_AtmP)

implicit none
type(fv_atmos_pert_type), intent(inout) :: FV_AtmP

!Prognostic
FV_AtmP%up = 0.0
FV_AtmP%vp = 0.0
FV_AtmP%ptp = 0.0
FV_AtmP%delpp = 0.0
FV_AtmP%qp = 0.0
FV_AtmP%wp = 0.0
FV_AtmP%delzP = 0.0

!Outputs
FV_AtmP%ze0p = 0.0
FV_AtmP%q_conp = 0.0
FV_AtmP%psp = 0.0
FV_AtmP%pep = 0.0
FV_AtmP%pkp = 0.0
FV_AtmP%pelnp = 0.0
FV_AtmP%pkzp = 0.0
FV_AtmP%omgap = 0.0
FV_AtmP%uap = 0.0
FV_AtmP%vap = 0.0
FV_AtmP%ucp = 0.0
FV_AtmP%vcp = 0.0
FV_AtmP%mfxp = 0.0
FV_AtmP%mfyp = 0.0
FV_AtmP%cxp = 0.0
FV_AtmP%cyp = 0.0

end subroutine zero_pert_vars

! ------------------------------------------------------------------------------

end module fv3jedi_tlm_mod
