! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_model_mod

use iso_c_binding
use config_mod
use duration_mod
use fv3jedi_geom_mod
use fv3jedi_fields_mod
use fv3jedi_trajectories
use fv3jedi_constants
use kinds

use fv_arrays_mod,  only: fv_atmos_type
use mpp_mod,        only: mpp_pe, mpp_root_pe 
use mpp_domains_mod, only: mpp_update_domains, mpp_get_boundary, DGRID_NE
use variable_transforms, only: compute_fv3_pressures, compute_fv3_pressures_tlm, compute_fv3_pressures_bwd
use fms_mod, only: set_domain, nullify_domain


#ifdef TLADPRES
use fv_arrays_nlm_mod, only: fv_atmos_pert_type
#endif

implicit none
private
public :: fv3jedi_model, & 
        & model_setup, model_delete, &
        & model_prepare_integration, model_prepare_integration_tl, model_prepare_integration_ad, &
        & model_propagate, model_propagate_tl, model_propagate_ad, &
        & model_prop_traj, model_wipe_traj, &
        & fv3jedi_model_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: fv3jedi_model
  real(kind=kind_real)                         :: DT                  !<Model big timestep
  real(kind_real), allocatable, dimension(:,:) :: ebuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: nbufferx            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: wbuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: sbufferx            !<Halo holder
  type(fv_atmos_type), allocatable             :: FV_Atm(:)           !<Main FV3 construct 
#ifdef TLADPRES
  type(fv_atmos_pert_type), allocatable        :: FV_AtmP(:)          !<Main FV3 construct perturbation variables
#endif
  logical, allocatable                         :: grids_on_this_pe(:) !<FV3 record
  integer                                      :: p_split = 1         !<FV3 record
  integer                                      :: isc,iec,jsc,jec     !<Convenience pointer to geom
  integer                                      :: isd,ied,jsd,jed     !<Convenience pointer to geom
  integer                                      :: npz                 !<Convenience
  integer                                      :: cp_dyn_ind          !<Module index for checkpointing
  integer                                      :: update_dgridwind=1  !<Update the fv3 pressures each time step
  integer                                      :: update_pressures=1  !<Update the fv3 pressures each time step
  integer                                      :: init_tlmadm         !<Is the TLM/ADM needed in this instance
end type fv3jedi_model

#define LISTED_TYPE fv3jedi_model

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine model_setup(model, geom, c_conf)

use fv_control_mod, only: fv_init, pelist_all

#ifdef TLADPRES
use fv_control_nlm_mod, only: fv_init_pert
use tapenade_iter, only: cp_iter, cp_iter_controls, initialize_cp_iter
#endif

implicit none
type(c_ptr), intent(in)    :: c_conf !< pointer to object of class Config
type(fv3jedi_model), target :: model  ! should I put intent on these?
type(fv3jedi_geom)          :: geom

character(len=20) :: ststep
type(duration) :: dtstep

integer :: i,j
integer :: tmp, fail

!Model time step
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
model%DT = real(duration_seconds(dtstep),kind_real)

!Call to fv_init
call fv_init(model%FV_Atm, model%DT, model%grids_on_this_pe, model%p_split)
deallocate(pelist_all)

!Compute grid must be same as geometry
if ( (geom%bd%isc .ne. model%FV_Atm(1)%bd%isc) .or. (geom%bd%iec .ne. model%FV_Atm(1)%bd%iec) .or. &
     (geom%bd%jsc .ne. model%FV_Atm(1)%bd%jsc) .or. (geom%bd%jec .ne. model%FV_Atm(1)%bd%jec) .or. &
     (geom%npz .ne. model%FV_Atm(1)%npz) ) then
   call abor1_ftn("fv3jedi model: compute areas for geometry and model do not agree")
endif

!Copy of geometry grid info, halos could be different between model and geometry
model%isc = model%FV_Atm(1)%bd%isc
model%iec = model%FV_Atm(1)%bd%iec
model%jsc = model%FV_Atm(1)%bd%jsc
model%jec = model%FV_Atm(1)%bd%jec
model%isd = model%FV_Atm(1)%bd%isd
model%ied = model%FV_Atm(1)%bd%ied
model%jsd = model%FV_Atm(1)%bd%jsd
model%jed = model%FV_Atm(1)%bd%jed
model%npz = model%FV_Atm(1)%npz

!Halo holders for domain grid
allocate(model%ebuffery(model%FV_Atm(1)%bd%jsd:model%FV_Atm(1)%bd%jed,model%FV_Atm(1)%npz))
allocate(model%wbuffery(model%FV_Atm(1)%bd%jsd:model%FV_Atm(1)%bd%jed,model%FV_Atm(1)%npz))
allocate(model%nbufferx(model%FV_Atm(1)%bd%isd:model%FV_Atm(1)%bd%ied,model%FV_Atm(1)%npz))
allocate(model%sbufferx(model%FV_Atm(1)%bd%isd:model%FV_Atm(1)%bd%ied,model%FV_Atm(1)%npz))

!Set ptop, ak, bk to be same as geometry
model%FV_Atm(1)%ak = geom%ak
model%FV_Atm(1)%bk = geom%bk
model%FV_Atm(1)%ptop = geom%ptop

!Always allocate w, delz, q_con for now
deallocate(model%FV_Atm(1)%w)
deallocate(model%FV_Atm(1)%delz)
deallocate(model%FV_Atm(1)%q_con)
allocate  ( model%FV_Atm(1)%w (model%FV_Atm(1)%bd%isd:model%FV_Atm(1)%bd%ied,model%FV_Atm(1)%bd%jsd:model%FV_Atm(1)%bd%jed,&
                               model%FV_Atm(1)%flagstruct%npz) )
allocate  ( model%FV_Atm(1)%delz (model%FV_Atm(1)%bd%isd:model%FV_Atm(1)%bd%ied,model%FV_Atm(1)%bd%jsd:model%FV_Atm(1)%bd%jed,&
                               model%FV_Atm(1)%flagstruct%npz) )
allocate  ( model%FV_Atm(1)%q_con(model%FV_Atm(1)%bd%isd:model%FV_Atm(1)%bd%ied,model%FV_Atm(1)%bd%jsd:model%FV_Atm(1)%bd%jed,&
                               model%FV_Atm(1)%flagstruct%npz) )

!fC and f0
if (model%FV_Atm(1)%flagstruct%grid_type == 4) then
   model%FV_Atm(1)%gridstruct%fC(:,:) = 2.*omega*sin(model%FV_Atm(1)%flagstruct%deglat/180.*pi)
   model%FV_Atm(1)%gridstruct%f0(:,:) = 2.*omega*sin(model%FV_Atm(1)%flagstruct%deglat/180.*pi)
else
   if (f_coriolis_angle == -999) then
      model%FV_Atm(1)%gridstruct%fC(:,:) = 0.0
      model%FV_Atm(1)%gridstruct%f0(:,:) = 0.0
   else
      do j=model%FV_Atm(1)%bd%jsd,model%FV_Atm(1)%bd%jed+1
         do i=model%FV_Atm(1)%bd%isd,model%FV_Atm(1)%bd%ied+1
            model%FV_Atm(1)%gridstruct%fC(i,j) = 2.*omega*( -COS(model%FV_Atm(1)%gridstruct%grid(i,j,1))*&
                                           COS(model%FV_Atm(1)%gridstruct%grid(i,j,2))*SIN(f_coriolis_angle) + &
                                           SIN(model%FV_Atm(1)%gridstruct%grid(i,j,2))*COS(f_coriolis_angle) )
         enddo
      enddo
      do j=model%FV_Atm(1)%bd%jsd,model%FV_Atm(1)%bd%jed
         do i=model%FV_Atm(1)%bd%isd,model%FV_Atm(1)%bd%ied
            model%FV_Atm(1)%gridstruct%f0(i,j) = 2.*omega*( -COS(model%FV_Atm(1)%gridstruct%agrid(i,j,1))*&
                                           COS(model%FV_Atm(1)%gridstruct%agrid(i,j,2))*SIN(f_coriolis_angle) + &
                                           SIN(model%FV_Atm(1)%gridstruct%agrid(i,j,2))*COS(f_coriolis_angle) )
         enddo
      enddo
   endif
endif

if (config_element_exists(c_conf,"update_dgridwind")) model%update_dgridwind = config_get_int(c_conf,"update_dgridwind")
if (config_element_exists(c_conf,"update_pressures")) model%update_pressures = config_get_int(c_conf,"update_pressures")

!Pointer to self when not nested
if (.not. model%FV_Atm(1)%gridstruct%nested) model%FV_Atm(1)%parent_grid => model%FV_Atm(1)

#ifdef TLADPRES

model%init_tlmadm = 0
if (config_element_exists(c_conf,"init_tlmadm")) &
  model%init_tlmadm = config_get_int(c_conf,"init_tlmadm")

if (model%init_tlmadm == 1) then

   !Initialize perturbation variables and read config
   call fv_init_pert(model%FV_Atm,model%FV_AtmP)
      
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
      model%cp_dyn_ind = config_get_int (c_conf,"CP_dyn_ind")
      cp_iter(model%cp_dyn_ind)%my_name(1:3) = 'dyn'
      
      cp_iter(model%cp_dyn_ind)%cp_test = .false.
      tmp = config_get_int (c_conf,"CP_dyn_test")
      if (tmp==1) cp_iter(model%cp_dyn_ind)%cp_test = .true.
      
      cp_iter(model%cp_dyn_ind)%cp_rep = .false.
      tmp = config_get_int (c_conf,"CP_dyn_rep")
      if (tmp==1) cp_iter(model%cp_dyn_ind)%cp_test = .true.
      
      !Hardwire these for now
      cp_iter(model%cp_dyn_ind)%check_st_control = .false.
      cp_iter(model%cp_dyn_ind)%check_st_integer = .false.
      cp_iter(model%cp_dyn_ind)%check_st_real_r4 = .false.
      cp_iter(model%cp_dyn_ind)%check_st_real_r8 = .false.
      
      cp_iter(model%cp_dyn_ind)%test_dim_st_control = 0
      cp_iter(model%cp_dyn_ind)%test_dim_st_integer = 0
      cp_iter(model%cp_dyn_ind)%test_dim_st_real_r4 = 0
      cp_iter(model%cp_dyn_ind)%test_dim_st_real_r8 = 0
      
      cp_iter(model%cp_dyn_ind)%test_dim_cp_control = 0
      cp_iter(model%cp_dyn_ind)%test_dim_cp_integer = 0
      cp_iter(model%cp_dyn_ind)%test_dim_cp_real_r4 = 0
      cp_iter(model%cp_dyn_ind)%test_dim_cp_real_r8 = 0
   
   endif

endif

#endif

end subroutine model_setup

! ------------------------------------------------------------------------------

subroutine model_delete(self)

use fv_arrays_mod, only: deallocate_fv_atmos_type

#ifdef TLADPRES
use fv_arrays_nlm_mod, only: deallocate_fv_atmos_pert_type
use tapenade_iter, only: cp_iter_controls, finalize_cp_iter
#endif

implicit none
type(fv3jedi_model) :: self

deallocate(self%ebuffery)
deallocate(self%wbuffery)
deallocate(self%nbufferx)
deallocate(self%sbufferx)

call deallocate_fv_atmos_type(self%FV_Atm(1))
deallocate(self%FV_Atm)

#ifdef TLADPRES
if (self%init_tlmadm == 1) then
   call deallocate_fv_atmos_pert_type(self%FV_AtmP(1))
   deallocate(self%FV_AtmP)
   
   if (cp_iter_controls%cp_i .ne. 0) call finalize_cp_iter
endif
#endif

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_prepare_integration(self, flds)

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds

!Get phis from fields, fixed for integration
!-------------------------------------------
call get_phi_from_flds(flds,self)

end subroutine model_prepare_integration

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_ad(self, flds)

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds

!Get phis from fields, fixed for integration
!-------------------------------------------
call get_phi_from_flds(flds,self)

end subroutine model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine model_prepare_integration_tl(self, flds)

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds

!Get phis from fields, fixed for integration
!-------------------------------------------
call get_phi_from_flds(flds,self)

end subroutine model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine model_propagate(self, flds)

use fv_dynamics_mod, only: fv_dynamics
use fv_sg_mod, only: fv_subgrid_z

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds

type(fv_atmos_type), pointer :: FV_Atm(:)
integer :: i,j,k

real(kind=kind_real), allocatable, dimension(:,:,:) :: u_dt, v_dt, t_dt

if (mpp_pe() == mpp_root_pe()) print*, 'Propagate nonlinear model'


!Convenience pointer to the main FV_Atm structure
!------------------------------------------------
FV_Atm => self%FV_Atm


!Copy to model variables
!-----------------------
call fields_to_model(flds,self)


!Update edges of d-grid winds
!----------------------------
if (self%update_dgridwind == 1) then

   !Zero the repeat edge points
   FV_Atm(1)%u(:,self%jec+1,:) = 0.0
   FV_Atm(1)%v(self%iec+1,:,:) = 0.0

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

endif

!Compute the other pressure variables needed by FV3
!--------------------------------------------------
if (self%update_pressures == 1) then
   call compute_fv3_pressures( self%isc, self%iec, self%jsc, self%jec, self%isd, self%ied, self%jsd, self%jed, &
                               self%npz, kappa, FV_Atm(1)%ptop, &
                               FV_Atm(1)%delp, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%pkz, FV_Atm(1)%peln )
endif

call set_domain(FV_Atm(1)%domain)

!Propagate FV3 one time step
!---------------------------
call fv_dynamics( FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,  &
                  self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,           &
                  FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                   &
                  cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,          &
                  FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                  &
                  FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,                       &
                  FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, &
                  FV_Atm(1)%ps, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,     &
                  FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                             &
                  FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,                      &
                  FV_Atm(1)%ak, FV_Atm(1)%bk,                                                  &
                  FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,     &
                  FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,   &
                  FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,  &
                  FV_Atm(1)%domain )

!Apply subgrid mixing (interface conflict between GMAO/NCEP so not straighforward to compile right now)
!--------------------
!if ( FV_Atm(1)%flagstruct%fv_sg_adj > 0 ) then
!     allocate ( u_dt(self%isd:self%ied,self%jsd:self%jed,self%npz) )
!     allocate ( v_dt(self%isd:self%ied,self%jsd:self%jed,self%npz) )
!     allocate ( t_dt(self%isc:self%iec,self%jsc:self%jec,self%npz) )
!     u_dt(:,:,:) = 0.0
!     v_dt(:,:,:) = 0.0
!     t_dt(:,:,:) = 0.0
!     call fv_subgrid_z(self%isd, self%ied, self%jsd, self%jed, self%isc, self%iec, self%jsc, self%jec, FV_Atm(1)%npz, &
!                       FV_Atm(1)%ncnst, self%DT, FV_Atm(1)%flagstruct%fv_sg_adj,      &
!                       FV_Atm(1)%flagstruct%nwat, FV_Atm(1)%delp, FV_Atm(1)%pe,     &
!                       FV_Atm(1)%peln, FV_Atm(1)%pkz, FV_Atm(1)%pt, FV_Atm(1)%q,       &
!                       FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%flagstruct%hydrostatic,&
!                       FV_Atm(1)%w, FV_Atm(1)%delz, u_dt, v_dt, t_dt, FV_Atm(1)%flagstruct%n_zfilter)
!     deallocate ( u_dt )
!     deallocate ( v_dt )
!     deallocate ( t_dt )
!endif

call nullify_domain()

!Copy back to fields
!-------------------
call model_to_fields(self,flds)

end subroutine model_propagate

! ------------------------------------------------------------------------------

subroutine model_propagate_ad(self, flds, traj)

#ifdef TLADPRES
use fv_dynamics_adm_mod, only: fv_dynamics_fwd, fv_dynamics_bwd
use mpp_domains_mod, only: mpp_get_boundary_ad
use tapenade_iter,   only: cp_iter_controls, cp_mod_ini, cp_mod_mid, cp_mod_end, pushrealarray, poprealarray
#endif

implicit none

type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds
type(fv3jedi_trajectory)    :: traj

#ifdef TLADPRES
type(fv_atmos_type), pointer :: FV_Atm(:)
type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
integer :: i,j,k


if (mpp_pe() == mpp_root_pe()) print*, 'Propagate model adjoint'

! Convenience pointer to the main FV_Atm structure
! ------------------------------------------------
FV_Atm => self%FV_Atm
FV_AtmP => self%FV_AtmP


! Get up the trajectory for this time step 
! ----------------------------------------
call get_traj( traj,flds%geom%bd%isd,flds%geom%bd%ied,flds%geom%bd%jsd,flds%geom%bd%jed,&
               flds%geom%npz,flds%geom%hydrostatic,FV_Atm(1)%ncnst,&
               FV_Atm(1)%u,FV_Atm(1)%v,FV_Atm(1)%pt,FV_Atm(1)%delp,FV_Atm(1)%q,FV_Atm(1)%w,FV_Atm(1)%delz)


! Update edges of d-grid winds
! ----------------------------
if (self%update_dgridwind == 1) then

   !Zero the repeat edge points
   FV_Atm(1)%u(:,self%jec+1,:) = 0.0
   FV_Atm(1)%v(self%iec+1,:,:) = 0.0

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

endif


! Copy to model variables
! -----------------------
call fields_to_model_pert(flds,self)

!Compute the other pressure variables needed by FV3
!--------------------------------------------------
if (self%update_pressures == 1) then
   call compute_fv3_pressures( self%isc, self%iec, self%jsc, self%jec, self%isd, self%ied, self%jsd, self%jed, &
                               self%npz, kappa, FV_Atm(1)%ptop, &
                               FV_Atm(1)%delp, FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%pkz, FV_Atm(1)%peln )
endif

call set_domain(FV_Atm(1)%domain)

! Initilize the module level checkpointing
! ----------------------------------------
if (cp_iter_controls%cp_i .ne. 0) then
   call cp_mod_ini(self%cp_dyn_ind)
endif


! Forward sweep of the dynamics with saving of checkpoints for use in backward sweep
! ----------------------------------------------------------------------------------
if (cp_iter_controls%cp_i <= 3) then

!   call fv_dynamics_fwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
!                        self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
!                        FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
!                        cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
!                        FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
!                        FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w,                                                           &
!                        FV_Atm(1)%delz, FV_Atm(1)%flagstruct%hydrostatic,                                                &
!                        FV_Atm(1)%pt, FV_Atm(1)%delp,                                                                    &
!                        FV_Atm(1)%q, FV_Atm(1)%ps, FV_Atm(1)%pe,                                                         &
!                        FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                                                     &
!                        FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga,                                                 &
!                        FV_Atm(1)%ua, FV_Atm(1)%va,                                                                      &
!                        FV_Atm(1)%uc, FV_Atm(1)%vc,                                                                      &
!                        FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
!                        FV_Atm(1)%mfx, FV_Atm(1)%mfy,                                                                    &
!                        FV_Atm(1)%cx, FV_Atm(1)%cy, FV_Atm(1)%ze0,                                                       &
!                        FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
!                        FV_AtmP(1)%flagstruct,                                                                           &
!                        FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )

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

call zero_pert_vars(FV_AtmP(1))

! Backward adjoint sweep of the dynamics
! --------------------------------------
!call fv_dynamics_bwd(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
!                     self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
!                     FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
!                     cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
!                     FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
!                     FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
!                     FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
!                     FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
!                     FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
!                     FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
!                     FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
!                     FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
!                     FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
!                     FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
!                     FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
!                     FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                      &
!                     FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct,                       &
!                     FV_AtmP(1)%flagstruct,                                                                           &
!                     FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )


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

call zero_pert_vars(FV_AtmP(1))

! Finalize iterative step and get ready for next iteration
! --------------------------------------------------------
if (cp_iter_controls%cp_i .ne. 0) then
   call cp_mod_end
endif


! Update edges of d-grid winds
! ----------------------------
if (self%update_dgridwind == 1) then

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

   !Zero the repeat edge points
   !FV_AtmP(1)%up(:,self%jec+1,:) = 0.0 !Might break things right now.
   !FV_AtmP(1)%vp(self%iec+1,:,:) = 0.0

endif

call nullify_domain()

! Copy back to fields 
! -------------------
call model_to_fields_pert(self,flds)

#endif

end subroutine model_propagate_ad

! ------------------------------------------------------------------------------

subroutine model_propagate_tl(self, flds, traj)

#ifdef TLADPRES
use fv_dynamics_tlm_mod, only: fv_dynamics_tlm
#endif

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_field)         :: flds
type(fv3jedi_trajectory)    :: traj

#ifdef TLADPRES
type(fv_atmos_type), pointer :: FV_Atm(:)
type(fv_atmos_pert_type), pointer :: FV_AtmP(:)
integer :: i,j,k


if (mpp_pe() == mpp_root_pe()) print*, 'Propagate tangent linear model'


!Convenience pointer to the main FV_Atm structure
!------------------------------------------------
FV_Atm => self%FV_Atm
FV_AtmP => self%FV_AtmP


!Get up the trajectory for this time step 
!----------------------------------------
call get_traj( traj,flds%geom%bd%isd,flds%geom%bd%ied,flds%geom%bd%jsd,flds%geom%bd%jed,&
               flds%geom%npz,flds%geom%hydrostatic,FV_Atm(1)%ncnst,&
               FV_Atm(1)%u,FV_Atm(1)%v,FV_Atm(1)%pt,FV_Atm(1)%delp,FV_Atm(1)%q,FV_Atm(1)%w,FV_Atm(1)%delz)


!Copy to model variables
!-----------------------
call fields_to_model_pert(flds,self)


!Update edges of d-grid winds, probably not needed ultimately
!------------------------------------------------------------
if (self%update_dgridwind == 1) then

   !Zero the repeat edge points
   FV_Atm(1)%u(:,self%jec+1,:) = 0.0
   FV_Atm(1)%v(self%iec+1,:,:) = 0.0

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
   
   !Zero the repeat edge points
   FV_AtmP(1)%up(:,self%jec+1,:) = 0.0
   FV_AtmP(1)%vp(self%iec+1,:,:) = 0.0

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

endif

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

call set_domain(FV_Atm(1)%domain)

call zero_pert_vars(FV_AtmP(1))

!Propagate TLM one time step
!---------------------------
!call fv_dynamics_tlm(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,                      &
!                     self%DT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill,                               &
!                     FV_Atm(1)%flagstruct%reproduce_sum, kappa,                                                       &
!                     cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst,                              &
!                     FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split,                                      &
!                     FV_Atm(1)%u, FV_AtmP(1)%up, FV_Atm(1)%v, FV_AtmP(1)%vp, FV_Atm(1)%w, FV_AtmP(1)%wp,                 &
!                     FV_Atm(1)%delz, FV_AtmP(1)%delzp, FV_Atm(1)%flagstruct%hydrostatic,                               &
!                     FV_Atm(1)%pt, FV_AtmP(1)%ptp, FV_Atm(1)%delp, FV_AtmP(1)%delpp,                                    &
!                     FV_Atm(1)%q, FV_AtmP(1)%qp, FV_Atm(1)%ps, FV_AtmP(1)%psp, FV_Atm(1)%pe, FV_AtmP(1)%pep,             &
!                     FV_Atm(1)%pk, FV_AtmP(1)%pkp, FV_Atm(1)%peln, FV_AtmP(1)%pelnp, FV_Atm(1)%pkz, FV_AtmP(1)%pkzp,     &
!                     FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_AtmP(1)%omgap,                                &
!                     FV_Atm(1)%ua, FV_AtmP(1)%uap, FV_Atm(1)%va, FV_AtmP(1)%vap,                                        &
!                     FV_Atm(1)%uc, FV_AtmP(1)%ucp, FV_Atm(1)%vc, FV_AtmP(1)%vcp,                                        &
!                     FV_Atm(1)%ak, FV_Atm(1)%bk,                                                                      &
!                     FV_Atm(1)%mfx, FV_AtmP(1)%mfxp, FV_Atm(1)%mfy, FV_AtmP(1)%mfyp,                                    &
!                     FV_Atm(1)%cx, FV_AtmP(1)%cxp, FV_Atm(1)%cy,  FV_AtmP(1)%cyp, FV_Atm(1)%ze0,                       &
!                     FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, FV_AtmP(1)%flagstruct, &
!                     FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid,FV_Atm(1)%domain      )
!
call zero_pert_vars(FV_AtmP(1))

call nullify_domain()

!Copy back to fields
!-------------------
call model_to_fields_pert(self,flds)

#endif

end subroutine model_propagate_tl

! ------------------------------------------------------------------------------

subroutine model_prop_traj(self, flds, traj)
implicit none
type(fv3jedi_model)      :: self
type(fv3jedi_field)      :: flds
type(fv3jedi_trajectory) :: traj

call fields_to_model(flds,self)

call set_traj( traj,flds%geom%bd%isd,flds%geom%bd%ied,flds%geom%bd%jsd,flds%geom%bd%jed,flds%geom%npz, &
               self%FV_Atm(1)%ncnst,flds%geom%hydrostatic, &
               self%FV_Atm(1)%u,self%FV_Atm(1)%v,self%FV_Atm(1)%pt,self%FV_Atm(1)%delp,self%FV_Atm(1)%q,self%FV_Atm(1)%w,self%FV_Atm(1)%delz)

end subroutine model_prop_traj

! ------------------------------------------------------------------------------

subroutine model_wipe_traj(traj)
implicit none
type(fv3jedi_trajectory), pointer :: traj

call delete_traj(traj)

end subroutine model_wipe_traj

! ------------------------------------------------------------------------------

subroutine fields_to_model(flds,self)

implicit none
type(fv3jedi_field), intent(in)    :: flds
type(fv3jedi_model), intent(inout) :: self

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

!To zero the halos
self%FV_Atm(1)%u    = 0.0
self%FV_Atm(1)%v    = 0.0
self%FV_Atm(1)%pt   = 0.0
self%FV_Atm(1)%delp = 0.0
self%FV_Atm(1)%q    = 0.0
if (.not. flds%geom%hydrostatic) then
   self%FV_Atm(1)%w   = 0.0
   self%FV_Atm(1)%delz = 0.0
endif

!Only copy compute grid incase model versus geometry halos are different
self%FV_Atm(1)%u   (isc:iec  ,jsc:jec+1,:  ) = flds%Atm%u   (isc:iec  ,jsc:jec+1,:  )
self%FV_Atm(1)%v   (isc:iec+1,jsc:jec  ,:  ) = flds%Atm%v   (isc:iec+1,jsc:jec  ,:  )
self%FV_Atm(1)%pt  (isc:iec  ,jsc:jec  ,:  ) = flds%Atm%pt  (isc:iec  ,jsc:jec  ,:  )
self%FV_Atm(1)%delp(isc:iec  ,jsc:jec  ,:  ) = flds%Atm%delp(isc:iec  ,jsc:jec  ,:  )
self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,:) = flds%Atm%q   (isc:iec  ,jsc:jec  ,:,:)
if (.not. flds%geom%hydrostatic) then
   self%FV_Atm(1)%delz(isc:iec  ,jsc:jec  ,:  ) = flds%Atm%delz(isc:iec  ,jsc:jec  ,:  )
   self%FV_Atm(1)%w   (isc:iec  ,jsc:jec  ,:  ) = flds%Atm%w   (isc:iec  ,jsc:jec  ,:  )
endif

end subroutine fields_to_model

! ------------------------------------------------------------------------------

subroutine model_to_fields(self,flds)

implicit none
type(fv3jedi_model), intent(in   ) :: self
type(fv3jedi_field), intent(inout) :: flds

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

!To zero the halos
flds%Atm%u    = 0.0
flds%Atm%v    = 0.0
flds%Atm%pt   = 0.0
flds%Atm%delp = 0.0
flds%Atm%q    = 0.0
if (.not. flds%geom%hydrostatic) then
   flds%Atm%delz = 0.0
   flds%Atm%w    = 0.0
endif

!Only copy compute grid incase model versus geometry halos are different
flds%Atm%u   (isc:iec  ,jsc:jec+1,:  ) = self%FV_Atm(1)%u   (isc:iec  ,jsc:jec+1,:  )
flds%Atm%v   (isc:iec+1,jsc:jec  ,:  ) = self%FV_Atm(1)%v   (isc:iec+1,jsc:jec  ,:  )
flds%Atm%pt  (isc:iec  ,jsc:jec  ,:  ) = self%FV_Atm(1)%pt  (isc:iec  ,jsc:jec  ,:  )
flds%Atm%delp(isc:iec  ,jsc:jec  ,:  ) = self%FV_Atm(1)%delp(isc:iec  ,jsc:jec  ,:  )
flds%Atm%q   (isc:iec  ,jsc:jec  ,:,:) = self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,:)
if (.not. flds%geom%hydrostatic) then
  flds%Atm%delz(isc:iec  ,jsc:jec  ,:  ) = self%FV_Atm(1)%delz(isc:iec  ,jsc:jec  ,:  )
  flds%Atm%w   (isc:iec  ,jsc:jec  ,:  ) = self%FV_Atm(1)%w   (isc:iec  ,jsc:jec  ,:  )
endif

end subroutine model_to_fields

! ------------------------------------------------------------------------------

subroutine fields_to_model_pert(flds,self)

implicit none
type(fv3jedi_field), intent(in)    :: flds
type(fv3jedi_model), intent(inout) :: self

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

!To zero the halos
self%FV_AtmP(1)%up    = 0.0
self%FV_AtmP(1)%vp    = 0.0
self%FV_AtmP(1)%ptp   = 0.0
self%FV_AtmP(1)%delpp = 0.0
self%FV_AtmP(1)%qp    = 0.0
if (.not. flds%geom%hydrostatic) then
   self%FV_AtmP(1)%wp   = 0.0
   self%FV_AtmP(1)%delzp = 0.0
endif

!Only copy compute grid incase model versus geometry halos are different
self%FV_AtmP(1)%up   (isc:iec  ,jsc:jec+1,:  ) = flds%Atm%u   (isc:iec  ,jsc:jec+1,:  )
self%FV_AtmP(1)%vp   (isc:iec+1,jsc:jec  ,:  ) = flds%Atm%v   (isc:iec+1,jsc:jec  ,:  )
self%FV_AtmP(1)%ptp  (isc:iec  ,jsc:jec  ,:  ) = flds%Atm%pt  (isc:iec  ,jsc:jec  ,:  )
self%FV_AtmP(1)%delpp(isc:iec  ,jsc:jec  ,:  ) = flds%Atm%delp(isc:iec  ,jsc:jec  ,:  )
self%FV_AtmP(1)%qp   (isc:iec  ,jsc:jec  ,:,:) = flds%Atm%q   (isc:iec  ,jsc:jec  ,:,:)
if (.not. flds%geom%hydrostatic) then
   self%FV_AtmP(1)%delzp(isc:iec  ,jsc:jec  ,:  ) = flds%Atm%delz(isc:iec  ,jsc:jec  ,:  )
   self%FV_AtmP(1)%wp   (isc:iec  ,jsc:jec  ,:  ) = flds%Atm%w   (isc:iec  ,jsc:jec  ,:  )
endif

end subroutine fields_to_model_pert

! ------------------------------------------------------------------------------

subroutine model_to_fields_pert(self,flds)

implicit none
type(fv3jedi_model), intent(in   ) :: self
type(fv3jedi_field), intent(inout) :: flds

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

!To zero the halos
flds%Atm%u    = 0.0
flds%Atm%v    = 0.0
flds%Atm%pt   = 0.0
flds%Atm%delp = 0.0
flds%Atm%q    = 0.0
if (.not. flds%geom%hydrostatic) then
   flds%Atm%delz = 0.0
   flds%Atm%w    = 0.0
endif

!Only copy compute grid incase model versus geometry halos are different
flds%Atm%u   (isc:iec  ,jsc:jec+1,:  ) = self%FV_AtmP(1)%up   (isc:iec  ,jsc:jec+1,:  )
flds%Atm%v   (isc:iec+1,jsc:jec  ,:  ) = self%FV_AtmP(1)%vp   (isc:iec+1,jsc:jec  ,:  )
flds%Atm%pt  (isc:iec  ,jsc:jec  ,:  ) = self%FV_AtmP(1)%ptp  (isc:iec  ,jsc:jec  ,:  )
flds%Atm%delp(isc:iec  ,jsc:jec  ,:  ) = self%FV_AtmP(1)%delpp(isc:iec  ,jsc:jec  ,:  )
flds%Atm%q   (isc:iec  ,jsc:jec  ,:,:) = self%FV_AtmP(1)%qp   (isc:iec  ,jsc:jec  ,:,:)
if (.not. flds%geom%hydrostatic) then
  flds%Atm%delz(isc:iec  ,jsc:jec  ,:  ) = self%FV_AtmP(1)%delzp(isc:iec  ,jsc:jec  ,:  )
  flds%Atm%w   (isc:iec  ,jsc:jec  ,:  ) = self%FV_AtmP(1)%wp   (isc:iec  ,jsc:jec  ,:  )
endif

end subroutine model_to_fields_pert

! ------------------------------------------------------------------------------

subroutine get_phi_from_flds(flds, self)

implicit none
type(fv3jedi_field) :: flds
type(fv3jedi_model) :: self

!To zero halo
self%FV_Atm(1)%phis = 0.0

!Get compute domain from fields
self%FV_Atm(1)%phis(self%isc:self%iec,self%jsc:self%jec) = flds%Atm%phis(self%isc:self%iec,self%jsc:self%jec)

!Fill halos
call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)

end subroutine get_phi_from_flds

! ------------------------------------------------------------------------------

subroutine zero_pert_vars(FV_AtmP)

implicit none
type(fv_atmos_pert_type), intent(inout) :: FV_AtmP

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

end module fv3jedi_model_mod
