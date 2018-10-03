! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_model_mod

use iso_c_binding
use config_mod
use duration_mod
use fv3jedi_geom_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 
use fv3jedi_constants
use kinds

use fv_arrays_mod,  only: fv_atmos_type
use mpp_mod,        only: mpp_pe, mpp_root_pe 
use mpp_domains_mod, only: mpp_update_domains, mpp_get_boundary, DGRID_NE
use pressure_vt_mod, only: compute_fv3_pressures, compute_fv3_pressures_tlm, compute_fv3_pressures_bwd
use fms_mod, only: set_domain, nullify_domain

use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index

implicit none
private

public :: fv3jedi_model
public :: model_create
public :: model_delete
public :: model_initialize
public :: model_step
public :: model_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: fv3jedi_model
  real(kind=kind_real)                         :: DT                  !<Model big timestep
  real(kind_real), allocatable, dimension(:,:) :: ebuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: nbufferx            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: wbuffery            !<Halo holder
  real(kind_real), allocatable, dimension(:,:) :: sbufferx            !<Halo holder
  type(fv_atmos_type), allocatable             :: FV_Atm(:)           !<Main FV3 construct 
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
end type fv3jedi_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine model_create(model, geom, c_conf)

use fv_control_mod, only: fv_init, pelist_all

implicit none
type(c_ptr), intent(in)     :: c_conf !< pointer to object of class Config
type(fv3jedi_model), target :: model  ! should I put intent on these?
type(fv3jedi_geom)          :: geom

character(len=20) :: ststep
type(duration) :: dtstep

integer :: i,j
integer :: tmp

!Model time step
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
model%DT = real(duration_seconds(dtstep),kind_real)

!Call to fv_init
call fv_init(model%FV_Atm, model%DT, model%grids_on_this_pe, model%p_split)
deallocate(pelist_all)

!Compute grid must be same as geometry
if ( (geom%isc .ne. model%FV_Atm(1)%bd%isc) .or. (geom%iec .ne. model%FV_Atm(1)%bd%iec) .or. &
     (geom%jsc .ne. model%FV_Atm(1)%bd%jsc) .or. (geom%jec .ne. model%FV_Atm(1)%bd%jec) .or. &
     (geom%npz .ne. model%FV_Atm(1)%npz) ) then
   call abor1_ftn("fv3jedi model: compute areas for geometry and model do not agree")
endif

!Copy of grid info for convenience
model%isc = model%FV_Atm(1)%bd%isc
model%iec = model%FV_Atm(1)%bd%iec
model%jsc = model%FV_Atm(1)%bd%jsc
model%jec = model%FV_Atm(1)%bd%jec
model%isd = model%FV_Atm(1)%bd%isd
model%ied = model%FV_Atm(1)%bd%ied
model%jsd = model%FV_Atm(1)%bd%jsd
model%jed = model%FV_Atm(1)%bd%jed
model%npz = model%FV_Atm(1)%npz
model%hydrostatic = model%FV_Atm(1)%flagstruct%hydrostatic

!Halo holders for domain grid
allocate(model%wbuffery(model%FV_Atm(1)%bd%jsc:model%FV_Atm(1)%bd%jec,model%FV_Atm(1)%npz))
allocate(model%sbufferx(model%FV_Atm(1)%bd%isc:model%FV_Atm(1)%bd%iec,model%FV_Atm(1)%npz))
allocate(model%ebuffery(model%FV_Atm(1)%bd%jsc:model%FV_Atm(1)%bd%jec,model%FV_Atm(1)%npz))
allocate(model%nbufferx(model%FV_Atm(1)%bd%isc:model%FV_Atm(1)%bd%iec,model%FV_Atm(1)%npz))

!Set ptop, ak, bk to be same as geometry
model%FV_Atm(1)%ak = geom%ak
model%FV_Atm(1)%bk = geom%bk
model%FV_Atm(1)%ptop = geom%ptop

!Tracer indexes
model%ti_q  = 1 !get_tracer_index (MODEL_ATMOS, 'sphum')
model%ti_ql = 2 !get_tracer_index (MODEL_ATMOS, 'liq_wat')
model%ti_qi = 3 !get_tracer_index (MODEL_ATMOS, 'ice_wat')
model%ti_o3 = 4 !get_tracer_index (MODEL_ATMOS, 'o3mr')

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
model%FV_Atm(1)%w = 0.0
model%FV_Atm(1)%delz = 0.0
model%FV_Atm(1)%q_con = 0.0

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

!Harwire some flags
model%FV_Atm(1)%flagstruct%reproduce_sum = .false.
model%FV_Atm(1)%flagstruct%fill = .false.
model%FV_Atm(1)%flagstruct%fv_debug = .false.
model%FV_Atm(1)%flagstruct%adiabatic = .false.
model%FV_Atm(1)%flagstruct%do_sat_adj = .false.
model%FV_Atm(1)%flagstruct%breed_vortex_inline = .false.

end subroutine model_create

! ------------------------------------------------------------------------------

subroutine model_delete(self)

use fv_arrays_mod, only: deallocate_fv_atmos_type

implicit none
type(fv3jedi_model) :: self

deallocate(self%ebuffery)
deallocate(self%wbuffery)
deallocate(self%nbufferx)
deallocate(self%sbufferx)

call deallocate_fv_atmos_type(self%FV_Atm(1))
deallocate(self%FV_Atm)

end subroutine model_delete

! ------------------------------------------------------------------------------

subroutine model_initialize(self, state)

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_state)         :: state

end subroutine model_initialize

! ------------------------------------------------------------------------------

subroutine model_step(geom, self, state)

use fv_dynamics_mod, only: fv_dynamics
use fv_sg_mod, only: fv_subgrid_z

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_state)         :: state
type(fv3jedi_geom)          :: geom

type(fv_atmos_type), pointer :: FV_Atm(:)
integer :: i,j,k


if (mpp_pe() == mpp_root_pe()) print*, 'Propagate nonlinear model'


!Convenience pointer to the main FV_Atm structure
!------------------------------------------------
FV_Atm => self%FV_Atm


!Copy to model variables
!-----------------------
call state_to_model(state,self)


!Get phis from state, fixed for integration
!-------------------------------------------
call get_phi_from_state(state,self)


! Zero local variables
! --------------------
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
if (self%update_dgridwind == 1) then

   call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, geom%domain, &
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


! MPP set domain
! --------------
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

! MPP nulify
! ----------
call nullify_domain()


!Copy back to state
!------------------
call model_to_state(self,state)


end subroutine model_step

! ------------------------------------------------------------------------------

subroutine model_finalize(self, state)

implicit none
type(fv3jedi_model), target :: self
type(fv3jedi_state)         :: state

end subroutine model_finalize

! ------------------------------------------------------------------------------

subroutine state_to_model(state,self)

implicit none
type(fv3jedi_state), intent(in)    :: state
type(fv3jedi_model), intent(inout) :: self

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

!NOTE: while the variable name is pt, FV3 expects dry temperature

!To zero the halos
self%FV_Atm(1)%u    = 0.0
self%FV_Atm(1)%v    = 0.0
self%FV_Atm(1)%pt   = 0.0
self%FV_Atm(1)%delp = 0.0
self%FV_Atm(1)%q    = 0.0
self%FV_Atm(1)%w    = 0.0
self%FV_Atm(1)%delz = 0.0

!Only copy compute grid incase of halo differences
self%FV_Atm(1)%u   (isc:iec  ,jsc:jec+1,:)            = state%ud  (isc:iec  ,jsc:jec+1,:)
self%FV_Atm(1)%v   (isc:iec+1,jsc:jec  ,:)            = state%vd  (isc:iec+1,jsc:jec  ,:)
self%FV_Atm(1)%pt  (isc:iec  ,jsc:jec  ,:)            = state%t   (isc:iec  ,jsc:jec  ,:)
self%FV_Atm(1)%delp(isc:iec  ,jsc:jec  ,:)            = state%delp(isc:iec  ,jsc:jec  ,:)
self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_q ) = state%q   (isc:iec  ,jsc:jec  ,:)
self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_qi) = state%qi  (isc:iec  ,jsc:jec  ,:)
self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_ql) = state%ql  (isc:iec  ,jsc:jec  ,:)
self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_o3) = state%o3  (isc:iec  ,jsc:jec  ,:)
if (.not. state%hydrostatic) then
   self%FV_Atm(1)%delz(isc:iec,jsc:jec,:) = state%delz(isc:iec,jsc:jec,:)
   self%FV_Atm(1)%w   (isc:iec,jsc:jec,:) = state%w   (isc:iec,jsc:jec,:)
endif

end subroutine state_to_model

! ------------------------------------------------------------------------------

subroutine model_to_state(self,state)

implicit none
type(fv3jedi_model), intent(in   ) :: self
type(fv3jedi_state), intent(inout) :: state

integer :: isc,iec,jsc,jec

isc = self%FV_Atm(1)%bd%isc
iec = self%FV_Atm(1)%bd%iec
jsc = self%FV_Atm(1)%bd%jsc
jec = self%FV_Atm(1)%bd%jec

state%ud   = 0.0
state%vd   = 0.0
state%t    = 0.0
state%delp = 0.0
state%q    = 0.0
if (.not. state%hydrostatic) then
   state%delz = 0.0
   state%w    = 0.0
endif
state%ua    = 0.0
state%va    = 0.0

state%ud  (isc:iec  ,jsc:jec+1,:) = self%FV_Atm(1)%u   (isc:iec  ,jsc:jec+1,:  )
state%vd  (isc:iec+1,jsc:jec  ,:) = self%FV_Atm(1)%v   (isc:iec+1,jsc:jec  ,:  )
state%t   (isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%pt  (isc:iec  ,jsc:jec  ,:  )
state%delp(isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%delp(isc:iec  ,jsc:jec  ,:  )
state%q   (isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_q)
state%qi  (isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_qi)
state%ql  (isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_ql)
state%o3  (isc:iec  ,jsc:jec  ,:) = self%FV_Atm(1)%q   (isc:iec  ,jsc:jec  ,:,self%ti_o3)
if (.not. state%hydrostatic) then
  state%delz(isc:iec,jsc:jec,:) = self%FV_Atm(1)%delz(isc:iec,jsc:jec,:)
  state%w   (isc:iec,jsc:jec,:) = self%FV_Atm(1)%w   (isc:iec,jsc:jec,:)
endif
state%ua(isc:iec,jsc:jec,:) = self%FV_Atm(1)%ua(isc:iec,jsc:jec,:)
state%va(isc:iec,jsc:jec,:) = self%FV_Atm(1)%va(isc:iec,jsc:jec,:)

end subroutine model_to_state

! ------------------------------------------------------------------------------

subroutine get_phi_from_state(state, self)

implicit none
type(fv3jedi_state) :: state
type(fv3jedi_model) :: self

!To zero halo
self%FV_Atm(1)%phis = 0.0

!Get compute domain from state
self%FV_Atm(1)%phis(self%isc:self%iec,self%jsc:self%jec) = state%phis(self%isc:self%iec,self%jsc:self%jec)

!Fill halos
call mpp_update_domains(self%FV_Atm(1)%phis, self%FV_Atm(1)%domain, complete=.true.)

end subroutine get_phi_from_state

! ------------------------------------------------------------------------------

end module fv3jedi_model_mod
