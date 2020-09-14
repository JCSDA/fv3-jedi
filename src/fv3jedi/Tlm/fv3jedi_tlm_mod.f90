! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_tlm_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use duration_mod
use oops_variables_mod

! fv3-jedi-linearmodel
use fv3jedi_lm_mod,        only: fv3jedi_lm_type

! fv3-jedi
use fv3jedi_field_mod,     only: copy_subset, has_field, pointer_field_array
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_increment_mod, only: fv3jedi_increment, create_inc=>create
use fv3jedi_state_mod,     only: fv3jedi_state
use fv3jedi_traj_mod,      only: fv3jedi_traj

implicit none
private
public :: fv3jedi_tlm

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold tlm definition
type:: fv3jedi_tlm
  type(fv3jedi_lm_type) :: fv3jedi_lm  !<Linearized model object
  contains
    procedure :: create
    procedure :: delete
    procedure :: initialize_tl
    procedure :: initialize_ad
    procedure :: step_tl
    procedure :: step_ad
    procedure :: finalize_tl
    procedure :: finalize_ad
end type fv3jedi_tlm

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
class(fv3jedi_tlm),        intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt
character(len=:), allocatable :: str

! Model time step
! ---------------
call conf%get_or_die("tstep",str)
ststep = str
deallocate(str)
dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)

! Model configuration and creation
! --------------------------------
call conf%get_or_die("lm_do_dyn",self%fv3jedi_lm%conf%do_dyn)
call conf%get_or_die("lm_do_trb",self%fv3jedi_lm%conf%do_phy_trb)
call conf%get_or_die("lm_do_mst",self%fv3jedi_lm%conf%do_phy_mst)

call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3jedi_tlm), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize_ad(self, inc)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine initialize_ad

! --------------------------------------------------------------------------------------------------

subroutine initialize_tl(self, inc)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%init_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine initialize_tl

! --------------------------------------------------------------------------------------------------

subroutine step_ad(self, inc, traj)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%step_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine step_ad

! --------------------------------------------------------------------------------------------------

subroutine step_tl(self, inc, traj)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc
type(fv3jedi_traj),      intent(in)    :: traj

call traj_to_traj(traj,self%fv3jedi_lm)

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%step_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine step_tl

! --------------------------------------------------------------------------------------------------

subroutine finalize_ad(self, inc)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%final_ad()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine finalize_ad

! --------------------------------------------------------------------------------------------------

subroutine finalize_tl(self, inc)

implicit none
class(fv3jedi_tlm),      intent(inout) :: self
type(fv3jedi_increment), intent(inout) :: inc

call inc_to_lm(inc,self%fv3jedi_lm)
call self%fv3jedi_lm%final_tl()
call lm_to_inc(self%fv3jedi_lm,inc)

end subroutine finalize_tl

! --------------------------------------------------------------------------------------------------

subroutine inc_to_lm(inc, lm)

implicit none
type(fv3jedi_increment), intent(in)    :: inc
type(fv3jedi_lm_type),   intent(inout) :: lm

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: delp
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: qi
real(kind=kind_real), pointer, dimension(:,:,:) :: ql
real(kind=kind_real), pointer, dimension(:,:,:) :: o3
real(kind=kind_real), pointer, dimension(:,:,:) :: w
real(kind=kind_real), pointer, dimension(:,:,:) :: delz

call pointer_field_array(inc%fields, 'ud'     , ud  )
call pointer_field_array(inc%fields, 'vd'     , vd  )
call pointer_field_array(inc%fields, 't'      , t   )
call pointer_field_array(inc%fields, 'delp'   , delp)
call pointer_field_array(inc%fields, 'sphum'  , q   )
call pointer_field_array(inc%fields, 'ice_wat', qi  )
call pointer_field_array(inc%fields, 'liq_wat', ql  )
if (has_field(inc%fields, 'o3mr') ) call pointer_field_array(inc%fields, 'o3mr'   , o3  )
if (has_field(inc%fields, 'o3ppmv') ) call pointer_field_array(inc%fields, 'o3ppmv'   , o3  )
lm%pert%u    = ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
lm%pert%v    = vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)
lm%pert%ua   = 0.0_kind_real
lm%pert%va   = 0.0_kind_real
lm%pert%t    = t
lm%pert%delp = delp
lm%pert%qv   = q
lm%pert%qi   = qi
lm%pert%ql   = ql
lm%pert%o3   = o3
if (.not. inc%hydrostatic) then
  call pointer_field_array(inc%fields, 'delz', delz)
  call pointer_field_array(inc%fields, 'delz', w)
   lm%pert%delz = delz
   lm%pert%w    = w
endif

end subroutine inc_to_lm

! --------------------------------------------------------------------------------------------------

subroutine lm_to_inc(lm, inc)

implicit none
type(fv3jedi_lm_type),   intent(in)    :: lm
type(fv3jedi_increment), intent(inout) :: inc

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: delp
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: qi
real(kind=kind_real), pointer, dimension(:,:,:) :: ql
real(kind=kind_real), pointer, dimension(:,:,:) :: o3
real(kind=kind_real), pointer, dimension(:,:,:) :: w
real(kind=kind_real), pointer, dimension(:,:,:) :: delz

call pointer_field_array(inc%fields, 'ud'     , ud  )
call pointer_field_array(inc%fields, 'vd'     , vd  )
call pointer_field_array(inc%fields, 't'      , t   )
call pointer_field_array(inc%fields, 'delp'   , delp)
call pointer_field_array(inc%fields, 'sphum'  , q   )
call pointer_field_array(inc%fields, 'ice_wat', qi  )
call pointer_field_array(inc%fields, 'liq_wat', ql  )
if (has_field(inc%fields, 'o3mr') ) call pointer_field_array(inc%fields, 'o3mr'   , o3  )
if (has_field(inc%fields, 'o3ppmv') ) call pointer_field_array(inc%fields, 'o3ppmv'   , o3  )
ud(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%u
vd(inc%isc:inc%iec,inc%jsc:inc%jec,1:inc%npz)   = lm%pert%v
t    = lm%pert%t
delp = lm%pert%delp
q    = lm%pert%qv
qi   = lm%pert%qi
ql   = lm%pert%ql
o3   = lm%pert%o3
if (.not. inc%hydrostatic) then
  call pointer_field_array(inc%fields, 'delz', delz)
  call pointer_field_array(inc%fields, 'delz', w)
  delz = lm%pert%delz
  w    = lm%pert%w
endif

end subroutine lm_to_inc

! --------------------------------------------------------------------------------------------------

subroutine traj_to_traj( traj, lm )

implicit none
type(fv3jedi_traj),    intent(in)    :: traj
type(fv3jedi_lm_type), intent(inout) :: lm

lm%traj%u    = traj%u
lm%traj%v    = traj%v
lm%traj%ua   = traj%ua
lm%traj%va   = traj%va
lm%traj%t    = traj%t
lm%traj%delp = traj%delp
lm%traj%qv   = traj%qv
lm%traj%ql   = traj%ql
lm%traj%qi   = traj%qi
lm%traj%o3   = traj%o3

if (.not. lm%conf%hydrostatic) then
  lm%traj%w    = traj%w
  lm%traj%delz = traj%delz
endif

if (lm%conf%do_phy_mst .ne. 0) then
  lm%traj%qls  = traj%qls
  lm%traj%qcn  = traj%qcn
  lm%traj%cfcn = traj%cfcn
endif

!> Rank two
lm%traj%phis    = traj%phis
lm%traj%frocean = traj%frocean
lm%traj%frland  = traj%frland
lm%traj%varflt  = traj%varflt
lm%traj%ustar   = traj%ustar
lm%traj%bstar   = traj%bstar
lm%traj%zpbl    = traj%zpbl
lm%traj%cm      = traj%cm
lm%traj%ct      = traj%ct
lm%traj%cq      = traj%cq
lm%traj%kcbl    = traj%kcbl
lm%traj%ts      = traj%ts
lm%traj%khl     = traj%khl
lm%traj%khu     = traj%khu

end subroutine traj_to_traj

! --------------------------------------------------------------------------------------------------

end module fv3jedi_tlm_mod
