! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3lm_mod

! oops uses
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! Linear model
use fv3jedi_lm_mod,        only: fv3jedi_lm_type

! fv3-jedi uses
use fv3jedi_fmsnamelist_mod, only: fv3jedi_fmsnamelist
use fv3jedi_kinds_mod,       only: kind_real
use fv3jedi_geom_mod,        only: fv3jedi_geom
use fv3jedi_state_mod,       only: fv3jedi_state

implicit none
private
public :: fv3lm_model

! --------------------------------------------------------------------------------------------------

type :: fv3lm_model
  type(fv3jedi_lm_type) :: fv3jedi_lm  !<Linearized model object
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: initialize
    procedure, public :: step
    procedure, public :: finalize
end type fv3lm_model

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
class(fv3lm_model),        intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt
character(len=:), allocatable :: str
type(fv3jedi_fmsnamelist) :: fmsnamelist

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

! Prepare namelist (nonlinear part)
! ---------------------------------
call fmsnamelist%replace_namelist(conf)

! Call model constructor
! ----------------------
call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

! Safety checks
! -------------

!The full trajecotory of the tlm/adm is not output by this simplified model
!so if being used to generate the trajectry with physics the traj must be read
!from file or obtained by running GEOS or GFS.
if ((self%fv3jedi_lm%conf%do_phy_trb .ne. 0) .or. &
    (self%fv3jedi_lm%conf%do_phy_mst .ne. 0) ) then
   call abor1_ftn("fv3lm_model | FV3LM : unless reading the trajecotory physics should be off")
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3lm_model), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state)

implicit none
class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

call self%fv3jedi_lm%init_nl()

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, geom)

implicit none
class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_geom),  intent(inout) :: geom

call state_to_lm(state,self%fv3jedi_lm)
call self%fv3jedi_lm%step_nl()
call lm_to_state(self%fv3jedi_lm,state)

end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine finalize(self, state)

implicit none
class(fv3lm_model),  intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state

call self%fv3jedi_lm%final_nl()

end subroutine finalize

! --------------------------------------------------------------------------------------------------

subroutine state_to_lm( state, lm )

implicit none
type(fv3jedi_state),   intent(in)    :: state
type(fv3jedi_lm_type), intent(inout) :: lm

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: delp
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: qi
real(kind=kind_real), pointer, dimension(:,:,:) :: ql
real(kind=kind_real), pointer, dimension(:,:,:) :: o3
real(kind=kind_real), pointer, dimension(:,:,:) :: w
real(kind=kind_real), pointer, dimension(:,:,:) :: delz
real(kind=kind_real), pointer, dimension(:,:,:) :: phis

call state%get_field('ud'     , ud  )
call state%get_field('vd'     , vd  )
call state%get_field('t'      , t   )
call state%get_field('delp'   , delp)
call state%get_field('sphum'  , q   )
call state%get_field('ice_wat', qi  )
call state%get_field('liq_wat', ql  )
if (state%has_field('o3mr'  )) call state%get_field('o3mr'  , o3)
if (state%has_field('o3ppmv')) call state%get_field('o3ppmv', o3)
lm%traj%ua = 0.0_kind_real
lm%traj%va = 0.0_kind_real

lm%traj%u       = ud(state%isc:state%iec,state%jsc:state%jec,:)
lm%traj%v       = vd(state%isc:state%iec,state%jsc:state%jec,:)
lm%traj%t       = t
lm%traj%delp    = delp
lm%traj%qv      = q
lm%traj%qi      = qi
lm%traj%ql      = ql
lm%traj%o3      = o3

if (state%has_field('ua')) then
  call state%get_field('ua',   ua  )
  call state%get_field('va',   va  )
  lm%traj%ua = ua
  lm%traj%va = va
endif

if (.not. lm%conf%hydrostatic) then
  call state%get_field('w   ', w   )
  call state%get_field('delz', delz)
  lm%traj%w       = w
  lm%traj%delz    = delz
endif

call state%get_field('phis', phis )
lm%traj%phis = phis(:,:,1)

end subroutine state_to_lm

! --------------------------------------------------------------------------------------------------

subroutine lm_to_state( lm, state )

implicit none
type(fv3jedi_lm_type), intent(in)    :: lm
type(fv3jedi_state),   intent(inout) :: state

real(kind=kind_real), pointer, dimension(:,:,:) :: ud
real(kind=kind_real), pointer, dimension(:,:,:) :: vd
real(kind=kind_real), pointer, dimension(:,:,:) :: ua
real(kind=kind_real), pointer, dimension(:,:,:) :: va
real(kind=kind_real), pointer, dimension(:,:,:) :: t
real(kind=kind_real), pointer, dimension(:,:,:) :: delp
real(kind=kind_real), pointer, dimension(:,:,:) :: q
real(kind=kind_real), pointer, dimension(:,:,:) :: qi
real(kind=kind_real), pointer, dimension(:,:,:) :: ql
real(kind=kind_real), pointer, dimension(:,:,:) :: o3
real(kind=kind_real), pointer, dimension(:,:,:) :: w
real(kind=kind_real), pointer, dimension(:,:,:) :: delz
real(kind=kind_real), pointer, dimension(:,:,:) :: phis

call state%get_field('ud'     , ud  )
call state%get_field('vd'     , vd  )
call state%get_field('t'      , t   )
call state%get_field('delp'   , delp)
call state%get_field('sphum'  , q   )
call state%get_field('ice_wat', qi  )
call state%get_field('liq_wat', ql  )
if ( state%has_field('o3mr' ) )call state%get_field('o3mr'   , o3  )
if ( state%has_field('o3ppmv' ) )call state%get_field('o3ppmv'   , o3  )
ud(state%isc:state%iec,state%jsc:state%jec,:)      = lm%traj%u
vd(state%isc:state%iec,state%jsc:state%jec,:)      = lm%traj%v
t       = lm%traj%t
delp    = lm%traj%delp
q       = lm%traj%qv
qi      = lm%traj%qi
ql      = lm%traj%ql
o3      = lm%traj%o3

if (state%has_field('ua')) then
  call state%get_field('ua',   ua  )
  call state%get_field('va',   va  )
  ua = lm%traj%ua
  va = lm%traj%va
endif

if (.not. lm%conf%hydrostatic) then
  call state%get_field('w   ', w   )
  call state%get_field('delz', delz)
  w       = lm%traj%w
  delz    = lm%traj%delz
endif

call state%get_field('phis', phis )
phis(:,:,1)    = lm%traj%phis

end subroutine lm_to_state

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fv3lm_mod
