! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#define I_AM_MAIN
#include "MAPL_Generic.h"
#include "unused_dummy.H"

module fv3jedi_geos_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment

use fckit_mpi_module, only: fckit_mpi_comm

use MPI
use ESMF
use MAPL_Mod
use MAPL_CapMod, only: MAPL_Cap
use GEOS_GcsGridCompMod, only: GEOS_GcsSS => SetServices

use FV_StateMod, only: fv_hydrostatic

implicit none
private

public :: geos_model
public :: geos_create
public :: geos_delete
public :: geos_initialize
public :: geos_step
public :: geos_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: geos_model
  type(MAPL_Cap) :: cap
  integer :: GEOSsubsteps
end type geos_model

! ------------------------------------------------------------------------------

  type :: Field_Attributes
     type(ESMF_Field) :: field(1)
     character(len=ESMF_MAXSTR) :: short_name, long_name, units
  end type Field_Attributes

contains

! ------------------------------------------------------------------------------

subroutine geos_create(self, geom, c_conf)

implicit none
type(c_ptr),        intent(in)    :: c_conf
type(geos_model),   intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom

integer :: rc
integer :: subcommunicator, mapl_communicator
integer :: geos_dt

character(len=20) :: ststep
type(duration) :: dtstep
integer :: jedi_dt, i

type (ESMF_VM) :: vm
type(MAPL_Communicators) :: mapl_comm
type(fckit_mpi_comm) :: f_comm

type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)

! FCKIT MPI wrapper for communicator
! ----------------------------------
f_comm = fckit_mpi_comm()


! Duplicate the communicator
! --------------------------
call MPI_Comm_dup(f_comm%communicator(), mapl_communicator, rc);
if (rc.ne.0) call abor1_ftn("geos_mod: MPI_Comm_dup failed")


! Create the MAPL_Cap object
! --------------------------
self%Cap = MAPL_Cap( name='GCM', set_services=GEOS_GcsSS, &
                     comm=mapl_communicator, cap_rc_file='CAP.rc')


! MPI
! ---
call self%cap%initialize_mpi() !This is only needed to set cap%rank


! Cap default values
! ------------------
call initialize_cap_default_values(self%cap)


! IO Server commincator (not on by default)
! -----------------------------------------
subcommunicator = self%cap%create_member_subcommunicator(self%cap%get_comm_world(), rc=rc);
if (rc.ne.0) call abor1_ftn("geos_mod: create_member_subcommunicator failed")
if (subcommunicator /= MPI_COMM_NULL) then
  call self%cap%initialize_io_servers(subcommunicator, rc = rc);
  if (rc.ne.0) call abor1_ftn("geos_mod: initialize_io_servers failed")
endif


! Initialize ESMF
! ---------------
mapl_comm = self%cap%get_mapl_comm()
call ESMF_Initialize (vm=vm, logKindFlag=ESMF_LOGKIND_NONE, mpiCommunicator=mapl_comm%esmf%comm, rc=rc)


! Intialize the Cap GridComp
! --------------------------
call self%cap%initialize_cap_gc(mapl_comm)


! GEOS SetServices
! ----------------
call self%cap%cap_gc%set_services(rc = rc);
if (rc.ne.0) call abor1_ftn("geos_mod: set_services failed")


! GEOS Initialize
! ---------------
call self%cap%cap_gc%initialize(rc = rc);
if (rc.ne.0) call abor1_ftn("geos_mod: initialize failed")


! Make sure required exports are allocated
! ----------------------------------------
call allocate_exports(self)


! Time step checks and get number of GEOSsubsteps if any
! ------------------------------------------------------
call f_conf%get_or_die("tstep",str)
ststep = str
deallocate(str)
dtstep = trim(ststep)
jedi_dt = int(duration_seconds(dtstep))

geos_dt = self%cap%cap_gc%get_heartbeat_dt(rc = rc)
if (rc.ne.0) call abor1_ftn("geos_mod: error retrieving heartbeat")

if (jedi_dt < geos_dt) then
  call abor1_ftn("geos_mod: JEDI model time step should not be less than GEOS time step")
elseif (mod(jedi_dt,geos_dt) .ne. 0) then
  call abor1_ftn("geos_mod: JEDI time step needs to be divisible by GEOS time step")
endif

self%GEOSsubsteps = jedi_dt/geos_dt

if (f_comm%rank() == 0) then
   print*, "There are ", int(self%GEOSsubsteps,2), " time steps of GEOS for each time step of JEDI"
endif

end subroutine geos_create

! ------------------------------------------------------------------------------

subroutine geos_delete(self)

implicit none
type(geos_model), intent(inout) :: self

integer :: rc

! Finalize GEOS
!call self%cap%cap_gc%finalize(rc = rc);
!if (rc.ne.0) call abor1_ftn("geos_mod: finalize failed")

! Finalize ESMF
!call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

end subroutine geos_delete

! ------------------------------------------------------------------------------

subroutine geos_initialize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

!There is no GEOS equivalent to JEDI initialize

!Could be where replay mode is activated for outer loops

end subroutine geos_initialize

! ------------------------------------------------------------------------------

subroutine geos_step(self, state, vdate)

implicit none
type(geos_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate !< Valid datetime after step

integer :: n, rc

!Convert JEDI state to GEOS state
call state_to_geos( state, self )

!Cycle GEOS through this time step of JEDI
do n = 1,self%GEOSsubsteps
  rc = ESMF_SUCCESS
  call self%cap%step_model(rc = rc)
  if (rc.ne.0) call abor1_ftn("geos_mod: step_model failed")
enddo

!Retrieve GEOS state and put into JEDI state
call geos_to_state( self, state )

end subroutine geos_step

! ------------------------------------------------------------------------------

subroutine geos_finalize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

!There is no GEOS equivalent to JEDI finalize

end subroutine geos_finalize

! ------------------------------------------------------------------------------

subroutine state_to_geos( state, self )

implicit none
type(fv3jedi_state), intent(in)    :: state
type(geos_model),    intent(inout) :: self

!TODO

! 1. Overwrite the MKIAU exports?
! 2. Dont compile with MKIAU and intercept EXT data?
! 3. Cinderella grid comp for data?

end subroutine state_to_geos

! ------------------------------------------------------------------------------

subroutine geos_to_state( self, state )

implicit none
type(geos_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state

integer :: num_items, i, rc, rank
type(ESMF_Field) :: field
character(len=ESMF_MAXSTR), allocatable :: item_names(:)

real(kind=ESMF_KIND_R4), pointer :: farrayPtr2(:,:)
real(kind=ESMF_KIND_R4), pointer :: farrayPtr3(:,:,:)

integer :: isc,iec,jsc,jec,npz
integer :: lb2(2), ub2(2)
integer :: lb3(3), ub3(3)


! Convenience
! -----------
isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
npz = state%npz


! Get number of items
! -------------------
call ESMF_StateGet(self%cap%cap_gc%export_state, itemcount = num_items, rc = rc)
if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_StateGet itemcount failed")


! Get names of the items
! ----------------------
allocate(item_names(num_items))
call ESMF_StateGet(self%cap%cap_gc%export_state, itemnamelist = item_names, rc = rc)
if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_StateGet itemnamelist failed")


! Loop over states coming from GEOS and convert to JEDI state
! -----------------------------------------------------------
do i = 1, num_items

   !Get field from the state
   call ESMF_StateGet(self%cap%cap_gc%export_state, item_names(i), field, rc = rc)
   if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_StateGet field failed")

   !Validate the field
   call ESMF_FieldValidate(field, rc = rc)
   if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldValidate failed")

   !Get the field rank
   call ESMF_FieldGet(field, rank = rank, rc = rc)
   if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldGet rank failed")


   !Convert field to pointer and pointer bounds
   if (rank == 2) then

     call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb2, totalUBound = ub2, rc = rc )
     if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldGet 2D failed")

    if (((ub2(1)-lb2(1)+1) .ne. (iec-isc+1)) .or. &
         ((ub2(2)-lb2(2)+1) .ne. (jec-jsc+1)) ) then
       call abor1_ftn("geos_to_state: rank 2 dimension mismatch between JEDI and GEOS grid")
     endif

   elseif (rank == 3) then

     call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb3, totalUBound = ub3, rc = rc )
     if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldGet 3D failed")

     if (((ub3(1)-lb3(1)+1) .ne. (iec-isc+1)) .or. &
         ((ub3(2)-lb3(2)+1) .ne. (jec-jsc+1)) .or. &
         ((ub3(3)-lb3(3)+1) .ne. (npz)) ) then
       call abor1_ftn("geos_to_state: rank 3 dimension mismatch between JEDI and GEOS grid")
     endif

   else

     call abor1_ftn("geos_mod: rank error in getting the state")

   endif
!TODO precision change
   !Fill up JEDI state
   select case (trim(item_names(i)))
   case ("U_DGRID")
     state%ud  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("V_DGRID")
     state%vd  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("U")
     state%ua  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("V")
     state%va  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("T")
     state%t   (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("DELP")
     state%delp(isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("DZ")
     state%delz(isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("W")
     state%w   (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("Q")
     state%q   (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("QITOT")
     state%qi  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("QLTOT")
     state%ql  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("O3")
     state%o3  (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("QCLSX0")
     state%qls (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("QCCNX0")
     state%qcn (isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("CLCNX0")
     state%cfcn(isc:iec,jsc:jec,1:npz) = farrayPtr3(lb3(1):ub3(1),lb3(2):ub3(2),lb3(3):ub3(3))
   case ("PHIS")
     state%phis   (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("FRLAND")
     state%frland (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("FROCEAN")
     state%frocean(isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("KCBL_moist")
     state%kcbl   (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("TS_moist")
     state%ts     (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("KHl_moist")
     state%khl    (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("KHu_moist")
     state%khu    (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("VARFLT")
     state%varflt (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("USTAR")
     state%ustar  (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("BSTAR")
     state%bstar  (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("ZPBL")
     state%zpbl   (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("CM")
     state%cm     (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("CT")
     state%ct     (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case ("CQ")
     state%cq     (isc:iec,jsc:jec,1) = farrayPtr2(lb2(1):ub2(1),lb2(2):ub2(2))
   case default
     call abor1_ftn("geos_to_state unknown variable")
   end select

   if (associated(farrayPtr3)) nullify(farrayPtr3)
   if (associated(farrayPtr2)) nullify(farrayPtr2)

end do

deallocate(item_names)

end subroutine geos_to_state

! ------------------------------------------------------------------------------

subroutine initialize_cap_default_values(cap)

implicit none
class(MAPL_Cap), intent(inout) :: cap
integer :: ierror, npes_model

call cap%set_n_members(1)
call cap%set_npes_input_server(0)
call cap%set_npes_output_server(0)
call cap%set_nodes_input_server(0)
call cap%set_nodes_output_server(0)
call MPI_Comm_size(cap%get_comm_world(), npes_model, ierror)
call cap%set_npes_model(npes_model)
call cap%set_npes_member(cap%get_npes_model() / cap%get_n_members())
call cap%set_ensemble_subdir_prefix("mem")

end subroutine initialize_cap_default_values

! ------------------------------------------------------------------------------

subroutine allocate_exports(self)

implicit none
type(geos_model),    intent(inout) :: self

integer :: num_items, i, rc, rank
type(ESMF_Field) :: field
character(len=ESMF_MAXSTR), allocatable :: item_names(:)


 ! MAPL exports are only allocated if a connection is found between the
 ! export and something else, e.g. HISTORY. This routine creates this
 ! link between exports and MAPL_CapGridComp

 ! Get number of items
 ! -------------------
 call ESMF_StateGet(self%cap%cap_gc%export_state, itemcount = num_items, rc = rc)
 if (rc.ne.0) call abor1_ftn("allocate_exports: ESMF_StateGet itemcount failed")


 ! Get names of the items
 ! ----------------------
 allocate(item_names(num_items))
 call ESMF_StateGet(self%cap%cap_gc%export_state, itemnamelist = item_names, rc = rc)
 if (rc.ne.0) call abor1_ftn("allocate_exports: ESMF_StateGet items failed")


 ! Loop over states and allocate coupling
 ! --------------------------------------
 do i = 1, num_items

   !Get field from the state
   call ESMF_StateGet(self%cap%cap_gc%export_state, item_names(i), field, rc = rc)
   if (rc.ne.0) call abor1_ftn("allocate_exports: ESMF_StateGet field failed")

   !Validate the field
   call ESMF_FieldValidate(field, rc = rc)
   if (rc.ne.0) call abor1_ftn("allocate_exports: ESMF_FieldValidate failed")

   call MAPL_AllocateCoupling(field, rc)
   if (rc.ne.0) call abor1_ftn("allocate_exports: MAPL_AllocateCoupling failed")

 end do

 deallocate(item_names)

end subroutine allocate_exports

! ------------------------------------------------------------------------------

end module fv3jedi_geos_mod
