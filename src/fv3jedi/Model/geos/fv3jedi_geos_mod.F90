! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#define I_AM_MAIN
#include "MAPL_Generic.h"
#include "unused_dummy.H"

module fv3jedi_geos_mod

use iso_c_binding
use datetime_mod
use duration_mod
use netcdf

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_field_mod,     only: fv3jedi_field, field_clen
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state
use fckit_mpi_module,      only: fckit_mpi_comm

use MPI
use ESMF
use MAPL_CapMod, only: MAPL_Cap
use MAPL_ApplicationSupport

use MAPL_Mod
use MAPL_Profiler, only: get_global_time_profiler, BaseProfiler, TimeProfiler
use GEOS_GcmGridCompMod, only: GEOS_Gcm => SetServices

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
  type(fckit_mpi_comm) :: f_comm
  integer :: isc, iec, jsc, jec, npz
  type(ESMF_TIme) :: start_Time
  logical :: saved_state
  logical :: reforecast
end type geos_model

! --------------------------------------------------------------------------------------------------

  type :: Field_Attributes
     type(ESMF_Field) :: field(1)
     character(len=ESMF_MAXSTR) :: short_name, long_name, units
  end type Field_Attributes

contains

! --------------------------------------------------------------------------------------------------

subroutine geos_create(self, geom, conf)

implicit none
type(geos_model),          intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

type(ESMF_VM)            :: vm
type(MAPL_CapOptions)    :: cap_options
type(MAPL_Communicators) :: mapl_comm

integer :: rc, subcommunicator

type(duration) :: dtstep
integer :: geos_dt, jedi_dt
character(len=20) :: ststep
character(len=:), allocatable :: str
logical :: esmf_logging
type(ESMF_LOGKIND_FLAG) :: esmf_logkind

! Grid dimensions
! ---------------
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz
self%f_comm = geom%f_comm

! Set up options for the cap
! --------------------------
cap_options = MAPL_CapOptions(cap_rc_file='CAP.rc')
cap_options%use_comm_world = .false.
cap_options%comm = geom%f_comm%communicator()
cap_options%npes_model = geom%f_comm%size()

! Create the MAPL_Cap object
! --------------------------
self%Cap = MAPL_Cap( name='GCM', set_services=GEOS_Gcm, cap_options = cap_options)

! MPI
! ---
call self%cap%initialize_mpi()

! CapCom
! ------
mapl_comm%mapl%comm = geom%f_comm%communicator()
mapl_comm%global%comm = geom%f_comm%communicator()
mapl_comm%esmf%comm = geom%f_comm%communicator()
mapl_comm%io%comm = MPI_COMM_NULL

call fill_comm(mapl_comm%mapl)
call fill_comm(mapl_comm%global)
call fill_comm(mapl_comm%esmf)
call fill_comm(mapl_comm%io)

! IO Server communicator (not on by default)
! -----------------------------------------
subcommunicator = self%cap%create_member_subcommunicator(self%cap%get_comm_world(), rc=rc)
call self%cap%initialize_io_clients_servers(subcommunicator, rc = rc)

if (conf%has('ESMF_Logging')) call conf%get_or_die('ESMF_Logging',esmf_logging)
! Initialize ESMF
! ---------------
if (esmf_logging) then
   esmf_logkind = ESMF_LOGKIND_MULTI
else
   esmf_logkind = ESMF_LOGKIND_NONE
end if
call ESMF_Initialize (vm=vm, logKindFlag=esmf_logkind, mpiCommunicator=mapl_comm%esmf%comm, rc=rc)

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
call conf%get_or_die("tstep",str)
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

if (geom%f_comm%rank() == 0) then
   print*, "There are ", int(self%GEOSsubsteps,2), " time steps of GEOS for each time step of JEDI"
endif

self%start_time = self%cap%cap_gc%get_current_time(rc=rc)
self%saved_state = .false.
if (conf%has('reforecast')) call conf%get_or_die('reforecast',self%reforecast)

end subroutine geos_create

! --------------------------------------------------------------------------------------------------

subroutine geos_delete(self)

implicit none
type(geos_model), intent(inout) :: self

integer :: rc

! Finalize GEOS
! -------------
call self%cap%cap_gc%destroy_state(rc = rc)
call self%cap%cap_gc%finalize(rc = rc);
if (rc.ne.0) call abor1_ftn("geos_mod: finalize failed")
call self%cap%finalize_io_clients_servers()
call MAPL_Finalize(comm=self%f_comm%communicator())

! Finalize ESMF
! -------------
call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

end subroutine geos_delete

! --------------------------------------------------------------------------------------------------

subroutine geos_initialize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state
integer :: rc
!There is no GEOS equivalent to JEDI initialize

if (self%saved_state) then
   call self%cap%rewind_model(self%start_time,rc=rc)
   call self%cap%cap_gc%refresh_state(rc=rc)
end if

if (self%reforecast .and. (.not.self%saved_state) ) then
   call self%cap%cap_gc%record_state(rc=rc)
   self%saved_state = .true.
end if



end subroutine geos_initialize

! --------------------------------------------------------------------------------------------------

subroutine geos_step(self, state)

implicit none
type(geos_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state

integer :: n, rc

!Convert JEDI state to GEOS state
!call state_to_geos( state, self )

!Cycle GEOS through this time step of JEDI
do n = 1,self%GEOSsubsteps
  rc = ESMF_SUCCESS
  call self%cap%step_model(rc = rc)
  if (rc.ne.0) call abor1_ftn("geos_mod: step_model failed")
enddo

!Retrieve GEOS state and put into JEDI state
call geos_to_state( self, state )

end subroutine geos_step

! --------------------------------------------------------------------------------------------------

subroutine geos_finalize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

!There is no GEOS equivalent to JEDI finalize

end subroutine geos_finalize

! --------------------------------------------------------------------------------------------------

! subroutine state_to_geos( state, self )
!
! implicit none
! type(fv3jedi_state), intent(in)    :: state
! type(geos_model),    intent(inout) :: self
!
! !TODO
!
! ! 1. Overwrite the MKIAU exports?
! ! 2. Dont compile with MKIAU and intercept EXT data?
! ! 3. Cinderella grid comp for data?
!
! end subroutine state_to_geos

! --------------------------------------------------------------------------------------------------

subroutine geos_to_state( self, state )

implicit none
type(geos_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state

integer :: num_items, i, rc, rank, lb(3), ub(3), fnpz
type(ESMF_Field) :: field
character(len=ESMF_MAXSTR), allocatable :: item_names(:)
real(kind=ESMF_KIND_R4), pointer :: farrayPtr2(:,:)
real(kind=ESMF_KIND_R4), pointer :: farrayPtr3(:,:,:)
character(len=field_clen) :: fv3jedi_name
type(fv3jedi_field), pointer :: field_ptr
real(kind=kind_real), allocatable, dimension(:,:,:) :: field_geos


! Array to hold output from GEOS in JEDI precision
! ------------------------------------------------
allocate(field_geos(self%isc:self%iec, self%jsc:self%jec, self%npz+1))


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

  ! Create map between GEOS name and fv3-jedi name
  ! ----------------------------------------------
  select case (trim(item_names(i)))

    ! The below need to be listed in Cap.rc as e.g.
    ! CAP_EXPORTS:
    ! U_DGRID,DYN
    ! Q,MOIST
    ! ::

    ! DYN
    case ("U_DGRID")
      fv3jedi_name = 'ud'
    case ("V_DGRID")
      fv3jedi_name = 'vd'
    case ("PT")
      fv3jedi_name = 'pt'
    case ("PE")
      fv3jedi_name = 'pe'
    case ("T")
      fv3jedi_name = 't'
    case ("DELP")
      fv3jedi_name = 'delp'
    case ("PS")
      fv3jedi_name = 'ps'
    case ("TA")
      fv3jedi_name = 'ts'
    case ("DZ")
      fv3jedi_name = 'delz'
    case ("W")
      fv3jedi_name = 'w'

    ! AGCM
    case ("PHIS")
      fv3jedi_name = 'phis'
    case ("QITOT")
      fv3jedi_name = 'ice_wat'
    case ("QLTOT")
      fv3jedi_name = 'liq_wat'
    case ("VARFLT")
      fv3jedi_name = 'varflt'

    ! MOIST
    case ("Q")
      fv3jedi_name = 'sphum'
    case ("QICN")
      fv3jedi_name = "qicn"
    case ("QLCN")
      fv3jedi_name = "qlcn"
    case ("QILS")
      fv3jedi_name = "qils"
    case ("QLLS")
      fv3jedi_name = "qlls"
    case ("QSTOT")
      fv3jedi_name = "qs"
    case ("QRTOT")
      fv3jedi_name = "qr"
    case ("QCLSX0")
      fv3jedi_name = 'qls'
    case ("QCCNX0")
      fv3jedi_name = 'qcn'
    case ("CLCNX0")
      fv3jedi_name = 'cfcn'
    case ("KCBL_moist")
      fv3jedi_name = 'kcbl'
    case ("TS_moist")
      fv3jedi_name = 'tsm'
    case ("KHl_moist")
      fv3jedi_name = 'khl'
    case ("KHu_moist")
      fv3jedi_name = 'khu'

    ! CHEMISTRY
    case ("O3")
      fv3jedi_name = 'o3ppmv'

    ! SURFACE
    case ("FRLAND")
      fv3jedi_name = 'frland'
    case ("FRLANDICE")
      fv3jedi_name = 'frlandice'
    case ("FRLAKE")
      fv3jedi_name = 'frlake'
    case ("FROCEAN")
      fv3jedi_name = 'frocean'
    case ("FRACI")
      fv3jedi_name = 'frseaice'
    case ("USTAR")
      fv3jedi_name = 'ustar'
    case ("BSTAR")
      fv3jedi_name = 'bstar'
    case ("CM")
      fv3jedi_name = 'cm'
    case ("CT")
      fv3jedi_name = 'ct'
    case ("CQ")
      fv3jedi_name = 'cq'
    case ("U10N")
      fv3jedi_name = 'u_srf'
    case ("V10N")
      fv3jedi_name = 'v_srf'
    case ('SNOMAS')
      fv3jedi_name = 'sheleg'
    case ('TSOIL1')
      fv3jedi_name = 'soilt'
    case ('WET1')
      fv3jedi_name = 'soilm'

    !TURBULENCE
    case ('ZPBL')
      fv3jedi_name = 'zpbl'

    ! NO MAP
    case default
      call abor1_ftn("geos_to_state no map between GEOS name and JEDI name "//trim(item_names(i)))

  end select


  ! Only need to extract field from GEOS if fv3-jedi needs it
  ! ---------------------------------------------------------
  if (state%has_field(trim(fv3jedi_name))) then

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
    field_geos = 0.0_kind_real
    if (rank == 2) then

      call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb(1:2), totalUBound = ub(1:2), rc = rc )
      if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldGet 2D failed")

      fnpz = 1
      field_geos(self%isc:self%iec,self%jsc:self%jec,1) = farrayPtr2(lb(1):ub(1),lb(2):ub(2))
      nullify(farrayPtr2)

    elseif (rank == 3) then

      call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )
      if (rc.ne.0) call abor1_ftn("geos_to_state: ESMF_FieldGet 3D failed")

      fnpz = ub(3)-lb(3)+1
      field_geos(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))
      nullify(farrayPtr3)

    else

      call abor1_ftn("geos_mod: can only handle rank 2 or rank 3 fields from GEOS")

    endif

    ! Check that dimensions match
    if ((ub(1)-lb(1)+1 .ne. self%iec-self%isc+1) .or. (ub(2)-lb(2)+1 .ne. self%jec-self%jsc+1) ) then
      call abor1_ftn("geos_to_state: dimension mismatch between JEDI and GEOS horizontal grid")
    endif

    ! Get pointer to fv3-jedi side field
    call state%get_field(trim(fv3jedi_name), field_ptr)

    if (field_ptr%npz .ne. fnpz) &
      call abor1_ftn("geos_to_state: dimension mismatch between JEDI and GEOS vertical grid")

    ! Copy from GEOS to fv3-jedi
    field_ptr%array(:,:,1:fnpz) = field_geos(:,:,1:fnpz)
  endif

end do

deallocate(item_names)

end subroutine geos_to_state

! --------------------------------------------------------------------------------------------------

subroutine allocate_exports(self)

implicit none
type(geos_model),    intent(inout) :: self

integer :: num_items, i, rc
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

! --------------------------------------------------------------------------------------------------

subroutine fill_comm(mcomm)

  type (MAPL_Communicator), intent(inout) :: mcomm

  integer :: comm
  integer :: status

  comm = mcomm%comm
  if (comm == MPI_COMM_NULL) then
     mcomm%size = 0
     mcomm%rank = MPI_UNDEFINED
     mcomm%root = MPI_UNDEFINED
  else
     call MPI_Comm_Size(comm, mcomm%size,status)
     _VERIFY(status)
     call MPI_Comm_Rank(comm, mcomm%rank,status)
     _VERIFY(status)
     mcomm%root = 0
  end if
end subroutine fill_comm

! --------------------------------------------------------------------------------------------------


end module fv3jedi_geos_mod
