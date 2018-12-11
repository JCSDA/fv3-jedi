! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

#include "MAPL_Generic.h"

module fv3jedi_geos_mod

use iso_c_binding
use config_mod
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

use MPI
use ESMF
use MAPL_Mod
use MAPL_CapMod, only: MAPL_Cap
use GEOS_GcsGridCompMod, only: GEOS_GcsSS => SetServices

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
  type(MAPL_Cap), pointer :: cap
  integer :: GEOSsubsteps
end type geos_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine geos_create(self, geom, c_conf)

implicit none
type(c_ptr),        intent(in)    :: c_conf
type(geos_model),   intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom

type(fckit_mpi_comm) :: f_comm
integer :: rc
integer :: subcommunicator
integer :: geos_dt

character(len=20) :: ststep
type(duration) :: dtstep
integer :: jedi_dt


! FCKIT MPI wrapper for communicator
! ----------------------------------
f_comm = fckit_mpi_comm()


! Initialize ESMF
! ---------------
call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
  defaultCalkind=ESMF_CALKIND_GREGORIAN, mpiCommunicator=f_comm%communicator(), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)


! Create the MAPL_Cap object
! --------------------------
allocate(self%Cap)
self%Cap = MAPL_Cap(name='GEOS', set_services=GEOS_GcsSS, &
                    comm=f_comm%communicator(), cap_rc_file='CAP.rc')


! Cap default values
! ------------------
call self%cap%set_n_members(1)
call self%cap%set_npes_input_server(0)
call self%cap%set_npes_output_server(0)
call self%cap%set_npes_model(f_comm%size())
call self%cap%set_npes_member(self%cap%get_npes_model() / self%cap%get_n_members())
call self%cap%set_ensemble_subdir_prefix("mem")


! MPI
! ---
call self%cap%set_comm_world(f_comm%communicator())
call self%cap%set_rank(f_comm%rank())


! IO Server commincator (not on by default)
! -----------------------------------------
subcommunicator = self%cap%create_member_subcommunicator(self%cap%get_comm_world(), rc=rc); _VERIFY(rc)
call self%cap%initialize_io_servers(subcommunicator, rc = rc); _VERIFY(rc)


! Intialize the Cap GridComp
! --------------------------
call self%cap%initialize_cap_gc(self%cap%get_mapl_comm())


! GEOS SetServices
! ----------------
call self%cap%cap_gc%set_services(rc = rc); _VERIFY(rc)


! GEOS Initialize
! ---------------
call self%cap%cap_gc%initialize(rc = rc); _VERIFY(rc)


! Time step checks and get number of GEOSsubsteps if any
! ------------------------------------------------------
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
jedi_dt = int(duration_seconds(dtstep))

geos_dt = self%cap%cap_gc%get_heartbeat_dt(rc = rc)

if (jedi_dt < geos_dt) then
  call abor1_ftn("JEDI model time step should not be less than GEOS time step")
elseif (mod(jedi_dt,geos_dt) .ne. 0) then
  call abor1_ftn("JEDI time step needs to be divisible by GEOS time step")
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
call self%cap%cap_gc%finalize(rc = rc); _VERIFY(rc)
deallocate(self%cap)

! Finalize ESMF
call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

end subroutine geos_delete

! ------------------------------------------------------------------------------

subroutine geos_initialize(self, state)

implicit none
type(geos_model), target :: self
type(fv3jedi_state)      :: state

!There is no GEOS equivalent to JEDI initialize

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
 

end subroutine state_to_geos


!CAP_EXPORTS:
!	U,DYN
!	V,DYN
!::
!
!CAP_IMPORTS:
!	U
!	V
!::

! ------------------------------------------------------------------------------

subroutine geos_to_state( self, state )

implicit none
type(geos_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state
 
end subroutine geos_to_state

! ------------------------------------------------------------------------------

end module fv3jedi_geos_mod
