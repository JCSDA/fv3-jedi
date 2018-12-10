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

! FCKIT MPI wrapper for global communicator
f_comm = fckit_mpi_comm()

! Initialize ESMF
call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
  defaultCalkind=ESMF_CALKIND_GREGORIAN, mpiCommunicator=f_comm%communicator(), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Create the MAPL_Cap object
allocate(self%Cap)

!Create CAP
self%Cap = MAPL_Cap(name='GEOS', set_services=GEOS_GcsSS, &
                    comm=f_comm%communicator(), cap_rc_file='CAP.rc')

!MPI
call self%cap%initialize_mpi(rc = rc); _VERIFY(rc)

!IO Server commincator (not on by default)
!subcommunicator = self%cap%create_member_subcommunicator(self%cap%get_comm_world(), rc=rc); _VERIFY(rc)
!call self%cap%initialize_io_servers(subcommunicator, rc = rc); _VERIFY(rc)

call self%cap%initialize_cap_gc(self%cap%get_mapl_comm())

print*, 'fv3jedi_geos_mod: CAP created'

!GEOS GCS SetServices
call self%cap%cap_gc%set_services(rc = rc); _VERIFY(rc)

print*, 'fv3jedi_geos_mod: SetServices complete'

call self%cap%cap_gc%initialize(rc = rc); _VERIFY(rc)

print*, 'fv3jedi_geos_mod: Initialize complete'

end subroutine geos_create

! ------------------------------------------------------------------------------

subroutine geos_delete(self)

implicit none
type(geos_model), intent(inout) :: self

integer :: rc

call self%cap%cap_gc%finalize(rc = rc); _VERIFY(rc)
deallocate(self%cap)

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

integer :: rc

call state_to_geos( state, self )

call self%cap%step_model(rc = rc)

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
