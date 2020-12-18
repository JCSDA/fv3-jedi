! (C) Copyright 2020 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#define esmf_err_abort(rc) if (esmf_LogFoundError(rc, msg="Aborting ESMF", line=__LINE__, file=__FILE__)) call esmf_Finalize(endflag=esmf_END_ABORT)

module fv3jedi_ufs_mod

! oops
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! fv3jedi
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state

use ESMF
!
!-----------------------------------------------------------------------
!

implicit none
private


public :: model_ufs

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: model_ufs
  type(esmf_GridComp) :: esmComp !ufs driver
  type(esmf_State) :: importState, exportState 
  type(esmf_time) :: starttime                                         !<-- the esmf start time.
  type(esmf_Clock) :: clock 
  type(esmf_vm) :: vm
  type(esmf_config) :: cf_main                                         !<-- the configure object
  integer :: dt
  integer :: mype
  contains
    procedure :: create
    procedure :: delete
    procedure :: initialize
    procedure :: step
    procedure :: finalize
end type model_ufs

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, geom)

class(model_ufs),          intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fv3jedi_geom),        intent(in)    :: geom


end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(model_ufs), intent(inout) :: self

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state, vdate)

class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state
type(datetime),      intent(in)    :: vdate


end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, vdate_start, vdate_final)


class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate_start
type(datetime),      intent(in)    :: vdate_final


end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine destructor(self)
class(model_ufs),    intent(inout) :: self

end subroutine destructor

subroutine finalize(self, state, vdate)

class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state
type(datetime),intent(in)    :: vdate


end subroutine finalize

end module fv3jedi_ufs_mod
