! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_nuopc_mod

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
use NUOPC_Comp, only: NUOPC_CompSetClock
use NUOPC_Driver, only: NUOPC_DriverGetComp
use NUOPC_Model, only: NUOPC_ModelGet
use ESM, only: esmSS => SetServices

implicit none
private

public :: model_nuopc_type
public :: model_nuopc_create
public :: model_nuopc_delete
public :: model_nuopc_initialize
public :: model_nuopc_step
public :: model_nuopc_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: model_nuopc_type
  integer             :: urc
  integer             :: dt
  type(ESMF_GridComp) :: esmComp !NUOPC driver
end type model_nuopc_type

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine model_nuopc_create(self, geom, c_conf)

implicit none
type(c_ptr),            intent(in)    :: c_conf
type(model_nuopc_type), intent(inout) :: self
type(fv3jedi_geom),     intent(in)    :: geom

integer :: rc

character(len=20) :: ststep
type(duration) :: dtstep
character(len=20) :: cdate_start
character(len=20) :: cdate_final

type(fckit_mpi_comm) :: f_comm

! FCKIT MPI wrapper for global communicator
f_comm = fckit_mpi_comm()

! Initialize ESMF
call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
  defaultCalkind=ESMF_CALKIND_GREGORIAN, mpiCommunicator=f_comm%communicator(), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_LogWrite("esm-JEDI STARTING", ESMF_LOGMSG_INFO, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Create the Earth system Component
self%esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
! SetServices for the Earth system Component
call ESMF_GridCompSetServices(self%esmComp, esmSS, userRc=self%urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=self%urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)


! Reset the clock based on what JEDI provides
! -------------------------------------------

ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
self%dt = int(duration_seconds(dtstep))

cdate_start = '2018-04-14T21:00:00Z'
cdate_final = '2018-04-15T03:00:00Z'

call nuopc_reset_clock(self%esmComp, self%dt, cdate_start, cdate_final)

end subroutine model_nuopc_create

! ------------------------------------------------------------------------------

subroutine model_nuopc_delete(self)

implicit none
type(model_nuopc_type), intent(inout) :: self

integer :: rc

! Destroy the Earth system Component
call ESMF_GridCompDestroy(self%esmComp, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_LogWrite("esm-JEDI FINISHED", ESMF_LOGMSG_INFO, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Finalize ESMF
call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

end subroutine model_nuopc_delete

! ------------------------------------------------------------------------------

subroutine model_nuopc_initialize(self, state, vdate)

implicit none
type(model_nuopc_type), intent(inout) :: self
type(fv3jedi_state),    intent(in)    :: state
type(datetime),         intent(in)    :: vdate !< Valid datetime after step

integer :: rc

! Call Initialize for the Earth system Component
call ESMF_GridCompInitialize(self%esmComp, userRc=self%urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=self%urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine model_nuopc_initialize

! ------------------------------------------------------------------------------

subroutine model_nuopc_step(self, state, vdate1, vdate2)

implicit none
type(model_nuopc_type), intent(inout) :: self
type(fv3jedi_state),    intent(inout) :: state
type(datetime),         intent(in)    :: vdate1 !< Valid datetime before advance
type(datetime),         intent(in)    :: vdate2 !< Valid datetime after advance

integer :: rc
character(len=20) :: vdatec1, vdatec2

!Convert datetimes
call datetime_to_string(vdate1, vdatec1)
call datetime_to_string(vdate2, vdatec2)

!Reset the GridComp clock for this advance step
call nuopc_reset_clock(self%esmComp, self%dt, vdatec1, vdatec2)

! JEDI state to NUOPC model state
call state_to_nuopc( state, self )

! Call Run  for Earth the system Component
call ESMF_GridCompRun(self%esmComp, userRc=self%urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=self%urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! JEDI state to NUOPC model state to JEDI state
call nuopc_to_state( self, state )

end subroutine model_nuopc_step

! ------------------------------------------------------------------------------

subroutine model_nuopc_finalize(self, state, vdate)

implicit none
type(model_nuopc_type), intent(inout) :: self
type(fv3jedi_state),    intent(in)    :: state
type(datetime),         intent(in)    :: vdate !< Valid datetime after step

integer :: rc

! Call Finalize for the Earth system Component
call ESMF_GridCompFinalize(self%esmComp, userRc=self%urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=self%urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine model_nuopc_finalize

! ------------------------------------------------------------------------------

subroutine state_to_nuopc( state, self )

implicit none
type(fv3jedi_state),    intent(in)    :: state
type(model_nuopc_type), intent(inout) :: self

type(ESMF_GridComp) :: atmComp
type(ESMF_Clock)    :: clock
type(ESMF_State)    :: importState, exportState
integer :: rc
type (ESMF_FieldBundle) :: bundle
type (ESMF_Field) :: field

real(kind_real), pointer :: sst(:,:), pmsl(:,:)
type(ESMF_FieldStatus_Flag) :: fieldStatus

integer :: LB(2), UB(2)

integer, save :: init = 0

!Get the atmosphere child component
call NUOPC_DriverGetComp(self%esmComp, "ATM", comp=atmComp, rc=rc)

!Get atmosphere import state
call NUOPC_ModelGet(atmComp, modelClock=clock, importState=importState, rc=rc)

!Get field from import state
call ESMF_StateGet(importState, "sst", field, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

!Get pointer to the air_pressure_at_sea_level field
call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_FieldGet(field, 0, sst, computationalLBound=LB, computationalUBound=UB, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

print*, "JED-DATA sst", sst(LB(1),LB(2))

if (associated(sst)) nullify(sst)
if (associated(pmsl)) nullify(pmsl)


end subroutine state_to_nuopc

! ------------------------------------------------------------------------------

subroutine nuopc_to_state( self, state )

implicit none
type(model_nuopc_type), intent(in)    :: self
type(fv3jedi_state),    intent(inout) :: state

type(ESMF_GridComp) :: atmComp
type(ESMF_Clock)    :: clock
type(ESMF_State)    :: importState, exportState
integer :: rc

!Get the atmosphere child component
call NUOPC_DriverGetComp(self%esmComp, "ATM", comp=atmComp, rc=rc)

end subroutine nuopc_to_state

! ------------------------------------------------------------------------------

subroutine nuopc_reset_clock( driver, dt, cdate_start, cdate_final)

implicit none
type(ESMF_GridComp), intent(inout)  :: driver
integer,             intent(in)     :: dt
character(len=20),   intent(in)     :: cdate_start
character(len=20),   intent(in)     :: cdate_final

type(ESMF_Time)               :: startTime
type(ESMF_Time)               :: stopTime
type(ESMF_TimeInterval)       :: timeStep
type(ESMF_Clock)              :: externalClock

integer :: rc

integer yy1,mm1,dd1,hh1,mn1
integer yy2,mm2,dd2,hh2,mn2

!Convert character dates to integers
read(cdate_start(1:4),'(i)') yy1
read(cdate_final(1:4),'(i)') yy2
read(cdate_start(6:7),'(i)') mm1
read(cdate_final(6:7),'(i)') mm2
read(cdate_start(9:10),'(i)') dd1
read(cdate_final(9:10),'(i)') dd2
read(cdate_start(12:13),'(i)') hh1
read(cdate_final(12:13),'(i)') hh2
read(cdate_start(15:16),'(i)') mn1
read(cdate_final(15:16),'(i)') mn2

call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_TimeSet(startTime, yy=yy1, mm=mm1, dd=dd1, h=hh1, m=mn1, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_TimeSet(stopTime, yy=yy2, mm=mm2, dd=dd2, h=hh2, m=mn2, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

externalClock = ESMF_ClockCreate(name="Application Clock", &
  timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call NUOPC_CompSetClock(driver, externalClock=externalClock, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine nuopc_reset_clock

! ------------------------------------------------------------------------------

end module fv3jedi_nuopc_mod
