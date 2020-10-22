! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_nems_mod

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
use NUOPC, only: NUOPC_Advertise, NUOPC_Write, NUOPC_SetTimestamp
use NUOPC_Comp, only: NUOPC_CompSetClock, NUOPC_CompSearchPhaseMap, &
                      NUOPC_CompAttributeSet
use NUOPC_Driver, only: NUOPC_DriverGetComp
use module_EARTH_GRID_COMP, only: esmSS => EARTH_REGISTER

implicit none
private

public :: nems_model
public :: nems_create
public :: nems_delete
public :: nems_initialize
public :: nems_step
public :: nems_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: nems_model
 type(ESMF_GridComp) :: esmComp !NUOPC driver
 type(ESMF_State)    :: importState, exportState
 type(ESMF_Clock)    :: clock
 integer :: dt
 character(len=20) :: startTime
end type nems_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine nems_create(self, geom, c_conf)

implicit none
type(nems_model),   intent(inout) :: self
type(fv3jedi_geom), intent(in)    :: geom
type(c_ptr),        intent(in)    :: c_conf

integer :: rc, urc, phase
type(duration) :: dtstep
character(len=20) :: ststep
character(len=20) :: cdate_start, cdate_final
type(fckit_mpi_comm) :: f_comm

! FCKIT MPI wrapper for communicator
f_comm = fckit_mpi_comm()

! Initialize ESMF
call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
  defaultCalkind=ESMF_CALKIND_GREGORIAN, &
  mpiCommunicator=f_comm%communicator(), &
  rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_LogWrite("JEDI control of ESM STARTING", ESMF_LOGMSG_INFO, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Create the ESM component
self%esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! SetServices for the ESM component
call ESMF_GridCompSetServices(self%esmComp, esmSS, userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
! Set ESM's Verbosity
call NUOPC_CompAttributeSet(self%esmComp, name="Verbosity", value="low", rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Reset the clock based on what JEDI provides
! -------------------------------------------
#ifdef JEDIDUMMYAPP

self%dt=8*3600  ! 8h steps
! Hard-coded settings for the protos
cdate_start = "2018-04-15T00:00:00Z"
cdate_final = "2019-04-15T00:00:00Z" ! end date not critical, just be in future

#else

ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
self%dt = int(duration_seconds(dtstep))
! User provides model start and end clocks
cdate_start = config_get_string(c_conf,len(cdate_start),"ModelStartTime")
cdate_final = config_get_string(c_conf,len(cdate_final),"ModelFinalTime")

#endif


call construct_clock(self%dt, cdate_start, cdate_final, clock=self%clock)

! create import- and exportState to be used as conduits in/out NUOPC ESM
self%importState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_IMPORT, &
  rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
self%exportState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_EXPORT, &
  rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

#if 1
! Advertise fields on the exportState, for data coming out of ESM component
call NUOPC_Advertise(self%exportState, &
#if 0
  ! Test with FV3 3D fields
  StandardNames=(/"inst_temp_levels              ", &
                  "inst_zonal_wind_levels_Dgrid  ", &
                  "inst_merid_wind_levels_Dgrid  ", &
                  "inst_spec_humid_levels        ", &
                  "inst_pres_thickness_levels    ", &
                  "inst_liq_water_levels         ", &
                  "inst_ice_water_levels         ", &
                  "inst_ozone_levels             "/), &
#else
  ! Test with a couple of 2D surface fields
  StandardNames=(/"inst_down_lw_flx              ", &
                  "inst_down_sw_flx              "/), &
#endif
  TransferOfferGeomObject="cannot provide", rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

#if 0
! Advertise fields on the importState, for data going into of ESM component
call NUOPC_Advertise(self%importState, &
  ! Test with FV3 3D fields
  StandardNames=(/"inst_temp_levels              ", &
                  "inst_zonal_wind_levels_Dgrid  ", &
                  "inst_merid_wind_levels_Dgrid  ", &
                  "inst_spec_humid_levels        ", &
                  "inst_pres_thickness_levels    ", &
                  "inst_liq_water_levels         ", &
                  "inst_ice_water_levels         ", &
                  "inst_ozone_levels             "/), &
  TransferOfferGeomObject="cannot provide", rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

! call ExternalAdvertise phase
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalAdvertise", phaseIndex=phase, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! call ExternalRealize phase
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalRealize", phaseIndex=phase, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine nems_create

! ------------------------------------------------------------------------------

subroutine nems_delete(self)

implicit none
type(nems_model), intent(inout) :: self

integer :: rc, urc

! completely take down the ESM system
call ESMF_GridCompFinalize(self%esmComp, &
  importState=self%importState, exportState=self%exportState, &
  userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Destroy the Earth system Component
call ESMF_GridCompDestroy(self%esmComp, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

call ESMF_LogWrite("JEDI control of ESM FINISHED", ESMF_LOGMSG_INFO, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Finalize ESMF
call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

end subroutine nems_delete

! ------------------------------------------------------------------------------

subroutine nems_initialize(self, state, vdate)

implicit none
type(nems_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate

type(ESMF_Time)   :: time
integer           :: rc, urc, phase
character(len=20) :: vdatec
integer, save     :: step=1
character(len=80) :: fileName

!Convert datetimes
call datetime_to_string(vdate, vdatec)

self%startTime = vdatec 

! Set the currTime of the self%clock to new vdate before passing to ESM
call construct_time(vdatec, time=time)
call ESMF_ClockSet(self%clock, currTime=time, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! call ExternalDataInitialize phase
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalDataInitialize", phaseIndex=phase, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

#if 1
! for testing, write out the fields in the exportState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_export_init",I2.2,".nc")') step
call FV3_StateWrite(self%exportState, fileName=trim(fileName), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

#if 1
! for testing, write out the fields in the importState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_import_init",I2.2,".nc")') step
call FV3_StateWrite(self%importState, fileName=trim(fileName), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

step = step+1

end subroutine nems_initialize

! ------------------------------------------------------------------------------

subroutine nems_step(self, state, vdate1, vdate2)

implicit none
type(nems_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate1
type(datetime),      intent(in)    :: vdate2

integer :: rc, urc
character(len=20) :: vdatec1, vdatec2
integer, save     :: step=1
character(len=80) :: fileName

!Convert JEDI state to model state
call state_to_nems( state, self )

!Convert datetimes
call datetime_to_string(vdate1, vdatec1)
call datetime_to_string(vdate2, vdatec2)

!Reset the GridComp clock for this advance step
call construct_clock(self%dt, vdatec1, vdatec2, clock=self%clock)

! timestamp the data going into the ESM or else NUOPC will flag incompatible
call NUOPC_SetTimestamp(self%importState, self%clock, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

! Step the ESM forward from vdate1 -> vdate2
call ESMF_GridCompRun(self%esmComp, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

#if 1
! for testing, write out the fields in the exportState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_export_step",I2.2,".nc")') step
call FV3_StateWrite(self%exportState, fileName=trim(fileName), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

#if 1
! for testing, write out the fields in the importState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_import_step",I2.2,".nc")') step
call FV3_StateWrite(self%importState, fileName=trim(fileName), rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
#endif

step = step+1

!Convert model state to JEDI state
call nems_to_state( self, state )

end subroutine nems_step

! ------------------------------------------------------------------------------

subroutine nems_finalize(self, state, vdate)

implicit none
type(nems_model),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate

integer           :: rc, urc
integer           :: phase
type(ESMF_Time)   :: time

! Finalize the ESM forward stepping by resetting the startTime on self%clock

! Set the currTime of the self%clock to new vdate before passing to ESM
call construct_time(self%startTime, time=time)
call ESMF_ClockSet(self%clock, startTime=time, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
! FinalizeReset the ESM
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_FINALIZE, &
  phaseLabel="ExternalFinalizeReset", phaseIndex=phase, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
call ESMF_GridCompFinalize(self%esmComp, phase=phase, &
  importState=self%importState, exportState=self%exportState, &
  userRc=urc, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)
if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine nems_finalize

! ------------------------------------------------------------------------------

subroutine state_to_nems( state, self )

implicit none
type(fv3jedi_state), intent(in)    :: state
type(nems_model),    intent(inout) :: self

integer :: rc

! import- and exportState, as well as a clock are stored in self
! Do not access them from the model that is under NUOPC ESM, since it may not
! exist across all PETs.


!Get the atmosphere child component
!call NUOPC_DriverGetComp(self%esmComp, "ATM", comp=atmComp, rc=rc)

!Get atmosphere import state
!call NUOPC_ModelGet(atmComp, modelClock=clock, importState=importState, rc=rc)

end subroutine state_to_nems

! ------------------------------------------------------------------------------

subroutine nems_to_state( self, state )

implicit none
type(nems_model),    intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state

integer :: rc

! import- and exportState, as well as a clock are stored in self
! Do not access them from the model that is under NUOPC ESM, since it may not
! exist across all PETs.

!!Get the atmosphere child component
!call NUOPC_DriverGetComp(self%esmComp, "ATM", comp=atmComp, rc=rc)

!Get atmosphere import state
!call NUOPC_ModelGet(atmComp, modelClock=clock, importState=importState, rc=rc)

end subroutine nems_to_state

! ------------------------------------------------------------------------------

subroutine construct_clock(dt, cdate_start, cdate_final, clock)

implicit none
integer,            intent(in)  :: dt
character(len=20),  intent(in)  :: cdate_start
character(len=20),  intent(in)  :: cdate_final
type(ESMF_Clock),   intent(out) :: clock

type(ESMF_Time)               :: startTime
type(ESMF_Time)               :: stopTime
type(ESMF_TimeInterval)       :: timeStep
integer :: rc

call construct_time(cdate_start, time=startTime)
call construct_time(cdate_final, time=stopTime)

call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

clock = ESMF_ClockCreate(name="Application Clock", &
  timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine construct_clock

! ------------------------------------------------------------------------------

subroutine construct_time(cdate, time)

implicit none
character(len=20),  intent(in)  :: cdate
type(ESMF_Time),    intent(out) :: time

integer :: rc
integer :: yy,mm,dd,hh,mn

!Convert character dates to integers
read(cdate(1:4),'(i4)') yy
read(cdate(6:7),'(i2)') mm
read(cdate(9:10),'(i2)') dd
read(cdate(12:13),'(i2)') hh
read(cdate(15:16),'(i2)') mn

call ESMF_TimeSet(time, yy=yy, mm=mm, dd=dd, h=hh, m=mn, rc=rc)
if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
  line=__LINE__, &
  file=__FILE__)) &
  call ESMF_Finalize(endflag=ESMF_END_ABORT)

end subroutine construct_time

! ------------------------------------------------------------------------------

subroutine FV3_StateWrite(state, fileName, rc)
  type(ESMF_State)      :: state
  character(len=*)      :: fileName
  integer, intent(out)  :: rc
  
  integer               :: itemCount, i
  type(ESMF_Field), allocatable  :: fieldList(:)
  character(len=80), allocatable :: itemNameList(:)
  type(ESMF_Grid)     :: grid
  type(ESMF_GridComp) :: ioComp
  
  rc=ESMF_SUCCESS
  
  call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) return
    
  if (itemCount==0) return
  
  allocate(fieldList(itemCount), itemNameList(itemCount))
  call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) return
  do i=1, itemCount
    call ESMF_StateGet(state, itemName=itemNameList(i), field=fieldList(i), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) return
  enddo
  
  call ESMF_FieldGet(fieldList(1), grid=grid, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) return
  
  ioComp = ESMFIO_Create(grid, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out
  call ESMFIO_Write(ioComp, fileName, fieldList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    return  ! bail out

  deallocate(fieldList, itemNameList)
 
  
end subroutine FV3_StateWrite

! ------------------------------------------------------------------------------

end module fv3jedi_nems_mod
