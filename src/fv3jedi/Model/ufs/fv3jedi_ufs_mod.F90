! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#define esmf_err_abort(rc) if (esmf_LogFoundError(rc, msg="Aborting NEMS", line=__LINE__, file=__FILE__)) call esmf_Finalize(endflag=esmf_END_ABORT)

module fv3jedi_ufs_mod

! oops
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! fv3jedi
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state

! ufs
use ESMF
use NUOPC, only: NUOPC_Advertise, NUOPC_Write, NUOPC_SetTimestamp
use NUOPC_Comp, only: NUOPC_CompSetClock, NUOPC_CompSearchPhaseMap, &
                      NUOPC_CompAttributeSet
use NUOPC_Driver, only: NUOPC_DriverGetComp
use module_EARTH_GRID_COMP, only: esmSS => EARTH_REGISTER
use module_nems_utils, only: message_check
use mpp_mod,            only: read_input_nml,mpp_pe
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
  character(len=20) :: startT                                         !<-- the esmf start time.
  type(esmf_Clock) :: clock 
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

implicit none
class(model_ufs),          intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf
type(fv3jedi_geom),        intent(in)    :: geom

integer :: rc,urc,phase
type(duration) :: dtstep
character(len=20) :: ststep
character(len=20) :: cdate_start, cdate_final
character(240)        :: msgString

! Initialize ESMF
call esmf_initialize(logkindflag=esmf_LOGKIND_MULTI, defaultCalkind=esmf_CALKIND_GREGORIAN, &
                     mpiCommunicator=geom%f_comm%communicator(), rc=rc)
esmf_err_abort(rc)

call ESMF_LogWrite("JEDI control of ESM STARTING", ESMF_LOGMSG_INFO, rc=rc)
esmf_err_abort(rc)

call esmf_logwrite("getting config", esmf_LOGMSG_INFO, rc=rc)
self%cf_main=esmf_configcreate(rc=rc)
!!
call esmf_configloadfile(config  =self%cf_main                         &  !<-- the configure object
                       ,filename='model_configure'               &  !<-- the name of the configure file
                       ,rc      =rc)
call esmf_logwrite("done getting config", esmf_LOGMSG_INFO, rc=rc)
esmf_err_abort(rc)

call read_input_nml
call esmf_logwrite("done reading input nml", esmf_LOGMSG_INFO, rc=rc)


! Create the ESM component
call esmf_logwrite("creating grid comp", esmf_LOGMSG_INFO, rc=rc)
self%esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
call esmf_logwrite("done creating grid comp", esmf_LOGMSG_INFO, rc=rc)
esmf_err_abort(rc)

! SetServices for the ESM component
call esmf_logwrite("calling setservices", esmf_LOGMSG_INFO, rc=rc)
call ESMF_GridCompSetServices(self%esmComp, esmSS, userRc=urc, rc=rc)
esmf_err_abort(rc)

  
! Set ESM's Verbosity
call esmf_logwrite("calling att set", esmf_LOGMSG_INFO, rc=rc)
call NUOPC_CompAttributeSet(self%esmComp, name="Verbosity", value="low", rc=rc)
esmf_err_abort(rc)

! Reset the clock based on what JEDI provides
! -------------------------------------------

! Hard-coded settings for the protos
self%dt=3600  ! 1h steps
cdate_start = "2019-08-29T00:00:00Z"
cdate_final = "2019-08-29T04:00:00Z" ! end date not critical, just be in future



call esmf_logwrite("constructing clock", esmf_LOGMSG_INFO, rc=rc)
call construct_clock(self%dt, cdate_start, cdate_final, clock=self%clock)
! create import- and exportState to be used as conduits in/out NUOPC ESM
call esmf_logwrite("creating states", esmf_LOGMSG_INFO, rc=rc)
self%importState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_IMPORT, &
  rc=rc)
esmf_err_abort(rc)

self%exportState = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_EXPORT, &
  rc=rc)
esmf_err_abort(rc)

! Advertise fields on the exportState, for data coming out of ESM component
call esmf_logwrite("calling advertise", esmf_LOGMSG_INFO, rc=rc)
call NUOPC_Advertise(self%exportState, &
  ! Test with a couple of 2D surface fields
  StandardNames=(/"inst_down_lw_flx              ", &
                  "inst_down_sw_flx              "/), &
  TransferOfferGeomObject="cannot provide", rc=rc)
  esmf_err_abort(rc)

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
  esmf_err_abort(rc)
#endif

! call ExternalAdvertise phase
call esmf_logwrite("calling search phase map", esmf_LOGMSG_INFO, rc=rc)
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalAdvertise", phaseIndex=phase, rc=rc)
esmf_err_abort(rc)

call esmf_logwrite("calling gc init", esmf_LOGMSG_INFO, rc=rc)
call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
call ESMF_GridCompSet(self%esmComp, clock=self%clock, rc=rc)

call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
esmf_err_abort(rc)

! call ExternalRealize phase
call esmf_logwrite("calling ext phase search", esmf_LOGMSG_INFO, rc=rc)
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalRealize", phaseIndex=phase, rc=rc)
esmf_err_abort(rc)

call esmf_logwrite("done calling ext phase search", esmf_LOGMSG_INFO, rc=rc)
! call ExternalDataInitialize phase
call NUOPC_CompSearchPhaseMap(self%esmComp, &
  methodflag=ESMF_METHOD_INITIALIZE, &
  phaseLabel="ExternalDataInitialize", phaseIndex=phase, rc=rc)
esmf_err_abort(rc)
call esmf_logwrite("done calling extdata phase search", esmf_LOGMSG_INFO, rc=rc)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(model_ufs), intent(inout) :: self

integer :: rc

! Finalize ESMF
! -------------
call esmf_Finalize(endflag=esmf_END_KEEPMPI, rc=rc)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state, vdate)

implicit none
include 'mpif.h'

class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state
type(datetime),      intent(in)    :: vdate

type(esmf_timeinterval) :: runduration                            &  !<-- the esmf time. the total forecast hours.
                         ,timestep                                  !<-- the esmf timestep length (we only need a dummy here)
integer :: rc, rc_user                                               !<-- the running error signal
!
character(len=20) :: vdatec
integer, save     :: tstep=1
character(len=80) :: fileName

!Convert datetimes
call datetime_to_string(vdate, vdatec)
self%mype = mpp_pe()
self%startT = vdatec 
call esmf_logwrite("in initialize", esmf_LOGMSG_INFO, rc=rc)

! for testing, write out the fields in the exportState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_export_init",I2.2,".nc")') tstep
call FV3_StateWrite(self%exportState, fileName=trim(fileName), rc=rc)
esmf_err_abort(rc)


! for testing, write out the fields in the importState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_import_init",I2.2,".nc")') tstep
call FV3_StateWrite(self%importState, fileName=trim(fileName), rc=rc)
esmf_err_abort(rc)


tstep = tstep+1
call esmf_logwrite("leaving initialize", esmf_LOGMSG_INFO, rc=rc)

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, vdate_start, vdate_final)

implicit none
include 'mpif.h'

class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate_start
type(datetime),      intent(in)    :: vdate_final
integer                                :: rc
    
! local variables
type(ESMF_Time)                        :: currTime
type(ESMF_TimeInterval)                :: timeStep
type(ESMF_Time)                        :: startTime, stopTime
type(ESMF_TimeInterval)                :: time_elapsed,RunDuration
integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
!
integer :: na, i, urc
logical :: isAlarmEnabled, isAlarmRinging, lalarm, reconcileFlag
character(len=20)  :: cdate_start, cdate_final
character(240)              :: msgString
character(240)              :: import_timestr, export_timestr
character(len=160) :: nuopcMsg
integer :: mype,date(6), fieldcount, fcst_nfld
real(kind=8)   :: timeri, timewri, timewr, timerhi, timerh
character(len=80) :: fileName
character(len=20) :: vdatec1, vdatec2
integer, save     :: tstep=1


!-----------------------------------------------------------------------------

mype = mpp_pe()
rc = ESMF_SUCCESS
!Convert JEDI state to model state
call state_to_nems( state, self )

!Convert datetimes
call datetime_to_string(vdate_start, vdatec1)
call datetime_to_string(vdate_final, vdatec2)

call esmf_logwrite("starting step", esmf_LOGMSG_INFO, rc=rc)
!Reset the GridComp clock for this advance step
call construct_clock(self%dt, vdatec1, vdatec2, clock=self%clock)

! timestamp the data going into the ESM or else NUOPC will flag incompatible
call NUOPC_SetTimestamp(self%importState, self%clock, rc=rc)
esmf_err_abort(rc)

! Step the ESM forward from vdate1 -> vdate2
call ESMF_GridCompRun(self%esmComp, &
  importState=self%importState, exportState=self%exportState, &
  clock=self%clock, userRc=urc, rc=rc)
esmf_err_abort(rc)


! for testing, write out the fields in the exportState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_export_step",I2.2,".nc")') tstep
call FV3_StateWrite(self%exportState, fileName=trim(fileName), rc=rc)
esmf_err_abort(rc)

! for testing, write out the fields in the importState to file
!TODO: only works for 2D surface fields
write(fileName, '("fields_in_esm_import_step",I2.2,".nc")') tstep
call FV3_StateWrite(self%importState, fileName=trim(fileName), rc=rc)
esmf_err_abort(rc)

tstep = tstep+1

!Convert model state to JEDI state
call nems_to_state( self, state )
!
end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine destructor(self)
implicit none
class(model_ufs),    intent(inout) :: self
integer :: rc, phase
integer :: rc_user                                               !<-- the running error signal
message_check="execute the nems component finalize step"
call esmf_logwrite(message_check,esmf_logmsg_info,rc=rc)
!!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!
call esmf_gridcompfinalize(gridcomp   =self%esmComp&  !<-- the nems component
 ,importstate=self%importState &  !<-- the nems component import state
 ,exportstate=self%exportState &  !<-- the nems component export state
 ,clock=self%clock     &  !<-- the main esmf clock
 ,phase=1  &
 ,userrc     =rc_user  &
 ,rc   =rc)
esmf_err_abort(rc)
esmf_err_abort(rc_user)
end subroutine destructor

subroutine finalize(self, state, vdate)

implicit none
class(model_ufs),    intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state
type(datetime),intent(in)    :: vdate


end subroutine finalize

! --------------------------------------------------------------------------------------------------
subroutine construct_clock(dt, cdate_start, cdate_final, clock)

implicit none
integer,            intent(in)  :: dt
character(len=20),  intent(in)  :: cdate_start
character(len=20),  intent(in)  :: cdate_final
type(ESMF_Clock),   intent(out) :: clock

type(ESMF_Time)               :: startTime
type(ESMF_Time)               :: stopTime
type(ESMF_TimeInterval)       :: timeStep
type(ESMF_TimeInterval)       :: RunDuration
integer :: rc

call construct_time(cdate_start, time=startTime)
call construct_time(cdate_final, time=stopTime)
runduration = stopTime - startTime 
esmf_err_abort(rc)

clock = ESMF_ClockCreate(name="main_clock", &
  timeStep=runduration, startTime=startTime, stopTime=stopTime, rc=rc)
esmf_err_abort(rc)

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
esmf_err_abort(rc)

end subroutine construct_time

subroutine FV3_StateWrite(state, fileName, rc)
  type(ESMF_State)      :: state
  character(len=*)      :: fileName
  integer, intent(out)  :: rc
  
  integer               :: itemCount, i, mype
  type(ESMF_Field), allocatable  :: fieldList(:)
  character(len=80), allocatable :: itemNameList(:)
  type(ESMF_Grid)     :: grid
  type(ESMF_GridComp) :: ioComp
  integer:: mpic
  
  rc=ESMF_SUCCESS
  
  call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
  esmf_err_abort(rc)
  if (itemCount==0) return
  
# if 0
  allocate(fieldList(itemCount), itemNameList(itemCount))
  call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
  esmf_err_abort(rc)
    call ESMF_StateGet(state, itemName=itemNameList(i), field=fieldList(i), &
      rc=rc)
    esmf_err_abort(rc)
#endif
   call ESMF_StateWrite(state, 'fcstState', rc=rc)
   esmf_err_abort(rc)

! deallocate(fieldList, itemNameList)
 
  
end subroutine FV3_StateWrite
subroutine state_to_nems( state, self )

implicit none
type(fv3jedi_state), intent(in)    :: state
type(model_ufs),    intent(inout) :: self

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
type(model_ufs),    intent(in)    :: self
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
end module fv3jedi_ufs_mod
