! (C) Copyright 2020 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#define esmf_err_abort(rc) if (esmf_LogFoundError(rc, msg="Aborting UFS", line=__LINE__, file=__FILE__)) call esmf_Finalize(endflag=esmf_END_ABORT)

module fv3jedi_ufs_mod

  ! oops
  use datetime_mod
  use duration_mod

  ! fckit
  use fckit_configuration_module, only: fckit_configuration

  ! fv3jedi
  use fv3jedi_geom_mod,      only: fv3jedi_geom
  use fv3jedi_state_mod,     only: fv3jedi_state
  use fv3jedi_field_mod,     only: fv3jedi_field, field_clen

  ! ufs
  use ESMF
  use NUOPC
  use NUOPC_Driver
  use module_EARTH_GRID_COMP, only: esmSS => EARTH_REGISTER
  use mpp_mod,            only: read_input_nml,mpp_pe


  implicit none
  private

  public :: model_ufs

  !> Fortran derived type to hold model definition
  type :: model_ufs
     type(ESMF_GridComp) :: esmComp
     type(ESMF_State) :: toJedi, fromJedi
     integer :: isc, iec, jsc, jec, npz
     type(esmf_Clock) :: clock
     type(esmf_config) :: cf_main                                         !<-- the configure object
   contains
     procedure :: create
     procedure :: delete
     procedure :: initialize
     procedure :: step
     procedure :: finalize
  end type model_ufs

  character(len=*), parameter :: modname='fv3jedi_ufs_mod'

  ! --------------------------------------------------------------------------------------------------

contains

  ! --------------------------------------------------------------------------------------------------

  subroutine create(self, conf, geom)

    implicit none
    class(model_ufs),          intent(inout) :: self
    type(fckit_configuration), intent(in)    :: conf
    type(fv3jedi_geom),        intent(in)    :: geom

    integer :: rc, urc, phase, i, cnt
    character(len=20) :: cdate_start, cdate_stop

    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep

    character(len=*),parameter :: subname = modname//' (create)'
    type(ESMF_CplComp),  pointer       :: connectors(:)
    character(len=128) :: name, msg

    ! Initialize ESMF
    call ESMF_Initialize(logkindflag=esmf_LOGKIND_MULTI, &
         defaultCalkind=esmf_CALKIND_GREGORIAN, &
         mpiCommunicator=geom%f_comm%communicator(), rc=rc)
    esmf_err_abort(rc)
    ! Flush log output while debugging
    call ESMF_LogSet(flush=.true., rc=rc)
    esmf_err_abort(rc)
    call ESMF_LogWrite("ESMF Initialized in "//subname, ESMF_LOGMSG_INFO)

    self%isc = geom%isc
    self%iec = geom%iec
    self%jsc = geom%jsc
    self%jec = geom%jec
    self%npz = geom%npz

    self%cf_main=esmf_configcreate(rc=rc)
    call ESMF_ConfigLoadFile(config=self%cf_main, &
         filename='model_configure', &
         rc=rc)

    ! This call to read_input_nml() seems to be required
    ! for CCPP.  However, it does not belong at this level
    ! but should be handled inside the model itself
    call read_input_nml()
    call ESMF_LogWrite("done reading input nml", ESMF_LOGMSG_INFO)

    ! Create the ESM component
    self%esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
    esmf_err_abort(rc)

    ! SetServices for the ESM component
    call ESMF_GridCompSetServices(self%esmComp, esmSS, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)


    ! Set ESM's Verbosity (High)  - 32513
    call NUOPC_CompAttributeSet(self%esmComp, name="Verbosity", &
         value="32513", rc=rc)
    esmf_err_abort(rc)



    ! Initialize the clock based on contents of model_configure
    ! -------------------------------------------
    call setUFSClock(self,startTime,stopTime)
    call ESMF_GridCompSet(self%esmComp, clock=self%clock, rc=rc)

    ! Create import and export states from
    ! perspective of the exernal system:
    !   toJedi is an IMPORT into Jedi and an EXPORT from ESM
    !   fromJedi is an EXPORT from Jedi and an IMPORT into ESM

    self%toJedi = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_IMPORT, &
         rc=rc)
    esmf_err_abort(rc)


    self%fromJedi = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_EXPORT, &
         rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("Advertising export from ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the exportState, for data coming out of ESM component
    ! Note--only certain fields are available. Check in GFS_surface_generic to see if they are filled
    call NUOPC_Advertise(self%toJedi, &
         StandardNames=(/ &
                        "u                                    ", &   ! Example fields
                        "v                                    ", &   ! Example fields
                        "ua                                   ", &   ! Example fields
                        "va                                   ", &   ! Example fields
                        "t                                    ", &   ! Example fields
                        "delp                                 ", &   ! Example fields
                        "sphum                                ", &   ! Example fields
                        "ice_wat                              ", &   ! Example fields
                        "liq_wat                              ", &   ! Example fields
                        "o3mr                                 ", &   ! Example fields
                        "phis                                 ", &   ! Example fields
                        "slmsk                                ", &   ! Example fields
                        "weasd                                ", &   ! Example fields
                        "tsea                                 ", &   ! Example fields
                        "vtype                                ", &   ! Example fields
                        "stype                                ", &   ! Example fields
                        "vfrac                                ", &   ! Example fields
                        "stc                                  ", &   ! Example fields
                        "smc                                  ", &   ! Example fields
                        "snwdph                               ", &   ! Example fields
                        "u_srf                                ", &   ! Example fields
                        "v_srf                                ", &   ! Example fields
                        "f10m                                 "/), &   ! Example fields
         SharePolicyField="share", &
         TransferOfferGeomObject="cannot provide", rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("Advertising imports to ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the importState, for data going into ESM componenta

! imports are not yet implemented
#if 0
    call NUOPC_Advertise(self%fromJedi, &
         StandardNames=(/ &
                        "u                                    ", &   ! Example fields
                        "v                                    ", &   ! Example fields
                        "ua                                   ", &   ! Example fields
                        "va                                   ", &   ! Example fields
                        "t                                    ", &   ! Example fields
                        "delp                                 ", &   ! Example fields
                        "sphum                                ", &   ! Example fields
                        "ice_wat                              ", &   ! Example fields
                        "liq_wat                              ", &   ! Example fields
                        "o3mr                                 ", &   ! Example fields
                        "phis                                 ", &   ! Example fields
                        "slmsk                                ", &   ! Example fields
                        "weasd                                ", &   ! Example fields
                        "tsea                                 ", &   ! Example fields
                        "vtype                                ", &   ! Example fields
                        "stype                                ", &   ! Example fields
                        "vfrac                                ", &   ! Example fields
                        "stc                                  ", &   ! Example fields
                        "smc                                  ", &   ! Example fields
                        "snwdph                               ", &   ! Example fields
                        "u_srf                                ", &   ! Example fields
                        "v_srf                                ", &   ! Example fields
                        "f10m                                 "/), &   ! Example fields
         TransferOfferGeomObject="cannot provide", rc=rc)

    esmf_err_abort(rc)
#endif

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After filling advertise toJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)


    ! call ExternalAdvertise phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalAdvertise, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)


    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)


    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After calling advertise toJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)


    ! Set verbosity flag on connectors
    nullify(connectors);
    call NUOPC_DriverGetComp(self%esmComp, &
         compList=connectors, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to set connector verbosity", ESMF_LOGMSG_INFO)
    do i=lbound(connectors,1), ubound(connectors,1)
       call ESMF_CplCompGet(connectors(i), name=name, rc=rc)
       esmf_err_abort(rc)

       call NUOPC_CompAttributeSet(connectors(i), name="Verbosity", &
            value="max", rc=rc)
       esmf_err_abort(rc)

       call ESMF_LogWrite(" --> Set verbosity on connector: "//trim(name), &
            ESMF_LOGMSG_INFO)
    enddo

    deallocate(connectors)

    ! call ExternalRealize phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalRealize, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)

    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt

    call ESMF_LogWrite("Dumping toJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)

    ! call ExternalDataInit phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalDataInit, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)

    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine create

! --------------------------------------------------------------------------------------------------

  subroutine initialize(self, state)

    implicit none

    class(model_ufs),    intent(inout) :: self
    type(fv3jedi_state), intent(in)    :: state

    character(len=*),parameter :: subname = modname//' (initialize)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine initialize

! --------------------------------------------------------------------------------------------------

  subroutine step(self, state, vdate_start, vdate_final)

    implicit none

    class(model_ufs),    intent(inout) :: self
    type(fv3jedi_state), intent(inout) :: state
    type(datetime),      intent(in)    :: vdate_start
    type(datetime),      intent(in)    :: vdate_final

    ! local variables
    integer :: rc, urc, cnt
    character(len=20) :: strStartTime, strStopTime
    character(len=128) :: name, msg
    type(ESMF_Time) :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer, save     :: tstep=1
    character(len=80) :: fileName

!-----------------------------------------------------------------------------

    character(len=*),parameter :: subname = modname//' (step)'
    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call datetime_to_string(vdate_start, strStartTime)
    call datetime_to_string(vdate_final, strStopTime)

    call ESMF_LogWrite(" --> REQUESTED START TIME:"//trim(strStartTime), ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(" --> REQUESTED STOP  TIME:"//trim(strStopTime), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(self%clock, startTime=startTime, &
         stopTime=stopTime, rc=rc)
    esmf_err_abort(rc)


    call ESMF_TimeSet(startTime, timeString=strStartTime, rc=rc)
    esmf_err_abort(rc)


    call ESMF_TimeSet(stopTime, timeString=strStopTime, rc=rc)
    esmf_err_abort(rc)


    timeStep = stopTime - startTime

    call ESMF_ClockSet(self%clock, startTime=startTime, &
         stopTime=stopTime, currTime=startTime, timeStep=timeStep, rc=rc)
     esmf_err_abort(rc)

    ! step the model forward
    call ESMF_GridCompRun(self%esmComp, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)
    write(msg, "(I2)") cnt
    call ESMF_LogWrite("after step toJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)
    write(fileName, '("fields_in_esm_import_step",I2.2,".nc")') tstep
    call fv3_to_state(self, state)
    call ESMF_LogWrite("after state write "//trim(msg)//" rc", &
         ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine step

! --------------------------------------------------------------------------------------------------

  subroutine delete(self)

    implicit none
    class(model_ufs), intent(inout) :: self
    integer :: rc
    character(len=*),parameter :: subname = modname//' (delete)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_GridCompDestroy(self%esmComp, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to destroy toJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%toJedi, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to destroy fromJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%fromJedi, rc=rc)
    esmf_err_abort(rc)

    ! Finalize ESMF
    ! -------------
    call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=rc)
    if (rc /= ESMF_SUCCESS) then
       call ESMF_LogWrite("ERROR FINALIZING ESMF "//subname, ESMF_LOGMSG_INFO)
    else
       call ESMF_LogWrite("SUCCESSFULLY FINALIZED ESMF"//subname, ESMF_LOGMSG_INFO)
    endif
  end subroutine delete

  ! --------------------------------------------------------------------------------------------------

  subroutine finalize(self, state)

    implicit none
    class(model_ufs),    intent(inout) :: self
    type(fv3jedi_state), intent(in)    :: state

    character(len=*),parameter :: subname = modname//' (finalize)'
    ! Clean up is being done in the delete method
    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine finalize

  subroutine fv3_to_state( self, state )

  implicit none
  type(model_ufs),    intent(in)    :: self
  type(fv3jedi_state), intent(inout) :: state

  integer :: num_items, i, rc, rank, lb(3), ub(3), fnpz
  type(ESMF_Field) :: field
  character(len=ESMF_MAXSTR), allocatable :: item_names(:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr2(:,:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr3(:,:,:)
  character(len=field_clen) :: short_name
  type(fv3jedi_field), pointer :: field_ptr

  real(kind=ESMF_KIND_R8),allocatable,dimension(:,:,:)      :: field_fv3


  ! Array to hold output from UFS in JEDI precision
  ! ------------------------------------------------
  allocate(field_fv3(self%isc:self%iec, self%jsc:self%jec, self%npz+1))


  ! Get number of items
  ! -------------------
  call ESMF_StateGet(self%toJedi, itemcount = num_items, rc = rc)
  if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_StateGet itemcount failed")


  ! Get names of the items
  ! ----------------------
  allocate(item_names(num_items))
  call ESMF_StateGet(self%toJedi, itemnamelist = item_names, rc = rc)
  if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_StateGet itemnamelist failed")


  ! Loop over states coming from UFS and convert to JEDI state
  ! -----------------------------------------------------------
  do i = 1, num_items
    ! Create map between UFS name and fv3-jedi name
    ! ----------------------------------------------
    short_name = trim(item_names(i))
    call ESMF_LogWrite("item name is "//short_name, ESMF_LOGMSG_INFO)
    if(trim(item_names(i)) == 'u') short_name = 'ud'
    if(trim(item_names(i)) == 'v') short_name = 'vd'
    if(trim(item_names(i)) == 'weasd') short_name = 'sheleg'

    ! Only need to extract field from UFS if fv3-jedi needs it
    ! ---------------------------------------------------------
    if (state%has_field(trim(short_name))) then

      !Get field from the state
      call ESMF_StateGet(self%toJedi, item_names(i), field, rc = rc)
      if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_StateGet field failed")

      !Validate the field
      call ESMF_FieldValidate(field, rc = rc)
      if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldValidate failed")

      !Get the field rank
      call ESMF_FieldGet(field, rank = rank, rc = rc)

      if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldGet rank failed")

      !Convert field to pointer and pointer bounds
      field_fv3 = 0.0_ESMF_KIND_R8
      if (rank == 2) then

        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb(1:2), totalUBound = ub(1:2), rc = rc )
        if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldGet 2D failed")

        fnpz = 1
        field_fv3(self%isc:self%iec,self%jsc:self%jec,1) = farrayPtr2(lb(1):ub(1),lb(2):ub(2))
        nullify(farrayPtr2)

      elseif (rank == 3) then
        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )
        if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldGet 3D failed",rc)

        fnpz = ub(3)-lb(3)+1
        field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))
        nullify(farrayPtr3)

      else

        call abor1_ftn("fv3_mod: can only handle rank 2 or rank 3 fields from UFS")

      endif

      ! Check that dimensions match
      if ((ub(1)-lb(1)+1 .ne. self%iec-self%isc+1) .or. (ub(2)-lb(2)+1 .ne. self%jec-self%jsc+1) ) then
        call abor1_ftn("fv3_to_state: dimension mismatch between JEDI and UFS horizontal grid")
      endif

      ! Get pointer to fv3-jedi side field
      call state%get_field(trim(short_name), field_ptr)

      if (field_ptr%npz .ne. fnpz) &
        call abor1_ftn("fv3_to_state: dimension mismatch between JEDI and UFS vertical grid")

      ! Copy from UFS to fv3-jedi
      field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz)
    else
      call ESMF_LogWrite("Not needed by JEDI is "//short_name, ESMF_LOGMSG_INFO)
    endif

  end do

  deallocate(item_names)

  end subroutine fv3_to_state


  subroutine setUFSClock(self,date_start,date_final)

    class(model_ufs),    intent(inout) :: self
    type(esmf_time),     intent(inout) :: date_start, date_final

    type(ESMF_TimeInterval)            :: runDuration, timestep
    integer :: yy,mm,dd,hh,mns,sec,rc
    !-----------------------------------------------------------------------
    !***  extract the start time from the configuration
    !-----------------------------------------------------------------------
    !
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =YY                            &
                            ,label ='start_year:'                 &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =MM                            &
                            ,label ='start_month:'                &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =DD                            &
                            ,label ='start_day:'                  &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =HH                            &
                            ,label ='start_hour:'                 &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =MNS                           &
                            ,label ='start_minute:'               &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =SEC                           &
                            ,label ='start_second:'               &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    ! Set date_start
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    call ESMF_TimeSet(time=date_start                                 &  !<-- The start time of the forecast (ESMF)
                 ,yy  =YY                                         &  !<-- Year from config file
                 ,mm  =MM                                         &  !<-- Month from config file
                 ,dd  =DD                                         &  !<-- Day from config file
                 ,h   =HH                                         &  !<-- Hour from config file
                 ,m   =MNS                                        &  !<-- Minute from config file
                 ,s   =SEC                                        &  !<-- Second from config file
                ,rc  =RC)
    esmf_err_abort(rc)

    ! Set any time interval here It will be overridden later
    call ESMF_timeintervalset(runduration, d=1, rc=rc)
    call ESMF_timeintervalset(timestep, h=1, rc=rc)
    date_final = date_start + runduration

    self%clock = ESMF_ClockCreate(name="main_clock", &
         timeStep=timestep, startTime=date_start, stopTime=date_final, rc=rc)

    esmf_err_abort(rc)

    call ESMF_ClockPrint(self%clock, options="startTime", &
    preString="Printing startTime to stdout: ", rc=rc)


    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  end subroutine setUFSClock

end module fv3jedi_ufs_mod
