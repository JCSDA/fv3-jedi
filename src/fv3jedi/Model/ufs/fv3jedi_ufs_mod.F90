! (C) Copyright 2020 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#define esmf_err_abort(rc) if (esmf_LogFoundError(rc, msg="Aborting UFS", line=__LINE__, file=__FILE__)) call esmf_Finalize(endflag=esmf_END_ABORT)

module fv3jedi_ufs_mod

  ! oops
  use datetime_mod
  use duration_mod
  use oops_variables_mod

  ! fckit
  use fckit_configuration_module, only: fckit_configuration

  ! fv3jedi
  use fv3jedi_geom_mod,      only: fv3jedi_geom
  use fv3jedi_state_mod,     only: fv3jedi_state
  use fv3jedi_field_mod,     only: fv3jedi_field, field_clen
  use fckit_mpi_module,           only: fckit_mpi_comm
  ! ufs
  use ESMF
  use NUOPC
  use NUOPC_Driver
  use UFSDriver, only: esmSS => UFSDriver_SS
  use mpp_mod,            only: read_input_nml,mpp_pe,mpp_set_current_pelist

  implicit none
  private

  public :: model_ufs

  !> Fortran derived type to hold model definition
  type :: model_ufs
     type(ESMF_GridComp) :: esmComp
     type(ESMF_State) :: toJedi, fromJedi
     integer :: isc, iec, jsc, jec, npz
     type(esmf_Clock) :: clock
     type(esmf_config) :: cf_main !<-- the configure object
     type(fckit_mpi_comm) :: comm
     logical :: initialized
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

    character(len=*),parameter :: subname = modname//' (create)'
    character(len=128) :: msg
    integer :: rc

    self%comm = geom%f_comm
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

    write(msg, "(a,5i7)") 'My dimensions are:', self%isc, self%iec, self%jsc, self%jec, self%npz
    call ESMF_LogWrite(trim(msg))

    self%initialized = .false.

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine create

! --------------------------------------------------------------------------------------------------

  subroutine initialize(self, state, vars, vdate_start, vdate_final)

    implicit none

    class(model_ufs),    intent(inout) :: self
    type(fv3jedi_state), intent(in)    :: state
    type(oops_variables),  intent(in)    :: vars

    type(datetime),      intent(in)    :: vdate_start
    type(datetime),      intent(in)    :: vdate_final

    integer :: rc, urc, phase, cnt, var

    type(ESMF_Time)         :: currTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    character(len=ESMF_MAXSTR), allocatable :: stdnames(:)

    character(len=20) :: strCurrTime, strStopTime

    type(ESMF_CplComp),  pointer       :: connectors(:)
    character(len=128) :: name, msg

    character(len=*),parameter :: subname = modname//' (initialize)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    if (self%initialized) then
        ! May need to reset/adjust the clock here, check later what is needed
        call ESMF_LogWrite("Model already initialized, do nothing", ESMF_LOGMSG_INFO)

        call ESMF_ClockPrint(self%clock, options="startTime", &
        preString="Printing startTime to stdout: ", rc=rc)

        call ESMF_ClockPrint(self%clock, options="currTime", &
        preString="Printing currTime to stdout: ", rc=rc)

        call ESMF_ClockPrint(self%clock, options="stopTime", &
        preString="Printing stopTime to stdout: ", rc=rc)

        call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)
        return
    end if

    ! Do only for very first initialization

    ! This may seem confusing. In the UFS, the model cold start time (ESMF lingo: startTime)
    ! never changes. Instead, the model warmstart/restart time (ESMF lingo: currTime) is
    ! updated as the model is integrated towards the stop time (ESMF lingo: stopTime).
    ! The latter can be updated (pushed out) as needed.
    call datetime_to_string(vdate_start, strCurrTime)
    call datetime_to_string(vdate_final, strStopTime)

    call ESMF_LogWrite(" --> REQUESTED START TIME:"//trim(strCurrTime), ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(" --> REQUESTED STOP  TIME:"//trim(strStopTime), ESMF_LOGMSG_INFO)

    call ESMF_TimeSet(currTime, timeString=strCurrTime, rc=rc)
    esmf_err_abort(rc)

    call ESMF_TimeSet(stopTime, timeString=strStopTime, rc=rc)
    esmf_err_abort(rc)

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
    call setUFSClock(self, currTime, stopTime)
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

    allocate(stdnames(vars%nvars()))
    do var = 1, vars%nvars()
       stdnames(var) = trim(vars%variable(var))
    enddo
    call ESMF_LogWrite("Advertising export from ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the exportState, for data coming out of ESM component
    ! Note--only certain fields are available. Check in GFS_surface_generic to see if they are filled
    ! Do only for very first initialization
    call NUOPC_Advertise(self%toJedi, &
         StandardNames=stdnames, &
         SharePolicyField="share", &
         TransferOfferGeomObject="cannot provide", rc=rc)
    esmf_err_abort(rc)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After filling advertise fromJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Advertising imports to ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the importState, for data going into ESM component
    ! Note--only certain fields are available. Check ???
    call NUOPC_Advertise(self%fromJedi, &
         StandardNames=stdnames, &
         SharePolicyField="share", &
         TransferOfferGeomObject="cannot provide", rc=rc)
    esmf_err_abort(rc)

    call ESMF_StateGet(self%fromJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After filling advertise fromJedi state has "//trim(msg)//" items.", &
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

    call ESMF_StateGet(self%fromJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After calling advertise fromJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)

    ! Set verbosity flag on connectors
    nullify(connectors);
    call NUOPC_DriverGetComp(self%esmComp, &
         compList=connectors, rc=rc)
    esmf_err_abort(rc)

    deallocate(connectors)
    deallocate(stdnames)

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

    call ESMF_StateGet(self%fromJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt

    call ESMF_LogWrite("Dumping fromJedi state with "//trim(msg)//" items", &
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

    call ESMF_LogWrite("Setting model as initialized", ESMF_LOGMSG_INFO)
    self%initialized = .true.

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
    character(len=20) :: strCurrTime, strStopTime
    character(len=128) :: name, msg
    type(ESMF_Time) :: currTime, stopTime
    type(ESMF_TimeInterval) :: timeStep

!-----------------------------------------------------------------------------

    character(len=*),parameter :: subname = modname//' (step)'
    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call datetime_to_string(vdate_start, strCurrTime)
    call datetime_to_string(vdate_final, strStopTime)

    call ESMF_LogWrite(" --> REQUESTED MODEL START TIME:"//trim(strCurrTime), ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(" --> REQUESTED MODEL STOP  TIME:"//trim(strStopTime), ESMF_LOGMSG_INFO)

    !call ESMF_ClockGet(self%clock, startTime=startTime, &
    !     currTime=currTime, stopTime=stopTime, rc=rc)
    !esmf_err_abort(rc)

    call ESMF_TimeSet(currTime, timeString=strCurrTime, rc=rc)
    esmf_err_abort(rc)

    call ESMF_TimeSet(stopTime, timeString=strStopTime, rc=rc)
    esmf_err_abort(rc)
!   call ESMF_TimeIntervalSet(timeStep, s=3600, rc=rc)
    timeStep = stopTime - currTime

    call ESMF_ClockSet(self%clock, &
         stopTime=stopTime, currTime=currTime, timeStep=timeStep, rc=rc)
    esmf_err_abort(rc)

    call ESMF_ClockPrint(self%clock, options="startTime", &
    preString="Printing startTime to stdout: ", rc=rc)

    call ESMF_ClockPrint(self%clock, options="currTime", &
    preString="Printing currTime to stdout: ", rc=rc)

    call ESMF_ClockPrint(self%clock, options="stopTime", &
    preString="Printing stopTime to stdout: ", rc=rc)

    call NUOPC_SetTimestamp(self%fromJedi, self%clock, rc=rc)
    esmf_err_abort(rc)
    write(msg, "(I3)") rc
    call ESMF_LogWrite("after NUOPC_SetTimestamp (self%fromJedi, self%clock, rc=rc): "//trim(msg)//" rc", &
         ESMF_LOGMSG_INFO)

    ! Update UFS state from JEDI
    call ESMF_StateGet(self%fromJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)
    write(msg, "(I2)") cnt
    call ESMF_LogWrite("before step fromJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)
    call state_to_fv3(self, state, strCurrTime)
    call ESMF_LogWrite("after UFS state write "//trim(msg)//" rc", &
         ESMF_LOGMSG_INFO)

    ! Step the model forward
    call ESMF_GridCompRun(self%esmComp, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    ! Update JEDI state from UFS
    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)
    write(msg, "(I2)") cnt
    call ESMF_LogWrite("after step toJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)
    call fv3_to_state(self, state, strCurrTime)
    call ESMF_LogWrite("after JEDI state write "//trim(msg)//" rc", &
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

    call ESMF_LogWrite("Setting model as uninitialized", ESMF_LOGMSG_INFO)
    self%initialized = .false.

    ! Finalize ESMF
    ! -------------
    call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=rc)
    if (rc /= ESMF_SUCCESS) then
       call ESMF_LogWrite("ERROR FINALIZING ESMF "//subname, ESMF_LOGMSG_INFO)
    endif

  end subroutine delete

  ! --------------------------------------------------------------------------------------------------

  subroutine finalize(self, state)

    implicit none
    class(model_ufs),    intent(inout) :: self
    type(fv3jedi_state), intent(in)    :: state
    integer :: rc
    character(len=*),parameter :: subname = modname//' (finalize)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_GridCompDestroy(self%esmComp, rc=rc)
    esmf_err_abort(rc)

    call ESMF_LogWrite("About to destroy toJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%toJedi, rc=rc)
    esmf_err_abort(rc)

    call ESMF_LogWrite("About to destroy fromJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%fromJedi, rc=rc)
    esmf_err_abort(rc)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)
    call mpp_set_current_pelist()

  end subroutine finalize

  subroutine fv3_to_state( self, state, strCurrTime )

  implicit none
  type(model_ufs),    intent(in)    :: self
  type(fv3jedi_state), intent(inout) :: state
  character(len=*), intent(in)      :: strCurrTime

  integer :: num_items, i, rc, rank, lb(3), ub(3), fnpz
  type(ESMF_Field) :: field
  character(len=ESMF_MAXSTR), allocatable :: item_names(:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr2(:,:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr3(:,:,:)
  character(len=field_clen) :: short_name
  type(fv3jedi_field), pointer :: field_ptr

  real(kind=ESMF_KIND_R8),allocatable,dimension(:,:,:)      :: field_fv3
  character(len=256) :: msg

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

        write(msg, "(4i6)") lb(1), ub(1), lb(2), ub(2)
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 2D dims: " // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(6i6)") self%isc, self%iec, self%jsc, self%jec, 1, fnpz
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", JEDI 2D dims: " // trim(msg),&
           ESMF_LOGMSG_INFO)

        write(msg, "(e16.7)") minval(farrayPtr2(lb(1):ub(1),lb(2):ub(2)))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 2D minval 1=" // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(e16.7)") maxval(farrayPtr2(lb(1):ub(1),lb(2):ub(2)))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 2D maxval 1=" // trim(msg),&
           ESMF_LOGMSG_INFO)

        write(msg, "(e16.7)") minval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 2D minval 2=" // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(e16.7)") maxval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 2D maxval 2=" // trim(msg),&
           ESMF_LOGMSG_INFO)

        nullify(farrayPtr2)

      elseif (rank == 3) then
        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )
        if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldGet 3D failed")

        fnpz = ub(3)-lb(3)+1
        field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))

        write(msg, "(6i6)") lb(1), ub(1), lb(2), ub(2), lb(3), ub(3)
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 3D dims: " // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(6i6)") self%isc, self%iec, self%jsc, self%jec, 1, fnpz
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", JEDI 3D dims: " // trim(msg),&
           ESMF_LOGMSG_INFO)

        write(msg, "(e16.7)") minval(farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 3D minval 1=" // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(e16.7)") maxval(farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 3D maxval 1=" // trim(msg),&
           ESMF_LOGMSG_INFO)

        write(msg, "(e16.7)") minval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 3D minval 2=" // trim(msg),&
           ESMF_LOGMSG_INFO)
        write(msg, "(e16.7)") maxval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", UFS 3D maxval 2=" // trim(msg),&
           ESMF_LOGMSG_INFO)

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

      write(msg, "(3i6)") lbound(field_ptr%array)
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", REAL JEDI 3D lbounds: " // trim(msg),&
           ESMF_LOGMSG_INFO)
      write(msg, "(3i6)") ubound(field_ptr%array)
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", REAL JEDI 3D ubounds: " // trim(msg),&
           ESMF_LOGMSG_INFO)

      if (field_ptr%npz .ne. fnpz) &
        call abor1_ftn("fv3_to_state: dimension mismatch between JEDI and UFS vertical grid")

      write(msg, "(e16.7)") minval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", JEDI minval=" // trim(msg),&
           ESMF_LOGMSG_INFO)
      write(msg, "(e16.7)") maxval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", JEDI maxval=" // trim(msg),&
           ESMF_LOGMSG_INFO)
      !
      ! Copy from UFS to fv3-jedi
      field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz)
      !
      write(msg, "(e16.7)") minval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("After updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", minval=" // trim(msg), &
           ESMF_LOGMSG_INFO)
      write(msg, "(e16.7)") maxval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("After updating field "// trim(short_name) // " in JEDI from UFS at " // trim(strCurrTime) // ", maxval=" // trim(msg), &
           ESMF_LOGMSG_INFO)
      !
      if (abs(minval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)))>1.0E30 .or. abs(maxval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)))>1.0E30) then
        write(msg, "(a,6i5,1x,e16.7,1x,3i6)") trim(short_name), self%isc, self%iec, self%jsc, self%jec, 1, fnpz, minval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)), minloc(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("DOM-ERROR MIN: " // trim(msg),ESMF_LOGMSG_INFO)
        write(msg, "(a,6i5,1x,e16.7,1x,3i6)") trim(short_name), self%isc, self%iec, self%jsc, self%jec, 1, fnpz, maxval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)), maxloc(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
        call ESMF_LogWrite("DOM-ERROR MAX: " // trim(msg),ESMF_LOGMSG_INFO)
      end if
    else
      call ESMF_LogWrite("Not needed by JEDI is "//short_name, ESMF_LOGMSG_INFO)
    endif

  end do

  deallocate(item_names)
  deallocate(field_fv3)

  end subroutine fv3_to_state


  subroutine state_to_fv3( self, state, strCurrTime )

  implicit none
  type(model_ufs),    intent(inout) :: self
  type(fv3jedi_state), intent(in)   :: state
  character(len=*), intent(in) :: strCurrTime

  integer :: num_items, i, rc, rank, lb(3), ub(3), fnpz
  type(ESMF_Field) :: field
  character(len=ESMF_MAXSTR), allocatable :: item_names(:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr2(:,:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr3(:,:,:)
  character(len=field_clen) :: short_name
  type(fv3jedi_field), pointer :: field_ptr

  real(kind=ESMF_KIND_R8),allocatable,dimension(:,:,:)      :: field_fv3
  character(len=256) :: msg

  ! Array to hold output from JEDI in UFS precision
  ! ------------------------------------------------
  allocate(field_fv3(self%isc:self%iec, self%jsc:self%jec, self%npz+1))


  ! Get number of items
  ! -------------------
  call ESMF_StateGet(self%fromJedi, itemcount = num_items, rc = rc)
  if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_StateGet itemcount failed")


  ! Get names of the items
  ! ----------------------
  allocate(item_names(num_items))
  call ESMF_StateGet(self%fromJedi, itemnamelist = item_names, rc = rc)
  if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_StateGet itemnamelist failed")


  ! Loop over states coming from UFS and convert to JEDI state
  ! -----------------------------------------------------------
  do i = 1, num_items
    ! Create map between UFS name and fv3-jedi name
    ! ----------------------------------------------
    short_name = trim(item_names(i))
    call ESMF_LogWrite("state_to_fv3: item name is "//short_name, ESMF_LOGMSG_INFO)
    ! DH*
    !if(trim(item_names(i)) == 't') short_name = 'air_temperature'
    !if(trim(item_names(i)) == 'T') short_name = 'air_temperature'
    call ESMF_LogWrite("Mapped UFS/ESMF field " // trim(item_names(i)) // " to JEDI field " // trim(short_name))
    ! *DH

    ! Only need to update field in UFS if fv3-jedi has it
    ! ---------------------------------------------------------
    if (state%has_field(trim(short_name))) then

      !Get field from the state
      call ESMF_StateGet(self%fromJedi, item_names(i), field, rc = rc)
      if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_StateGet field failed")

      call ESMF_LogWrite("Retrieved field "//short_name, ESMF_LOGMSG_INFO)

      !Validate the field
      call ESMF_FieldValidate(field, rc = rc)
      if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_FieldValidate failed")

      call ESMF_LogWrite("Validated field "//short_name, ESMF_LOGMSG_INFO)

      !Get the field rank
      call ESMF_FieldGet(field, rank = rank, rc = rc)
      if (rc.ne.0) call abor1_ftn("fv3_to_state: ESMF_FieldGet rank failed")

      call ESMF_LogWrite("Retrieved rank for field "//short_name, ESMF_LOGMSG_INFO)

      !Convert field to pointer and pointer bounds
      if (rank == 2) then
        call ESMF_LogWrite("Calling 'call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb(1:2), totalUBound = ub(1:2), rc = rc )'")
        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb(1:2), totalUBound = ub(1:2), rc = rc )
        if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_FieldGet 2D failed")
        fnpz = 1
      elseif (rank == 3) then
        call ESMF_LogWrite("Calling 'call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )'")
        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )
        if (rc.ne.0) call abor1_ftn("state_to_fv3: ESMF_FieldGet 3D failed")
        fnpz = ub(3)-lb(3)+1
      else
        call abor1_ftn("fv3_mod: can only handle rank 2 or rank 3 fields from UFS")
      endif

      ! Check that dimensions match
      if ((ub(1)-lb(1)+1 .ne. self%iec-self%isc+1) .or. (ub(2)-lb(2)+1 .ne. self%jec-self%jsc+1) ) then
        call abor1_ftn("state_to_fv3: dimension mismatch between JEDI and UFS horizontal grid")
      endif

      ! Get pointer to fv3-jedi side field
      call state%get_field(trim(short_name), field_ptr)

      call ESMF_LogWrite("Got field pointer for field "//short_name, ESMF_LOGMSG_INFO)
      write(msg, "(a,e16.7,a,e16.7)") "field_ptr for " // trim(short_name) // " has minval ", minval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)), " and maxval ", maxval(field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

      if (field_ptr%npz .ne. fnpz) &
        call abor1_ftn("state_to_fv3: dimension mismatch between JEDI and UFS vertical grid")

      write(msg, "(e16.7)") minval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in UFS from JEDI at " // trim(strCurrTime) // ", minval=" // trim(msg))
      write(msg, "(e16.7)") maxval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("Before updating field "// trim(short_name) // " in UFS from JEDI at " // trim(strCurrTime) // ", maxval=" // trim(msg))
      !
      ! Copy from fv3-jedi to UFS
      !
      field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = field_ptr%array(self%isc:self%iec,self%jsc:self%jec,1:fnpz)
      !
      write(msg, "(e16.7)") minval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("After updating field "// trim(short_name) // " in UFS from JEDI at " // trim(strCurrTime) // ", minval=" // trim(msg))
      write(msg, "(e16.7)") maxval(field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz))
      call ESMF_LogWrite("After updating field "// trim(short_name) // " in UFS from JEDI at " // trim(strCurrTime) // ", maxval=" // trim(msg))

      ! Update UFS state fieldpointer
      if (rank == 2) then
        farrayPtr2(lb(1):ub(1),lb(2):ub(2)) = field_fv3(self%isc:self%iec,self%jsc:self%jec,1)
        write(msg, "(a,e16.7,a,e16.7)") "farrayPtr2 for " // trim(short_name) // " has minval ", minval(farrayPtr2), " and maxval ", maxval(farrayPtr2)
        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
        nullify(farrayPtr2)
      elseif (rank == 3) then
        farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = field_fv3(self%isc:self%iec,self%jsc:self%jec,1:fnpz)
        write(msg, "(a,e16.7,a,e16.7)") "farrayPtr3 for " // trim(short_name) // " has minval ", minval(farrayPtr3), " and maxval ", maxval(farrayPtr3)
        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
        nullify(farrayPtr3)
      else
        call abor1_ftn("fv3_mod: can only handle rank 2 or rank 3 fields from UFS")
      endif

    else
      call ESMF_LogWrite("Not provided by JEDI is "//short_name, ESMF_LOGMSG_INFO)
    endif

  end do

  deallocate(item_names)
  deallocate(field_fv3)

  end subroutine state_to_fv3


  subroutine setUFSClock(self,date_current,date_final)

    class(model_ufs),    intent(inout) :: self
    type(esmf_time),     intent(in) :: date_current, date_final

    type(esmf_time)       :: date_start
    type(ESMF_TimeInterval) :: timestep
    !type(ESMF_TimeInterval)            :: runDuration, timestep
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

    ! Set any time interval here, it will be overwritten later

    timeStep = date_final - date_current
!   call ESMF_TimeIntervalSet(timeStep, s=3600, rc=rc)

    self%clock = ESMF_ClockCreate(name="main_clock", &
         timeStep=timestep, startTime=date_start, stopTime=date_final, rc=rc)
    esmf_err_abort(rc)

    call ESMF_ClockSet(self%clock, currTime=date_current, rc=rc)
    esmf_err_abort(rc)

    call ESMF_ClockPrint(self%clock, options="startTime", &
    preString="Printing startTime to stdout: ", rc=rc)

    call ESMF_ClockPrint(self%clock, options="currTime", &
    preString="Printing currTime to stdout: ", rc=rc)

    call ESMF_ClockPrint(self%clock, options="stopTime", &
    preString="Printing stopTime to stdout: ", rc=rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  end subroutine setUFSClock

end module fv3jedi_ufs_mod
