! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_fv3_mod

use iso_c_binding
use config_mod
use datetime_mod
use duration_mod
use netcdf

use fv3jedi_kinds_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment 

use fv3jedi_lm_mod, only: fv3jedi_lm_type

implicit none
private

public :: fv3_model
public :: fv3_create
public :: fv3_delete
public :: fv3_initialize
public :: fv3_step
public :: fv3_finalize

! ------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: fv3_model
  type(fv3jedi_lm_type)                        :: fv3jedi_lm          !<Linearized model object
  integer                                      :: readtraj            !<Read trajectory from file
  character(len=255)                           :: trajpath            !<User specified path to traj files
  character(len=255)                           :: trajfile            !<User specified path to traj files
end type fv3_model

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3_create(self, geom, c_conf)

implicit none
type(c_ptr),         intent(in)    :: c_conf
type(fv3_model), intent(inout) :: self
type(fv3jedi_geom),  intent(in)    :: geom

!Locals
character(len=20) :: ststep
type(duration) :: dtstep
real(kind=kind_real) :: dt


! Model time step
! ---------------
ststep = config_get_string(c_conf,len(ststep),"tstep")
dtstep = trim(ststep)
dt = real(duration_seconds(dtstep),kind_real)


! Option to read traj from file instead of propagating model
! ----------------------------------------------------------
self%readtraj = config_get_int(c_conf,"readtraj")
if (self%readtraj == 1) then
  self%trajpath = config_get_string(c_conf,len(self%trajpath),"trajpath")
  self%trajfile = config_get_string(c_conf,len(self%trajfile),"trajfile")
endif


! Model configuration and creation
! --------------------------------
self%fv3jedi_lm%conf%do_dyn     = config_get_int(c_conf,"lm_do_dyn")
self%fv3jedi_lm%conf%do_phy_trb = config_get_int(c_conf,"lm_do_trb")
self%fv3jedi_lm%conf%do_phy_mst = config_get_int(c_conf,"lm_do_mst")

call self%fv3jedi_lm%create(dt,geom%npx,geom%npy,geom%npz,geom%ptop,geom%ak,geom%bk)


! Safety checks
! -------------

!The full trajecotory of the tlm/adm is not output by this simplified model
!so if being used to generate the trajectry with physics the traj must be read
!from file or obtained by running GEOS or GFS. 
if ((self%fv3jedi_lm%conf%do_phy_trb .ne. 0 .and. self%readtraj == 0) .or. &
    (self%fv3jedi_lm%conf%do_phy_mst .ne. 0 .and. self%readtraj == 0) ) then
   call abor1_ftn("fv3_model | FV3 : unless reading the trajecotory physics should be off")
endif

end subroutine fv3_create

! ------------------------------------------------------------------------------

subroutine fv3_delete(self)

implicit none
type(fv3_model), intent(inout) :: self

!Delete the model
!----------------
call self%fv3jedi_lm%delete()

end subroutine fv3_delete

! ------------------------------------------------------------------------------

subroutine fv3_initialize(self, state)

implicit none
type(fv3_model), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

call self%fv3jedi_lm%init_nl()

end subroutine fv3_initialize

! ------------------------------------------------------------------------------

subroutine fv3_step(self, state, vdate)

implicit none
type(fv3_model), intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(datetime),      intent(in)    :: vdate !< Valid datetime after step

character(len=20) :: vdatec

if (self%readtraj == 0) then

  call state_to_lm(state,self%fv3jedi_lm)
  call self%fv3jedi_lm%step_nl()
  call lm_to_state(self%fv3jedi_lm,state)

else

  call datetime_to_string(vdate, vdatec)
  call read_state( self, state, vdatec)

endif

end subroutine fv3_step

! ------------------------------------------------------------------------------

subroutine fv3_finalize(self, state)

implicit none
type(fv3_model), target :: self
type(fv3jedi_state)         :: state

call self%fv3jedi_lm%final_nl()

end subroutine fv3_finalize

! ------------------------------------------------------------------------------

subroutine state_to_lm( state, lm )

implicit none
type(fv3jedi_state),   intent(in)    :: state
type(fv3jedi_lm_type), intent(inout) :: lm
 
lm%traj%u       = state%fields(state%ud)%field
lm%traj%v       = state%fields(state%vd)%field
lm%traj%ua      = state%fields(state%ua)%field
lm%traj%va      = state%fields(state%va)%field
lm%traj%t       = state%fields(state%t)%field
lm%traj%delp    = state%fields(state%delp)%field
lm%traj%qv      = state%fields(state%q)%field
lm%traj%qi      = state%fields(state%qi)%field
lm%traj%ql      = state%fields(state%ql)%field
lm%traj%o3      = state%fields(state%o3)%field

if (.not. lm%conf%hydrostatic) then
lm%traj%w       = state%fields(state%w)%field
lm%traj%delz    = state%fields(state%delz)%field
endif

lm%traj%phis = state%fields(state%phis)%field(:,:,1)

end subroutine state_to_lm

! ------------------------------------------------------------------------------

subroutine lm_to_state( lm, state )

implicit none
type(fv3jedi_lm_type), intent(in)    :: lm
type(fv3jedi_state),   intent(inout) :: state
 
state%fields(state%ud)%field      = lm%traj%u
state%fields(state%vd)%field      = lm%traj%v
state%fields(state%ua)%field      = lm%traj%ua
state%fields(state%va)%field      = lm%traj%va
state%fields(state%t)%field       = lm%traj%t
state%fields(state%delp)%field    = lm%traj%delp
state%fields(state%q)%field       = lm%traj%qv
state%fields(state%qi)%field      = lm%traj%qi
state%fields(state%ql)%field      = lm%traj%ql
state%fields(state%o3)%field      = lm%traj%o3

if (.not. lm%conf%hydrostatic) then
state%fields(state%w)%field       = lm%traj%w
state%fields(state%delz)%field    = lm%traj%delz
endif

state%fields(state%phis)%field(:,:,1)    = lm%traj%phis

end subroutine lm_to_state

! ------------------------------------------------------------------------------

subroutine read_state( self, state, vdatec)

implicit none
type(fv3_model), intent(in)    :: self
type(fv3jedi_state), intent(inout) :: state
character(len=20),   intent(in)    :: vdatec

character(len=255) :: date, path, fname1, fname2, filename
character(len=4)   :: yyyy,mm,dd,hh,mn
character(len=20)  :: var
integer, allocatable :: istart(:), icount(:)
integer :: ncid, ncstat, dimid, varid
integer :: im, jm

integer :: tileoff
logical :: tiledimension = .false.

integer :: isc,iec,jsc,jec

write(date,*) vdatec(1:4),vdatec(6:7),vdatec(9:10),'_',vdatec(12:13),vdatec(15:16),'z.nc4'

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec

path = self%trajpath
fname1 = self%trajfile

!> Build filename
yyyy = vdatec(1 :4 )
mm   = vdatec(6 :7 )
dd   = vdatec(9 :10)
hh   = vdatec(12:13)
mn   = vdatec(15:16)
fname2 = 'z.nc4'

filename = trim(path)//trim(fname1)//trim(yyyy)//trim(mm)//trim(dd)//"_"//trim(hh)//trim(mn)//trim(fname2)

if (self%fv3jedi_lm%conf%rpe) print*, ' '
if (self%fv3jedi_lm%conf%rpe) print*, 'Reading trajectory: ', trim(filename)
if (self%fv3jedi_lm%conf%rpe) print*, ' '

!> Open the file
ncstat = nf90_open(trim(filename), NF90_NOWRITE, ncid)
if(ncstat /= nf90_noerr) print *, "OPEN: "//trim(nf90_strerror(ncstat))

!> Get dimensions, lon,lat
ncstat = nf90_inq_dimid(ncid, "lon", dimid)
if(ncstat /= nf90_noerr) print *, "lon: "//trim(nf90_strerror(ncstat))
ncstat = nf90_inquire_dimension(ncid, dimid, len = im)
if(ncstat /= nf90_noerr) print *, "lon:"//trim(nf90_strerror(ncstat))

ncstat = nf90_inq_dimid(ncid, "lat", dimid)
if(ncstat /= nf90_noerr) print *, "lat: "//trim(nf90_strerror(ncstat))
ncstat = nf90_inquire_dimension(ncid, dimid, len = jm)
if(ncstat /= nf90_noerr) print *, "lat: "//trim(nf90_strerror(ncstat))

!> GEOS can use concatenated tiles or tile as a dimension
if ( (im == state%npx-1) .and. (jm == 6*(state%npy-1) ) ) then
  tiledimension = .false.
  tileoff = (state%ntile-1)*(jm/state%ntiles)
else
  tiledimension = .true.
  tileoff = 0
  call abor1_ftn("Trajectory GEOS: tile dimension in file not done yet")
endif

allocate(istart(4))
allocate(icount(4))
istart(1) = isc
istart(2) = tileoff + jsc
istart(3) = 1
istart(4) = 1

icount(1) = iec-isc+1
icount(2) = jec-jsc+1
icount(3) = 72
icount(4) = 1

var = 'ud'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ud)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'vd'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%vd)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'ua'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ua)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'va'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%va)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 't'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%t)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'delp'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%delp)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'q'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%q)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'qi'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%qi)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'ql'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ql)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'o3mr'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%o3)%field, istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

if (.not.self%fv3jedi_lm%conf%hydrostatic) then

  var = 'w'
  ncstat = nf90_inq_varid (ncid, trim(var), varid)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
  ncstat = nf90_get_var(ncid, varid, state%fields(state%w)%field, istart, icount)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 
  var = 'delz'
  ncstat = nf90_inq_varid (ncid, trim(var), varid)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
  ncstat = nf90_get_var(ncid, varid, state%fields(state%delz)%field, istart, icount)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

endif

if (self%fv3jedi_lm%conf%do_phy_mst .ne. 0) then

  var = 'qls'
  ncstat = nf90_inq_varid (ncid, trim(var), varid)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
  ncstat = nf90_get_var(ncid, varid, state%fields(state%qls)%field, istart, icount)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 
  var = 'qcn'
  ncstat = nf90_inq_varid (ncid, trim(var), varid)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
  ncstat = nf90_get_var(ncid, varid, state%fields(state%qcn)%field, istart, icount)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
 
  var = 'cfcn'
  ncstat = nf90_inq_varid (ncid, trim(var), varid)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
  ncstat = nf90_get_var(ncid, varid, state%fields(state%cfcn)%field, istart, icount)
  if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

endif

deallocate(istart,icount)

allocate(istart(3))
allocate(icount(3))
istart(1) = isc
istart(2) = tileoff + jsc
istart(3) = 1

icount(1) = iec-isc+1
icount(2) = jec-jsc+1
icount(3) = 1

var = 'phis'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%phis)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'frocean'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%frocean)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'frland'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%frland)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'varflt'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%varflt)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'ustar'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ustar)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'bstar'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%bstar)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'zpbl'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%zpbl)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'cm'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%cm)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'ct'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ct)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'cq'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%cq)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'kcbl'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%kcbl)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'ts'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%ts)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'khl'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%khl)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

var = 'khu'
ncstat = nf90_inq_varid (ncid, trim(var), varid)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, state%fields(state%khu)%field(:,:,1), istart, icount)
if(ncstat /= nf90_noerr) print *, trim(var)//trim(nf90_strerror(ncstat))

!Close this file
ncstat = nf90_close(ncid)
if(ncstat /= nf90_noerr) print *, "CLOSE: "//trim(nf90_strerror(ncstat))

deallocate(istart,icount)

end subroutine read_state

! ------------------------------------------------------------------------------

end module fv3jedi_fv3_mod
