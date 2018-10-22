! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Handle state for the FV3JEDI odel

module fv3jedi_state_mod

use iso_c_binding
use config_mod
use datetime_mod

use ioda_locs_mod
use ufo_vars_mod
use ufo_geovals_mod

use fv3jedi_constants_mod, only: rad2deg, constoz
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_state_io_mod 
use fv3jedi_state_utils_mod, only: fv3jedi_state
use fv3jedi_vars_mod, only: fv3jedi_vars
use fv3jedi_getvalues_mod, only: getvalues

implicit none

private
public :: create, delete, zeros, copy, axpy, add_incr, &
          read_file, write_file, gpnorm, staterms, &
          change_resol, getvalues, analytic_IC
public :: fv3jedi_state

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_geom), target,  intent(in)    :: geom
type(fv3jedi_vars),  intent(in)    :: vars

integer :: var

! Copy the variable names
self%vars%nv = vars%nv
allocate(self%vars%fldnames(self%vars%nv))
self%vars%fldnames = vars%fldnames

! Allocate variables based on names
do var = 1, self%vars%nv

   select case (trim(self%vars%fldnames(var)))

     case("ud")
       if (.not.allocated(  self%ud)) allocate (  self%ud(geom%isc:geom%iec,  geom%jsc:geom%jec+1, geom%npz))
     case("vd")
       if (.not.allocated(  self%vd)) allocate (  self%vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  , geom%npz))
     case("ua")
       if (.not.allocated(  self%ua)) allocate (  self%ua(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("va")
       if (.not.allocated(  self%va)) allocate (  self%va(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("t")
       if (.not.allocated(   self%t)) allocate (   self%t(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("delp")
       if (.not.allocated(self%delp)) allocate (self%delp(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("q")
       if (.not.allocated(   self%q)) allocate (   self%q(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("qi")
       if (.not.allocated(  self%qi)) allocate (  self%qi(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("ql")
       if (.not.allocated(  self%ql)) allocate (  self%ql(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("o3")
       if (.not.allocated(  self%o3)) allocate (  self%o3(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("w")
       if (.not.allocated(   self%w)) allocate (   self%w(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("delz")
       if (.not.allocated(self%delz)) allocate (self%delz(geom%isc:geom%iec,  geom%jsc:geom%jec  , geom%npz))
     case("ps")
     case default 
       call abor1_ftn("Create: unknown variable "//trim(self%vars%fldnames(var)))

   end select

enddo

self%hydrostatic = .true.
if (allocated(self%w).and.allocated(self%delz)) self%hydrostatic = .false.

if (.not.allocated(self%phis)) allocate(self%phis(geom%isc:geom%iec,geom%jsc:geom%jec))

!CRTM surface variables
if (.not.allocated(self%slmsk )) allocate(self%slmsk (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%sheleg)) allocate(self%sheleg(geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%tsea  )) allocate(self%tsea  (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%vtype )) allocate(self%vtype (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%stype )) allocate(self%stype (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%vfrac )) allocate(self%vfrac (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%stc   )) allocate(self%stc   (geom%isc:geom%iec,geom%jsc:geom%jec,4))
if (.not.allocated(self%smc   )) allocate(self%smc   (geom%isc:geom%iec,geom%jsc:geom%jec,4))
if (.not.allocated(self%snwdph)) allocate(self%snwdph(geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%u_srf )) allocate(self%u_srf (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%v_srf )) allocate(self%v_srf (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%f10m  )) allocate(self%f10m  (geom%isc:geom%iec,geom%jsc:geom%jec))

!Linearized model trajectory
if (.not.allocated(self%qls    )) allocate(self%qls    (geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
if (.not.allocated(self%qcn    )) allocate(self%qcn    (geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
if (.not.allocated(self%cfcn   )) allocate(self%cfcn   (geom%isc:geom%iec,geom%jsc:geom%jec,geom%npz))
if (.not.allocated(self%frocean)) allocate(self%frocean(geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%frland )) allocate(self%frland (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%varflt )) allocate(self%varflt (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%ustar  )) allocate(self%ustar  (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%bstar  )) allocate(self%bstar  (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%zpbl   )) allocate(self%zpbl   (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%cm     )) allocate(self%cm     (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%ct     )) allocate(self%ct     (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%cq     )) allocate(self%cq     (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%kcbl   )) allocate(self%kcbl   (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%ts     )) allocate(self%ts     (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%khl    )) allocate(self%khl    (geom%isc:geom%iec,geom%jsc:geom%jec))
if (.not.allocated(self%khu    )) allocate(self%khu    (geom%isc:geom%iec,geom%jsc:geom%jec))

! Initialize all domain arrays to zero
call zeros(self)
self%phis   = 0.0_kind_real

! For convenience
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%isd = geom%isd
self%ied = geom%ied
self%jsd = geom%jsd
self%jed = geom%jed
self%npx = geom%npx
self%npy = geom%npy
self%npz = geom%npz

self%ntile = geom%ntile
self%ntiles = geom%ntiles

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(fv3jedi_state), intent(inout) :: self

if (allocated(self%ud  )) deallocate (self%ud  )
if (allocated(self%vd  )) deallocate (self%vd  )
if (allocated(self%ua  )) deallocate (self%ua  )
if (allocated(self%va  )) deallocate (self%va  )
if (allocated(self%t   )) deallocate (self%t   )
if (allocated(self%delp)) deallocate (self%delp)
if (allocated(self%q   )) deallocate (self%q   )
if (allocated(self%qi  )) deallocate (self%qi  )
if (allocated(self%ql  )) deallocate (self%ql  )
if (allocated(self%o3  )) deallocate (self%o3  )
if (allocated(self%phis)) deallocate (self%phis)
if (allocated(self%w   )) deallocate (self%w   )
if (allocated(self%delz)) deallocate (self%delz)

if (allocated(self%slmsk )) deallocate(self%slmsk )
if (allocated(self%sheleg)) deallocate(self%sheleg)
if (allocated(self%tsea  )) deallocate(self%tsea  )
if (allocated(self%vtype )) deallocate(self%vtype )
if (allocated(self%stype )) deallocate(self%stype )
if (allocated(self%vfrac )) deallocate(self%vfrac )
if (allocated(self%stc   )) deallocate(self%stc   )
if (allocated(self%smc   )) deallocate(self%smc   )
if (allocated(self%snwdph)) deallocate(self%snwdph)
if (allocated(self%u_srf )) deallocate(self%u_srf )
if (allocated(self%v_srf )) deallocate(self%v_srf )
if (allocated(self%f10m  )) deallocate(self%f10m  )

if (allocated(self%qls    )) deallocate(self%qls    )
if (allocated(self%qcn    )) deallocate(self%qcn    )
if (allocated(self%cfcn   )) deallocate(self%cfcn   )
if (allocated(self%frocean)) deallocate(self%frocean)
if (allocated(self%frland )) deallocate(self%frland )
if (allocated(self%varflt )) deallocate(self%varflt )
if (allocated(self%ustar  )) deallocate(self%ustar  )
if (allocated(self%bstar  )) deallocate(self%bstar  )
if (allocated(self%zpbl   )) deallocate(self%zpbl   )
if (allocated(self%cm     )) deallocate(self%cm     )
if (allocated(self%ct     )) deallocate(self%ct     )
if (allocated(self%cq     )) deallocate(self%cq     )
if (allocated(self%kcbl   )) deallocate(self%kcbl   )
if (allocated(self%ts     )) deallocate(self%ts     )
if (allocated(self%khl    )) deallocate(self%khl    )
if (allocated(self%khu    )) deallocate(self%khu    )

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(fv3jedi_state), intent(inout) :: self

!Zero out the entire domain

!Model
if(allocated(self%ud  )) self%ud   = 0.0_kind_real
if(allocated(self%vd  )) self%vd   = 0.0_kind_real
if(allocated(self%ua  )) self%ua   = 0.0_kind_real
if(allocated(self%va  )) self%va   = 0.0_kind_real
if(allocated(self%t   )) self%t    = 0.0_kind_real
if(allocated(self%delp)) self%delp = 0.0_kind_real
if(allocated(self%q   )) self%q    = 0.0_kind_real
if(allocated(self%qi  )) self%qi   = 0.0_kind_real
if(allocated(self%ql  )) self%ql   = 0.0_kind_real
if(allocated(self%o3  )) self%o3   = 0.0_kind_real
if(allocated(self%w   )) self%w    = 0.0_kind_real
if(allocated(self%delz)) self%delz = 0.0_kind_real

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: rhs

self%isc            = rhs%isc           
self%iec            = rhs%iec           
self%jsc            = rhs%jsc           
self%jec            = rhs%jec           
self%isd            = rhs%isd           
self%ied            = rhs%ied           
self%jsd            = rhs%jsd           
self%jed            = rhs%jed           
self%npx            = rhs%npx           
self%npy            = rhs%npy           
self%npz            = rhs%npz           
self%havecrtmfields = rhs%havecrtmfields
self%hydrostatic    = rhs%hydrostatic   
self%calendar_type  = rhs%calendar_type 
self%date           = rhs%date          
self%date_init      = rhs%date_init     
self%ntile          = rhs%ntile
self%ntiles         = rhs%ntiles

if(allocated(self%ud  ).and.allocated(rhs%ud  )) self%ud   = rhs%ud  
if(allocated(self%vd  ).and.allocated(rhs%vd  )) self%vd   = rhs%vd  
if(allocated(self%ua  ).and.allocated(rhs%ua  )) self%ua   = rhs%ua  
if(allocated(self%va  ).and.allocated(rhs%va  )) self%va   = rhs%va  
if(allocated(self%t   ).and.allocated(rhs%t   )) self%t    = rhs%t   
if(allocated(self%delp).and.allocated(rhs%delp)) self%delp = rhs%delp
if(allocated(self%q   ).and.allocated(rhs%q   )) self%q    = rhs%q   
if(allocated(self%qi  ).and.allocated(rhs%qi  )) self%qi   = rhs%qi  
if(allocated(self%ql  ).and.allocated(rhs%ql  )) self%ql   = rhs%ql  
if(allocated(self%o3  ).and.allocated(rhs%o3  )) self%o3   = rhs%o3  
if(allocated(self%w   ).and.allocated(rhs%w   )) self%w    = rhs%w   
if(allocated(self%delz).and.allocated(rhs%delz)) self%delz = rhs%delz

self%phis   = rhs%phis
self%slmsk  = rhs%slmsk
self%sheleg = rhs%sheleg
self%tsea   = rhs%tsea
self%vtype  = rhs%vtype
self%stype  = rhs%stype
self%vfrac  = rhs%vfrac
self%stc    = rhs%stc
self%smc    = rhs%smc
self%snwdph = rhs%snwdph
self%u_srf  = rhs%u_srf
self%v_srf  = rhs%v_srf
self%f10m   = rhs%f10m

self%qls     = rhs%qls    
self%qcn     = rhs%qcn    
self%cfcn    = rhs%cfcn   
self%frocean = rhs%frocean
self%frland  = rhs%frland 
self%varflt  = rhs%varflt 
self%ustar   = rhs%ustar  
self%bstar   = rhs%bstar  
self%zpbl    = rhs%zpbl   
self%cm      = rhs%cm     
self%ct      = rhs%ct     
self%cq      = rhs%cq     
self%kcbl    = rhs%kcbl   
self%ts      = rhs%ts     
self%khl     = rhs%khl    
self%khu     = rhs%khu    

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(fv3jedi_state), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(fv3jedi_state), intent(in)    :: rhs

if(allocated(self%ud  ).and.allocated(rhs%ud  )) self%ud   = self%ud   + zz * rhs%ud  
if(allocated(self%vd  ).and.allocated(rhs%vd  )) self%vd   = self%vd   + zz * rhs%vd  
if(allocated(self%ua  ).and.allocated(rhs%ua  )) self%ua   = self%ua   + zz * rhs%ua  
if(allocated(self%va  ).and.allocated(rhs%va  )) self%va   = self%va   + zz * rhs%va  
if(allocated(self%t   ).and.allocated(rhs%t   )) self%t    = self%t    + zz * rhs%t   
if(allocated(self%delp).and.allocated(rhs%delp)) self%delp = self%delp + zz * rhs%delp
if(allocated(self%q   ).and.allocated(rhs%q   )) self%q    = self%q    + zz * rhs%q   
if(allocated(self%qi  ).and.allocated(rhs%qi  )) self%qi   = self%qi   + zz * rhs%qi  
if(allocated(self%ql  ).and.allocated(rhs%ql  )) self%ql   = self%ql   + zz * rhs%ql  
if(allocated(self%o3  ).and.allocated(rhs%o3  )) self%o3   = self%o3   + zz * rhs%o3  
if(allocated(self%w   ).and.allocated(rhs%w   )) self%w    = self%w    + zz * rhs%w   
if(allocated(self%delz).and.allocated(rhs%delz)) self%delz = self%delz + zz * rhs%delz

return
end subroutine axpy

! ------------------------------------------------------------------------------

subroutine add_incr(geom,self,rhs)

use wind_vt_mod, only: a2d

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_state),     intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,npz,k

real(kind=kind_real), allocatable, dimension(:,:,:) :: ud, vd

!Check for matching resolution between state and increment
if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1)==0) then

  isc = rhs%isc
  iec = rhs%iec
  jsc = rhs%jsc
  jec = rhs%jec
  isd = rhs%isd
  ied = rhs%ied
  jsd = rhs%jsd
  jed = rhs%jed
  npz = rhs%npz

  !Convert A-Grid increment to D-Grid
  allocate(ud(isc:iec  ,jsc:jec+1,1:npz))
  allocate(vd(isc:iec+1,jsc:jec  ,1:npz))
  ud = 0.0_kind_real
  vd = 0.0_kind_real

  call a2d(geom, rhs%ua(isc:iec,jsc:jec,1:npz), rhs%va(isc:iec,jsc:jec,1:npz), ud, vd)

  if(allocated(self%ud  )) self%ud(isc:iec  ,jsc:jec+1,:)   = self%ud(isc:iec  ,jsc:jec+1,:)   + ud  (isc:iec  ,jsc:jec+1,:)
  if(allocated(self%vd  )) self%vd(isc:iec+1,jsc:jec  ,:)   = self%vd(isc:iec+1,jsc:jec  ,:)   + vd  (isc:iec+1,jsc:jec  ,:)

  deallocate(ud,vd)

  if(allocated(self%ua  )) self%ua   = self%ua   + rhs%ua  
  if(allocated(self%va  )) self%va   = self%va   + rhs%va  
  if(allocated(self%t   )) self%t    = self%t    + rhs%t
  if(allocated(self%delp)) then
    if (allocated(rhs%ps)) then
      do k = 1,geom%npz
        self%delp(:,:,k) = self%delp(:,:,k) + (geom%bk(k+1)-geom%bk(k))*rhs%ps
      enddo
    elseif (allocated(rhs%delp)) then
      self%delp = self%delp + rhs%delp
    endif
  endif
  if(allocated(self%q   )) self%q    = self%q    + rhs%q   
  if(allocated(self%qi  )) self%qi   = self%qi   + rhs%qi  
  if(allocated(self%ql  )) self%ql   = self%ql   + rhs%ql  
  if(allocated(self%o3  )) self%o3   = self%o3   + rhs%o3  
  if(allocated(self%w   )) self%w    = self%w    + rhs%w   
  if(allocated(self%delz)) self%delz = self%delz + rhs%delz 
else
   call abor1_ftn("fv3jedi state:  add_incr not implemented for low res increment yet")
endif

return
end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine change_resol(state,rhs)
implicit none
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_state), intent(in)    :: rhs

integer :: check
check = (rhs%iec-rhs%isc+1) - (state%iec-state%isc+1)

if (check==0) then
   call copy(state, rhs)
else
   call abor1_ftn("fv3jedi_state: change_resol not implmeneted yet")
endif

return
end subroutine change_resol

! ------------------------------------------------------------------------------
!> Analytic Initialization for the FV3 Model
!!
!! \details **analytic_IC()** initializes the FV3JEDI state and State objects using one of
!! several alternative idealized analytic models.  This is intended to facilitate testing by
!! eliminating the need to read in the initial state from a file and by providing exact expressions
!! to test interpolations.  This function is activated by setting the "analytic_init" state in the
!! "initial" or "StateFile" section of the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * invent-state: Backward compatibility with original analytic init option
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!! * dcmip-test-3-1: Non-hydrostatic gravity wave
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \author M. Miesch (adapted from a pre-existing call to invent_state)
!! \date March, 2018: Created
!! \date May, 2018: Added dcmip-test-3-1
!! \date June, 2018: Added dcmip-test-4-0
!!
!! \warning This routine initializes the fv3jedi_state object.  However, since the fv_atmos_type
!! component of fv3jedi_state is a subset of the corresponding object in the fv3 model,
!! this initialization routine is not sufficient to comprehensively define the full fv3 state.
!! So, this intitialization can be used for interpolation and other tests within JEDI but it is
!! cannot currently be used to initiate a forecast with gfs.
!!
!! \warning This routine does not initialize the fv3jedi_interp member of the fv3jedi_state object
!!
!! \warning Though an input state file is not required for these analytic initialization routines,
!! some grid information (in particular the hybrid vertical grid coefficients ak and bk)
!! is still read in from an input file when creating the geometry object that is a required
!! member of fv3jedi_state; see c_fv3jedi_geo_setup() in fv3jedi_geom_mod.F90.
!!
!! \warning It's unclear whether the pt member of the fv_atmos_type structure is potential temperature
!! or temperature.  This routine assumes the latter.  If this is not correct, then we will need to
!! implement a conversion
!!
subroutine analytic_IC(state, geom, c_conf, vdate)

  use fv3jedi_kinds_mod
  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log
  use constants_mod, only: pi=>pi_8
  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

  !FV3 Test Cases
  use fv_arrays_mod,  only: fv_atmos_type, deallocate_fv_atmos_type
  use test_cases_mod, only: init_case, test_case
  use fv_control_mod, only: fv_init, pelist_all

  implicit none

  type(fv3jedi_state), intent(inout) :: state !< State
  type(fv3jedi_geom),  intent(inout) :: geom    !< Geometry 
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: i,j,k
  real(kind=kind_real) :: rlat, rlon, z
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  type(fv_atmos_type), allocatable :: FV_AtmIC(:)
  real(kind=kind_real)             :: DTdummy = 900.0
  logical, allocatable             :: grids_on_this_pe(:)
  integer                          :: p_split = 1

  If (config_element_exists(c_conf,"analytic_init")) Then
     IC = Trim(config_get_string(c_conf,len(IC),"analytic_init"))
  Else
     ! This default value is for backward compatibility
     IC = "invent-state"
  EndIf

  call log%warning("fv3jedi_state:analytic_init: "//IC)
  sdate = config_get_string(c_conf,len(sdate),"date")
  WRITE(buf,*) 'validity date is: '//sdate
  call log%info(buf)
  call datetime_set(sdate, vdate)

  !===================================================================
  int_option: Select Case (IC)

     Case("invent-state")

        call invent_state(state,c_conf,geom)

     Case("fv3_init_case")

        !Initialize temporary FV_Atm fv3 construct
        call fv_init(FV_AtmIC, DTdummy, grids_on_this_pe, p_split)
        deallocate(pelist_all)

        !Test case to run, see fv3: /tools/test_cases.F90 for possibilities
        test_case = config_get_int(c_conf,"fv3_test_case")

        call init_case( FV_AtmIC(1)%u,FV_AtmIC(1)%v,FV_AtmIC(1)%w,FV_AtmIC(1)%pt,FV_AtmIC(1)%delp,FV_AtmIC(1)%q, &
                        FV_AtmIC(1)%phis, FV_AtmIC(1)%ps,FV_AtmIC(1)%pe, FV_AtmIC(1)%peln,FV_AtmIC(1)%pk,FV_AtmIC(1)%pkz, &
                        FV_AtmIC(1)%uc,FV_AtmIC(1)%vc, FV_AtmIC(1)%ua,FV_AtmIC(1)%va,        & 
                        FV_AtmIC(1)%ak, FV_AtmIC(1)%bk, FV_AtmIC(1)%gridstruct, FV_AtmIC(1)%flagstruct,&
                        FV_AtmIC(1)%npx, FV_AtmIC(1)%npy, FV_AtmIC(1)%npz, FV_AtmIC(1)%ng, &
                        FV_AtmIC(1)%flagstruct%ncnst, FV_AtmIC(1)%flagstruct%nwat,  &
                        FV_AtmIC(1)%flagstruct%ndims, FV_AtmIC(1)%flagstruct%ntiles, &
                        FV_AtmIC(1)%flagstruct%dry_mass, &
                        FV_AtmIC(1)%flagstruct%mountain,       &
                        FV_AtmIC(1)%flagstruct%moist_phys, FV_AtmIC(1)%flagstruct%hydrostatic, &
                        FV_AtmIC(1)%flagstruct%hybrid_z, FV_AtmIC(1)%delz, FV_AtmIC(1)%ze0, &
                        FV_AtmIC(1)%flagstruct%adiabatic, FV_AtmIC(1)%ks, FV_AtmIC(1)%neststruct%npx_global, &
                        FV_AtmIC(1)%ptop, FV_AtmIC(1)%domain, FV_AtmIC(1)%tile, FV_AtmIC(1)%bd )

        !Copy from temporary structure into state
        state%ud = FV_AtmIC(1)%u
        state%vd = FV_AtmIC(1)%v
        state%t = FV_AtmIC(1)%pt
        state%delp = FV_AtmIC(1)%delp
        state%q = FV_AtmIC(1)%q(:,:,:,1)
        state%phis = FV_AtmIC(1)%phis
        geom%ak = FV_AtmIC(1)%ak
        geom%ak = FV_AtmIC(1)%ak
        geom%ptop = FV_AtmIC(1)%ptop
        if (.not. state%hydrostatic) then
           state%w = FV_AtmIC(1)%w
           state%delz = FV_AtmIC(1)%delz
        endif

        !Deallocate temporary FV_Atm fv3 structure
        call deallocate_fv_atmos_type(FV_AtmIC(1))
        deallocate(FV_AtmIC)
        deallocate(grids_on_this_pe)

     Case ("dcmip-test-1-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
                                               phis0,ps,rho0,hum0,q1,q2,q3,q4)

              state%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                                                  phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 state%ud(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 state%vd(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
                 If (.not.state%hydrostatic) state%w(i,j,k) = w0
                 state%t(i,j,k) = t0
                 state%delp(i,j,k) = pe2-pe1
                 state%q (i,j,k) = hum0
                 state%qi(i,j,k) = q1
                 state%ql(i,j,k) = q2
                 state%o3(i,j,k) = q3
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                          t0,phis0,ps,rho0,hum0,q1)

              state%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                             t0,phis0,ps,rho0,hum0,q1)

                 state%ud(i,j,k) = u0 !ATTN comment above
                 state%vd(i,j,k) = v0
                 If (.not.state%hydrostatic) state%w(i,j,k) = w0
                 state%t(i,j,k) = t0
                 state%delp(i,j,k) = pe2-pe1
                 state%q (i,j,k) = hum0
                 state%qi(i,j,k) = q1
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                      t0,phis0,ps,rho0,hum0)

              state%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0)

                 state%ud(i,j,k) = u0 !ATTN comment above
                 state%vd(i,j,k) = v0
                 If (.not.state%hydrostatic) state%w(i,j,k) = w0
                 state%t(i,j,k) = t0
                 state%delp(i,j,k) = pe2-pe1
                 state%q(i,j,k) = hum0
                 
              enddo
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

              state%phis(i,j) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

                 state%ud(i,j,k) = u0 !ATTN comment above
                 state%vd(i,j,k) = v0
                 If (.not.state%hydrostatic) state%w(i,j,k) = w0
                 state%t(i,j,k) = t0
                 state%delp(i,j,k) = pe2-pe1
                 state%q(i,j,k) = hum0
                 
              enddo
           enddo
        enddo

     Case Default

        call invent_state(state,c_conf,geom)

     End Select int_option
        
end subroutine analytic_IC
  
! ------------------------------------------------------------------------------
subroutine invent_state(state,config,geom)

use fv3jedi_kinds_mod

implicit none

type(fv3jedi_state), intent(inout) :: state    !< Model state
type(c_ptr), intent(in)            :: config  !< Configuration structure
type(fv3jedi_geom),  intent(in)    :: geom

integer :: i,j,k

!ud
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%ud(i,j,k) = cos(0.25*geom%grid_lon(i,j)) + cos(0.25*geom%grid_lat(i,j))
    enddo
  enddo
enddo

!vd
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%vd(i,j,k) = 1.0_kind_real
    enddo
  enddo
enddo

!ua
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%ua(i,j,k) = cos(0.25*geom%grid_lon(i,j)) + cos(0.25*geom%grid_lat(i,j))
    enddo
  enddo
enddo

!va
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%va(i,j,k) = 1.0_kind_real
    enddo
  enddo
enddo

!t
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%t(i,j,k) = cos(0.25*geom%grid_lon(i,j)) + cos(0.25*geom%grid_lat(i,j))
    enddo
  enddo
enddo

!delp
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%delp(i,j,k) = k
    enddo
  enddo
enddo

!q
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%q(i,j,k) = 0.0
    enddo
  enddo
enddo

!qi
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%q(i,j,k) = 0.0
    enddo
  enddo
enddo

!ql
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%q(i,j,k) = 0.0
    enddo
  enddo
enddo

!o3
do k = 1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec
      state%q(i,j,k) = 0.0
    enddo
  enddo
enddo

return
end subroutine invent_state

! ------------------------------------------------------------------------------

subroutine read_file(geom, state, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_state), intent(inout) :: state      !< State
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call read_fms_restart(geom, state, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call read_geos_restart(state, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_state read: restart type not supported")
  endif

  return

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, state, c_conf, vdate)

  implicit none

  type(fv3jedi_geom), intent(inout)  :: geom     !< Geometry
  type(fv3jedi_state), intent(in)    :: state      !< State
  type(c_ptr), intent(in)            :: c_conf   !< Configuration
  type(datetime), intent(inout)      :: vdate    !< DateTime

  character(len=10) :: restart_type

  restart_type = config_get_string(c_conf,len(restart_type),"restart_type")

  if (trim(restart_type) == 'gfs') then
     call write_fms_restart(geom, state, c_conf, vdate)
  elseif (trim(restart_type) == 'geos') then
     call write_geos_restart(geom, state, c_conf, vdate)
  else
     call abor1_ftn("fv3jedi_state write: restart type not supported")
  endif

  return

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(state, nf, pstat)
implicit none
type(fv3jedi_state), intent(in) :: state
integer, intent(in) :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)

integer :: isc, iec, jsc, jec, gs

!1. Min
!2. Max
!3. RMS

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec

gs = (iec-isc+1)*(jec-jsc+1)*state%npz

!ud
if (allocated(state%ud)) then
  pstat(1,1) = minval(state%ud(isc:iec,jsc:jec,:))
  pstat(2,1) = maxval(state%ud(isc:iec,jsc:jec,:))
  pstat(3,1) = sqrt((sum(state%ud(isc:iec,jsc:jec,:))/gs)**2)
endif

!vd
if (allocated(state%vd)) then
  pstat(1,2) = minval(state%vd(isc:iec,jsc:jec,:))
  pstat(2,2) = maxval(state%vd(isc:iec,jsc:jec,:))
  pstat(3,2) = sqrt((sum(state%vd(isc:iec,jsc:jec,:))/gs)**2)
endif

!ua
if (allocated(state%ua)) then
  pstat(1,1) = minval(state%ua(isc:iec,jsc:jec,:))
  pstat(2,1) = maxval(state%ua(isc:iec,jsc:jec,:))
  pstat(3,1) = sqrt((sum(state%ua(isc:iec,jsc:jec,:))/gs)**2)
endif

!va
if (allocated(state%va)) then
  pstat(1,2) = minval(state%va(isc:iec,jsc:jec,:))
  pstat(2,2) = maxval(state%va(isc:iec,jsc:jec,:))
  pstat(3,2) = sqrt((sum(state%va(isc:iec,jsc:jec,:))/gs)**2)
endif

!t
if (allocated(state%t)) then
  pstat(1,3) = minval(state%t(isc:iec,jsc:jec,:))
  pstat(2,3) = maxval(state%t(isc:iec,jsc:jec,:))
  pstat(3,3) = sqrt((sum(state%t(isc:iec,jsc:jec,:))/gs)**2)
endif

!delp
if (allocated(state%delp)) then
  pstat(1,4) = minval(state%delp(isc:iec,jsc:jec,:))
  pstat(2,4) = maxval(state%delp(isc:iec,jsc:jec,:))
  pstat(3,4) = sqrt((sum(state%delp(isc:iec,jsc:jec,:))/gs)**2)
endif

!q
if (allocated(state%q)) then
  pstat(1,5) = minval(state%q(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(state%q(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(state%q(isc:iec,jsc:jec,:))/gs)**2)
endif

!qi
if (allocated(state%qi)) then
  pstat(1,5) = minval(state%qi(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(state%qi(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(state%qi(isc:iec,jsc:jec,:))/gs)**2)
endif

!ql
if (allocated(state%ql)) then
  pstat(1,5) = minval(state%ql(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(state%ql(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(state%ql(isc:iec,jsc:jec,:))/gs)**2)
endif

!o3
if (allocated(state%o3)) then
  pstat(1,5) = minval(state%o3(isc:iec,jsc:jec,:))
  pstat(2,5) = maxval(state%o3(isc:iec,jsc:jec,:))
  pstat(3,5) = sqrt((sum(state%o3(isc:iec,jsc:jec,:))/gs)**2)
endif

return

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine staterms(state, prms)
use fckit_mpi_module, only : fckit_mpi_comm, fckit_mpi_sum
implicit none
type(fv3jedi_state), intent(in) :: state
real(kind=kind_real), intent(out) :: prms

real(kind=kind_real) :: zz
integer i,j,k,ii,nt,ierr,npes,iisum
integer :: isc,iec,jsc,jec,npz
type(fckit_mpi_comm) :: f_comm

isc = state%isc
iec = state%iec
jsc = state%jsc
jec = state%jec
npz = state%npz

f_comm = fckit_mpi_comm()

zz = 0.0_kind_real
prms = 0.0_kind_real
ii = 0

!ud
if (allocated(state%ud)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%ud(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!vd
if (allocated(state%vd)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%vd(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ua
if (allocated(state%ua)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%ua(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!va
if (allocated(state%va)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%va(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!t
if (allocated(state%t)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%t(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!delp
if (allocated(state%delp)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%delp(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!q
if (allocated(state%q)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%q(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!qi
if (allocated(state%qi)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%qi(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!ql
if (allocated(state%ql)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%ql(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!o3
if (allocated(state%o3)) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        zz = zz + state%o3(i,j,k)**2
        ii = ii + 1
      enddo
    enddo
  enddo
endif

!Get global values
call f_comm%allreduce(zz,prms,fckit_mpi_sum())
call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

!if (ierr .ne. 0) then
!   print *,'error in staterms/mpi_allreduce, error code=',ierr
!endif
prms = sqrt(prms/real(iisum,kind_real))

end subroutine staterms

! ------------------------------------------------------------------------------

end module fv3jedi_state_mod
