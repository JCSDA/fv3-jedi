! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_linvarcha_nmcbal_mod

use fckit_configuration_module, only: fckit_configuration

use netcdf
use fv3jedi_netcdf_utils_mod, only: nccheck

use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_geom_mod, only: fv3jedi_geom
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use fv3jedi_kinds_mod
use fv3jedi_constants_mod, only: rad2deg, deg2rad, pi
use mpp_domains_mod, only: center

use pressure_vt_mod
use temperature_vt_mod
use moisture_vt_mod
use wind_vt_mod

use type_bump, only: bump_type
use type_gaugrid, only: gaussian_grid
use fv3jedi_field_mod
!use fv3jedi_bump_mod, only: bump_init, bump_apply, bump_apply_ad

implicit none
private

integer,parameter :: r_single = 4

public :: fv3jedi_linvarcha_nmcbal
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_linvarcha_nmcbal

 integer :: npe  ! size of PEs
 integer :: mype ! my rank

 !
 ! Cubed-sphere grid
 !
 integer :: isc,iec,jsc,jec,npz
 integer :: ngrid_cs ! data size for grid interpolation
 type(fv3jedi_geom) :: geom_cs
 ! increments
 real(kind=kind_real), allocatable :: wnd1(:,:,:)    ! ua/psi
 real(kind=kind_real), allocatable :: wnd2(:,:,:)    ! va/chi
 real(kind=kind_real), allocatable :: temp(:,:,:)    ! t/tv
 real(kind=kind_real), allocatable :: psrf(:,:,:)    ! ps
 real(kind=kind_real), allocatable :: ugfld(:,:,:) ! unstrucrtured grid field

 ! trajectores for varchange
 real(kind=kind_real), allocatable :: ttraj(:,:,:)    ! tempature
 real(kind=kind_real), allocatable :: tvtraj(:,:,:)   ! virtual tempature
 real(kind=kind_real), allocatable :: qtraj(:,:,:)    ! specific humidity
 real(kind=kind_real), allocatable :: qsattraj(:,:,:) ! saturation specific hum.

 !
 ! Gaussian grid
 !
 type(gaussian_grid) :: glb ! global Gaussian grid
 type(gaussian_grid) :: gg  ! local Gaussian grd
 integer :: nx,ny,nz,nv ! gobal grid size
 integer :: x2,y2,z2    ! local grid size with halo
 integer :: hx,hy       ! halo size
 integer :: layout(2)   ! domain layout
 integer,allocatable :: istrx(:),jstry(:) ! start index at each subdomain
 integer :: ngrid_gg, ngrid_gg1 ! data size for grid interpolation
 ! regression coeffs for GSI unbalanced variables
 character(len=256) :: path_to_nmcbalance_coeffs  ! path to netcdf file
 integer :: read_latlon_from_nc                   ! 1: read rlats/rlons of
                                                  !    global gaugrid from
                                                  !    netcdf
 real(kind=kind_real), allocatable :: agvz(:,:,:) ! coeffs. for psi and t 
 real(kind=kind_real), allocatable :: wgvz(:,:)   ! coeffs. for psi and ps
 real(kind=kind_real), allocatable :: bvz(:,:)    ! coeffs. for psi and chi
 ! data pointer
 real(kind=kind_real), pointer :: psi(:,:,:) => null() ! stream function 
 real(kind=kind_real), pointer :: chi(:,:,:) => null() ! velocity potential 
 real(kind=kind_real), pointer :: tv (:,:,:) => null() ! virtual tempature 
 real(kind=kind_real), pointer :: ps (:,:)   => null() ! surface pressure 

 !
 ! bump for grid interpolations
 !
 type(bump_type) :: bump_c2g ! bump for cubed-sphere to Gaussian grid
 type(bump_type) :: bump_g2c ! bump for Gaussian to cubed-sphere grid

end type fv3jedi_linvarcha_nmcbal

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, conf)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
type(fv3jedi_geom), target,     intent(in)    :: geom
type(fv3jedi_state), target,    intent(in)    :: bg
type(fv3jedi_state), target,    intent(in)    :: fg
type(fckit_configuration),      intent(in)    :: conf

real(kind=kind_real), pointer :: t   (:,:,:)
real(kind=kind_real), pointer :: q   (:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)

integer :: hx,hy,layoutx,layouty,read_latlon_from_nc
character(len=:),allocatable :: str

! Cubed-sphere grid geometry
call self%geom_cs%clone(geom,geom%fields)
self%npe = geom%f_comm%size()
self%mype = geom%f_comm%rank()

! 1. Get trajectores
!> Pointers to the background state
call pointer_field_array(bg%fields, 't'   , t)
call pointer_field_array(bg%fields, 'q'   , q)
call pointer_field_array(bg%fields, 'delp', delp)

!> Virtual temperature trajectory
allocate(self%tvtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
call T_to_Tv(geom,t,q,self%tvtraj)

!> Temperature trajectory
allocate(self%ttraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%ttraj = t

!> Specific humidity trajecotory
allocate(self%qtraj   (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
self%qtraj = q

!> Compute saturation specific humidity for q to RH transform
allocate(self%qsattraj(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

!> Compute saturation specific humidity
call get_qsat(geom,delp,t,q,self%qsattraj)

!> temporary fields for grid interpolation
self%isc = geom%isc; self%iec = geom%iec
self%jsc = geom%jsc; self%jec = geom%jec
self%npz = geom%npz
! ua/psi
allocate(self%wnd1(self%isc:self%iec,self%jsc:self%jec,self%npz))
! va/chi
allocate(self%wnd2(self%isc:self%iec,self%jsc:self%jec,self%npz))
! t/tv
allocate(self%temp(self%isc:self%iec,self%jsc:self%jec,self%npz))
! ps
allocate(self%psrf(self%isc:self%iec,self%jsc:self%jec,       1))
self%wnd1 = 0.0_kind_real; self%wnd2 = 0.0_kind_real
self%temp = 0.0_kind_real; self%psrf = 0.0_kind_real

! 2. Read balance operator grid from NetCDF
if (conf%has("path_to_nmcbalance_coeffs")) then
  call conf%get_or_die("path_to_nmcbalance_coeffs",str)
  self%path_to_nmcbalance_coeffs = str
  deallocate(str)
else
  call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.create:&
& path_to_nmcbalance_coeffs must be set in configuration file")
end if
call read_nmcbalance_grid(self)

! 3. Create Gaussian grid
self%hx = 1; self%hy = 1 ! default halo size
if(conf%has("hx")) then
  call conf%get_or_die("hx",hx); self%hx = hx
end if
if(conf%has("hy")) then
  call conf%get_or_die("hy",hy); self%hy = hy
end if

allocate(self%istrx(self%npe))
allocate(self%jstry(self%npe))

! set subdomain

if(conf%has("layoutx").and.conf%has("layouty")) then
  call conf%get_or_die("layoutx",layoutx)
  call conf%get_or_die("layouty",layouty)
  if(layoutx*layouty/=self%npe) then
    call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.create:&
& invalid layout values") 
  end if
  self%layout(1) = layoutx
  self%layout(2) = layouty
else
  if(self%npe==24)then
    self%layout(1)=6
    self%layout(2)=4
  else if(self%npe==18)then
    self%layout(1)=6
    self%layout(2)=3
  else if(self%npe==12)then
    self%layout(1)=4
    self%layout(2)=3
  else if(self%npe==6)then
    self%layout(1)=3
    self%layout(2)=2
  else
    self%layout(1)=self%npe
    self%layout(2)=1
  end if
end if
call deter_subdomain(self)

! set Gaussian grid parameters
self%nv = 4 ! psi,chi,tv,ps
self%z2 = self%nz
!self%z2 = (nz*3+1)/self%npe
!if(mod(nz*3+1,self%npe)/=0) self%z2 = self%z2 + 1
self%glb%nlat = self%ny; self%gg%nlat = self%y2
self%glb%nlon = self%nx; self%gg%nlon = self%x2
self%glb%nlev = self%nz; self%gg%nlev = self%z2
self%glb%nvar = self%nv; self%gg%nvar = self%nv

! allocate and initialize
call self%glb%create(); call self%gg%create()

! set Gaussian field pointer
call self%gg%fld3d_pointer(1,'psi',self%psi)
call self%gg%fld3d_pointer(2,'chi',self%chi)
call self%gg%fld3d_pointer(3,'tv' ,self%tv)
call self%gg%fld2d_pointer(4,'ps' ,self%ps)

! calculate and set Gaussian lat-lon
self%read_latlon_from_nc = 0
if(conf%has("read_latlon_from_nc")) then
  call conf%get_or_die("read_latlon_from_nc",read_latlon_from_nc)
  self%read_latlon_from_nc = read_latlon_from_nc
end if
if(self%read_latlon_from_nc==1) then
  call read_nmcbalance_latlon(self)
else
  call self%glb%calc_glb_latlon()
end if
call set_gaugrid_latlon(self)

! 4. Create interplation weights (bump)
call bump_init_gaugrid(self,geom)

! 5. Read balance coeffs from fixed file
allocate(self%agvz(self%gg%nlat,self%gg%nlev,self%gg%nlev))
allocate(self%wgvz(self%gg%nlat,self%gg%nlev))
allocate(self%bvz (self%gg%nlat,self%gg%nlev))
call read_nmcbalance_coeffs(self)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

! Deallacote everthing
if (allocated(self%tvtraj))   deallocate(self%tvtraj)
if (allocated(self%ttraj))    deallocate(self%ttraj)
if (allocated(self%qtraj))    deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

if (allocated(self%istrx))    deallocate(self%istrx)
if (allocated(self%jstry))    deallocate(self%jstry)
if (allocated(self%agvz))     deallocate(self%agvz)
if (allocated(self%wgvz))     deallocate(self%wgvz)
if (allocated(self%bvz))      deallocate(self%bvz)

if (allocated(self%wnd1))     deallocate(self%wnd1)
if (allocated(self%wnd2))     deallocate(self%wnd2)
if (allocated(self%temp))     deallocate(self%temp)
if (allocated(self%psrf))     deallocate(self%psrf)
if (allocated(self%ugfld))    deallocate(self%ugfld)

! Nullify pointers
if (associated(self%psi))     nullify(self%psi)
if (associated(self%chi))     nullify(self%chi)
if (associated(self%tv))      nullify(self%tv)
if (associated(self%ps))      nullify(self%ps)

call self%glb%delete()
call self%gg%delete()

call self%bump_c2g%dealloc()
call self%bump_g2c%dealloc()

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xuba,xbal)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xuba
type(fv3jedi_increment),        intent(inout) :: xbal

real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_rh
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_ps

real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_q
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ps

real(kind_real), allocatable :: psi_d(:,:,:), chi_d(:,:,:) ! psi/chi on D-grid

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xuba%fields,xbal%fields)

!Tangent linear of control to analysis variables
call pointer_field_array(xuba%fields, 'psi', xuba_psi)
call pointer_field_array(xuba%fields, 'chi', xuba_chi)
call pointer_field_array(xuba%fields, 'tv' , xuba_tv)
call pointer_field_array(xuba%fields, 'rh' , xuba_rh)
call pointer_field_array(xuba%fields, 'ps' , xuba_ps)
call pointer_field_array(xbal%fields, 'ua' , xbal_ua)
call pointer_field_array(xbal%fields, 'va' , xbal_va)
call pointer_field_array(xbal%fields, 't'  , xbal_t)
call pointer_field_array(xbal%fields, 'q'  , xbal_q)
call pointer_field_array(xbal%fields, 'ps' , xbal_ps)

! Apply Balance operator to psi,chi,tv(,ps)

! input cubed-sphere fields
self%wnd1 = xuba_psi
self%wnd2 = xuba_chi
self%temp = xuba_tv
self%psrf = xuba_ps

! 1. Intep to Gaussian grid (using bump)
call field_interp_to_gaugrid(self)

! 2. Apply balance operator
call balance(self)

! 3. Intep back to cubed-sphere grid (using bump)
call field_interp_from_gaugrid(self)

! 4. Linear variable changes (Control to Analysis)
!psi/chi -> ua/va
allocate(psi_d(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
allocate(chi_d(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
psi_d = 0.0_kind_real; chi_d = 0.0_kind_real
psi_d(geom%isc:geom%iec,geom%jsc:geom%jec,:) = self%wnd1
chi_d(geom%isc:geom%iec,geom%jsc:geom%jec,:) = self%wnd2
call psichi_to_uava(geom,psi_d,chi_d,xbal_ua,xbal_va)
deallocate(psi_d,chi_d)

!rh -> q
call rh_to_q_tl(geom,self%qsattraj,xuba_rh,xbal_q)

!Tv -> T
call tv_to_t_tl(geom,self%tvtraj,self%temp,self%qtraj,xbal_q,xbal_t)

!Ps
xbal_ps = self%psrf

! Copy calendar infomation
xbal%calendar_type = xuba%calendar_type
xbal%date_init = xuba%date_init

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xbal,xuba)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(inout) :: xbal
type(fv3jedi_increment),        intent(inout) :: xuba

real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_q
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ps

real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_rh
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_ps

real(kind_real), allocatable :: psi_d(:,:,:), chi_d(:,:,:) ! psi/chi on D-grid
real(kind_real), allocatable :: ps_not_copied(:,:,:)

! Copy fields that are the same in both
! -------------------------------------
allocate(ps_not_copied(geom%isc:geom%iec,geom%jsc:geom%jec,1))
call copy_field_array(xuba%fields, 'ps', ps_not_copied) ! temporary copy
call copy_subset(xbal%fields,xuba%fields)

!Tangent linear of control to analysis variables
call pointer_field_array(xbal%fields, 'ua' , xbal_ua)
call pointer_field_array(xbal%fields, 'va' , xbal_va)
call pointer_field_array(xbal%fields, 't'  , xbal_t)
call pointer_field_array(xbal%fields, 'q'  , xbal_q)
call pointer_field_array(xbal%fields, 'ps' , xbal_ps)
call pointer_field_array(xuba%fields, 'psi', xuba_psi)
call pointer_field_array(xuba%fields, 'chi', xuba_chi)
call pointer_field_array(xuba%fields, 'tv' , xuba_tv)
call pointer_field_array(xuba%fields, 'rh' , xuba_rh)
call pointer_field_array(xuba%fields, 'ps' , xuba_ps)

xuba_ps(:,:,1) = ps_not_copied(:,:,1)
deallocate(ps_not_copied)

self%wnd1 = 0.0_kind_real
self%wnd2 = 0.0_kind_real
self%temp = 0.0_kind_real
self%psrf = 0.0_kind_real

! 1. Adjoint of linear variable changes (Control to Analysis)
!Ps
self%psrf = xbal_ps; xbal_ps = 0.0_kind_real

!Tv -> T
call tv_to_t_ad(geom,self%tvtraj,self%temp,self%qtraj,xbal_q,xbal_t)

!rh -> q
call rh_to_q_ad(geom,self%qsattraj,xuba_rh,xbal_q)

!psi/chi -> ua/va
allocate(psi_d(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
allocate(chi_d(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
psi_d = 0.0_kind_real; chi_d = 0.0_kind_real
call psichi_to_uava_adm(geom,psi_d,chi_d,xbal_ua,xbal_va)
self%wnd1 = psi_d(geom%isc:geom%iec,geom%jsc:geom%jec,:)
self%wnd2 = chi_d(geom%isc:geom%iec,geom%jsc:geom%jec,:)
deallocate(psi_d,chi_d)

! 2. Adjoint of intep. back to cubed-sphere grid
call field_interp_from_gaugrid_ad(self)

! 3. Adjoint of applying balance operator
call tbalance(self)

! 4. Adjoint of intep. to Gaussian grid
call field_interp_to_gaugrid_ad(self)

xuba_ps   = xuba_ps  + self%psrf; self%psrf = 0.0_kind_real
xuba_tv   = xuba_tv  + self%temp; self%temp = 0.0_kind_real
xuba_psi  = xuba_psi + self%wnd1; self%wnd1 = 0.0_kind_real
xuba_chi  = xuba_chi + self%wnd2; self%wnd2 = 0.0_kind_real

! Copy calendar infomation
xuba%calendar_type = xbal%calendar_type
xuba%date_init = xbal%date_init

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xbal,xuba)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xbal
type(fv3jedi_increment),        intent(inout) :: xuba

real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_q

real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_rh


! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xbal%fields,xuba%fields)

!Tangent linear of control to analysis variables
call pointer_field_array(xbal%fields, 'ua' , xbal_ua)
call pointer_field_array(xbal%fields, 'va' , xbal_va)
call pointer_field_array(xbal%fields, 't'  , xbal_t)
call pointer_field_array(xbal%fields, 'q'  , xbal_q)
call pointer_field_array(xuba%fields, 'psi', xuba_psi)
call pointer_field_array(xuba%fields, 'chi', xuba_chi)
call pointer_field_array(xuba%fields, 'tv' , xuba_tv)
call pointer_field_array(xuba%fields, 'rh' , xuba_rh)

xuba_psi = xbal_ua
xuba_chi = xbal_va
xuba_tv  = xbal_t
xuba_rh  = xbal_q

! Copy calendar infomation
xuba%calendar_type = xbal%calendar_type
xuba%date_init = xbal%date_init

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xuba,xbal)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(in)    :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_increment),        intent(in)    :: xuba
type(fv3jedi_increment),        intent(inout) :: xbal

real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_psi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_chi
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_tv
real(kind=kind_real), pointer, dimension(:,:,:) :: xuba_rh

real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_t
real(kind=kind_real), pointer, dimension(:,:,:) :: xbal_q

! Copy fields that are the same in both
! -------------------------------------
call copy_subset(xuba%fields,xbal%fields)

!Tangent linear of control to analysis variables
call pointer_field_array(xuba%fields, 'psi', xuba_psi)
call pointer_field_array(xuba%fields, 'chi', xuba_chi)
call pointer_field_array(xuba%fields, 'tv' , xuba_tv)
call pointer_field_array(xuba%fields, 'rh' , xuba_rh)
call pointer_field_array(xbal%fields, 'ua' , xbal_ua)
call pointer_field_array(xbal%fields, 'va' , xbal_va)
call pointer_field_array(xbal%fields, 't'  , xbal_t)
call pointer_field_array(xbal%fields, 'q'  , xbal_q)

xbal_ua = xuba_psi
xbal_va = xuba_chi
xbal_t  = xuba_tv
xbal_q  = xuba_rh

! Copy calendar infomation
xbal%calendar_type = xuba%calendar_type
xbal%date_init = xuba%date_init

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------
! Determine subdomain settings 
! This subroutine was migrated from  gsi/general_sub2grid_mod.F90
subroutine deter_subdomain(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

if(self%layout(2)==1)then
  call deter_subdomain_nolayout(self)
else
  call deter_subdomain_withlayout(self)
end if

end subroutine deter_subdomain

! -------------------------------------

subroutine deter_subdomain_withlayout(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout)    :: self

integer:: nxpe,nype
integer,allocatable:: imxy(:), jmxy(:)
integer:: i,j,k,istart0, jstart0
integer:: im,jm,mm1

integer,allocatable,dimension(:) :: jlon1,ilat1

nxpe=self%layout(1)
nype=self%layout(2)
allocate(imxy(nxpe))
allocate(jmxy(nype))

! start

im=self%nx
jm=self%ny
call get_local_dims_ ( im,imxy,nxpe )
call get_local_dims_ ( jm,jmxy,nype )

allocate(jlon1(self%npe))
allocate(ilat1(self%npe))

! compute subdomain boundaries  (axis indices)

! compute local subdomain (offset and sizes)

K=0
jstart0 = 1
self%istrx=1
self%jstry=1
do j=1,nype
  istart0 = 1
  if (j>1) then
    jstart0 = jstart0 + jmxy(J-1)
  end if
  do i=1,nxpe
    k = k + 1
    ilat1(k) = jmxy(j)
    self%jstry(k) = jstart0
    jlon1(k) = imxy(i)
    if (i>1) then
      istart0 = istart0 + imxy(i-1)
    end if
    self%istrx(k) = istart0
  end do
end do

! Set number of latitude and longitude for given subdomain
mm1=self%mype+1
self%y2=ilat1(mm1)+self%hy*2
self%x2=jlon1(mm1)+self%hx*2
deallocate(ilat1,jlon1)
deallocate(imxy,jmxy)

end subroutine deter_subdomain_withlayout

! -------------------------------------

subroutine deter_subdomain_nolayout(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

integer :: npts,nrnc,iinum,iileft,jrows,jleft,k,i,jjnum
integer :: j,mm1,iicnt,ipts,jjleft
integer,allocatable,dimension(:) :: jlon1,ilat1
integer,allocatable,dimension(:) :: iiend,jjend,iistart
real(kind_real) :: anperpe

! Compute number of points on full grid
! and target number of point per mpi task (pe)
npts=self%nx*self%ny
if(self%npe<=0) then
  call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.deter_subdomain_nolayout:&
& npe has invalid value")
end if
anperpe=float(npts)/float(self%npe)

! Start with square subdomains
nrnc=int(sqrt(anperpe))
iinum=self%nx/nrnc
if(iinum==0) iinum=1
iicnt=self%nx/iinum
iileft=self%nx-iicnt*iinum
jrows=self%npe/iinum
jleft=self%npe-jrows*iinum

! Adjust subdomain boundaries
k=0
self%istrx=1
self%jstry=1
allocate(iistart(self%npe+1))
allocate(iiend(self%npe+1))
allocate(jjend(self%npe+1))
allocate(jlon1(self%npe))
allocate(ilat1(self%npe))
iistart(1)=1
do i=1,iinum
  ipts = iicnt
  if(i <= iileft)ipts=ipts+1
  iiend(i)=iistart(i)+ipts-1
  iistart(i+1)=iiend(i)+1
  jjnum=jrows
  if(i <= jleft)jjnum=jrows+1
  do j=1,jjnum
    k=k+1
    jlon1(k)=ipts
    self%istrx(k)= iistart(i)
    ilat1(k)=self%ny/jjnum
    jjleft=self%ny-ilat1(k)*jjnum
    if(j <= jjleft)ilat1(k)=ilat1(k)+1
    if(j > 1)self%jstry(k)=jjend(j-1)+1
    jjend(j)=self%jstry(k)+ilat1(k)-1
  end do
end do
deallocate(iistart,iiend,jjend)

! Set number of latitude and longitude for given subdomain
mm1=self%mype+1
self%y2=ilat1(mm1)+self%hy*2
self%x2=jlon1(mm1)+self%hx*2
deallocate(ilat1,jlon1)

end subroutine deter_subdomain_nolayout

! -------------------------------------

subroutine get_local_dims_(dim_world,dim,ndes)

implicit  none
integer,intent(in   ) :: dim_world, ndes
integer,intent(inout) :: dim(0:ndes-1)
integer :: n,im,rm

im = dim_world/ndes
rm = dim_world-ndes*im
do n=0,ndes-1
  dim(n) = im
  if( n<=rm-1 ) dim(n) = im+1
end do

end subroutine get_local_dims_

! ------------------------------------------------------------------------------
! Set local Gausisan latitudes and longitudes
subroutine set_gaugrid_latlon(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
integer :: i,j,ix,jy

! Latitudes
jy=self%jstry(self%mype+1)-self%hy
do j=1,self%gg%nlat
  if(jy<=0) jy=1
  if(jy>=self%glb%nlat+1) jy=self%glb%nlat
  self%gg%rlats(j)=self%glb%rlats(jy)
  jy = self%jstry(self%mype+1)-self%hy+j
end do

! Longitudes
ix=self%istrx(self%mype+1)-self%hx
do i=1,self%gg%nlon
  if(ix<=0) ix=ix+self%glb%nlon
  if(ix>self%glb%nlon) ix=ix-self%glb%nlon
  self%gg%rlons(i) = self%glb%rlons(ix)
  ix=self%istrx(self%mype+1)-self%hx+i
end do

end subroutine set_gaugrid_latlon

! ------------------------------------------------------------------------------
! Setup bump for grid interpolation
subroutine bump_init_gaugrid(self,geom)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
type(fv3jedi_geom), target,     intent(in)    :: geom

real(kind=kind_real), allocatable :: lon_cs(:,:) ! Cubed-sphere longitudes
real(kind=kind_real), allocatable :: lat_cs(:,:) ! Cubed-sphere latitudes
real(kind=kind_real), allocatable :: lon_gg(:,:) ! Gaussian longitudes
real(kind=kind_real), allocatable :: lat_gg(:,:) ! Gaussian latitudes
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

integer :: i,j,ij

if(  (geom%isc/=self%isc).or.(geom%iec/=self%iec) &
&.or.(geom%jsc/=self%jsc).or.(geom%jec/=self%jec)) then
  call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.bump_init_gaugrid:&
& grid does not match that in the geometry")
end if

self%ngrid_cs = (geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)
allocate(lon_cs(self%ngrid_cs,1))
allocate(lat_cs(self%ngrid_cs,1))
ij=0
do j=geom%jsc,geom%jec
  do i=geom%isc,geom%iec
    ij = ij+1
    lon_cs(ij,1)=geom%grid_lon(i,j)*rad2deg
    lat_cs(ij,1)=geom%grid_lat(i,j)*rad2deg
  end do
end do

! bump setup for cubed-shere to Gaussian grid
self%ngrid_gg=self%gg%nlat*self%gg%nlon ! packed data size with halo
allocate(lon_gg(self%ngrid_gg,1))
allocate(lat_gg(self%ngrid_gg,1))
ij=0
do j=1,self%gg%nlat
  do i=1,self%gg%nlon
    ij = ij+1
    lon_gg(ij,1)=self%gg%rlons(i)*rad2deg
    lat_gg(ij,1)=self%gg%rlats(j)*rad2deg
  end do
end do

! allocate unstructured gird field
allocate(self%ugfld(self%ngrid_cs,self%npz,4))
self%ugfld = 0.0_kind_real

!call bump_init(geom,self%ngrid_gg,lat_gg,lon_gg,self%bump_c2g)
call self%bump_c2g%nam%init(self%npe)
self%bump_c2g%nam%prefix = trim("dummy")
self%bump_c2g%nam%default_seed = .true.
self%bump_c2g%nam%new_obsop = .true.
self%bump_c2g%nam%write_obsop = .false.
self%bump_c2g%nam%verbosity = "none"
allocate(area(self%ngrid_cs))          ; area  = 1.0_kind_real
allocate(vunit(self%ngrid_cs,self%npz)); vunit = 1.0_kind_real
allocate(lmask(self%ngrid_cs,self%npz)); lmask = .true.
call self%bump_c2g%setup_online(              &
&     geom%f_comm,self%ngrid_cs,self%npz,1,1, &
&     lon_cs,lat_cs,area,vunit,lmask,         &
&     nobs=self%ngrid_gg,                     &
&     lonobs=lon_gg(:,1),latobs=lat_gg(:,1))
call self%bump_c2g%run_drivers()
deallocate(area,vunit,lmask)
deallocate(lon_gg,lat_gg)

! bump setup for Gaussian to cubed-sphere grid
self%ngrid_gg1=(self%gg%nlat-self%hy*2)*(self%gg%nlon-self%hx*2) ! paked data size without halo
allocate(lon_gg(self%ngrid_gg1,1))
allocate(lat_gg(self%ngrid_gg1,1))
ij=0
do j=1+self%hy,self%gg%nlat-self%hy
  do i=1+self%hx,self%gg%nlon-self%hx
    ij = ij+1
    lon_gg(ij,1)=self%gg%rlons(i)*rad2deg
    lat_gg(ij,1)=self%gg%rlats(j)*rad2deg
  end do
end do

!call bump_init(geom_gg1,self%ngrid_cs,lat_cs,lon_cs,self%bump_g2c)
call self%bump_g2c%nam%init(self%npe)
self%bump_g2c%nam%prefix = trim("dummy")
self%bump_g2c%nam%default_seed = .true.
self%bump_g2c%nam%new_obsop = .true.
self%bump_g2c%nam%write_obsop = .false.
self%bump_g2c%nam%verbosity = "none"
allocate(area(self%ngrid_gg1))         ; area  = 1.0_kind_real
allocate(vunit(self%ngrid_gg1,self%nz)); vunit = 1.0_kind_real
allocate(lmask(self%ngrid_gg1,self%nz)); lmask = .true.
call self%bump_g2c%setup_online(              &
&     geom%f_comm,self%ngrid_gg1,self%nz,1,1, &
&     lon_gg,lat_gg,area,vunit,lmask,         &
&     nobs=self%ngrid_cs,                     &
&     lonobs=lon_cs(:,1),latobs=lat_cs(:,1))
call self%bump_g2c%run_drivers()
deallocate(area,vunit,lmask)
deallocate(lon_cs,lat_cs)
deallocate(lon_gg,lat_gg)

end subroutine bump_init_gaugrid

! ------------------------------------------------------------------------------
! read NMC balance operator girds from NetCDF
subroutine read_nmcbalance_grid(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

integer :: iunit,ierr,ncid,dimid
integer :: grid(3)
character(len=128) :: log_grid
character(len=16)  :: nlat_c,nlon_c,nlev_c

if(self%mype==0) then
! Open NetCDF
  call nccheck(nf90_open(trim(self%path_to_nmcbalance_coeffs),NF90_NOWRITE,ncid), &
 &                       "nf90_open "//trim(self%path_to_nmcbalance_coeffs))
! Get grid
  call nccheck(nf90_inq_dimid(ncid,"nlat",dimid), "nf90_inq_dimid nlat")
  call nccheck(nf90_inquire_dimension(ncid,dimid,len=grid(1)), &
 &            "nf90_inquire_dimension nlat" )
  call nccheck(nf90_inq_dimid(ncid,"nlon",dimid), "nf90_inq_dimid nlon")
  call nccheck(nf90_inquire_dimension(ncid,dimid,len=grid(2)), &
 &            "nf90_inquire_dimension nlon" )
  call nccheck(nf90_inq_dimid(ncid,"nlev",dimid), "nf90_inq_dimid nlev")
  call nccheck(nf90_inquire_dimension(ncid,dimid,len=grid(3)), &
&            "nf90_inquire_dimension nlev" )

! Close NetCDF
  call nccheck(nf90_close(ncid),"nf90_close")
end if
call self%geom_cs%f_comm%broadcast(grid,0)
if(grid(3)/=self%npz) then
  call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.read_nmcbalance_grid:&
& nlev must be the same as npz on FV3 geometry")
end if
self%ny = grid(1)
self%nx = grid(2)
self%nz = grid(3)
if(self%mype==0) then
  write(nlat_c,*) grid(1); write(nlon_c,*) grid(2); write(nlev_c,*) grid(3);
  write(log_grid,*) "NMC balance operalor grid size : nx=",&
 &                          trim(adjustl(nlon_c)), &
 &                   ",ny=",trim(adjustl(nlat_c)), &
 &                   ",nz=",trim(adjustl(nlev_c))
  write(6,'(A)') log_grid
end if

end subroutine read_nmcbalance_grid

! ------------------------------------------------------------------------------
! read latlon from NetCDF
subroutine read_nmcbalance_latlon(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

integer :: iunit,ierr,ncid,varid
integer :: mm1,jx,i,j,k

real(kind=kind_real),allocatable :: rlats(:)
real(kind=kind_real),allocatable :: rlons(:)

allocate(rlats(self%glb%nlat))
allocate(rlons(self%glb%nlon))

if(self%mype==0) then
! Open NetCDF
  call nccheck(nf90_open(trim(self%path_to_nmcbalance_coeffs),NF90_NOWRITE,ncid), &
  &                      "nf90_open "//trim(self%path_to_nmcbalance_coeffs))

! Get coefficients
  call nccheck(nf90_inq_varid(ncid,"lats",varid), "nf90_inq_varid lats")
  call nccheck(nf90_get_var(ncid,varid,rlats),"nf90_get_var lats" )
  call nccheck(nf90_inq_varid(ncid,"lons",varid), "nf90_inq_varid lons")
  call nccheck(nf90_get_var(ncid,varid,rlons),"nf90_get_var lons" )

! Close NetCDF
  call nccheck(nf90_close(ncid),"nf90_close")
end if
call self%geom_cs%f_comm%broadcast(rlats,0)
call self%geom_cs%f_comm%broadcast(rlons,0)

do i=1,self%glb%nlon
  self%glb%rlons(i) = rlons(i)*deg2rad
end do
do i=1,self%glb%nlat
  self%glb%rlats(i) = rlats(i)*deg2rad
end do

end subroutine read_nmcbalance_latlon

! ------------------------------------------------------------------------------
! read NMC balance cofficients from NetCDF
subroutine read_nmcbalance_coeffs(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

integer :: iunit,ierr,ncid,varid
integer :: mm1,jx,i,j,k

real(kind=r_single),allocatable :: agvin(:,:,:,:)
real(kind=r_single),allocatable :: wgvin(:,:,:)
real(kind=r_single),allocatable :: bvin(:,:,:)

allocate(agvin(self%glb%nlat,self%glb%nlev,self%glb%nlev,1))
allocate(wgvin(self%glb%nlat,self%glb%nlev,1))
allocate(bvin (self%glb%nlat,self%glb%nlev,1))

if(self%mype==0) then
! Open NetCDF
  call nccheck(nf90_open(trim(self%path_to_nmcbalance_coeffs),NF90_NOWRITE,ncid), &
 &                       "nf90_open"//trim(self%path_to_nmcbalance_coeffs))
! Get coefficients
  call nccheck(nf90_inq_varid(ncid,"agvz",varid), "nf90_inq_varid agvz")
  call nccheck(nf90_get_var(ncid,varid,agvin),"nf90_get_var agvz" )
  call nccheck(nf90_inq_varid(ncid,"wgvz",varid), "nf90_inq_varid wgvz")
  call nccheck(nf90_get_var(ncid,varid,wgvin),"nf90_get_var wgvz" )
  call nccheck(nf90_inq_varid(ncid,"bvz",varid), "nf90_inq_varid bvz")
  call nccheck(nf90_get_var(ncid,varid,bvin),"nf90_get_var bvz" )

! Close NetCDF
  call nccheck(nf90_close(ncid),"nf90_close")
end if
call self%geom_cs%f_comm%broadcast(agvin,0)
call self%geom_cs%f_comm%broadcast(wgvin,0)
call self%geom_cs%f_comm%broadcast(bvin,0)
self%agvz = 0.0_kind_real
self%bvz  = 0.0_kind_real
self%wgvz = 0.0_kind_real
mm1=self%mype+1
do k=1,self%gg%nlev
  do j=1,self%gg%nlev
    do i=1,self%gg%nlat
      jx=self%jstry(mm1)+i-2
      jx=max(jx,2)
      jx=min(self%glb%nlat-1,jx)
      self%agvz(i,j,k)=agvin(jx,j,k,1)
    end do
  end do
  do i=1,self%gg%nlat
    jx=self%jstry(mm1)+i-2
    jx=max(jx,2)
    jx=min(self%glb%nlat-1,jx)
    self%wgvz(i,k)=wgvin(jx,k,1)
    self%bvz(i,k)=bvin(jx,k,1)
  end do
end do
deallocate(agvin,wgvin,bvin)

end subroutine read_nmcbalance_coeffs

! ------------------------------------------------------------------------------
! Pack increments
subroutine pack_inc_all(self,input_to_zero)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
logical,                        intent(in)    :: input_to_zero

integer :: imga,jx,jy,jl

if(.not.allocated(self%ugfld)) allocate(self%ugfld(self%ngrid_cs,self%npz,4))
self%ugfld = 0.0_kind_real
imga = 0
do jy=self%jsc,self%jec
  do jx=self%isc,self%iec
    imga = imga+1
    do jl=1,self%npz
      self%ugfld(imga,jl,1) = self%wnd1(jx,jy,jl)
      self%ugfld(imga,jl,2) = self%wnd2(jx,jy,jl)
      self%ugfld(imga,jl,3) = self%temp(jx,jy,jl)
    end do
    self%ugfld(imga,1,4) = self%psrf(jx,jy,1)
  end do
end do

if(input_to_zero)then
  self%wnd1 = 0.0_kind_real
  self%wnd2 = 0.0_kind_real
  self%temp = 0.0_kind_real
  self%psrf = 0.0_kind_real
end if

end subroutine pack_inc_all

! ------------------------------------------------------------------------------
! Unpack increments
subroutine unpack_inc_all(self,input_to_zero,do_adj)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self
logical,                        intent(in)    :: input_to_zero
logical,                        intent(in)    :: do_adj

integer :: imga,jx,jy,jl

if(.not.allocated(self%ugfld)) then
  call abor1_ftn("fv3jedi_linvarcha_nmcbal_mod.unpack_inc_all:&
& ugfld is not allocated")
end if

if(.not.do_adj) then
  self%wnd1 = 0.0_kind_real
  self%wnd2 = 0.0_kind_real
  self%temp = 0.0_kind_real
  self%psrf = 0.0_kind_real
end if

imga = 0
do jy=self%jsc,self%jec
  do jx=self%isc,self%iec
    imga = imga+1
    do jl=1,self%npz
      self%wnd1(jx,jy,jl) = self%wnd1(jx,jy,jl) + self%ugfld(imga,jl,1)
      self%wnd2(jx,jy,jl) = self%wnd2(jx,jy,jl) + self%ugfld(imga,jl,2)
      self%temp(jx,jy,jl) = self%temp(jx,jy,jl) + self%ugfld(imga,jl,3)
    end do
    self%psrf(jx,jy,1) = self%psrf(jx,jy,1) + self%ugfld(imga,1,4)
  end do
end do

if(input_to_zero) self%ugfld = 0.0_kind_real

end subroutine unpack_inc_all

! ------------------------------------------------------------------------------
! Interpolate cubed-sphere grid field to Gaussian gird using bump
subroutine field_interp_to_gaugrid(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

real(kind_real),allocatable :: ugfld_g(:,:) ! packed Gaussian gird field
integer :: i,j,k,n,ij

! Pack cubed-sphere field to unstrucured grid field
call pack_inc_all(self,.false.)

allocate(ugfld_g(self%ngrid_gg,self%gg%nlev))
do n=1,4
  ugfld_g=0.0_kind_real
  call self%bump_c2g%apply_obsop(self%ugfld(:,:,n),ugfld_g(:,:))
  do k=1,self%gg%nlev
    ij=0
    do j=1,self%gg%nlat
      do i=1,self%gg%nlon
        ij = ij+1
        self%gg%fld(j,i,k,n)=ugfld_g(ij,k)
      end do
    end do
  end do
end do

deallocate(ugfld_g)

end subroutine field_interp_to_gaugrid

! ------------------------------------------------------------------------------
! Interpolate Gaussian gird filed to cubed-sphere grid using bump
subroutine field_interp_from_gaugrid(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

real(kind_real),allocatable :: ugfld_g(:,:) ! packed Gaussian grid field
integer :: i,j,k,n,ij

allocate(ugfld_g(self%ngrid_gg1,self%gg%nlev))
do n=1,4
  ugfld_g=0.0_kind_real
  do k=1,self%gg%nlev
    ij=0
    do j=1+self%hy,self%gg%nlat-self%hy
      do i=1+self%hx,self%gg%nlon-self%hx
        ij = ij+1
        ugfld_g(ij,k)=self%gg%fld(j,i,k,n)
      end do
    end do
  end do
  call self%bump_g2c%apply_obsop(ugfld_g(:,:),self%ugfld(:,:,n))
end do

! Unack unstrucured grid field to cubed-sphere field
call unpack_inc_all(self,.false.,.false.)

deallocate(ugfld_g)

end subroutine field_interp_from_gaugrid

! ------------------------------------------------------------------------------

subroutine field_interp_to_gaugrid_ad(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

real(kind_real),allocatable :: ugfld_g(:,:) ! packed Gaussian grid field
integer :: i,j,k,n,ij

allocate(ugfld_g(self%ngrid_gg,self%gg%nlev))
do n=1,4
  ugfld_g=0.0_kind_real
  do k=1,self%gg%nlev
    ij=0
    do j=1,self%gg%nlat
      do i=1,self%gg%nlon
        ij = ij+1
        ugfld_g(ij,k)=self%gg%fld(j,i,k,n)
      end do
    end do
  end do
  call self%bump_c2g%apply_obsop_ad(ugfld_g(:,:),self%ugfld(:,:,n))
end do

! Adjoint of paking
call unpack_inc_all(self,.true.,.true.)

deallocate(ugfld_g)

end subroutine field_interp_to_gaugrid_ad

! ------------------------------------------------------------------------------

subroutine field_interp_from_gaugrid_ad(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout) :: self

real(kind_real),allocatable :: ugfld_g(:,:) ! paked Gaussian grid field
integer :: i,j,k,n,ij

! Adjoint of unpaking
call pack_inc_all(self,.true.)

allocate(ugfld_g(self%ngrid_gg1,self%gg%nlev))
self%gg%fld=0.0_kind_real
do n=1,4
  ugfld_g=0.0_kind_real
  call self%bump_g2c%apply_obsop_ad(self%ugfld(:,:,n),ugfld_g(:,:))
  do k=1,self%gg%nlev
    ij=0
    do j=1+self%hy,self%gg%nlat-self%hy
      do i=1+self%hx,self%gg%nlon-self%hx
        ij = ij+1
!        self%gg%fld(j,i,k,n)=self%gg%fld(j,i,k,n)+ugfld_g(ij,k)
        self%gg%fld(j,i,k,n)=ugfld_g(ij,k)
      end do
    end do
  end do
end do

deallocate(ugfld_g)

end subroutine field_interp_from_gaugrid_ad

! ------------------------------------------------------------------------------
! Apply balance operator
! This subroutine was migrated from gsi/balmod.F90
subroutine balance(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout)    :: self
integer :: i,j,k,l

! Add contribution from streamfunction and unbalanced vp
! to surface pressure.
do k=1,self%gg%nlev
  do i=1,self%gg%nlon
    do j=1,self%gg%nlat
      self%ps(j,i)=self%ps(j,i)+self%wgvz(j,k)*self%psi(j,i,k)
    end do
  end do
end do

do k=1,self%gg%nlev
! Add contribution from streamfunction to veloc. potential
  do i=1,self%gg%nlon
    do j=1,self%gg%nlat
      self%chi(j,i,k)=self%chi(j,i,k)+self%bvz(j,k)*self%psi(j,i,k)
    end do
  end do

! Add contribution from streamfunction to unbalanced temperature.
  do l=1,self%gg%nlev
    do i=1,self%gg%nlon
      do j=1,self%gg%nlat
        self%tv(j,i,k)=self%tv(j,i,k)+self%agvz(j,k,l)*self%psi(j,i,l)
      end do
    end do
  end do
end do

end subroutine balance

! ------------------------------------------------------------------------------
! Adjoint of balance operator
subroutine tbalance(self)

implicit none
type(fv3jedi_linvarcha_nmcbal), intent(inout)    :: self
integer :: i,j,k,l

do k=1,self%gg%nlev
! Adjoint of contribution to temperature from streamfunction.
  do l=1,self%gg%nlev
    do i=1,self%gg%nlon
      do j=1,self%gg%nlat
        self%psi(j,i,k)=self%psi(j,i,k)+self%agvz(j,l,k)*self%tv(j,i,l)
      end do
    end do
  end do

! Adjoint of contribution to velocity potential from streamfunction.
  do i=1,self%gg%nlon
    do j=1,self%gg%nlat
      self%psi(j,i,k)=self%psi(j,i,k)+self%bvz(j,k)*self%chi(j,i,k)
    end do
  end do
end do

! Adjoint of streamfunction and unbalanced velocity potential
! contribution to surface pressure.
do k=1,self%gg%nlev
  do i=1,self%gg%nlon
    do j=1,self%gg%nlat
      self%psi(j,i,k)=self%psi(j,i,k)+self%wgvz(j,k)*self%ps(j,i)
    end do
  end do
end do

end subroutine tbalance

! ------------------------------------------------------------------------------
! Use to file format change from binary to NetCDF for more general reading
subroutine bi2nc(path_to_gsi_bal_coeffs)
implicit none

character(len=*),               intent(in   ) :: path_to_gsi_bal_coeffs

integer :: iunit,ierr
integer :: nlev,nlat,nlon
integer :: fileopts, ncid, vc, varid(10)
integer :: x_dimid,y_dimid,z_dimid,t_dimid
integer,allocatable :: v_dimids(:)
real(kind=kind_real),allocatable :: lats(:),lons(:)
integer,allocatable :: zk(:)
integer :: i,j,k

character(len=64)  :: nlat_c,nlev_c
character(len=128) :: ncname
real(kind=r_single),allocatable :: agvin(:,:,:,:)
real(kind=r_single),allocatable :: wgvin(:,:,:)
real(kind=r_single),allocatable :: bvin(:,:,:)
type(gaussian_grid) :: gauss

iunit=get_fileunit()
open(iunit,file=trim(path_to_gsi_bal_coeffs),form='unformatted', &
&     convert='big_endian',status='old',iostat=ierr)
call check_iostat(ierr,"open balance coeffs header")
rewind(iunit)
read(iunit,iostat=ierr) nlev,nlat,nlon

gauss%nlat = nlat
gauss%nlon = nlon
gauss%nlev = nlev
gauss%nvar = 1
call gauss%create()
call gauss%calc_glb_latlon()

allocate(lats(nlat))
do j=1,nlat
  lats(j)=gauss%rlats(j)*rad2deg
end do
allocate(lons(nlon))
do i=1,nlon
  lons(i)=gauss%rlons(i)*rad2deg
end do
allocate(zk(nlev))
do k=1,nlev
  zk(k)=k
end do

allocate(agvin(nlat,nlev,nlev,1))
allocate(wgvin(nlat,nlev,1))
allocate(bvin (nlat,nlev,1))
read(iunit,iostat=ierr) agvin,bvin,wgvin
close(iunit)

fileopts = ior(NF90_NETCDF4, NF90_MPIIO)

ncid = 1
write(nlev_c,*) nlev; write(nlat_c,*) nlat
write(ncname,*) "global_berror.l",trim(adjustl(nlev_c)),"y",trim(adjustl(nlat_c)),".nc"
call nccheck(nf90_create(trim(adjustl(ncname)),fileopts,ncid),'nf90_create')

call nccheck(nf90_def_dim(ncid,'nlat',nlat,y_dimid),'nf90_def_dim nlat')
call nccheck(nf90_def_dim(ncid,'nlon',nlon,x_dimid),'nf90_def_dim nlon')
call nccheck(nf90_def_dim(ncid,'nlev',nlev,z_dimid),'nf90_def_dim nlev')
call nccheck(nf90_def_dim(ncid,'time',   1,t_dimid),'nf90_def_dim time')

vc = 1
call nccheck(nf90_def_var(ncid,"lats",NF90_DOUBLE,y_dimid,varid(vc)),"nf90_def_var lats")
call nccheck(nf90_put_att(ncid,varid(vc),"long_name","latitude"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","degrees_north"))
vc = 2
call nccheck(nf90_def_var(ncid,"lons",NF90_DOUBLE,x_dimid,varid(vc)),"nf90_def_var lons")
call nccheck(nf90_put_att(ncid,varid(vc),"long_name","longitude"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","degrees_west"))
vc = 3
call nccheck(nf90_def_var(ncid,"level",NF90_INT,z_dimid,varid(vc)),"nf90_def_var level")
call nccheck(nf90_put_att(ncid,varid(vc),"long_name","vertical level"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","none"))
vc = 4
call nccheck(nf90_def_var(ncid,"time",NF90_INT,t_dimid,varid(vc)),"nf90_def_var time")
call nccheck(nf90_put_att(ncid,varid(vc),"long_name","time"),"nf90_def_var time long_name")

vc = 5
allocate(v_dimids(4))
v_dimids(:) = (/y_dimid,z_dimid,z_dimid,t_dimid/)
call nccheck(nf90_def_var(ncid,'agvz',NF90_FLOAT,v_dimids,varid(vc)),'nf90_def_var agvz')
call nccheck(nf90_put_att(ncid,varid(vc), &
& "long_name","projection_of_streamfunction_onto_balanced_temperature"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","K s m-2"))
deallocate(v_dimids)
vc = 6
allocate(v_dimids(3))
v_dimids(:) = (/y_dimid,z_dimid,t_dimid/)
call nccheck(nf90_def_var(ncid,'bvz',NF90_FLOAT,v_dimids,varid(vc)),'nf90_def_var bvz')
call nccheck(nf90_put_att(ncid,varid(vc), &
& "long_name","projection_of_streamfunction_onto_balanced_velocity_potential"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","none"))
vc = 7
call nccheck(nf90_def_var(ncid,'wgvz',NF90_FLOAT,v_dimids,varid(vc)),'nf90_def_var wgvz')
call nccheck(nf90_put_att(ncid,varid(vc), &
& "long_name","projection_of_streamfunction_onto_balanced_surface_pressure"))
call nccheck(nf90_put_att(ncid,varid(vc),"units","Pa s m-2"))
deallocate(v_dimids)

call nccheck(nf90_enddef(ncid),"nf90_enddef")

vc = 1
call nccheck(nf90_put_var(ncid,varid(vc),lats),"nf90_put_var lats")
vc = 2
call nccheck(nf90_put_var(ncid,varid(vc),lons),"nf90_put_var lons")
vc = 3
call nccheck(nf90_put_var(ncid,varid(vc),zk),"nf90_put_var level")
vc = 4
call nccheck(nf90_put_var(ncid,varid(vc),1),"nf90_put_var time")

vc = 5
call nccheck(nf90_put_var(ncid,varid(vc),agvin),"nf90_put_var agvz")
vc = 6
call nccheck(nf90_put_var(ncid,varid(vc),bvin),"nf90_put_var bvz")
vc = 7
call nccheck(nf90_put_var(ncid,varid(vc),wgvin),"nf90_put_var wgvz")

call nccheck(nf90_close(ncid),'nf90_close')

deallocate(lats)
deallocate(zk)
deallocate(agvin,wgvin,bvin)

end subroutine bi2nc

! ------------------------------------------------------------------------------
! get unused file unit
integer function get_fileunit(iunit_in) result(iunit)

implicit none
integer,intent(in),optional :: iunit_in

logical :: lopened, lexist
integer,parameter :: IUNIT_STR=11, IUNIT_END=99
integer :: i

if(present(iunit_in)) then
  inquire(unit=iunit_in,opened=lopened,exist=lexist)
  if(lexist.and.(.not.lopened)) iunit=iunit_in
  return
end if

do i=IUNIT_STR,IUNIT_END
  inquire(unit=i,opened=lopened,exist=lexist)
  if(lexist.and.(.not.lopened)) then
    iunit=i
    exit
  end if
end do

end function get_fileunit

! ------------------------------------------------------------------------------
! check I/O return code
subroutine check_iostat(ierr,mess)

implicit none
integer,intent(in) :: ierr
character(len=*),intent(in) :: mess

character(len=256) :: abort_log

write(abort_log,'(3A,I3)') "IO Error : ",trim(mess)," : iostat=",ierr
if(ierr/=0) then
  call abor1_ftn(trim(adjustl(abort_log)))
end if

end subroutine check_iostat

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_nmcbal_mod
