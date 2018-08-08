
!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

!General JEDI uses
use kinds
use iso_c_binding
use config_mod
use netcdf

!Uses for setting up tracer numbers from field_table
use tracer_manager_mod, only: get_number_tracers
use field_manager_mod,  only: MODEL_ATMOS

!FMS/MPP uses
use mpp_domains_mod,    only: domain2D, mpp_deallocate_domain
use mpp_domains_mod,    only: mpp_define_layout, mpp_define_mosaic, mpp_define_io_domain
use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, NOTE

!Uses for generating geometry using FV3 routines
use fv_arrays_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
use fv_control_mod,     only: fv_init, pelist_all

implicit none
private
public :: fv3jedi_geom
public :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

!> Skinny version of fv_grid_bounds_type
type fv_grid_bounds_type
    integer :: isd, ied, jsd, jed ! data domain
    integer :: isc, iec, jsc, jec ! compute domain
end type fv_grid_bounds_type

!> Fortran derived type to hold geometry data for the FV3JEDI model
type :: fv3jedi_geom
  !From user, maybe via input.nml
  integer :: npx                                          !x-dir grid edge points per tile
  integer :: npy                                          !y-dir grid edge points per tile
  integer :: npz                                          !z-dir grid points global
  logical :: hydrostatic                                  !Are fields on this geometry hydrostatic
  integer :: layout(2)                                    !Processor layout for computation
  integer :: io_layout(2)                                 !Processor layout for read/write
  integer :: halo                                         !Number of halo points, normally 3
  character(len=255) :: nml_file                          !FV3 nml file associated with this geom
  character(len=255) :: trc_file                          !FV3 field_table associated with this geom
  character(len=255) :: wind_type                         !A-grid or D-grid in the state vector
  !Hardwired or determined
  logical :: am_i_root_pe = .false.                       !Is this the root process 
  integer :: size_cubic_grid                              !Size of cubed sphere grid (cell center)
  integer :: ntracers                                     !Number of tracers
  type(domain2D) :: domain                                !MPP domain
  type(fv_grid_bounds_type) :: bd                         !FV grid bounds
  integer :: ntile                                        !Tile ID
  integer :: ntiles = 6                                   !Number of tiles, always 6
  integer :: stackmax                                     !Stackmax
  real(kind=kind_real), allocatable :: grid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: grid_lat(:,:)      !Latitude at cell center
  real(kind=kind_real), allocatable :: egrid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: egrid_lat(:,:)      !Latitude at cell center
  real(kind=kind_real), allocatable :: area(:,:)          !Grid area
  real(kind=kind_real), allocatable :: ak(:),bk(:)        !Model level coefficients
  real(kind=kind_real) :: ptop                            !Pressure at top of domain
  real(kind=kind_real), allocatable :: sin_sg(:,:,:)
  real(kind=kind_real), allocatable :: cos_sg(:,:,:)
  real(kind=kind_real), allocatable :: cosa_u(:,:)
  real(kind=kind_real), allocatable :: cosa_v(:,:)
  real(kind=kind_real), allocatable :: cosa_s(:,:)
  real(kind=kind_real), allocatable :: rsin_u(:,:)
  real(kind=kind_real), allocatable :: rsin_v(:,:)
  real(kind=kind_real), allocatable :: rsin2(:,:)
  real(kind=kind_real), allocatable :: dxa(:,:)
  real(kind=kind_real), allocatable :: dya(:,:)
  real(kind=kind_real), allocatable :: dx(:,:)
  real(kind=kind_real), allocatable :: dy(:,:)
  real(kind=kind_real), allocatable :: dxc(:,:)
  real(kind=kind_real), allocatable :: dyc(:,:)
  real(kind=kind_real), allocatable :: rarea(:,:)
  real(kind=kind_real), allocatable :: rarea_c(:,:)
  real(kind=kind_real), allocatable :: edge_w(:)
  real(kind=kind_real), allocatable :: edge_e(:)
  real(kind=kind_real), allocatable :: edge_s(:)
  real(kind=kind_real), allocatable :: edge_n(:)
  real(kind=kind_real), allocatable :: grid(:,:,:)
  real(kind=kind_real), allocatable :: agrid(:,:,:)

  logical :: sw_corner, se_corner, ne_corner, nw_corner
end type fv3jedi_geom

#define LISTED_TYPE fv3jedi_geom

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_setup(c_key_self, c_conf) bind(c,name='fv3jedi_geo_setup_f90')

implicit none

!Arguments
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)       :: c_conf

!Locals
type(fv3jedi_geom), pointer       :: self
character(len=256)                :: filename_akbk
character(len=256)                :: filepath_akbk
character(len=256)                :: ak_var
character(len=256)                :: bk_var
type(fv_atmos_type), allocatable  :: FV_Atm(:)
logical, allocatable              :: grids_on_this_pe(:)
integer                           :: p_split = 1, fail
integer                           :: ncstat, ncid, dimid, varid, i, readdim, dcount
integer, dimension(nf90_max_var_dims) :: dimIDs, dimLens
integer, allocatable              :: istart(:), icount(:)

! Init, add and get key
! ---------------------
call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

! Is this the root
if (mpp_pe() == mpp_root_pe()) self%am_i_root_pe = .true.

! User input for grid and layout
! ------------------------------

!FV3 Layout files (npx,npy,npz,layout,layout_io,hydrostatic)
self%nml_file = config_get_string(c_conf,len(self%nml_file),"nml_file")
self%trc_file = config_get_string(c_conf,len(self%trc_file),"trc_file")

!Halo
self%halo = config_get_int(c_conf,"halo")

!State wind D or A grid
self%wind_type = config_get_string(c_conf,len(self%nml_file),"wind_type")
if (trim(self%wind_type) /= 'A-grid' .and. trim(self%wind_type) /= 'D-grid') then
   call abor1_ftn("fv3-jedi geometry: wind_type must be either A-grid or D-grid")
endif


! Get number of tracers from field_table
! --------------------------------------
call get_number_tracers(MODEL_ATMOS, num_tracers=self%ntracers)

! Set filenames for ak and bk
! ---------------------------
filepath_akbk = "Data/"
filename_akbk = 'grid_spec.nc'

if (config_element_exists(c_conf,"filepath_akbk")) then
   filepath_akbk = config_get_string(c_conf,len(filepath_akbk),"filepath_akbk")
endif

if (config_element_exists(c_conf,"filename_akbk")) then
   filename_akbk = config_get_string(c_conf,len(filename_akbk),"filename_akbk")
endif

filename_akbk = trim(filepath_akbk)//"/"//trim(filename_akbk)

ak_var = "AK"
bk_var = "BK"
if (config_element_exists(c_conf,"ak_var")) then
   ak_var = config_get_string(c_conf,len(ak_var),"ak_var")
endif
if (config_element_exists(c_conf,"bk_var")) then
   bk_var = config_get_string(c_conf,len(bk_var),"bk_var")
endif

!Intialize using the model setup routine
call fv_init(FV_Atm, 300.0_kind_real, grids_on_this_pe, p_split)
deallocate(pelist_all)

self%bd%isd = FV_Atm(1)%bd%isd
self%bd%ied = FV_Atm(1)%bd%ied
self%bd%jsd = FV_Atm(1)%bd%jsd
self%bd%jed = FV_Atm(1)%bd%jed
self%bd%isc = FV_Atm(1)%bd%isc
self%bd%iec = FV_Atm(1)%bd%iec
self%bd%jsc = FV_Atm(1)%bd%jsc
self%bd%jec = FV_Atm(1)%bd%jec
self%ntile  = FV_Atm(1)%tile

self%npx = FV_Atm(1)%npx
self%npy = FV_Atm(1)%npy
self%npz = FV_Atm(1)%npz
self%layout(1) = FV_Atm(1)%layout(1)
self%layout(2) = FV_Atm(1)%layout(2)
self%io_layout(1) = FV_Atm(1)%io_layout(1)
self%io_layout(2) = FV_Atm(1)%io_layout(2)
self%hydrostatic = FV_Atm(1)%flagstruct%hydrostatic

!Lat,lon and area from
allocate ( self%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
allocate ( self%grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
allocate ( self%grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
allocate ( self%egrid_lon(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
allocate ( self%egrid_lat(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )

self%area = FV_Atm(1)%gridstruct%area_64
self%grid_lon = real(FV_Atm(1)%gridstruct%agrid_64(:,:,1),kind_real)
self%grid_lat = real(FV_Atm(1)%gridstruct%agrid_64(:,:,2),kind_real)
self%egrid_lon = real(FV_Atm(1)%gridstruct%grid_64(:,:,1),kind_real)
self%egrid_lat = real(FV_Atm(1)%gridstruct%grid_64(:,:,2),kind_real)

allocate( self%sin_sg(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9))
allocate( self%cos_sg(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9))
allocate( self%cosa_u(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate( self%cosa_v(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate( self%cosa_s(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate( self%rsin_u(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate( self%rsin_v(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(  self%rsin2(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(    self%dxa(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(    self%dya(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(     self%dx(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(     self%dy(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate(    self%dxc(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed ))
allocate(    self%dyc(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(  self%rarea(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(self%rarea_c(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed+1))
allocate(self%edge_s(self%npx))
allocate(self%edge_n(self%npx))
allocate(self%edge_w(self%npy))
allocate(self%edge_e(self%npy))
allocate(self%grid (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed+1,1:2))
allocate(self%agrid(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,1:2))

self%sin_sg = Fv_Atm(1)%gridstruct%sin_sg
self%cos_sg = Fv_Atm(1)%gridstruct%cos_sg
self%cosa_u = Fv_Atm(1)%gridstruct%cosa_u
self%cosa_v = Fv_Atm(1)%gridstruct%cosa_v
self%cosa_s = Fv_Atm(1)%gridstruct%cosa_s
self%rsin_u = Fv_Atm(1)%gridstruct%rsin_u
self%rsin_v = Fv_Atm(1)%gridstruct%rsin_v
self%rsin2 = Fv_Atm(1)%gridstruct%rsin2
self%dxa = Fv_Atm(1)%gridstruct%dxa
self%dya = Fv_Atm(1)%gridstruct%dya
self%dx = Fv_Atm(1)%gridstruct%dx
self%dy = Fv_Atm(1)%gridstruct%dy
self%dxc = Fv_Atm(1)%gridstruct%dxc
self%dyc = Fv_Atm(1)%gridstruct%dyc
self%rarea = Fv_Atm(1)%gridstruct%rarea
self%rarea_c = Fv_Atm(1)%gridstruct%rarea_c
self%sw_corner = FV_Atm(1)%gridstruct%sw_corner
self%se_corner = FV_Atm(1)%gridstruct%se_corner
self%ne_corner = FV_Atm(1)%gridstruct%ne_corner
self%nw_corner = FV_Atm(1)%gridstruct%nw_corner
self%edge_s = FV_Atm(1)%gridstruct%edge_s
self%edge_n = FV_Atm(1)%gridstruct%edge_n
self%edge_w = FV_Atm(1)%gridstruct%edge_w
self%edge_e = FV_Atm(1)%gridstruct%edge_e
self%grid = FV_Atm(1)%gridstruct%grid
self%agrid = FV_Atm(1)%gridstruct%agrid

!ak and bk are read from file
allocate ( self%ak(self%npz+1) )
allocate ( self%bk(self%npz+1) )

ncstat = nf90_open(filename_akbk, nf90_nowrite, ncid)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

ncstat = nf90_inq_varid(ncid, ak_var, varid)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

dimids = 0
ncstat = nf90_inquire_variable(ncid, varid, dimids = dimids)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

readdim = -1
dcount = 0
do i = 1,nf90_max_var_dims
  if (dimIDs(i) > 0) then
     ncstat = nf90_inquire_dimension(ncid, dimIDs(i), len = dimlens(i))
     if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
     if (dimlens(i) == self%npz+1) then
        readdim = i
     endif
     dcount = dcount + 1
  endif
enddo

if (readdim == -1) call abor1_ftn("fv3-jedi geometry: ak/bk in file does not match dimension of npz from input.nml")

allocate(istart(dcount))
allocate(icount(dcount))

istart = 1
icount = 1
icount(readdim) = self%npz+1

ncstat = nf90_get_var(ncid, varid, self%ak)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

ncstat = nf90_inq_varid(ncid, bk_var, varid)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))
ncstat = nf90_get_var(ncid, varid, self%bk)
if(ncstat /= nf90_noerr) print *, trim(nf90_strerror(ncstat))

self%ptop = self%ak(1)

!Done with the FV_Atm stucture here
call deallocate_fv_atmos_type(FV_Atm(1))
deallocate(FV_Atm)  
deallocate(grids_on_this_pe)

!Misc
self%stackmax = 4000000
self%size_cubic_grid = self%npx-1

!Resetup domain to avoid risk of copied pointers
call setup_domain( self%domain, self%size_cubic_grid, self%size_cubic_grid, &
                   self%ntiles, self%layout, self%io_layout, self%halo)

end subroutine c_fv3jedi_geo_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_clone(c_key_self, c_key_other) bind(c,name='fv3jedi_geo_clone_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(fv3jedi_geom), pointer :: self, other

!add, get, get key
call fv3jedi_geom_registry%add(c_key_other)
call fv3jedi_geom_registry%get(c_key_other, other)
call fv3jedi_geom_registry%get(c_key_self, self)

allocate(other%grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed))
allocate(other%grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed))
allocate(other%egrid_lon(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
allocate(other%egrid_lat(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
allocate(other%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed))
allocate(other%ak(self%npz+1))
allocate(other%bk(self%npz+1))

allocate( other%sin_sg(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9))
allocate( other%cos_sg(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9))
allocate( other%cosa_u(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate( other%cosa_v(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate( other%cosa_s(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate( other%rsin_u(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate( other%rsin_v(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(  other%rsin2(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(    other%dxa(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(    other%dya(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(     other%dx(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(     other%dy(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate(    other%dxc(self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed  ))
allocate(    other%dyc(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1))
allocate(  other%rarea(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(other%rarea_c(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ))
allocate(other%edge_s(self%npx))
allocate(other%edge_n(self%npx))
allocate(other%edge_w(self%npy))
allocate(other%edge_e(self%npy))
allocate(other%grid (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed+1,1:2))
allocate(other%agrid(self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,1:2))

other%npx             = self%npx
other%npy             = self%npy
other%npz             = self%npz
other%hydrostatic     = self%hydrostatic
other%layout          = self%layout
other%io_layout       = self%io_layout
other%halo            = self%halo
other%ntracers        = self%ntracers
other%nml_file        = self%nml_file
other%trc_file        = self%trc_file
other%wind_type       = self%wind_type
other%size_cubic_grid = self%size_cubic_grid
other%ntracers        = self%ntracers
other%bd%isc          = self%bd%isc
other%bd%isd          = self%bd%isd
other%bd%iec          = self%bd%iec
other%bd%ied          = self%bd%ied
other%bd%jsc          = self%bd%jsc
other%bd%jsd          = self%bd%jsd
other%bd%jec          = self%bd%jec
other%bd%jed          = self%bd%jed
other%ntile           = self%ntile
other%ntiles          = self%ntiles
other%stackmax        = self%stackmax
other%grid_lon        = self%grid_lon
other%grid_lat        = self%grid_lat
other%egrid_lon       = self%egrid_lon
other%egrid_lat       = self%egrid_lat
other%area            = self%area
other%ak              = self%ak
other%bk              = self%bk
other%ptop            = self%ptop
other%am_i_root_pe    = self%am_i_root_pe

other%sin_sg = self%sin_sg
other%cos_sg = self%cos_sg
other%cosa_u = self%cosa_u
other%cosa_v = self%cosa_v
other%cosa_s = self%cosa_s
other%rsin_u = self%rsin_u
other%rsin_v = self%rsin_v
other%rsin2 = self%rsin2
other%dxa = self%dxa
other%dya = self%dya
other%dx = self%dx
other%dy = self%dy
other%dxc = self%dxc
other%dyc = self%dyc
other%rarea = self%rarea
other%rarea_c = self%rarea_c
other%sw_corner = self%sw_corner
other%se_corner = self%se_corner
other%ne_corner = self%ne_corner
other%nw_corner = self%nw_corner
other%edge_s = self%edge_s
other%edge_n = self%edge_n
other%edge_w = self%edge_w
other%edge_e = self%edge_e
other%grid   = self%grid
other%agrid  = self%agrid

call setup_domain( other%domain, other%size_cubic_grid, other%size_cubic_grid, &
                   other%ntiles, other%layout, other%io_layout, other%halo)

end subroutine c_fv3jedi_geo_clone

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_delete(c_key_self) bind(c,name='fv3jedi_geo_delete_f90')

implicit none

integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_geom), pointer :: self

! Get key
call fv3jedi_geom_registry%get(c_key_self, self)

! Deallocate
deallocate(self%grid_lon)
deallocate(self%grid_lat)
deallocate(self%egrid_lon)
deallocate(self%egrid_lat)
deallocate(self%area)
deallocate(self%ak)
deallocate(self%bk)

deallocate(self%sin_sg)
deallocate(self%cos_sg)
deallocate(self%cosa_u)
deallocate(self%cosa_v)
deallocate(self%cosa_s)
deallocate(self%rsin_u)
deallocate(self%rsin_v)
deallocate( self%rsin2)
deallocate(   self%dxa)
deallocate(   self%dya)
deallocate(    self%dx)
deallocate(    self%dy)
deallocate(    self%dxc)
deallocate(    self%dyc)
deallocate( self%rarea)
deallocate( self%rarea_c)
deallocate(self%edge_s)
deallocate(self%edge_n)
deallocate(self%edge_w)
deallocate(self%edge_e)
deallocate(self%grid)
deallocate(self%agrid)

call mpp_deallocate_domain(self%domain)

! Remove key
call fv3jedi_geom_registry%remove(c_key_self)

end subroutine c_fv3jedi_geo_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_info(c_key_self) bind(c,name='fv3jedi_geo_info_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
type(fv3jedi_geom), pointer :: self

call fv3jedi_geom_registry%get(c_key_self, self)

end subroutine c_fv3jedi_geo_info

! ------------------------------------------------------------------------------

subroutine setup_domain(domain, nx, ny, ntiles, layout_in, io_layout, halo)

 implicit none

 type(domain2D),   intent(inout) :: domain
 integer,          intent(in)    :: nx, ny, ntiles
 integer,          intent(in)    :: layout_in(:), io_layout(:)
 integer,          intent(in)    :: halo

 integer                              :: pe, npes, npes_per_tile, tile
 integer                              :: num_contact
 integer                              :: n, layout(2)
 integer, allocatable, dimension(:,:) :: global_indices, layout2D
 integer, allocatable, dimension(:)   :: pe_start, pe_end
 integer, allocatable, dimension(:)   :: tile1, tile2
 integer, allocatable, dimension(:)   :: istart1, iend1, jstart1, jend1
 integer, allocatable, dimension(:)   :: istart2, iend2, jstart2, jend2
 integer, allocatable :: tile_id(:)
 logical :: is_symmetry

  pe = mpp_pe()
  npes = mpp_npes()

  if (mod(npes,ntiles) /= 0) then
     call mpp_error(NOTE, "setup_domain: npes can not be divided by ntiles")
     return
  endif
  npes_per_tile = npes/ntiles
  tile = pe/npes_per_tile + 1

  if (layout_in(1)*layout_in(2) == npes_per_tile) then
     layout = layout_in
  else
     call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
  endif

  if (io_layout(1) <1 .or. io_layout(2) <1) call mpp_error(FATAL, &
          "setup_domain: both elements of variable io_layout must be positive integer")
  if (mod(layout(1), io_layout(1)) /= 0 ) call mpp_error(FATAL, &
       "setup_domain: layout(1) must be divided by io_layout(1)")
  if (mod(layout(2), io_layout(2)) /= 0 ) call mpp_error(FATAL, &
       "setup_domain: layout(2) must be divided by io_layout(2)")

  allocate(global_indices(4,ntiles), layout2D(2,ntiles), pe_start(ntiles), pe_end(ntiles) )
  do n = 1, ntiles
     global_indices(:,n) = (/1,nx,1,ny/)
     layout2D(:,n)       = layout
     pe_start(n)         = (n-1)*npes_per_tile
     pe_end(n)           = n*npes_per_tile-1
  enddo

  ! this code copied from domain_decomp in fv_mp_mod.f90
  num_contact = 12
  allocate(tile1(num_contact), tile2(num_contact) )
  allocate(tile_id(ntiles))
  allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
  allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )
  !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
  tile1(1) = 1; tile2(1) = 2
  istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
  istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
  !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
  tile1(2) = 1; tile2(2) = 3
  istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
  istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
  !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
  tile1(3) = 1; tile2(3) = 5
  istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
  istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
  !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
  tile1(4) = 1; tile2(4) = 6
  istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
  istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
  !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
  tile1(5) = 2; tile2(5) = 3
  istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
  istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
  !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
  tile1(6) = 2; tile2(6) = 4
  istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
  istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
  !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
  tile1(7) = 2; tile2(7) = 6
  istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
  istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
  !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
  tile1(8) = 3; tile2(8) = 4
  istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
  istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
  !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
  tile1(9) = 3; tile2(9) = 5
  istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
  istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
  !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
  tile1(10) = 4; tile2(10) = 5
  istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
  istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
  !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
  tile1(11) = 4; tile2(11) = 6
  istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
  istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
  !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
  tile1(12) = 5; tile2(12) = 6
  istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
  istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
  is_symmetry = .true.
  do n = 1, ntiles
     tile_id(n) = n
  enddo

  call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                         pe_start, pe_end, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo,    &
                         symmetry=is_symmetry, tile_id=tile_id, &
                         name='cubic_grid')

  if (io_layout(1) /= 1 .or. io_layout(2) /= 1) call mpp_define_io_domain(domain, io_layout)

  deallocate(pe_start, pe_end)
  deallocate(layout2D, global_indices)
  deallocate(tile1, tile2, tile_id)
  deallocate(istart1, iend1, jstart1, jend1)
  deallocate(istart2, iend2, jstart2, jend2)

end subroutine setup_domain

end module fv3jedi_geom_mod
