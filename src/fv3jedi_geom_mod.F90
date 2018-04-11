
!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

!General JEDI uses
use kinds
use iso_c_binding
use config_mod

!Uses for setting up tracer numbers from field_table
use tracer_manager_mod, only: get_number_tracers
use field_manager_mod,  only: MODEL_ATMOS

!Derived types within geometry
use fv3jedi_mod,         only: fv_grid_bounds_type, setup_domain

!FMS/MPP uses
use mpp_domains_mod,    only: domain2D, mpp_get_tile_id, mpp_deallocate_domain
use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_mod,            only: mpp_pe
use fms_mod,            only: get_mosaic_tile_grid
use fms_io_mod,         only: restart_file_type, register_restart_field, &
                              free_restart_type, restore_state, read_data

!Uses for generating geometry using FV3 routines
use fv_arrays_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
use fv_control_mod,     only: fv_init, pelist_all
use fv3jedi_constants,   only: rad2deg

implicit none
private
public :: fv3jedi_geom
public :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the FV3JEDI model

type :: fv3jedi_geom
  !From user
  integer :: npx                                          !x-dir grid edge points per tile
  integer :: npy                                          !y-dir grid edge points per tile
  integer :: npz                                          !z-dir grid points global
  logical :: hydrostatic                                  !Are fields on this geometry hydrostatic
  integer :: layout(2)                                    !Processor layout for computation
  integer :: io_layout(2)                                 !Processor layout for read/write
  integer :: halo                                         !Number of halo points, normally 3
  character(len=255) :: nml_file                          !Datapath for generating grid from file
  character(len=255) :: trc_file                          !Datapath for generating grid from file
  !Hardwired or determined
  integer :: size_cubic_grid                              !Size of cubed sphere grid (cell center)
  integer :: ntracers                                     !Number of tracers
  type(domain2D) :: domain                                !MPP domain
  type(fv_grid_bounds_type) :: bd                         !FV grid bounds
  integer :: ntile                                        !Tile ID
  integer :: ntiles = 6                                   !Number of tiles, always 6
  integer :: stackmax                                     !Stackmax
  real(kind=kind_real), allocatable :: grid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: grid_lat(:,:)      !Latitude at cell center
  real(kind=kind_real), allocatable :: area(:,:)          !Grid area
  real(kind=kind_real), allocatable :: ak(:),bk(:)        !Model level coefficients
  real(kind=kind_real) :: ptop                            !Pressure at top of domain
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
type(c_ptr), intent(in)    :: c_conf

!Things needed in general
type(fv3jedi_geom), pointer :: self
character(len=20) :: init_type
integer :: hydro_int

!Things needed for generting from file read
character(len=256)                :: datapath_in
character(len=256)                :: atm_hgrid, atm_hgridf
integer                           :: start(4), nread(4), id_restart
type(restart_file_type)           :: Fv_restart
integer                           :: ntile(1), i, j
real(kind=kind_real), allocatable :: tmpx(:,:), tmpy(:,:)
integer                           :: isc2, iec2, jsc2, jec2
character(len=256)                :: filename_spec 
character(len=256)                :: filename_core
character(len=256)                :: fullpath_spec
character(len=256)                :: fullpath_core

!Things needed for generating inline
type(fv_atmos_type), allocatable  :: FV_Atm(:)
logical, allocatable              :: grids_on_this_pe(:)
integer                           :: p_split = 1, fail

! Init, add and get key
! ---------------------
call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

! Setup from file or inline using fv_init
! ---------------------------------------
init_type = 'inline'
if (config_element_exists(c_conf,"init_type")) then
  init_type = config_get_string(c_conf,len(init_type),"init_type")
endif

! User input about grid and layout
! --------------------------------
self%npx = config_get_int(c_conf,"npx")
self%npy = config_get_int(c_conf,"npy")
self%npz = config_get_int(c_conf,"npz")
hydro_int = config_get_int(c_conf,"hydrostatic")
if (hydro_int == 1) then
  self%hydrostatic = .true.
else
  self%hydrostatic = .false.
endif
self%layout(1) = config_get_int(c_conf,"layout_1")
self%layout(2) = config_get_int(c_conf,"layout_2")
self%io_layout(1) = config_get_int(c_conf,"io_layout_1")
self%io_layout(2) = config_get_int(c_conf,"io_layout_2")
self%halo = config_get_int(c_conf,"halo")

self%nml_file = config_get_string(c_conf,len(self%nml_file),"nml_file")
self%trc_file = config_get_string(c_conf,len(self%trc_file),"trc_file")

if (self%npx /= self%npy) then
   call abor1_ftn("fv3-jedi geometry: cube faces not square (npx /= npy)")
endif

self%stackmax = 4000000
self%size_cubic_grid = self%npx-1

! Get number of tracers from field_table
! --------------------------------------
call get_number_tracers(MODEL_ATMOS, num_tracers=self%ntracers)

! Set filenames
! -------------
filename_spec = 'grid_spec.nc'
filename_core = 'fv_core.res.nc'

if (config_element_exists(c_conf,"filename_spec")) then
   filename_spec = config_get_string(c_conf,len(filename_spec),"filename_spec")
endif
if (config_element_exists(c_conf,"filename_core")) then
   filename_core = config_get_string(c_conf,len(filename_core),"filename_core")
endif

!Main data path for restarts
datapath_in = config_get_string(c_conf,len(datapath_in), "datapath_geom")

!Full paths to restarts
fullpath_spec = trim(adjustl(datapath_in))//trim(adjustl(filename_spec))
fullpath_core = trim(adjustl(datapath_in))//trim(adjustl(filename_core))


if (trim(init_type) .ne. "inline") then

   if (mpp_pe() == 0) print*, 'Grid generation method: read from file'
  
   ! Create geometry based on grid type and layout
   call setup_domain(self%domain, self%size_cubic_grid, self%size_cubic_grid, self%ntiles, self%layout, self%io_layout, self%halo)
   
   ! Get compute and data domain information
   call mpp_get_compute_domain(self%domain, self%bd%isc, self%bd%iec, self%bd%jsc, self%bd%jec)
   call mpp_get_data_domain(self%domain, self%bd%isd, self%bd%ied, self%bd%jsd, self%bd%jed)
   ntile = mpp_get_tile_id(self%domain)
         
   ! get grid_lon/grid_lat from superset grid in grid tile files.
   ! easy to add lon/lat values for d-grid u and v from this.
   call get_mosaic_tile_grid(atm_hgrid, fullpath_spec, self%domain)

   atm_hgridf = trim(adjustl(datapath_in))//trim(adjustl(atm_hgrid))

   isc2 = 2*self%bd%isc-1; iec2 = 2*self%bd%iec+1
   jsc2 = 2*self%bd%jsc-1; jec2 = 2*self%bd%jec+1
   allocate(tmpx(isc2:iec2, jsc2:jec2) )
   allocate(tmpy(isc2:iec2, jsc2:jec2) )
   allocate ( self%grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   allocate ( self%grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   start = 1; nread = 1
   start(1) = isc2; nread(1) = iec2 - isc2 + 1
   start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
   call read_data(atm_hgridf, 'x', tmpx, start, nread, no_domain=.true.)
   call read_data(atm_hgridf, 'y', tmpy, start, nread, no_domain=.true.)
   self%grid_lon = 0; self%grid_lat = 0
   do j = self%bd%jsc, self%bd%jec
      do i = self%bd%isc, self%bd%iec
         self%grid_lon(i,j) = tmpx(2*i,2*j)
         self%grid_lat(i,j) = tmpy(2*i,2*j)
      enddo
   enddo
   deallocate(tmpx,tmpy)
   
   ! get grid areas from superset grid in grid tile files.
   allocate ( self%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   isc2 = 2*self%bd%isc-1; iec2 = 2*self%bd%iec
   jsc2 = 2*self%bd%jsc-1; jec2 = 2*self%bd%jec
   start = 1; nread = 1
   start(1) = isc2; nread(1) = iec2 - isc2 + 1
   start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
   allocate(tmpx(isc2:iec2, jsc2:jec2) )
   call read_data(atm_hgridf, 'area', tmpx, start, nread, no_domain=.true.)
   self%area = 0.
   do j = self%bd%jsc, self%bd%jec
      do i = self%bd%isc, self%bd%iec
         self%area(i,j) = tmpx(2*i,2*j)  +tmpx(2*i-1,2*j-1)+&
                          tmpx(2*i-1,2*j)+tmpx(2*i,2*j-1)
      enddo
   enddo
   deallocate(tmpx)
   
   ! get ak, bk from fv_core.res.nc
   allocate ( self%ak(self%npz+1) )
   allocate ( self%bk(self%npz+1) )

   ! register and read ak,bk
   id_restart = register_restart_field(Fv_restart, fullpath_core, 'ak', self%ak, &
                 no_domain=.true.)
   id_restart = register_restart_field(Fv_restart, fullpath_core, 'bk', self%bk, &
                 no_domain=.true.)
   call restore_state(Fv_restart, directory='')
   call free_restart_type(Fv_restart)
   
   self%ptop = self%ak(1)

else
   
   if (mpp_pe() == 0) print*, 'Grid generation method: inline'

   !Intialize using the model setup routine
   call fv_init(FV_Atm, 300.0_kind_real, grids_on_this_pe, p_split)
   deallocate(pelist_all)

   !Perform some checks
   fail = 0
   if (self%npx /= FV_Atm(1)%npx) fail = fail + 1
   if (self%npy /= FV_Atm(1)%npy) fail = fail + 1
   if (self%npz /= FV_Atm(1)%npz) fail = fail + 1
   if (fail /= 0) then
      print*, 'Consistency check failure ', fail
      call abor1_ftn("fv3-jedi geometry: input.nml inconsistent with regular input")
   endif

   self%domain = FV_Atm(1)%domain
   self%bd%isd = FV_Atm(1)%bd%isd
   self%bd%ied = FV_Atm(1)%bd%ied
   self%bd%jsd = FV_Atm(1)%bd%jsd
   self%bd%jed = FV_Atm(1)%bd%jed
   self%bd%isc = FV_Atm(1)%bd%isc
   self%bd%iec = FV_Atm(1)%bd%iec
   self%bd%jsc = FV_Atm(1)%bd%jsc
   self%bd%jec = FV_Atm(1)%bd%jec
   self%ntile  = FV_Atm(1)%tile
   
   !Lat,lon and area from
   allocate ( self%grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   allocate ( self%grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   allocate ( self%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   self%grid_lon = rad2deg*real(FV_Atm(1)%gridstruct%agrid_64(:,:,1),kind_real)
   self%grid_lat = rad2deg*real(FV_Atm(1)%gridstruct%agrid_64(:,:,2),kind_real)
   self%area = FV_Atm(1)%gridstruct%area_64
   
   !ak and bk are still read from file
   allocate ( self%ak(self%npz+1) )
   allocate ( self%bk(self%npz+1) )
   id_restart = register_restart_field(Fv_restart, fullpath_core, 'ak', self%ak, no_domain=.true.)
   id_restart = register_restart_field(Fv_restart, fullpath_core, 'bk', self%bk, no_domain=.true.)
   call restore_state(Fv_restart, directory='')
   call free_restart_type(Fv_restart)

   self%ptop = self%ak(1)

   !Done with the FV_Atm stucture here
   call deallocate_fv_atmos_type(FV_Atm(1))
   deallocate(FV_Atm)  
   deallocate(grids_on_this_pe)

endif

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
allocate(other%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed))
allocate(other%ak(self%npz+1))
allocate(other%bk(self%npz+1))

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
other%area            = self%area
other%ak              = self%ak
other%bk              = self%bk
other%ptop            = self%ptop

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
deallocate(self%area)
deallocate(self%ak)
deallocate(self%bk)
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

end module fv3jedi_geom_mod
