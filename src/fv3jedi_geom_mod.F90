
!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

!General JEDI uses
use kinds
use iso_c_binding
use config_mod

!Uses for setting up tracer numbers from field_table
use tracer_manager_mod, only: get_number_tracers
use field_manager_mod,  only: MODEL_ATMOS

!FMS/MPP uses
use mpp_domains_mod,    only: domain2D, mpp_get_tile_id, mpp_deallocate_domain
use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod,    only: mpp_define_layout, mpp_define_mosaic, mpp_define_io_domain
use mpp_mod,            only: mpp_pe, mpp_npes, mpp_error, FATAL, NOTE
use fms_mod,            only: get_mosaic_tile_grid
use fms_io_mod,         only: restart_file_type, register_restart_field, &
                              free_restart_type, restore_state, read_data

!Uses for generating geometry using FV3 routines
use fv_arrays_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
use fv_control_mod,     only: fv_init, pelist_all
use fv3jedi_constants,  only: rad2deg

implicit none
private
public :: fv3jedi_geom
public :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

!> Skinny version of fv_grid_type
!type :: fv_grid_type
!  real(kind=kind_real), allocatable, dimension(:,:,:) :: sin_sg
!  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_u
!  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_v
!  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_s
!  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin_u
!  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin_v
!  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin2
!  real(kind=kind_real), allocatable, dimension(:,:)   :: dxa, dya
!  logical :: sw_corner, se_corner, ne_corner, nw_corner
!endtype

!> Skinny version of fv_grid_bounds_type
type fv_grid_bounds_type
    integer :: isd, ied, jsd, jed ! data domain
    integer :: isc, iec, jsc, jec ! compute domain
end type fv_grid_bounds_type

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
  character(len=255) :: wind_type                         !A-grid or D-grid in the control vector
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
  real(kind=kind_real), allocatable :: egrid_lon(:,:)      !Longitude at cell center
  real(kind=kind_real), allocatable :: egrid_lat(:,:)      !Latitude at cell center
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

self%wind_type = config_get_string(c_conf,len(self%nml_file),"wind_type")

if (trim(self%wind_type) /= 'A-grid' .and. trim(self%wind_type) /= 'D-grid') then
   call abor1_ftn("fv3-jedi geometry: wind_type must be either A-grid or D-grid")
endif

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
   
   allocate ( self%egrid_lon(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
   allocate ( self%egrid_lat(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
   self%egrid_lon = 0.0 !Not in the file
   self%egrid_lat = 0.0

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
   
   allocate ( self%egrid_lon(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
   allocate ( self%egrid_lat(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
   self%egrid_lon = rad2deg*real(FV_Atm(1)%gridstruct%grid_64(:,:,1),kind_real)
   self%egrid_lat = rad2deg*real(FV_Atm(1)%gridstruct%grid_64(:,:,2),kind_real)

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

   !Resetup domain to avoid risk of copied pointers
   call setup_domain( self%domain, self%size_cubic_grid, self%size_cubic_grid, &
                      self%ntiles, self%layout, self%io_layout, self%halo)


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
allocate(other%egrid_lon(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
allocate(other%egrid_lat(self%bd%isd:self%bd%ied+1, self%bd%jsd:self%bd%jed+1) )
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
other%egrid_lon        = self%egrid_lon
other%egrid_lat        = self%egrid_lat
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
deallocate(self%egrid_lon)
deallocate(self%egrid_lat)
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

! Setup domain for generating grid without using fv_init
! ------------------------------------------------------
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
