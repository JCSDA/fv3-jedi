
!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

use kinds
use iso_c_binding
use config_mod

use mpp_domains_mod,    only: domain2D, mpp_get_layout,mpp_get_tile_id
use mpp_domains_mod,    only: mpp_copy_domain, mpp_deallocate_domain
use mpp_domains_mod,    only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_io_mod,         only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
use mpp_mod,            only: mpp_pe, mpp_npes, mpp_init, mpp_exit, mpp_error, &
                              FATAL
use fms_io_mod,         only: file_exist, read_data
use fms_mod,            only: get_mosaic_tile_grid, fms_init, fms_end
use fms_io_mod,         only: restart_file_type, register_restart_field, &
                              free_restart_type, restore_state, save_restart, file_exist
use tracer_manager_mod, only: get_number_tracers
use field_manager_mod,  only: MODEL_ATMOS

use fv3jedi_mod,         only: fv_grid_bounds_type, fv_grid_type
use fv3jedi_mod,         only: setup_domain

use fv3jedi_constants,   only: rad2deg

use fv_arrays_mod,      only: fv_atmos_type, deallocate_fv_atmos_type
use fv_control_mod,     only: fv_init, pelist_all

implicit none
private
public :: fv3jedi_geom
public :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the FV3JEDI model
type :: fv3jedi_geom
  type(domain2D) :: domain
  type(fv_grid_bounds_type) :: bd
  integer :: ntile
  integer :: size_cubic_grid, nlevs, ntracers, halo, ntiles
  integer :: stackmax
  logical :: hydrostatic, agrid_vel_rst
  integer :: layout(2)
  integer :: io_layout(2)
  real(kind=kind_real), allocatable, dimension(:,:) :: grid_lon
  real(kind=kind_real), allocatable, dimension(:,:) :: grid_lat
  real(kind=kind_real), allocatable, dimension(:,:) :: area
  real(kind=kind_real), allocatable, dimension(:)   :: ak,bk
  character(len=255) :: datapath_in
!Needed for d to a grid interpolation
  type(fv_grid_type) :: gridstruct
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

integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(fv3jedi_geom), pointer :: self

character(len=256)                 :: atm_hgrid, atm_hgridf
integer                            :: start(4), nread(4), id_restart
type(restart_file_type) :: Fv_restart

integer :: ntile(1), pe, i, j
character(len=1) :: tileid
character(len=4) :: resid
real(kind=kind_real), allocatable, dimension(:,:) :: tmpx, tmpy
!real(kind=kind_real), allocatable, dimension(:,:) :: grid_lon, grid_lat
type(domain2D), pointer :: domain
character(len=32) :: gtype = 'cubic_grid'
integer :: isc2, iec2, jsc2, jec2
integer :: io_status, nmlunit
integer :: size_cubic_grid = 48
integer :: nlevs = 10
integer :: ntracers = 3
logical :: hydrostatic = .true.  ! nonhydrostatic fields in restart?
logical :: agrid_vel_rst = .false. ! agrid winds in restart?
integer :: halo = 0
integer :: ntiles = 6
integer :: stackmax = 4000000
integer :: layout(2) = (/1,1/) ! set ndivs_x and ndivs_y to divide each tile into io_layout(1)*io_layout(2)
                                    ! group and write out data from the root pe of each group.
integer :: io_layout(2) = (/1,1/) 
character(len=256) :: grid_spec_file = 'INPUT/grid_spec.nc'
character(len=256) :: fv_core_res_file = 'INPUT/fv_core.res.nc'
character(len=256) :: grid_spec_path
character(len=256) :: fv_core_res_path
 
type(fv_atmos_type), allocatable :: FV_Atm(:)
logical, allocatable             :: grids_on_this_pe(:)
integer                          :: p_split = 1

character(len=20) :: init_type

integer :: print_info = 0

namelist /fv3jedi_geom_nml/ ntiles, ntracers, nlevs, size_cubic_grid, &
                           hydrostatic, halo, stackmax, layout, io_layout, &
                           agrid_vel_rst,grid_spec_file,fv_core_res_file

pe = mpp_pe()

!Seems reasonable to do this here, since fms and mpp may be used thoughout
call fms_init()
call mpp_init()

call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

!Two methods exist for initializing the FV3 grid, either it can
!be read from file or initialzed using the model setup, fv_itit.
!fv_init does a lot more than generate what is required by the 
!jedi geom so may be slower. On the other hand no files are 
!required in advance of running. Files containing ak and bk
!are always required from fv3jedi_geom_nml/fv_core_res_file

init_type = config_get_string(c_conf,len(init_type),"init_type")

!Main data path for restarts
self%datapath_in = config_get_string(c_conf,len(self%datapath_in), "datapath_in")

!Full paths to restarts
grid_spec_path = trim(adjustl(self%datapath_in))//trim(adjustl(grid_spec_file))
fv_core_res_path = trim(adjustl(self%datapath_in))//trim(adjustl(fv_core_res_file))

!Read namelist.
if (file_exist('input.nml') )then
   nmlunit = 101
   call mpp_open(nmlunit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY)
   read(nmlunit,fv3jedi_geom_nml,iostat=io_status)
   call mpp_close(nmlunit)
endif
if (io_status > 0) then
   call mpp_error(FATAL, '=> fv3jedi_geom: Error reading fv3jedi_geom_nml')
endif

if (trim(init_type) .ne. "inline" .and. .not.agrid_vel_rst) then

  !In order to interpolate from D-grid to A-grid will require things that 
  !are obtained by generating the grid inline, if reading from file these
  !quantities are not present
  call mpp_error(FATAL, 'If A-grid winds not present in file, grid must be generated inline')

endif

!Get number of tracers from field_table
call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

if (trim(init_type) .ne. "inline") then

   if (pe == 0) print*, 'Grid generation method: read from file'
  
   ! Create geometry based on grid type and layout
   call setup_domain(self%domain, size_cubic_grid, size_cubic_grid, ntiles, layout, io_layout, halo)
   
   ! Get compute and data domain information
   call mpp_get_compute_domain(self%domain, self%bd%isc, self%bd%iec, self%bd%jsc, self%bd%jec)
   call mpp_get_data_domain(self%domain, self%bd%isd, self%bd%ied, self%bd%jsd, self%bd%jed)
   ntile = mpp_get_tile_id(self%domain)
   
   self%ntile = ntile(1)
   self%size_cubic_grid = size_cubic_grid
   self%nlevs = nlevs
   self%ntracers = ntracers
   self%hydrostatic = hydrostatic
   self%agrid_vel_rst = agrid_vel_rst
   self%halo = halo
   self%ntiles = ntiles
   self%layout = layout
   self%io_layout = io_layout
   self%stackmax = stackmax
   
   ! get grid_lon/grid_lat directly from grid_spec files.
   if (print_info == 1) print *,'pe,ntile,isc,iec,jsc,jec',pe,self%ntile,self%bd%isc,self%bd%iec,self%bd%jsc,self%bd%jec
   if (print_info == 1) print *,'pe,ntile,isd,ied,jsd,jed',pe,self%ntile,self%bd%isd,self%bd%ied,self%bd%jsd,self%bd%jed
   !write(tileid,'(i1)') ntile
   !write(resid,'(i4)') size_cubic_grid
   !atm_hgrid = 'INPUT/C'//trim(adjustl(resid))//'_grid_spec.tile'//tileid//'.nc'
   !print *,'1: pe,atm_hgrid ',pe,trim(atm_hgrid)
   !allocate ( grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   !allocate ( grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   !start = 1; nread = 1
   !start(1) = self%bd%isd; nread(1) = self%bd%ied - self%bd%isd + 1
   !start(2) = self%bd%jsd; nread(2) = self%bd%jed - self%bd%jsd + 1
   !call read_data(atm_hgrid, 'grid_lont', grid_lon, start, nread, domain=self%domain)
   !call read_data(atm_hgrid, 'grid_latt', grid_lat, start, nread, domain=self%domain)
   !print *,'1: pe,grid_lon ',pe,minval(grid_lon),maxval(grid_lon)
   !print *,'1: pe,grid_lat ',pe,minval(grid_lat),maxval(grid_lat)
   
   ! get grid_lon/grid_lat from superset grid in grid tile files.
   ! easy to add lon/lat values for d-grid u and v from this.
   call get_mosaic_tile_grid(atm_hgrid, grid_spec_path, self%domain)

   atm_hgridf = trim(adjustl(self%datapath_in))//trim(adjustl(atm_hgrid))

   !print *,'2: pe,atm_hgrid ',pe,trim(atm_hgrid)
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
   !print *,'2: pe,grid_lon ',pe,minval(self%grid_lon),maxval(self%grid_lon)
   !print *,'2: pe,grid_lat ',pe,minval(self%grid_lat),maxval(self%grid_lat)
   !print *,'2: pe,grid_lon diff ',pe,minval(self%grid_lon-grid_lon),maxval(self%grid_lon-grid_lon)
   !print *,'2: pe,grid_lat diff ',pe,minval(self%grid_lat-grid_lat),maxval(self%grid_lat-grid_lat)
   !deallocate(grid_lon, grid_lat)
   
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
   !print *,'pe,area ',pe,sqrt(minval(self%area))/1000.,sqrt(maxval(self%area))/1000.
   deallocate(tmpx)
   
   ! get ak, bk from fv_core.res.nc
   allocate ( self%ak(self%nlevs+1) )
   allocate ( self%bk(self%nlevs+1) )
   ! register and read ak,bk
   id_restart = register_restart_field(Fv_restart, fv_core_res_path, 'ak', self%ak, &
                 no_domain=.true.)
   id_restart = register_restart_field(Fv_restart, fv_core_res_path, 'bk', self%bk, &
                 no_domain=.true.)
   call restore_state(Fv_restart, directory='')
   call free_restart_type(Fv_restart)
   
else
   
   if (pe == 0) print*, 'Grid generation method: inline'

   !Intialize using the model setup routine
   call fv_init(FV_Atm, 300.0_kind_real, grids_on_this_pe, p_split)
                        !^some dummy value

   self%domain = FV_Atm(1)%domain
   self%bd%isd = FV_Atm(1)%bd%isd
   self%bd%ied = FV_Atm(1)%bd%ied
   self%bd%jsd = FV_Atm(1)%bd%jsd
   self%bd%jed = FV_Atm(1)%bd%jed
   self%bd%isc = FV_Atm(1)%bd%isc
   self%bd%iec = FV_Atm(1)%bd%iec
   self%bd%jsc = FV_Atm(1)%bd%jsc
   self%bd%jec = FV_Atm(1)%bd%jec
   self%ntile = FV_Atm(1)%tile
   self%size_cubic_grid = size_cubic_grid
   self%nlevs = nlevs
   self%ntracers = ntracers
   self%halo = halo
   self%ntiles = ntiles
   self%stackmax = stackmax
   self%hydrostatic = hydrostatic
   self%agrid_vel_rst = agrid_vel_rst
   self%layout = layout
   self%io_layout = io_layout
   
   !Lat,lon and area from
   allocate ( self%grid_lon(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   allocate ( self%grid_lat(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   allocate ( self%area(self%bd%isd:self%bd%ied, self%bd%jsd:self%bd%jed) )
   self%grid_lon = rad2deg*real(FV_Atm(1)%gridstruct%agrid_64(:,:,1),kind_real)
   self%grid_lat = rad2deg*real(FV_Atm(1)%gridstruct%agrid_64(:,:,2),kind_real)
   self%area = FV_Atm(1)%gridstruct%area_64
   
   !ak and bk are still read from file
   allocate ( self%ak(self%nlevs+1) )
   allocate ( self%bk(self%nlevs+1) )
   id_restart = register_restart_field(Fv_restart, fv_core_res_path, 'ak', self%ak, no_domain=.true.)
   id_restart = register_restart_field(Fv_restart, fv_core_res_path, 'bk', self%bk, no_domain=.true.)
   call restore_state(Fv_restart, directory='')
   call free_restart_type(Fv_restart)
   
   !Gridstruct needed to do interpolation from D to A grid
   if (.not. self%agrid_vel_rst) then
      allocate ( self%gridstruct%sin_sg (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9) )
      allocate ( self%gridstruct%cosa_u (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed    ) )
      allocate ( self%gridstruct%cosa_v (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1  ) )
      allocate ( self%gridstruct%cosa_s (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
      allocate ( self%gridstruct%rsin_u (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed    ) )
      allocate ( self%gridstruct%rsin_v (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1  ) )
      allocate ( self%gridstruct% rsin2 (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
      allocate ( self%gridstruct%   dxa (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
      allocate ( self%gridstruct%   dya (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )

      self%gridstruct%sin_sg = FV_Atm(1)%gridstruct%sin_sg
      self%gridstruct%cosa_u = FV_Atm(1)%gridstruct%cosa_u
      self%gridstruct%cosa_v = FV_Atm(1)%gridstruct%cosa_v
      self%gridstruct%cosa_s = FV_Atm(1)%gridstruct%cosa_s
      self%gridstruct%rsin_u = FV_Atm(1)%gridstruct%rsin_u
      self%gridstruct%rsin_v = FV_Atm(1)%gridstruct%rsin_v
      self%gridstruct% rsin2 = FV_Atm(1)%gridstruct% rsin2
      self%gridstruct%   dxa = FV_Atm(1)%gridstruct%   dxa
      self%gridstruct%   dya = FV_Atm(1)%gridstruct%   dya

      self%gridstruct%sw_corner = FV_Atm(1)%gridstruct%sw_corner
      self%gridstruct%se_corner = FV_Atm(1)%gridstruct%se_corner
      self%gridstruct%ne_corner = FV_Atm(1)%gridstruct%ne_corner
      self%gridstruct%nw_corner = FV_Atm(1)%gridstruct%nw_corner
   endif

   !Done with the FV_Atm stucture here
   call deallocate_fv_atmos_type(FV_Atm(1))
   deallocate(FV_Atm)
   deallocate(pelist_all) !From fv_control
   deallocate(grids_on_this_pe)

endif

end subroutine c_fv3jedi_geo_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_clone(c_key_self, c_key_other) bind(c,name='fv3jedi_geo_clone_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_key_other
integer layout(2)
integer io_layout(2)

type(fv3jedi_geom), pointer :: self, other

call fv3jedi_geom_registry%add(c_key_other)
call fv3jedi_geom_registry%get(c_key_other, other)
call fv3jedi_geom_registry%get(c_key_self, self)

! copy self to other
other%bd%isc = self%bd%isc
other%bd%isd = self%bd%isd
other%bd%iec = self%bd%iec
other%bd%ied = self%bd%ied
other%bd%jsc = self%bd%jsc
other%bd%jsd = self%bd%jsd
other%bd%jec = self%bd%jec
other%bd%jed = self%bd%jed
other%ntile = self%ntile
other%size_cubic_grid = self%size_cubic_grid
other%nlevs = self%nlevs
other%halo = self%halo
other%ntiles = self%ntiles
other%ntracers = self%ntracers
other%hydrostatic = self%hydrostatic
other%agrid_vel_rst = self%agrid_vel_rst
other%layout = self%layout
other%io_layout = self%io_layout
other%stackmax = self%stackmax
!call mpp_copy_domain(self%domain, other%domain) ! doesn't work
call setup_domain(other%domain, other%size_cubic_grid, other%size_cubic_grid, &
                other%ntiles, other%layout, other%io_layout, other%halo)
!call mpp_get_layout(self%domain, layout)
!print *,'clone get domain layout 1',layout
!call mpp_get_layout(other%domain, layout)
!print *,'clone get domain layout 2',layout
allocate ( other%grid_lon(other%bd%isd:other%bd%ied, other%bd%jsd:other%bd%jed) )
allocate ( other%grid_lat(other%bd%isd:other%bd%ied, other%bd%jsd:other%bd%jed) )
allocate ( other%area(other%bd%isd:other%bd%ied, other%bd%jsd:other%bd%jed) )
other%grid_lon = self%grid_lon
other%grid_lat = self%grid_lat
other%area = self%area
allocate ( other%ak(other%nlevs+1) )
allocate ( other%bk(other%nlevs+1) )
other%ak = self%ak
other%bk = self%bk
other%datapath_in = self%datapath_in
if (.not. self%agrid_vel_rst) then
   allocate ( other%gridstruct%sin_sg (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed  ,9) )
   allocate ( other%gridstruct%cosa_u (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed    ) )
   allocate ( other%gridstruct%cosa_v (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1  ) )
   allocate ( other%gridstruct%cosa_s (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
   allocate ( other%gridstruct%rsin_u (self%bd%isd:self%bd%ied+1,self%bd%jsd:self%bd%jed    ) )
   allocate ( other%gridstruct%rsin_v (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed+1  ) )
   allocate ( other%gridstruct% rsin2 (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
   allocate ( other%gridstruct%   dxa (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
   allocate ( other%gridstruct%   dya (self%bd%isd:self%bd%ied  ,self%bd%jsd:self%bd%jed    ) )
   other%gridstruct%sin_sg = self%gridstruct%sin_sg
   other%gridstruct%cosa_u = self%gridstruct%cosa_u
   other%gridstruct%cosa_v = self%gridstruct%cosa_v
   other%gridstruct%cosa_s = self%gridstruct%cosa_s
   other%gridstruct%rsin_u = self%gridstruct%rsin_u
   other%gridstruct%rsin_v = self%gridstruct%rsin_v
   other%gridstruct% rsin2 = self%gridstruct% rsin2
   other%gridstruct%   dxa = self%gridstruct%   dxa
   other%gridstruct%   dya = self%gridstruct%   dya
   other%gridstruct%sw_corner = self%gridstruct%sw_corner
   other%gridstruct%se_corner = self%gridstruct%se_corner
   other%gridstruct%ne_corner = self%gridstruct%ne_corner
   other%gridstruct%nw_corner = self%gridstruct%nw_corner
endif

end subroutine c_fv3jedi_geo_clone

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_delete(c_key_self) bind(c,name='fv3jedi_geo_delete_f90')

implicit none

integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_geom), pointer :: self

call fv3jedi_geom_registry%get(c_key_self, self)

! deallocate whatever is in self
call mpp_deallocate_domain(self%domain)
deallocate(self%grid_lon)
deallocate(self%grid_lat)
deallocate(self%area)
deallocate(self%ak)
deallocate(self%bk)
if (.not. self%agrid_vel_rst) then
   deallocate ( self%gridstruct%sin_sg )
   deallocate ( self%gridstruct%cosa_u )
   deallocate ( self%gridstruct%cosa_v )
   deallocate ( self%gridstruct%cosa_s )
   deallocate ( self%gridstruct%rsin_u )
   deallocate ( self%gridstruct%rsin_v )
   deallocate ( self%gridstruct% rsin2 )
   deallocate ( self%gridstruct%   dxa )
   deallocate ( self%gridstruct%   dya )
endif

!call fms_end()

call fv3jedi_geom_registry%remove(c_key_self)

end subroutine c_fv3jedi_geo_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_info(c_key_self) bind(c,name='fv3jedi_geo_info_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
type(fv3jedi_geom), pointer :: self

call fv3jedi_geom_registry%get(c_key_self, self)
! get a few numbers back to C++ to print so that one can have a quick idea
! about the definition of the FV3JEDI geom
!c_n = self%...

end subroutine c_fv3jedi_geo_info

! ------------------------------------------------------------------------------

end module fv3jedi_geom_mod
