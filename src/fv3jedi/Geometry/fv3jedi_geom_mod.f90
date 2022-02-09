! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_mod

use netcdf
use mpi
use string_f_c_mod

! atlas uses
use atlas_module, only: atlas_field, atlas_fieldset, atlas_real, atlas_functionspace

! fckit uses
use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

! fms uses
use fms_io_mod,                 only: set_domain, nullify_domain
use fms_mod,                    only: fms_init
use mpp_mod,                    only: mpp_exit, mpp_pe, mpp_npes, mpp_error, FATAL, NOTE
use mpp_domains_mod,            only: domain2D, mpp_deallocate_domain, mpp_define_layout, &
                                      mpp_define_mosaic, mpp_define_io_domain, mpp_domains_exit, &
                                      mpp_domains_set_stack_size
use field_manager_mod,          only: fm_string_len, field_manager_init

! fv3 uses
use fv_arrays_mod,              only: fv_atmos_type, deallocate_fv_atmos_type

! fv3jedi uses
use fields_metadata_mod,         only: fields_metadata, field_metadata
use fv3jedi_constants_mod,       only: ps, rad2deg, kap1, kapr
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_netcdf_utils_mod,    only: nccheck
use fv_init_mod,                 only: fv_init
use fv3jedi_fmsnamelist_mod,     only: fv3jedi_fmsnamelist
use fv3jedi_io_fms_mod,          only: fv3jedi_io_fms, read_fields
use fv3jedi_field_mod,           only: fv3jedi_field

implicit none
private
public :: fv3jedi_geom, getVerticalCoord, getVerticalCoordLogP, initialize, pedges2pmidlayer

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the FV3JEDI model
type :: fv3jedi_geom
  integer :: isd, ied, jsd, jed                                                     !data domain
  integer :: isc, iec, jsc, jec, kec                                                !compute domain
  integer :: npx,npy,npz,ngrid                                                      !x/y/z-dir grid edge points per tile
  integer :: layout(2), io_layout(2)                                                !Processor layouts
  integer :: ntile, ntiles                                                          !Tile number and total
  integer :: iterator_dimension                                                     !iterator dimension
  real(kind=kind_real) :: ptop                                                      !Pressure at top of domain
  type(domain2D) :: domain_fix                                                      !MPP domain
  type(domain2D), pointer :: domain                                                 !MPP domain
  character(len=10) :: interp_method                                                !Interpolation type
  real(kind=kind_real), allocatable, dimension(:)       :: ak, bk                   !Model level coefficients
  real(kind=kind_real), allocatable, dimension(:,:)     :: grid_lon, grid_lat       !Lat/lon centers
  real(kind=kind_real), allocatable, dimension(:,:)     :: egrid_lon, egrid_lat     !Lat/lon edges
  real(kind=kind_real), allocatable, dimension(:)       :: lon_us, lat_us           !Lat/lon centers unstructured
  real(kind=kind_real), allocatable, dimension(:,:)     :: area                     !Grid area
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: orography                !Grid surface elevation
  real(kind=kind_real), allocatable, dimension(:,:)     :: dx, dy                   !dx/dy at edges
  real(kind=kind_real), allocatable, dimension(:,:)     :: dxc, dyc                 !dx/dy c grid
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: grid, vlon, vlat
  real(kind=kind_real), allocatable, dimension(:)       :: edge_vect_n, edge_vect_e
  real(kind=kind_real), allocatable, dimension(:)       :: edge_vect_s, edge_vect_w
  real(kind=kind_real), allocatable, dimension(:,:,:,:) :: es, ew
  real(kind=kind_real), allocatable, dimension(:,:)     :: a11, a12, a21, a22
  type(fckit_mpi_comm) :: f_comm
  type(fields_metadata) :: fields
  ! Vertical Coordinate
  real(kind=kind_real), allocatable, dimension(:)       :: vCoord                   !Model vertical coordinate
  real(kind=kind_real), allocatable, dimension(:,:)     :: surface_pressure         !Grid surface pressure
  ! For D to (A to) C grid
  real(kind=kind_real), allocatable, dimension(:,:)     :: rarea
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: sin_sg
  real(kind=kind_real), allocatable, dimension(:,:)     :: cosa_u
  real(kind=kind_real), allocatable, dimension(:,:)     :: cosa_v
  real(kind=kind_real), allocatable, dimension(:,:)     :: cosa_s
  real(kind=kind_real), allocatable, dimension(:,:)     :: rsin_u
  real(kind=kind_real), allocatable, dimension(:,:)     :: rsin_v
  real(kind=kind_real), allocatable, dimension(:,:)     :: rsin2
  real(kind=kind_real), allocatable, dimension(:,:)     :: dxa, dya
  logical :: ne_corner, se_corner, sw_corner, nw_corner
  logical :: nested = .false.
  logical :: bounded_domain = .false.
  logical :: logp = .false.

  integer :: grid_type = 0
  logical :: dord4 = .true.
  type(atlas_functionspace) :: afunctionspace
  contains
    procedure, public :: create
    procedure, public :: clone
    procedure, public :: delete
    procedure, public :: set_atlas_lonlat
    procedure, public :: fill_atlas_fieldset
end type fv3jedi_geom

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine initialize(conf, comm)

type(fckit_configuration), intent(in) :: conf
type(fckit_mpi_comm),      intent(in) :: comm

integer :: stackmax
character(len=1024) :: nml_filename, field_table_filename
character(len=:), allocatable :: str

! Path for input.nml file
call conf%get_or_die("namelist filename",str)
if (len(str) > 1024) call abor1_ftn("Length of fms namelist filename too long")
nml_filename = str
deallocate(str)

! Field table file
call conf%get_or_die("field table filename",str)
if (len(str) > fm_string_len) call abor1_ftn("Length of fms field table filename too long")
field_table_filename = str
deallocate(str)

! Initialize fms, mpp, etc.
call fms_init(localcomm=comm%communicator(), alt_input_nml_path = nml_filename)

! Set max stacksize
call conf%get_or_die("stackmax", stackmax)
call mpp_domains_set_stack_size(stackmax)

! Initialize the tracers
call field_manager_init(table_name = field_table_filename)

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine create(self, conf, comm, nlevs)

!Arguments
class(fv3jedi_geom), target, intent(inout) :: self
type(fckit_configuration),   intent(in)    :: conf
type(fckit_mpi_comm),        intent(in)    :: comm
integer,                     intent(out)   :: nlevs

!Locals
character(len=256)                    :: file_akbk
type(fv_atmos_type), allocatable      :: Atm(:)
logical, allocatable                  :: grids_on_this_pe(:)
integer                               :: i, j, jj, gtile
integer                               :: p_split = 1
integer                               :: ncstat, ncid, akvarid, bkvarid, readdim, dcount
integer, dimension(nf90_max_var_dims) :: dimids, dimlens

character(len=:), allocatable :: str
logical :: do_write_geom = .false.
logical :: logp = .false.
integer :: iterator_dimension = 2

type(fv3jedi_fmsnamelist) :: fmsnamelist

! Add the communicator to the geometry
! ------------------------------------
self%f_comm = comm

! Interpolation type
! ------------------
call conf%get_or_die("interpolation method",str)
self%interp_method = str
deallocate(str)

call conf%get_or_die("iterator dimension", iterator_dimension)
self%iterator_dimension = iterator_dimension

! Update the fms namelist with this Geometry
! ------------------------------------------
call fmsnamelist%replace_namelist(conf)

!Intialize using the model setup routine
! --------------------------------------
call fv_init(Atm, 300.0_kind_real, grids_on_this_pe, p_split, gtile, .true.)

! Copy relevant contents of Atm
! -----------------------------
self%isd = Atm(1)%bd%isd
self%ied = Atm(1)%bd%ied
self%jsd = Atm(1)%bd%jsd
self%jed = Atm(1)%bd%jed

self%isc = Atm(1)%bd%isc
self%iec = Atm(1)%bd%iec
self%jsc = Atm(1)%bd%jsc
self%jec = Atm(1)%bd%jec
self%kec = Atm(1)%npz

self%ntile  = gtile
self%ntiles = Atm(1)%flagstruct%ntiles

self%npx = Atm(1)%npx
self%npy = Atm(1)%npy
self%npz = Atm(1)%npz

nlevs = self%npz

self%layout(1) = Atm(1)%layout(1)
self%layout(2) = Atm(1)%layout(2)
self%io_layout(1) = Atm(1)%io_layout(1)
self%io_layout(2) = Atm(1)%io_layout(2)

!Allocatable arrays
allocate(self%ak(self%npz+1) )
allocate(self%bk(self%npz+1) )

allocate(self%grid_lon   (self%isd  :self%ied,  self%jsd  :self%jed  ))
allocate(self%grid_lat   (self%isd  :self%ied,  self%jsd  :self%jed  ))
allocate(self%egrid_lon  (self%isd  :self%ied+1,self%jsd  :self%jed+1))
allocate(self%egrid_lat  (self%isd  :self%ied+1,self%jsd  :self%jed+1))
allocate(self%area       (self%isd  :self%ied,  self%jsd  :self%jed  ))
allocate(self%dx         (self%isd  :self%ied  ,self%jsd  :self%jed+1))
allocate(self%dy         (self%isd  :self%ied+1,self%jsd  :self%jed  ))
allocate(self%dxc        (self%isd  :self%ied+1,self%jsd  :self%jed  ))
allocate(self%dyc        (self%isd  :self%ied  ,self%jsd  :self%jed+1))

allocate(self%grid       (self%isd  :self%ied+1,self%jsd  :self%jed+1,2))
allocate(self%vlon       (self%isc-2:self%iec+2,self%jsc-2:self%jec+2,3))
allocate(self%vlat       (self%isc-2:self%iec+2,self%jsc-2:self%jec+2,3))

allocate(self%edge_vect_n(self%isd:self%ied))
allocate(self%edge_vect_e(self%jsd:self%jed))
allocate(self%edge_vect_s(self%isd:self%ied))
allocate(self%edge_vect_w(self%jsd:self%jed))

allocate(self%es(3,self%isd:self%ied  ,self%jsd:self%jed+1,2))
allocate(self%ew(3,self%isd:self%ied+1,self%jsd:self%jed,  2))

allocate(self%a11(self%isc-1:self%iec+1,self%jsc-1:self%jec+1) )
allocate(self%a12(self%isc-1:self%iec+1,self%jsc-1:self%jec+1) )
allocate(self%a21(self%isc-1:self%iec+1,self%jsc-1:self%jec+1) )
allocate(self%a22(self%isc-1:self%iec+1,self%jsc-1:self%jec+1) )

allocate(self%rarea (self%isd:self%ied  ,self%jsd:self%jed  ))
allocate(self%sin_sg(self%isd:self%ied  ,self%jsd:self%jed  ,9))
allocate(self%cosa_u(self%isd:self%ied+1,self%jsd:self%jed  ))
allocate(self%cosa_v(self%isd:self%ied  ,self%jsd:self%jed+1))
allocate(self%cosa_s(self%isd:self%ied  ,self%jsd:self%jed  ))
allocate(self%rsin_u(self%isd:self%ied+1,self%jsd:self%jed  ))
allocate(self%rsin_v(self%isd:self%ied  ,self%jsd:self%jed+1))
allocate(self%rsin2 (self%isd:self%ied  ,self%jsd:self%jed  ))
allocate(self%dxa   (self%isd:self%ied  ,self%jsd:self%jed  ))
allocate(self%dya   (self%isd:self%ied  ,self%jsd:self%jed  ))

! ak and bk hybrid coordinate coefficients
! ----------------------------------------
if (self%npz > 1) then

  ! Set path/filename for ak and bk file
  call conf%get_or_die("akbk",str)
  file_akbk = str

  !Open file
  call nccheck ( nf90_open(file_akbk, nf90_nowrite, ncid), "nf90_open "//file_akbk )

  !Search for ak in the file
  ncstat = nf90_inq_varid(ncid, "ak", akvarid)
  if(ncstat /= nf90_noerr) call abor1_ftn("Failed to find ak in file "//file_akbk)

  !Search for bk in the file
  ncstat = nf90_inq_varid(ncid, "bk", bkvarid)
  if(ncstat /= nf90_noerr) call abor1_ftn("Failed to find bk in file "//file_akbk)

  ! Check that dimension of ak/bk in the file match vertical levels of model
  dimids = 0
  call nccheck ( nf90_inquire_variable(ncid, akvarid, dimids = dimids), "nf90_inq_var ak" )
  readdim = -1
  dcount = 0
  do i = 1,nf90_max_var_dims
    if (dimids(i) > 0) then
       call nccheck( nf90_inquire_dimension(ncid, dimids(i), len = dimlens(i)), &
                     "nf90_inquire_dimension" )
       if (dimlens(i) == self%npz+1) then
          readdim = i
       endif
       dcount = dcount + 1
    endif
  enddo
  if (readdim == -1) call abor1_ftn("ak/bk in file does not match dimension of npz from input.nml")

  !Read ak and bk from the file
  call nccheck( nf90_get_var(ncid, akvarid, self%ak), "fv3jedi_geom, nf90_get_var ak" )
  call nccheck( nf90_get_var(ncid, bkvarid, self%bk), "fv3jedi_geom, nf90_get_var bk" )
else
  self%ak = 0.0_kind_real
  self%bk = 0.0_kind_real
endif

! Arrays from the Atm Structure
! --------------------------------

self%grid_lon  = real(Atm(1)%gridstruct%agrid_64(:,:,1),kind_real)
self%grid_lat  = real(Atm(1)%gridstruct%agrid_64(:,:,2),kind_real)
self%egrid_lon = real(Atm(1)%gridstruct%grid_64(:,:,1),kind_real)
self%egrid_lat = real(Atm(1)%gridstruct%grid_64(:,:,2),kind_real)
self%area      = real(Atm(1)%gridstruct%area_64,kind_real)
self%dx        = real(Atm(1)%gridstruct%dx ,kind_real)
self%dy        = real(Atm(1)%gridstruct%dy ,kind_real)
self%dxc       = real(Atm(1)%gridstruct%dxc,kind_real)
self%dyc       = real(Atm(1)%gridstruct%dyc,kind_real)

self%grid      = real(Atm(1)%gridstruct%grid,kind_real)
self%vlon      = real(Atm(1)%gridstruct%vlon,kind_real)
self%vlat      = real(Atm(1)%gridstruct%vlat,kind_real)

self%edge_vect_n = real(Atm(1)%gridstruct%edge_vect_n,kind_real)
self%edge_vect_e = real(Atm(1)%gridstruct%edge_vect_e,kind_real)
self%edge_vect_s = real(Atm(1)%gridstruct%edge_vect_s,kind_real)
self%edge_vect_w = real(Atm(1)%gridstruct%edge_vect_w,kind_real)

self%es = real(Atm(1)%gridstruct%es,kind_real)
self%ew = real(Atm(1)%gridstruct%ew,kind_real)

self%a11 = real(Atm(1)%gridstruct%a11,kind_real)
self%a12 = real(Atm(1)%gridstruct%a12,kind_real)
self%a21 = real(Atm(1)%gridstruct%a21,kind_real)
self%a22 = real(Atm(1)%gridstruct%a22,kind_real)

self%rarea     = real(Atm(1)%gridstruct%rarea ,kind_real)
self%sin_sg    = real(Atm(1)%gridstruct%sin_sg,kind_real)
self%cosa_u    = real(Atm(1)%gridstruct%cosa_u,kind_real)
self%cosa_v    = real(Atm(1)%gridstruct%cosa_v,kind_real)
self%cosa_s    = real(Atm(1)%gridstruct%cosa_s,kind_real)
self%rsin_u    = real(Atm(1)%gridstruct%rsin_u,kind_real)
self%rsin_v    = real(Atm(1)%gridstruct%rsin_v,kind_real)
self%rsin2     = real(Atm(1)%gridstruct%rsin2 ,kind_real)
self%dxa       = real(Atm(1)%gridstruct%dxa   ,kind_real)
self%dya       = real(Atm(1)%gridstruct%dya   ,kind_real)
self%ne_corner = Atm(1)%gridstruct%ne_corner
self%se_corner = Atm(1)%gridstruct%se_corner
self%sw_corner = Atm(1)%gridstruct%sw_corner
self%nw_corner = Atm(1)%gridstruct%nw_corner
self%nested    = Atm(1)%gridstruct%nested
self%bounded_domain =  Atm(1)%gridstruct%bounded_domain

allocate(self%vCoord(self%npz))
allocate(self%surface_pressure(self%isd:self%ied, self%jsd:self%jed))

self%surface_pressure = real(Atm(1)%ps ,kind_real)

call conf%get_or_die("logp",logp)
self%logp = logp


!Unstructured lat/lon
self%ngrid = (self%iec-self%isc+1)*(self%jec-self%jsc+1)
allocate(self%lat_us(self%ngrid))
allocate(self%lon_us(self%ngrid))

jj = 0
do j = self%jsc,self%jec
  do i = self%isc,self%iec
     jj = jj + 1
     self%lat_us(jj) = self%grid_lat(i,j)
     self%lon_us(jj) = self%grid_lon(i,j)
  enddo
enddo

!Set Ptop
self%ptop = self%ak(1)

!Done with the Atm stucture here
call deallocate_fv_atmos_type(Atm(1))
deallocate(Atm)
deallocate(grids_on_this_pe)

!Resetup domain to avoid risk of copied pointers
call setup_domain( self%domain_fix, self%npx-1, self%npy-1, &
                   self%ntiles, self%layout, self%io_layout, 3)

self%domain => self%domain_fix
call nullify_domain()

! Optionally write the geometry to file
! -------------------------------------
call conf%get_or_die("write geom",do_write_geom)

if (do_write_geom) then
  call write_geom(self)
endif

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine clone(self, other, fields)

class(fv3jedi_geom),        intent(inout) :: self
type(fv3jedi_geom), target, intent(in)    :: other
type(fields_metadata),      intent(in)    :: fields

allocate(self%ak(other%npz+1) )
allocate(self%bk(other%npz+1) )

allocate(self%grid_lon   (other%isd  :other%ied,  other%jsd  :other%jed  ))
allocate(self%grid_lat   (other%isd  :other%ied,  other%jsd  :other%jed  ))
allocate(self%egrid_lon  (other%isd  :other%ied+1,other%jsd  :other%jed+1))
allocate(self%egrid_lat  (other%isd  :other%ied+1,other%jsd  :other%jed+1))
allocate(self%area       (other%isd  :other%ied,  other%jsd  :other%jed  ))
allocate(self%dx         (other%isd  :other%ied  ,other%jsd  :other%jed+1))
allocate(self%dy         (other%isd  :other%ied+1,other%jsd  :other%jed  ))
allocate(self%dxc        (other%isd  :other%ied+1,other%jsd  :other%jed  ))
allocate(self%dyc        (other%isd  :other%ied  ,other%jsd  :other%jed+1))

allocate(self%grid       (other%isd  :other%ied+1,other%jsd  :other%jed+1,2))
allocate(self%vlon       (other%isc-2:other%iec+2,other%jsc-2:other%jec+2,3))
allocate(self%vlat       (other%isc-2:other%iec+2,other%jsc-2:other%jec+2,3))

allocate(self%edge_vect_n(other%isd:other%ied))
allocate(self%edge_vect_e(other%jsd:other%jed))
allocate(self%edge_vect_s(other%isd:other%ied))
allocate(self%edge_vect_w(other%jsd:other%jed))

allocate(self%es(3,other%isd:other%ied  ,other%jsd:other%jed+1,2))
allocate(self%ew(3,other%isd:other%ied+1,other%jsd:other%jed,  2))

allocate(self%a11(other%isc-1:other%iec+1,other%jsc-1:other%jec+1) )
allocate(self%a12(other%isc-1:other%iec+1,other%jsc-1:other%jec+1) )
allocate(self%a21(other%isc-1:other%iec+1,other%jsc-1:other%jec+1) )
allocate(self%a22(other%isc-1:other%iec+1,other%jsc-1:other%jec+1) )

allocate(self%rarea (other%isd:other%ied  ,other%jsd:other%jed  ))
allocate(self%sin_sg(other%isd:other%ied  ,other%jsd:other%jed  ,9))
allocate(self%cosa_u(other%isd:other%ied+1,other%jsd:other%jed  ))
allocate(self%cosa_v(other%isd:other%ied  ,other%jsd:other%jed+1))
allocate(self%cosa_s(other%isd:other%ied  ,other%jsd:other%jed  ))
allocate(self%rsin_u(other%isd:other%ied+1,other%jsd:other%jed  ))
allocate(self%rsin_v(other%isd:other%ied  ,other%jsd:other%jed+1))
allocate(self%rsin2 (other%isd:other%ied  ,other%jsd:other%jed  ))
allocate(self%dxa   (other%isd:other%ied  ,other%jsd:other%jed  ))
allocate(self%dya   (other%isd:other%ied  ,other%jsd:other%jed  ))

allocate(self%lat_us(other%ngrid))
allocate(self%lon_us(other%ngrid))

self%npx             = other%npx
self%npy             = other%npy
self%npz             = other%npz
self%ngrid           = other%ngrid
self%layout          = other%layout
self%io_layout       = other%io_layout
self%isc             = other%isc
self%isd             = other%isd
self%iec             = other%iec
self%ied             = other%ied
self%jsc             = other%jsc
self%jsd             = other%jsd
self%jec             = other%jec
self%jed             = other%jed
self%ntile           = other%ntile
self%ntiles          = other%ntiles
self%iterator_dimension = other%iterator_dimension

self%ptop            = other%ptop
self%ak              = other%ak
self%bk              = other%bk
self%grid_lon        = other%grid_lon
self%grid_lat        = other%grid_lat
self%egrid_lon       = other%egrid_lon
self%egrid_lat       = other%egrid_lat
self%area            = other%area
self%dx              = other%dx
self%dy              = other%dy
self%dxc             = other%dxc
self%dyc             = other%dyc
self%grid            = other%grid
self%vlon            = other%vlon
self%vlat            = other%vlat
self%edge_vect_n     = other%edge_vect_n
self%edge_vect_e     = other%edge_vect_e
self%edge_vect_s     = other%edge_vect_s
self%edge_vect_w     = other%edge_vect_w
self%es              = other%es
self%ew              = other%ew
self%a11             = other%a11
self%a12             = other%a12
self%a21             = other%a21
self%a22             = other%a22
self%f_comm          = other%f_comm
self%interp_method   = other%interp_method

self%rarea     = other%rarea
self%sin_sg    = other%sin_sg
self%cosa_u    = other%cosa_u
self%cosa_v    = other%cosa_v
self%cosa_s    = other%cosa_s
self%rsin_u    = other%rsin_u
self%rsin_v    = other%rsin_v
self%rsin2     = other%rsin2
self%dxa       = other%dxa
self%dya       = other%dya
self%ne_corner = other%ne_corner
self%se_corner = other%se_corner
self%sw_corner = other%sw_corner
self%nw_corner = other%nw_corner

self%domain => other%domain

self%afunctionspace = atlas_functionspace(other%afunctionspace%c_ptr())

self%fields = fields

self%lat_us = other%lat_us
self%lon_us = other%lon_us

self%nested = other%nested
self%bounded_domain = other%bounded_domain

self%logp = other%logp

if (allocated(other%orography)) then
  allocate(self%orography(other%isc:other%iec,other%jsc:other%jec,1))
  self%orography = other%orography
endif

end subroutine clone

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_geom), intent(inout) :: self

! Deallocate
deallocate(self%ak)
deallocate(self%bk)
deallocate(self%grid_lon)
deallocate(self%grid_lat)
deallocate(self%egrid_lon)
deallocate(self%egrid_lat)
deallocate(self%area)
deallocate(self%dx)
deallocate(self%dy)
deallocate(self%dxc)
deallocate(self%dyc)
deallocate(self%grid)
deallocate(self%vlon)
deallocate(self%vlat)
deallocate(self%edge_vect_n)
deallocate(self%edge_vect_e)
deallocate(self%edge_vect_s)
deallocate(self%edge_vect_w)
deallocate(self%es)
deallocate(self%ew)
deallocate(self%a11)
deallocate(self%a12)
deallocate(self%a21)
deallocate(self%a22)

deallocate(self%rarea)
deallocate(self%sin_sg)
deallocate(self%cosa_u)
deallocate(self%cosa_v)
deallocate(self%cosa_s)
deallocate(self%rsin_u)
deallocate(self%rsin_v)
deallocate(self%rsin2 )
deallocate(self%dxa   )
deallocate(self%dya   )

deallocate(self%lat_us)
deallocate(self%lon_us)

if (allocated(self%orography)) deallocate(self%orography)

! Required memory leak, since copying this causes problems
!call mpp_deallocate_domain(self%domain_fix)

call self%afunctionspace%final()

! Could finalize the fms routines. Possibly needs to be done only when key = 0
!call fms_io_exit
!call mpp_domains_exit
!call mpp_exit

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine set_atlas_lonlat(self, afieldset)

!Arguments
class(fv3jedi_geom),  intent(inout) :: self
type(atlas_fieldset), intent(inout) :: afieldset

!Locals
real(kind_real), pointer :: real_ptr(:,:)
type(atlas_field) :: afield

! Create lon/lat field
afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,(self%iec-self%isc+1)*(self%jec-self%jsc+1)/))
call afield%data(real_ptr)
real_ptr(1,:) = rad2deg*reshape(self%grid_lon(self%isc:self%iec,self%jsc:self%jec),(/self%ngrid/))
real_ptr(2,:) = rad2deg*reshape(self%grid_lat(self%isc:self%iec,self%jsc:self%jec),(/self%ngrid/))
call afieldset%add(afield)

end subroutine set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

subroutine fill_atlas_fieldset(self, afieldset)

!Arguments
class(fv3jedi_geom),  intent(inout) :: self
type(atlas_fieldset), intent(inout) :: afieldset

!Locals
integer :: jl
real(kind=kind_real) :: sigmaup, sigmadn
real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
type(atlas_field) :: afield
real(kind=kind_real) :: plevli(self%npz+1),logp(self%npz)

! Add area
afield = self%afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=0)
call afield%data(real_ptr_1)
real_ptr_1 = reshape(self%area(self%isc:self%iec,self%jsc:self%jec),(/self%ngrid/))
call afieldset%add(afield)
call afield%final()

! Add vertical unit
afield = self%afunctionspace%create_field(name='vunit', kind=atlas_real(kind_real), levels=self%npz)
call afield%data(real_ptr_2)

if (.not. self%logp) then
   do jl=1,self%npz
      sigmaup = self%ak(jl+1)/ps+self%bk(jl+1) ! si are now sigmas
      sigmadn = self%ak(jl  )/ps+self%bk(jl  )
      real_ptr_2(jl,:) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
   enddo
else
   call getVerticalCoordLogP(self,logp,self%npz,ps)
   do jl=1,self%npz
      real_ptr_2(jl,:) = logp(jl)
   enddo
endif

call afieldset%add(afield)
call afield%final()

end subroutine fill_atlas_fieldset

! --------------------------------------------------------------------------------------------------

subroutine setup_domain(domain, nx, ny, ntiles, layout_in, io_layout, halo)

 type(domain2D),   intent(inout) :: domain
 integer,          intent(in)    :: nx, ny, ntiles
 integer,          intent(in)    :: layout_in(:), io_layout(:)
 integer,          intent(in)    :: halo

 integer                              :: pe, npes, npes_per_tile, tile
 integer                              :: num_contact, num_alloc
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

  ! select case based off of 1 or 6 tiles
  select case(ntiles)
  case ( 1 ) ! FV3-SAR
    num_contact = 0
  case ( 6 ) ! FV3 global
    num_contact = 12
  case default
    call mpp_error(FATAL, "setup_domain: ntiles != 1 or 6")
  end select

  do n = 1, ntiles
     global_indices(:,n) = (/1,nx,1,ny/)
     layout2D(:,n)       = layout
     pe_start(n)         = (n-1)*npes_per_tile
     pe_end(n)           = n*npes_per_tile-1
  enddo
  num_alloc = max(1, num_contact)
  ! this code copied from domain_decomp in fv_mp_mod.f90
  allocate(tile1(num_alloc), tile2(num_alloc) )
  allocate(tile_id(ntiles))
  allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
  allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )
  ! select case based off of 1 or 6 tiles
  select case(ntiles)
  case ( 1 ) ! FV3-SAR
    ! No contacts, do nothing
  case ( 6 ) ! FV3 global
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
  end select
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

! --------------------------------------------------------------------------------------------------

subroutine write_geom(self)

  type(fv3jedi_geom), intent(in) :: self

  type(fckit_mpi_comm) :: f_comm
  character(len=255) :: filename
  integer :: ncid, xf_dimid, yf_dimid, xv_dimid, yv_dimid, ti_dimid, pe_dimid
  integer :: mydims(3,3), ijdims(1), ijdimf(1), tmpij(1)
  integer :: varid(8)


  ! Pointer to fv3jedi geom communicator
  f_comm = self%f_comm

  write(filename,"(A9,I0.4,A4)") 'fv3grid_c', self%npx-1, '.nc4'

  ! Create and open the file for parallel write
  call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                             comm = f_comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

  !Dimensions
  call nccheck ( nf90_def_dim(ncid, 'fxdim', self%npx-1   , xf_dimid), "nf90_def_dim fxdim" )
  call nccheck ( nf90_def_dim(ncid, 'fydim', self%npy-1   , yf_dimid), "nf90_def_dim fydim" )
  call nccheck ( nf90_def_dim(ncid, 'vxdim', self%npx     , xv_dimid), "nf90_def_dim vxdim" )
  call nccheck ( nf90_def_dim(ncid, 'vydim', self%npy     , yv_dimid), "nf90_def_dim vydim" )
  call nccheck ( nf90_def_dim(ncid, 'ntile', 6            , ti_dimid), "nf90_def_dim ntile" )
  call nccheck ( nf90_def_dim(ncid, 'nproc', f_comm%size(), pe_dimid), "nf90_def_dim ntile" )

  !Define variables
  call nccheck( nf90_def_var(ncid, "flons", NF90_DOUBLE, (/ xf_dimid, yf_dimid, ti_dimid /), varid(1)), "nf90_def_var flons" )
  call nccheck( nf90_put_att(ncid, varid(1), "long_name", "longitude of faces") )
  call nccheck( nf90_put_att(ncid, varid(1), "units", "degrees_east") )

  call nccheck( nf90_def_var(ncid, "flats", NF90_DOUBLE, (/ xf_dimid, yf_dimid, ti_dimid /), varid(2)), "nf90_def_var flats" )
  call nccheck( nf90_put_att(ncid, varid(2), "long_name", "latitude of faces") )
  call nccheck( nf90_put_att(ncid, varid(2), "units", "degrees_north") )

  call nccheck( nf90_def_var(ncid, "vlons", NF90_DOUBLE, (/ xv_dimid, yv_dimid, ti_dimid /), varid(3)), "nf90_def_var vlons" )
  call nccheck( nf90_put_att(ncid, varid(3), "long_name", "longitude of vertices") )
  call nccheck( nf90_put_att(ncid, varid(3), "units", "degrees_east") )

  call nccheck( nf90_def_var(ncid, "vlats", NF90_DOUBLE, (/ xv_dimid, yv_dimid, ti_dimid /), varid(4)), "nf90_def_var vlats" )
  call nccheck( nf90_put_att(ncid, varid(4), "long_name", "latitude of vertices") )
  call nccheck( nf90_put_att(ncid, varid(4), "units", "degrees_north") )

  call nccheck( nf90_def_var(ncid, "isc", NF90_INT, (/ pe_dimid /), varid(5)), "nf90_def_var isc" )
  call nccheck( nf90_put_att(ncid, varid(5), "long_name", "starting index i direction") )
  call nccheck( nf90_put_att(ncid, varid(5), "units", "1") )

  call nccheck( nf90_def_var(ncid, "iec", NF90_INT, (/ pe_dimid /), varid(6)), "nf90_def_var iec" )
  call nccheck( nf90_put_att(ncid, varid(6), "long_name", "ending index i direction") )
  call nccheck( nf90_put_att(ncid, varid(6), "units", "1") )

  call nccheck( nf90_def_var(ncid, "jsc", NF90_INT, (/ pe_dimid /), varid(7)), "nf90_def_var jsc" )
  call nccheck( nf90_put_att(ncid, varid(7), "long_name", "starting index j direction") )
  call nccheck( nf90_put_att(ncid, varid(7), "units", "1") )

  call nccheck( nf90_def_var(ncid, "jec", NF90_INT, (/ pe_dimid /), varid(8)), "nf90_def_var jec" )
  call nccheck( nf90_put_att(ncid, varid(8), "long_name", "ending index j direction") )
  call nccheck( nf90_put_att(ncid, varid(8), "units", "1") )

  ! End define mode
  call nccheck( nf90_enddef(ncid), "nf90_enddef" )

  ! Write variables
  mydims(1,1) = 1;          mydims(2,1) = self%npx-1
  mydims(1,2) = 1;          mydims(2,2) = self%npy-1
  mydims(1,3) = self%ntile; mydims(2,3) = 1

  call nccheck( nf90_put_var( ncid, varid(1), self%grid_lon(self%isc:self%iec,self%jsc:self%jec), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var flons" )

  call nccheck( nf90_put_var( ncid, varid(2), self%grid_lat(self%isc:self%iec,self%jsc:self%jec), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var flats" )

  mydims(1,1) = 1;          mydims(2,1) = self%npx
  mydims(1,2) = 1;          mydims(2,2) = self%npy
  mydims(1,3) = self%ntile; mydims(2,3) = 1

  call nccheck( nf90_put_var( ncid, varid(3), self%egrid_lon(self%isc:self%iec+1,self%jsc:self%jec+1), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var vlons" )

  call nccheck( nf90_put_var( ncid, varid(4), self%egrid_lat(self%isc:self%iec+1,self%jsc:self%jec+1), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var vlats" )

  ijdims(1) = f_comm%rank()+1
  ijdimf(1) = 1

  tmpij = self%isc
  call nccheck( nf90_put_var( ncid, varid(5), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var isc" )

  tmpij = self%iec
  call nccheck( nf90_put_var( ncid, varid(6), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var iec" )

  tmpij = self%jsc
  call nccheck( nf90_put_var( ncid, varid(7), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var jsc" )

  tmpij = self%jec
  call nccheck( nf90_put_var( ncid, varid(8), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var jec" )

  ! Close the file
  call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine write_geom

!----------------------------------------------------------------------------
! 1d pressure_edge to pressure_mid
!----------------------------------------------------------------------------

subroutine pedges2pmidlayer(npz,ptype,pe1d,p1d)
 integer,              intent(in)  :: npz       !number of model layers
 character(len=*),     intent(in)  :: ptype     !midlayer pressure definition: 'average' or 'Philips'
 real(kind=kind_real), intent(in)  :: pe1d(npz+1) !pressure edge
 real(kind=kind_real), intent(out) :: p1d(npz)    !pressure mid

 select case (ptype)
   case('Philips')
     p1d = ((pe1d(2:npz+1)**kap1 - pe1d(1:npz)**kap1)/&
            (kap1*(pe1d(2:npz+1) - pe1d(1:npz))))**kapr
   case default
     p1d = 0.5*(pe1d(2:npz+1) + pe1d(1:npz))
 end select

end subroutine pedges2pmidlayer

!--------------------------------------------------------------------------------------------------
subroutine getVerticalCoord(self, vc, npz, psurf)
  ! returns log(pressure) at mid level of the vertical column with surface
  ! prsssure of psurf
  ! coded using an example from Jeff Whitaker used in GSI ENKF pacakge

  type(fv3jedi_geom),   intent(in) :: self
  integer,              intent(in) :: npz
  real(kind=kind_real), intent(in) :: psurf
  real(kind=kind_real), intent(out) :: vc(npz)

  real(kind=kind_real) :: plevli(npz+1), p(npz)
  integer :: k

  ! compute interface pressure
  do k=1,npz+1
    plevli(k) = self%ak(k) + self%bk(k)*psurf
  enddo

  ! compute presure at mid level and convert it to logp
  call pedges2pmidlayer(npz,'Philips',plevli,vc)

end subroutine getVerticalCoord

!--------------------------------------------------------------------------------------------------
subroutine getVerticalCoordLogP(self, vc, npz, psurf)
  ! returns log(pressure) at mid level of the vertical column with surface prsssure of psurf
  ! coded using an example from Jeff Whitaker used in GSI ENKF pacakge

  type(fv3jedi_geom),   intent(in) :: self
  integer,              intent(in) :: npz
  real(kind=kind_real), intent(in) :: psurf
  real(kind=kind_real), intent(out) :: vc(npz)

  real(kind=kind_real) :: p(npz)

  call getVerticalCoord(self, p, npz, psurf)
  vc = - log(p)

end subroutine getVerticalCoordLogP
! --------------------------------------------------------------------------------------------------

end module fv3jedi_geom_mod
