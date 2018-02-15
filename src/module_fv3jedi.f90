module fv3jedi_mod

!    #include <fms_platform.h>

    use kinds
    use mpp_mod,         only: mpp_pe, mpp_npes, mpp_init, mpp_exit
    use mpp_mod,         only: stdout, mpp_error, FATAL, NOTE
    use mpp_mod,         only: input_nml_file
    use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_mosaic
    use mpp_domains_mod, only: mpp_domains_init, mpp_domains_exit
    use mpp_domains_mod, only: mpp_domains_set_stack_size, mpp_define_io_domain
    use mpp_io_mod,      only: mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY
    use fms_io_mod,      only: fms_io_init, fms_io_exit
    use fms_io_mod,      only: file_exist

    !nicas
    use type_geom,  only: geomtype
    use type_odata, only: odatatype

    implicit none

    private

    public :: fv_atmos_type
    public :: fv_grid_bounds_type
    public :: allocate_fv_atmos_type
    public :: setup_domain
    public :: fv_grid_type
    public :: fv3jedi_interp_type

type fv_atmos_type
    ! A lean version from the model (fv_arrays.F90)

    logical :: allocated = .false.
    logical :: hydrostatic = .false.   ! nonhydrostatic fields in restart?
    logical :: agrid_vel_rst = .false. ! agrid winds in restart?

!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.
    real(kind=kind_real), allocatable, dimension(:,:,:) :: u      ! D grid zonal wind (m/s)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: v      ! D grid meridional wind (m/s)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: pt     ! temperature (K)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: delp   ! pressure thickness (pascal)
    real(kind=kind_real), allocatable, dimension(:,:,:,:) :: q    ! tracers (specific humidity and prognostic constituents)

!----------------------
! non-hydrostatic state:
!----------------------------------------------------------------------
    real(kind=kind_real), allocatable, dimension(:,:,:) :: w      ! cell center vertical wind (m/s)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: delz   ! layer thickness (meters)

!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real(kind=kind_real), allocatable, dimension(:,:) :: phis     ! Surface geopotential (g*Z_surf)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: ua     ! (ua, va) are mostly used as the A grid winds
    real(kind=kind_real), allocatable, dimension(:,:,:) :: va

    integer :: calendar_type
    integer, dimension(6) :: date
    integer, dimension(6) :: date_init

end type fv_atmos_type

type fv_grid_bounds_type
    ! A lean version from the model (fv_arrays.F90)

    integer :: isd, ied, jsd, jed ! data domain
    integer :: isc, iec, jsc, jec ! compute domain

end type fv_grid_bounds_type

type :: fv_grid_type
    ! A lean version from the model (fv_arrays.F90)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: sin_sg
  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_u
  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_v
  real(kind=kind_real), allocatable, dimension(:,:)   :: cosa_s
  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin_u
  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin_v
  real(kind=kind_real), allocatable, dimension(:,:)   :: rsin2
  real(kind=kind_real), allocatable, dimension(:,:)   :: dxa, dya
  logical :: sw_corner, se_corner, ne_corner, nw_corner
endtype

type :: fv3jedi_interp_type
  logical         :: interp_initialized = .FALSE.
  type(geomtype)  :: geom
  type(odatatype) :: odata
endtype

contains

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

    !print *,'layout_in,io_layout',layout_in,io_layout
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

    if (io_layout(1) /= 1 .or. io_layout(2) /= 1) &
        call mpp_define_io_domain(domain, io_layout)

    deallocate(pe_start, pe_end)
    deallocate(layout2D, global_indices)
    deallocate(tile1, tile2, tile_id)
    deallocate(istart1, iend1, jstart1, jend1)
    deallocate(istart2, iend2, jstart2, jend2)

end subroutine setup_domain

subroutine allocate_fv_atmos_type(Atm, &
                                  isd, ied, jsd, jed, &
                                  isc, iec, jsc, jec, &
                                  nz, nq, hydrostatic, agrid_vel_rst)

    !WARNING: Before calling this routine, be sure to have set up the
    ! proper domain parameters.

    implicit none

    type(fv_atmos_type), intent(inout), target :: Atm
    logical, intent(in) :: hydrostatic, agrid_vel_rst
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: isc, iec, jsc, jec
    integer, intent(in) :: nz, nq

    Atm%hydrostatic = hydrostatic
    Atm%agrid_vel_rst = agrid_vel_rst

    if (Atm%allocated) return

    ! what is the difference between compute and data domains (ied,jed vs iec,jec)?
    ! why is data domain used here?

    allocate (    Atm%u(isd:ied,  jsd:jed+1,nz) )
    allocate (    Atm%v(isd:ied+1,jsd:jed,  nz) )
    allocate (   Atm%pt(isd:ied,  jsd:jed,  nz) )
    allocate ( Atm%delp(isd:ied,  jsd:jed,  nz) )
    allocate (    Atm%q(isd:ied,  jsd:jed,  nz, nq) )
    allocate ( Atm%phis(isd:ied,  jsd:jed) )

    !--- include agrid winds in restarts for use in data assimilation
    !if (Atm%agrid_vel_rst) then
        allocate ( Atm%ua(isd:ied, jsd:jed, nz) )
        allocate ( Atm%va(isd:ied, jsd:jed, nz) )
    !endif

    !--------------------------
    ! Non-hydrostatic dynamics:
    !--------------------------
    !if (.not. Atm%hydrostatic ) then
        allocate (    Atm%w(isd:ied, jsd:jed, nz  ) )
        allocate ( Atm%delz(isd:ied, jsd:jed, nz) )
    !endif

end subroutine allocate_fv_atmos_type

end module fv3jedi_mod
