module fv3jedi_bump_init_mod

use fckit_mpi_module,  only: fckit_mpi_comm
use type_bump,         only: bump_type

use fv3jedi_kinds_mod, only: kind_real

implicit none
private
public bump_init

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine bump_init(f_comm, ngrid_in, lat_in, lon_in, ngrid_ou, lat_ou, lon_ou, bump, bumpid)

implicit none

!Arguments
type(fckit_mpi_comm), intent(in)    :: f_comm
integer,              intent(in)    :: ngrid_in         !Input grid dimensions
real(kind=kind_real), intent(in)    :: lat_in(ngrid_in) !Degrees -90 to 90
real(kind=kind_real), intent(in)    :: lon_in(ngrid_in) !Degrees 0 to 360
integer,              intent(in)    :: ngrid_ou         !Output grid dimensions
real(kind=kind_real), intent(in)    :: lat_ou(ngrid_ou) !Degrees -90 to 90
real(kind=kind_real), intent(in)    :: lon_ou(ngrid_ou) !Degrees 0 to 360
type(bump_type),      intent(inout) :: bump
integer,              intent(in)    :: bumpid

!Locals
real(kind=kind_real), allocatable :: lat_in_us(:,:), lon_in_us(:,:) !Two dimensions
real(kind=kind_real), allocatable :: area(:),vunit(:,:)
logical, allocatable :: lmask(:,:)

character(len=5)    :: cbumpcount
character(len=1024) :: bump_nam_prefix

! Each bump%nam%prefix must be distinct
! -------------------------------------
write(cbumpcount,"(I0.5)") bumpid
bump_nam_prefix = 'fv3jedi_bump_interp_data_'//cbumpcount

! Put latlon into rank 2 unstructured format
! ------------------------------------------
allocate(lat_in_us(ngrid_in,1))
allocate(lon_in_us(ngrid_in,1))

lat_in_us(:,1) = lat_in
lon_in_us(:,1) = lon_in


! Namelist options
! ----------------

!Important namelist options
call bump%nam%init

!Less important namelist options (should not be changed)
bump%nam%prefix = trim(bump_nam_prefix)   ! Prefix for files output
bump%nam%default_seed = .true.
bump%nam%new_obsop = .true.

bump%nam%write_obsop = .false.
bump%nam%verbosity = "none"

! Initialize geometry
! -------------------
allocate(area(ngrid_in))
allocate(vunit(ngrid_in,1))
allocate(lmask(ngrid_in,1))
area = 1.0_kind_real   ! Dummy area
vunit = 1.0_kind_real  ! Dummy vertical unit
lmask = .true.         ! Mask

! Initialize BUMP
! ---------------
call bump%setup_online( f_comm,ngrid_in,1,1,1,lon_in_us,lat_in_us,area,vunit,lmask, &
                        nobs=ngrid_ou,lonobs=lon_ou,latobs=lat_ou)

!Run BUMP drivers
call bump%run_drivers

!Partial deallocate option
call bump%partial_dealloc

! Release memory
! --------------
deallocate(area)
deallocate(vunit)
deallocate(lmask)
deallocate(lat_in_us)
deallocate(lon_in_us)

end subroutine bump_init

! --------------------------------------------------------------------------------------------------

end module fv3jedi_bump_init_mod
