module fv3jedi_bump_mod

use atlas_module, only: atlas_field, atlas_fieldset, atlas_functionspace, atlas_real

use fckit_mpi_module,  only: fckit_mpi_comm
use type_bump,         only: bump_type

use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_geom_mod,      only: fv3jedi_geom

implicit none
private
public bump_init, bump_apply, bump_apply_ad

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine bump_init(geom_in, ngrid_ou, lat_ou_us, lon_ou_us, bump, bumpid)

implicit none

!Arguments
type(fv3jedi_geom), target, intent(in)    :: geom_in
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(in)    :: lat_ou_us(ngrid_ou)
real(kind=kind_real),       intent(in)    :: lon_ou_us(ngrid_ou)
type(bump_type),            intent(inout) :: bump
integer, optional,          intent(in)    :: bumpid

!Locals
character(len=5)    :: cbumpcount
character(len=1024) :: bump_nam_prefix
type(atlas_fieldset) :: afieldset

! Each bump%nam%prefix must be distinct
! -------------------------------------
if (present(bumpid)) then
  write(cbumpcount,"(I0.5)") bumpid
else
  write(cbumpcount,"(I0.5)") 99999
endif
bump_nam_prefix = 'fv3jedi_bump_interp_data_'//cbumpcount

! Namelist options
! ----------------
call bump%nam%init(geom_in%f_comm%size())

bump%nam%prefix = trim(bump_nam_prefix)   ! Prefix for files output
bump%nam%default_seed = .true.
bump%nam%new_obsop = .true.

bump%nam%write_obsop = .false.
bump%nam%verbosity = "none"

bump%nam%nl = geom_in%npz+1
bump%nam%nv = 1
bump%nam%varname(1) = "var"
bump%nam%nts = 1
bump%nam%timeslot(1) = "0"

!Empty atlas field set
afieldset = atlas_fieldset()

! Initialize BUMP
! -------------------
call bump%setup( geom_in%f_comm,geom_in%afunctionspace,afieldset,nobs=ngrid_ou,lonobs=lon_ou_us,latobs=lat_ou_us)

!Run BUMP drivers
call bump%run_drivers

!Partial deallocate option
call bump%partial_dealloc

! Release memory
! --------------
call afieldset%final()

end subroutine bump_init

! --------------------------------------------------------------------------------------------------

subroutine bump_apply( npz, geom_in, field_in, ngrid_ou, field_ou, bump )

implicit none

!Arguments
integer,                    intent(in)    :: npz
type(fv3jedi_geom), target, intent(in)    :: geom_in
real(kind=kind_real),       intent(in)    :: field_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,npz)
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(inout) :: field_ou(ngrid_ou,npz)
type(bump_type),            intent(inout) :: bump

!Locals
integer :: jl
real(kind_real), pointer :: real_ptr_2(:,:)
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset

! Set number of levels
bump%geom%nl0 = npz

! Define ATLAS fieldset
afieldset = atlas_fieldset()
afield = geom_in%afunctionspace%create_field(name='var_0', kind=atlas_real(kind_real), levels=npz)
call afieldset%add(afield)

! Put input field into ATLAS fieldset
call afield%data(real_ptr_2)
do jl=1,npz
  real_ptr_2(jl,:) = pack(field_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,jl),.true.)
enddo

! Apply BUMP interpolation
call bump%apply_obsop(afieldset,field_ou)

! Release memory
call afieldset%final()
call afield%final()

end subroutine bump_apply

! --------------------------------------------------------------------------------------------------

subroutine bump_apply_ad( npz, geom_in, field_in, ngrid_ou, field_ou, bump )

implicit none

!Arguments
integer,                    intent(in)    :: npz
type(fv3jedi_geom), target, intent(in)    :: geom_in
real(kind=kind_real),       intent(inout) :: field_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,npz)
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(in)    :: field_ou(ngrid_ou,npz)
type(bump_type),            intent(inout) :: bump

!Locals
integer :: jl
real(kind_real), pointer :: real_ptr_2(:,:)
logical :: umask(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec)
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset

! Set number of levels
bump%geom%nl0 = npz

! Define ATLAS fieldset
afieldset = atlas_fieldset()
afield = geom_in%afunctionspace%create_field(name='var_0', kind=atlas_real(kind_real), levels=npz)
call afieldset%add(afield)

! Apply BUMP interpolation adjoint
call bump%apply_obsop_ad(field_ou,afieldset)

! Get output field from ATLAS fieldset
call afield%data(real_ptr_2)
umask = .true.
do jl=1,npz
  field_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,jl) = unpack(real_ptr_2(jl,:),umask, &
 & field_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,jl))
enddo

! Release memory
call afieldset%final()
call afield%final()

end subroutine bump_apply_ad

! --------------------------------------------------------------------------------------------------

end module fv3jedi_bump_mod
