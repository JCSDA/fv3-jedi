! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


module fv3jedi_bump_interp_mod

! atlas
use atlas_module,          only: atlas_field, atlas_fieldset, atlas_functionspace, atlas_real, &
                                 atlas_functionspace_pointcloud

! fckit
use fckit_mpi_module,      only: fckit_mpi_comm

! saber
use type_bump,             only: bump_type

! fv3jedi
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_constants_mod, only: rad2deg

! --------------------------------------------------------------------------------------------------

implicit none
private
public fv3jedi_bump_interp

! --------------------------------------------------------------------------------------------------

type fv3jedi_bump_interp
  type(bump_type) :: bump
  integer :: isc, iec, jsc, jec
  type(atlas_functionspace) :: afunctionspace
  contains
   procedure :: setup
   procedure :: delete
   procedure :: apply
   procedure :: apply_ad
   final     :: dummy_final
end type fv3jedi_bump_interp

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup(self, comm, isc, iec, jsc, jec, npz, lon_in, lat_in, &
                 ngrid_ou, lon_ou_us, lat_ou_us, bumpid)

!Arguments
class(fv3jedi_bump_interp), intent(inout) :: self
type(fckit_mpi_comm),       intent(in)    :: comm
integer,                    intent(in)    :: isc, iec, jsc, jec, npz
real(kind=kind_real),       intent(in)    :: lon_in(isc:iec,jsc:jec)
real(kind=kind_real),       intent(in)    :: lat_in(isc:iec,jsc:jec)
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(in)    :: lon_ou_us(ngrid_ou)
real(kind=kind_real),       intent(in)    :: lat_ou_us(ngrid_ou)
integer, optional,          intent(in)    :: bumpid

!Locals
integer :: ngrid, jnode, jx, jy
real(kind=kind_real), allocatable :: lon_in_us(:)
real(kind=kind_real), allocatable :: lat_in_us(:)
character(len=5)    :: cbumpcount
character(len=1024) :: bump_nam_prefix

real(kind=kind_real), pointer :: real_ptr(:,:)
type(atlas_functionspace) :: afunctionspace
type(atlas_fieldset) :: afieldset
type(atlas_field) :: afield

! Save domain
! -----------
self%isc = isc
self%iec = iec
self%jsc = jsc
self%jec = jec

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
call self%bump%nam%init(comm%size())

self%bump%nam%prefix = trim(bump_nam_prefix)   ! Prefix for files output
self%bump%nam%new_obsop = .true.
self%bump%nam%write_obsop = .false.
self%bump%nam%verbosity = "none"
self%bump%nam%nl = npz+1
self%bump%nam%nv = 1
self%bump%nam%variables(1) = "var"

! Pack lonlat into atlas function space
! -------------------------------------
ngrid = (iec-isc+1)*(jec-jsc+1)
allocate(lon_in_us(ngrid))
allocate(lat_in_us(ngrid))
jnode = 0
do jy = jsc,jec
  do jx = isc,iec
    jnode = jnode+1
    lon_in_us(jnode) = lon_in(jx,jy)*rad2deg
    lat_in_us(jnode) = lat_in(jx,jy)*rad2deg
  end do
end do
afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2, ngrid/))
call afield%data(real_ptr)
real_ptr(1,:) = lon_in_us
real_ptr(2,:) = lat_in_us
self%afunctionspace = atlas_functionspace_pointcloud(afield)
call afield%final()

! Initialize BUMP
! -------------------
call self%bump%setup( comm, self%afunctionspace, nobs=ngrid_ou, lonobs=lon_ou_us, latobs=lat_ou_us)
call self%bump%run_drivers()

! Release memory
! --------------
call self%bump%partial_dealloc()
call afunctionspace%final()

end subroutine setup

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_bump_interp), intent(inout) :: self

call self%afunctionspace%final()
call self%bump%dealloc()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine apply( self, npz, field_in, ngrid_ou, field_ou)

implicit none

!Arguments
class(fv3jedi_bump_interp), intent(inout) :: self
integer,                    intent(in)    :: npz
real(kind=kind_real),       intent(in)    :: field_in(self%isc:self%iec,self%jsc:self%jec,npz)
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(inout) :: field_ou(ngrid_ou,npz)

!Locals
integer :: jl
real(kind_real), pointer :: real_ptr_2(:,:)
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset

! Set number of levels
self%bump%geom%nl0 = npz

! Define ATLAS fieldset
afieldset = atlas_fieldset()
afield = self%afunctionspace%create_field(name='var', kind=atlas_real(kind_real), levels=npz)
call afieldset%add(afield)

! Put input field into ATLAS fieldset
call afield%data(real_ptr_2)
do jl=1,npz
  real_ptr_2(jl,:) = pack(field_in(self%isc:self%iec,self%jsc:self%jec,jl),.true.)
enddo

! Apply BUMP interpolation
call self%bump%apply_obsop(afieldset,field_ou)

! Release memory
call afieldset%final()
call afield%final()

end subroutine apply

! --------------------------------------------------------------------------------------------------

subroutine apply_ad( self, npz, field_in, ngrid_ou, field_ou )

implicit none

!Arguments
class(fv3jedi_bump_interp), intent(inout) :: self
integer,                    intent(in)    :: npz
real(kind=kind_real),       intent(inout) :: field_in(self%isc:self%iec,self%jsc:self%jec,npz)
integer,                    intent(in)    :: ngrid_ou
real(kind=kind_real),       intent(in)    :: field_ou(ngrid_ou,npz)

!Locals
integer :: jl
real(kind_real), pointer :: real_ptr_2(:,:)
logical :: umask(self%isc:self%iec,self%jsc:self%jec)
type(atlas_field) :: afield
type(atlas_fieldset) :: afieldset

! Set number of levels
self%bump%geom%nl0 = npz

! Define ATLAS fieldset
afieldset = atlas_fieldset()
afield = self%afunctionspace%create_field(name='var', kind=atlas_real(kind_real), levels=npz)
call afieldset%add(afield)

! Apply BUMP interpolation adjoint
call self%bump%apply_obsop_ad(field_ou,afieldset)

! Get output field from ATLAS fieldset
call afield%data(real_ptr_2)
umask = .true.
do jl=1,npz
  field_in(self%isc:self%iec,self%jsc:self%jec,jl) = unpack(real_ptr_2(jl,:),umask, &
 & field_in(self%isc:self%iec,self%jsc:self%jec,jl))
enddo

! Release memory
call afieldset%final()
call afield%final()

end subroutine apply_ad

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
type(fv3jedi_bump_interp), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

end module fv3jedi_bump_interp_mod
