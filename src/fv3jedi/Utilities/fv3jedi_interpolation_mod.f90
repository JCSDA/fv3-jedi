! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_interpolation_mod

use type_bump, only: bump_type

use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_field_mod,     only: fv3jedi_field
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_constants_mod, only: rad2deg

use mpp_domains_mod,  only: east, north, center

implicit none
private

public :: bilinear_bump_interp

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine bilinear_bump_interp(nf, geom_in, fields_in, geom_ou, fields_ou)

  implicit none

  !Arguments
  integer,                    intent(in)    :: nf             !Number of fields
  type(fv3jedi_geom), target, intent(in)    :: geom_in        !Geometry of input grid
  type(fv3jedi_field),        intent(in)    :: fields_in(nf)  !Input fields
  type(fv3jedi_geom), target, intent(in)    :: geom_ou        !Geometry of output grid
  type(fv3jedi_field),        intent(inout) :: fields_ou(nf)  !Output fields

  !Locals
  integer :: var
  integer :: isc_in, iec_in, jsc_in, jec_in, npz_in
  integer :: isc_ou, iec_ou, jsc_ou, jec_ou, npz_ou
  real(kind=kind_real), pointer :: lat_in(:,:), lon_in(:,:)
  real(kind=kind_real), pointer :: lat_ou(:,:), lon_ou(:,:)
  logical, target  :: bump_c_alloc, bump_n_alloc, bump_e_alloc
  logical, pointer :: bump_alloc
  type(bump_type), target  :: bump_c, bump_n, bump_e
  type(bump_type), pointer :: bump


  bump_c_alloc = .false.
  bump_n_alloc = .false.
  bump_e_alloc = .false.

  !Interpolate all states to this resolution
  do var = 1,nf

    isc_in = fields_in(var)%isc
    iec_in = fields_in(var)%iec
    jsc_in = fields_in(var)%jsc
    jec_in = fields_in(var)%jec
    npz_in = fields_in(var)%npz

    isc_ou = fields_ou(var)%isc
    iec_ou = fields_ou(var)%iec
    jsc_ou = fields_ou(var)%jsc
    jec_ou = fields_ou(var)%jec
    npz_ou = fields_ou(var)%npz

    if (.not.associated(lat_in)) allocate(lat_in(isc_in:iec_in,jsc_in:jec_in))
    if (.not.associated(lon_in)) allocate(lon_in(isc_in:iec_in,jsc_in:jec_in))
    if (.not.associated(lat_ou)) allocate(lat_ou(isc_ou:iec_ou,jsc_ou:jec_ou))
    if (.not.associated(lon_ou)) allocate(lon_ou(isc_ou:iec_ou,jsc_ou:jec_ou))

    ! Set pointers
    if (fields_ou(var)%staggerloc == center) then

      bump_alloc => bump_c_alloc
      bump => bump_c

      if (.not. bump_alloc) then
        lat_in => geom_in%grid_lat(isc_in:iec_in,jsc_in:jec_in)
        lon_in => geom_in%grid_lon(isc_in:iec_in,jsc_in:jec_in)
        lat_ou => geom_ou%grid_lat(isc_ou:iec_ou,jsc_ou:jec_ou)
        lon_ou => geom_ou%grid_lon(isc_ou:iec_ou,jsc_ou:jec_ou)
      endif

    elseif (fields_ou(var)%staggerloc == north) then

      bump_alloc => bump_n_alloc
      bump => bump_n

      if (.not. bump_alloc) then
        lat_in => geom_in%egrid_lat(isc_in:iec_in,jsc_in:jec_in)
        lon_in => geom_in%grid_lon(isc_in:iec_in,jsc_in:jec_in)
        lat_ou => geom_ou%egrid_lat(isc_ou:iec_ou,jsc_ou:jec_ou)
        lon_ou => geom_ou%grid_lon(isc_ou:iec_ou,jsc_ou:jec_ou)
      endif

    elseif (fields_ou(var)%staggerloc == east) then

      bump_alloc => bump_e_alloc
      bump => bump_e

      if (.not. bump_alloc) then
        lat_in => geom_in%grid_lat(isc_in:iec_in,jsc_in:jec_in)
        lon_in => geom_in%egrid_lon(isc_in:iec_in,jsc_in:jec_in)
        lat_ou => geom_ou%grid_lat(isc_ou:iec_ou,jsc_ou:jec_ou)
        lon_ou => geom_ou%egrid_lon(isc_ou:iec_ou,jsc_ou:jec_ou)
      endif

    endif

    if (.not. bump_alloc) then
      call bilinear_bump_init(isc_in, iec_in, jsc_in, jec_in, rad2deg*lat_in, rad2deg*lon_in-180.0_kind_real, &
                              isc_ou, iec_ou, jsc_ou, jec_ou, rad2deg*lat_ou, rad2deg*lon_ou-180.0_kind_real, bump, 99999)
      bump_alloc = .true.
    endif


    call bilinear_bump_apply( isc_in, iec_in, jsc_in, jec_in, npz_in, fields_in(var)%array, &
                              isc_ou, iec_ou, jsc_ou, jec_ou, npz_ou, fields_ou(var)%array, bump )

  enddo

  nullify(lat_in,lon_in,lat_ou,lon_ou)
  if (bump_c_alloc) call bump_c%dealloc()
  if (bump_n_alloc) call bump_n%dealloc()
  if (bump_e_alloc) call bump_e%dealloc()

end subroutine bilinear_bump_interp

! ------------------------------------------------------------------------------

subroutine bilinear_bump_init(isc_in, iec_in, jsc_in, jec_in, lat_in, lon_in, &
                              isc_ou, iec_ou, jsc_ou, jec_ou, lat_ou, lon_ou, bump, bumpid)

  implicit none

  !Arguments
  integer,              intent(in)    :: isc_in, iec_in, jsc_in, jec_in      !Input grid dimensions
  real(kind=kind_real), intent(in)    :: lat_in(isc_in:iec_in,jsc_in:jec_in) !Degrees -90 to 90
  real(kind=kind_real), intent(in)    :: lon_in(isc_in:iec_in,jsc_in:jec_in) !Degrees -180 to 180
  integer,              intent(in)    :: isc_ou, iec_ou, jsc_ou, jec_ou      !Output grid dimensions
  real(kind=kind_real), intent(in)    :: lat_ou(isc_ou:iec_ou,jsc_ou:jec_ou) !Degrees -90 to 90
  real(kind=kind_real), intent(in)    :: lon_ou(isc_ou:iec_ou,jsc_ou:jec_ou) !Degrees -180 to 180
  type(bump_type),      intent(inout) :: bump
  integer,              intent(in)    :: bumpid

  !Locals
  integer :: ngrid_in, ngrid_ou
  real(kind=kind_real), allocatable :: lat_in_us(:,:), lon_in_us(:,:) !Unstructured
  real(kind=kind_real), allocatable :: lat_ou_us(:,:), lon_ou_us(:,:) !Unstructured

  real(kind=kind_real), allocatable :: area(:),vunit(:,:)
  logical, allocatable :: lmask(:,:)

  character(len=5)    :: cbumpcount
  character(len=1024) :: bump_nam_prefix

  integer :: ii, ji, jj

  ! Each bump%nam%prefix must be distinct
  ! -------------------------------------
  write(cbumpcount,"(I0.5)") bumpid
  bump_nam_prefix = 'fv3jedi_bumpobsop_data_'//cbumpcount


  ! Put latlon into unstructured format
  ! -----------------------------------
  ngrid_in = (iec_in - isc_in + 1) * (jec_in - jsc_in + 1)
  allocate(lat_in_us(ngrid_in,1))
  allocate(lon_in_us(ngrid_in,1))

  ii = 0
  do jj = jsc_in, jec_in
    do ji = isc_in, iec_in
      ii = ii + 1
      lat_in_us(ii, 1) = lat_in(ji, jj)
      lon_in_us(ii, 1) = lon_in(ji, jj)
    enddo
  enddo

  ngrid_ou = (iec_ou - isc_ou + 1) * (jec_ou - jsc_ou + 1)
  allocate(lat_ou_us(ngrid_ou,1))
  allocate(lon_ou_us(ngrid_ou,1))

  ii = 0
  do jj = jsc_ou, jec_ou
    do ji = isc_ou, iec_ou
      ii = ii + 1
      lat_ou_us(ii, 1) = lat_ou(ji, jj)
      lon_ou_us(ii, 1) = lon_ou(ji, jj)
    enddo
  enddo


  ! Namelist options
  ! ----------------

  !Important namelist options
  call bump%nam%init
  bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)

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
  call bump%setup_online( ngrid_in,1,1,1,lon_in_us,lat_in_us,area,vunit,lmask, &
                          nobs=ngrid_ou,lonobs=lon_ou_us(:,1),latobs=lat_ou_us(:,1))

  !Run BUMP drivers
  call bump%run_drivers

  !Partial deallocate option
  !call bump%partial_dealloc

  ! Release memory
  ! --------------
  deallocate(area)
  deallocate(vunit)
  deallocate(lmask)
  deallocate(lat_in_us,lat_ou_us)
  deallocate(lon_in_us,lon_ou_us)

end subroutine bilinear_bump_init

! ------------------------------------------------------------------------------

subroutine bilinear_bump_apply( isc_in, iec_in, jsc_in, jec_in, npz_in, field_in, &
                                isc_ou, iec_ou, jsc_ou, jec_ou, npz_ou, field_ou, bump )

  implicit none

  !Arguments
  integer,              intent(in)    :: isc_in, iec_in, jsc_in, jec_in, npz_in
  real(kind=kind_real), intent(in)    :: field_in(isc_in:iec_in,jsc_in:jec_in,npz_in)
  integer,              intent(in)    :: isc_ou, iec_ou, jsc_ou, jec_ou, npz_ou
  real(kind=kind_real), intent(inout) :: field_ou(isc_ou:iec_ou,jsc_ou:jec_ou,npz_ou)
  type(bump_type),      intent(inout) :: bump

  integer :: ngrid_in, ngrid_ou
  integer :: ii, ji, jj, jk
  real(kind=kind_real), allocatable :: field_in_us(:,:) !Unstructured
  real(kind=kind_real), allocatable :: field_ou_us(:,:) !Unstructured

  ! Check for matching vertical number of levels
  if (npz_in .ne. npz_ou) &
    call abor1_ftn("fv3jedi_interpolation_mod.apply_bump does not support different vertical levels")


  ! Put field into unstructured format
  ! ----------------------------------
  ngrid_in = (iec_in - isc_in + 1) * (jec_in - jsc_in + 1)
  allocate(field_in_us(ngrid_in,1))

  ngrid_ou = (iec_ou - isc_ou + 1) * (jec_ou - jsc_ou + 1)
  allocate(field_ou_us(ngrid_ou,1))
  field_ou_us = 0.0_kind_real

  do jk = 1,npz_in

    ii = 0
    do jj = jsc_in, jec_in
      do ji = isc_in, iec_in
        ii = ii + 1
        field_in_us(ii, 1) = field_in(ji, jj, jk)
      enddo
    enddo

    call bump%apply_obsop(field_in_us,field_ou_us)

    ii = 0
    do jj = jsc_ou, jec_ou
      do ji = isc_ou, iec_ou
        ii = ii + 1
        field_ou(ji, jj, jk) = field_ou_us(ii, 1)
      enddo
    enddo

  enddo

  deallocate(field_in_us,field_ou_us)

end subroutine bilinear_bump_apply

! ------------------------------------------------------------------------------

end module fv3jedi_interpolation_mod
