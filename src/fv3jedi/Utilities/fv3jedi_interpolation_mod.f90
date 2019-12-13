! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_interpolation_mod

use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_field_mod,     only: fv3jedi_field, get_field
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_constants_mod, only: rad2deg
use wind_vt_mod,           only: d2a, a2d

use type_bump, only: bump_type

use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_max
use fckit_geometry_module, only: sphere_distance
use fckit_kdtree_module, only: kdtree,kdtree_create,kdtree_destroy,kdtree_k_nearest_neighbors

implicit none
private

public :: field2field_interp

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine field2field_interp(nf, geom_in, fields_in, geom_ou, fields_ou)

  implicit none

  !Arguments
  integer,                    intent(in)    :: nf             !Number of fields
  type(fv3jedi_geom), target, intent(inout) :: geom_in        !Geometry of input grid
  type(fv3jedi_field),        intent(in)    :: fields_in(nf)  !Input fields
  type(fv3jedi_geom), target, intent(inout) :: geom_ou        !Geometry of output grid
  type(fv3jedi_field),        intent(inout) :: fields_ou(nf)  !Output fields

  !Locals
  integer :: var
  logical :: dgrid_winds_done
  real(kind=kind_real), allocatable :: ua_in(:,:,:)
  real(kind=kind_real), allocatable :: va_in(:,:,:)
  real(kind=kind_real), allocatable :: ua_ou(:,:,:)
  real(kind=kind_real), allocatable :: va_ou(:,:,:)
  type(fv3jedi_field), pointer :: ud_in
  type(fv3jedi_field), pointer :: vd_in
  type(fv3jedi_field), pointer :: ud_ou
  type(fv3jedi_field), pointer :: vd_ou
  character(len=10) :: fname

  ! Interpolation objects bump
  logical  :: bump_alloc, bary_alloc
  type(bump_type)  :: bump

  ! Interpolation objects barycenric
  integer :: i, j, k, n, nnearest = 4, ngrid_ou, maxtype, mintype, disp
  real(kind=kind_real), allocatable :: lats_ou(:), lons_ou(:)
  real(kind=kind_real), allocatable :: interp_w(:,:)
  integer, allocatable :: interp_i(:,:)
  real(kind=kind_real), allocatable :: field_neighbours(:,:)
  real(kind=kind_real), allocatable :: field_ou(:), field_types(:)


  ! Intepolation weights for regular fields
  ! ---------------------------------------
  bump_alloc = .false.
  do var = 1,nf
    if (.not. fields_in(var)%integerfield .and. .not. bump_alloc) then
      call bump_init(geom_in%f_comm, &
                     geom_in%isc, geom_in%iec, geom_in%jsc, geom_in%jec, &
                     rad2deg*geom_in%grid_lat(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec), &
                     rad2deg*geom_in%grid_lon(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec), &
                     geom_ou%isc, geom_ou%iec, geom_ou%jsc, geom_ou%jec, &
                     rad2deg*geom_ou%grid_lat(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec), &
                     rad2deg*geom_ou%grid_lon(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec), &
                     bump, 99999)
      bump_alloc = .true.

    endif
  enddo

  ! Intepolation weights and factors for integer fields
  ! ---------------------------------------------------
  bary_alloc = .false.
  do var = 1,nf
    if (fields_in(var)%integerfield .and. .not. bary_alloc) then

      ngrid_ou = (geom_ou%iec-geom_ou%isc+1)*(geom_ou%jec-geom_ou%jsc+1)
      allocate(lats_ou(ngrid_ou))
      allocate(lons_ou(ngrid_ou))
      allocate(interp_w(ngrid_ou, nnearest))
      allocate(interp_i(ngrid_ou, nnearest))

      ! Unstructured lat/lon
      n = 0
      do j = geom_ou%jsc,geom_ou%jec
        do i = geom_ou%isc,geom_ou%iec
          n = n + 1
          lats_ou(n) = rad2deg*geom_ou%grid_lat(i,j)
          lons_ou(n) = rad2deg*geom_ou%grid_lon(i,j)
        enddo
      enddo

      call barycentric_init(geom_in, ngrid_ou, lats_ou, lons_ou, nnearest, interp_w, interp_i)

      deallocate(lats_ou,lons_ou)

      bary_alloc = .true.

    endif
  enddo


  ! Special case of D-grid winds
  ! ----------------------------

  ! Interpolation is only for A-Grid winds, otherwise the tangential properties of the wind
  ! vectors may not be properly maintained. Need a D to D grid interpolation

  dgrid_winds_done = .false.

  do var = 1,nf

    if (.not. dgrid_winds_done) then

      if (fields_in(var)%fv3jedi_name == 'ud' .or. fields_in(var)%fv3jedi_name == 'vd') then

        fname = 'ud'
        call get_field(nf,fields_in,fname,ud_in)
        call get_field(nf,fields_ou,fname,ud_ou)
        fname = 'vd'
        call get_field(nf,fields_in,fname,vd_in)
        call get_field(nf,fields_ou,fname,vd_ou)

        allocate(ua_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,ud_in%npz))
        allocate(va_in(geom_in%isc:geom_in%iec,geom_in%jsc:geom_in%jec,vd_in%npz))
        allocate(ua_ou(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,ud_ou%npz))
        allocate(va_ou(geom_ou%isc:geom_ou%iec,geom_ou%jsc:geom_ou%jec,vd_ou%npz))

        call d2a(geom_in, ud_in%array, vd_in%array, ua_in, va_in)

        call bump_apply( geom_in%isc, geom_in%iec, geom_in%jsc, geom_in%jec, ud_in%npz, &
                         ua_in, &
                         geom_ou%isc, geom_ou%iec, geom_ou%jsc, geom_ou%jec, ud_ou%npz, &
                         ua_ou, &
                         bump )
        call bump_apply( geom_in%isc, geom_in%iec, geom_in%jsc, geom_in%jec, vd_in%npz, &
                         va_in, &
                         geom_ou%isc, geom_ou%iec, geom_ou%jsc, geom_ou%jec, vd_ou%npz, &
                         va_ou, &
                         bump )

        call a2d(geom_ou, ua_ou, va_ou, ud_ou%array, vd_ou%array)

        deallocate(ua_in)
        deallocate(va_in)
        deallocate(ua_ou)
        deallocate(va_ou)
        nullify(ud_in)
        nullify(vd_in)
        nullify(ud_ou)
        nullify(vd_ou)

        dgrid_winds_done = .true.

      endif

    endif

  enddo


  !Interpolate all states to this resolution
  do var = 1,nf

    if (fields_in(var)%fv3jedi_name .ne. 'ud' .and. fields_in(var)%fv3jedi_name .ne. 'vd') then

      if (.not. fields_in(var)%integerfield ) then

        call bump_apply( geom_in%isc, geom_in%iec, geom_in%jsc, geom_in%jec, fields_in(var)%npz, &
                         fields_in(var)%array, &
                         geom_ou%isc, geom_ou%iec, geom_ou%jsc, geom_ou%jec, fields_ou(var)%npz, &
                         fields_ou(var)%array, &
                         bump )

      else

        ! Integer interpolation

        if (.not.allocated(field_neighbours)) allocate(field_neighbours(ngrid_ou,nnearest))
        if (.not.allocated(field_ou)) allocate(field_ou(ngrid_ou))

        ! Get field maximum and create array
        if (.not.allocated(field_types)) then
          call geom_in%f_comm%allreduce(int(maxval(fields_in(var)%array)),maxtype,fckit_mpi_max())
          call geom_in%f_comm%allreduce(int(minval(fields_in(var)%array)),mintype,fckit_mpi_max())
          allocate(field_types(mintype:maxtype))
          disp = mintype - 1
        endif

        ! Get neighbours
        do k = 1, fields_ou(var)%npz

          call barycentric_getneighbours(geom_in, ngrid_ou, nnearest, interp_i, fields_in(var)%array(:,:,k), field_neighbours)

          do i = 1,ngrid_ou
            field_types = 0.0
            do n = 1,nnearest
              field_types(int(field_neighbours(i,n))) = field_types(int(field_neighbours(i,n))) + interp_w(i,n)
            enddo
            field_ou(i) = real(maxloc(field_types,1)+disp,kind_real)
          enddo

          n = 0
          do j = geom_ou%jsc,geom_ou%jec
            do i = geom_ou%isc,geom_ou%iec
              n = n + 1
              fields_ou(var)%array(i,j,k) = field_ou(n)
            enddo
          enddo

        enddo

      endif

    endif

  enddo

  ! Deallocate (bump)
  if (bump_alloc) call bump%dealloc()

  ! Deallocate (barycentric)
  if (allocated(interp_w)) deallocate(interp_w)
  if (allocated(interp_i)) deallocate(interp_i)

  if (allocated(field_neighbours)) deallocate(field_neighbours)
  if (allocated(field_ou)) deallocate(field_ou)
  if (allocated(field_types)) deallocate(field_types)

end subroutine field2field_interp

! ------------------------------------------------------------------------------

subroutine bump_init(f_comm, isc_in, iec_in, jsc_in, jec_in, lat_in, lon_in, &
                     isc_ou, iec_ou, jsc_ou, jec_ou, lat_ou, lon_ou, bump, bumpid)

  implicit none

  !Arguments
  type(fckit_mpi_comm), intent(in)    :: f_comm
  integer,              intent(in)    :: isc_in, iec_in, jsc_in, jec_in      !Input grid dimensions
  real(kind=kind_real), intent(in)    :: lat_in(isc_in:iec_in,jsc_in:jec_in) !Degrees -90 to 90
  real(kind=kind_real), intent(in)    :: lon_in(isc_in:iec_in,jsc_in:jec_in) !Degrees 0 to 360
  integer,              intent(in)    :: isc_ou, iec_ou, jsc_ou, jec_ou      !Output grid dimensions
  real(kind=kind_real), intent(in)    :: lat_ou(isc_ou:iec_ou,jsc_ou:jec_ou) !Degrees -90 to 90
  real(kind=kind_real), intent(in)    :: lon_ou(isc_ou:iec_ou,jsc_ou:jec_ou) !Degrees 0 to 360
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

end subroutine bump_init

! ------------------------------------------------------------------------------

subroutine bump_apply( isc_in, iec_in, jsc_in, jec_in, npz_in, field_in, &
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

end subroutine bump_apply

!----------------------------------------------------------------------------

subroutine barycentric_init(geom, nobs, lats_ob, lons_ob, nn, interp_w, interp_i)

implicit none

!Arguments
type(fv3jedi_geom),   intent(in)  :: geom               !Model geometry
integer,              intent(in)  :: nobs               !Number of obs on this processor
real(kind=kind_real), intent(in)  :: lats_ob(nobs)      !Observation locations, lats
real(kind=kind_real), intent(in)  :: lons_ob(nobs)      !Observation locations, lons
integer,              intent(in)  :: nn                 !Number of neighbours to get back
real(kind=kind_real), intent(out) :: interp_w(nobs,nn)  !Interpolation weights
integer,              intent(out) :: interp_i(nobs,nn)  !Interpolation indices (global unstructured grid)

!Locals
type(kdtree) :: kd
integer :: i, j, ob, jj, kk, ngrid_loc, ngrid_all
integer, allocatable :: displs(:), rcvcnt(:), nn_index(:)
real(kind=kind_real) :: wprod, dist, bw(nn), bsw
real(kind=kind_real), allocatable :: lat_loc(:), lon_loc(:), lat_all(:), lon_all(:)
real(kind=kind_real), allocatable :: nn_dist(:)


! Gather the model grid to all processors
! ---------------------------------------

! Allocate arrays for unstructured lat/lon
ngrid_loc = (geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)
ngrid_all = 6 * (geom%npx-1) * (geom%npy-1)

allocate(lat_loc(ngrid_loc))
allocate(lon_loc(ngrid_loc))
allocate(lat_all(ngrid_all))
allocate(lon_all(ngrid_all))

! Fill local unstructured lat/lon
jj = 0
do j = geom%jsc,geom%jec
  do i = geom%isc,geom%iec
     jj = jj + 1
     lat_loc(jj) = rad2deg*geom%grid_lat(i,j)
     lon_loc(jj) = rad2deg*geom%grid_lon(i,j)
  enddo
enddo

! Gather receive count from each processor
allocate(rcvcnt(geom%f_comm%size()))
call geom%f_comm%allgather(jj, rcvcnt)

! Diplacement for each processor
allocate(displs(geom%f_comm%size()))
displs(1) = 0
do j = 2,geom%f_comm%size()
   displs(j) = displs(j-1) + rcvcnt(j-1)
enddo

! Gather the lat/lons to all
call geom%f_comm%allgather(lat_loc,lat_all,jj,rcvcnt,displs)
call geom%f_comm%allgather(lon_loc,lon_all,jj,rcvcnt,displs)

! Deallocate
deallocate(rcvcnt,displs)
deallocate(lat_loc,lon_loc)


! Create a KDTree for finding nearest neighbours
! ----------------------------------------------

! Create kdtree
kd = kdtree_create(ngrid_all,lon_all,lat_all)

! Loop over observations calling kdtree to generate barycentric interpolation weights
! -----------------------------------------------------------------------------------

allocate(nn_index(nn))
allocate(nn_dist(nn))

do ob = 1,nobs

  nn_index = 0
  nn_dist = 0.0_kind_real

  call kdtree_k_nearest_neighbors(kd,lons_ob(ob),lats_ob(ob),nn,nn_index)
  do kk=1,nn
    nn_dist(kk) = sphere_distance(lons_ob(ob),lats_ob(ob),lon_all(nn_index(kk)),lat_all(nn_index(kk)))
  enddo

  !if (geom%f_comm%rank()==0) then
    !!Check the kdtree returned close neighbours
    !print*, '==========================='
    !print*, lats_ob(ob)
    !print*, lat_all(nn_index)
    !print*, ' '
    !print*, lons_ob(ob)
    !print*, lon_all(nn_index)
    !print*, '==========================='
  !endif

  !Barycentric weights formula
  bw(:) = 0.0_kind_real
  do jj = 1,nn
    wprod = 1.0_kind_real
    do kk = 1,nn
      if (jj.ne.kk) then
        dist = sphere_distance(lon_all(nn_index(jj)),lat_all(nn_index(jj)),&
                               lon_all(nn_index(kk)),lat_all(nn_index(kk)))
        wprod = wprod * dist
      endif
    enddo
    bw(jj) = 1.0_kind_real / wprod
  enddo

  !Barycentric weights
  bsw = 0.0_kind_real
  do jj = 1,nn
    bsw = bsw + (bw(jj) / nn_dist(jj))
  enddo

  interp_w(ob,:) = 0.0_kind_real
  do jj = 1,nn
    interp_w(ob,jj) = ( bw(jj) / nn_dist(jj) ) / bsw
    interp_i(ob,jj) = nn_index(jj)
  enddo

enddo

!Deallocate
call kdtree_destroy(kd)
deallocate(nn_index,nn_dist)
deallocate(lon_all,lat_all)

end subroutine barycentric_init

!----------------------------------------------------------------------------

subroutine barycentric_getneighbours(geom, nobs, nn, interp_i, field_in, field_out)

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom               !Model geometry
integer             , intent(in)  :: nobs               !Number of obs on this processor
integer             , intent(in)  :: nn                 !Number of neighbours
integer             , intent(in)  :: interp_i(nobs,nn)  !Interpolation index
real(kind=kind_real), intent(in)  :: field_in(geom%isc:geom%iec,geom%jsc:geom%jec) !Fields in
real(kind=kind_real), intent(out) :: field_out(nobs,nn)                            !Field nearest neighbours

!Locals
integer :: i, j, jj, ob, ngrid_loc, ngrid_all
integer, allocatable :: displs(:), rcvcnt(:)
real(kind=kind_real), allocatable :: field_loc(:), field_all(:)


! Gather the field to all processors
! ---------------------------------------

! Allocate arrays for unstructured fields
ngrid_loc = (geom%iec-geom%isc+1)*(geom%jec-geom%jsc+1)
ngrid_all = 6 * (geom%npx-1) * (geom%npy-1)

allocate(field_loc(ngrid_loc))
allocate(field_all(ngrid_all))

! Fill local unstructured field
jj = 0
do j = geom%jsc,geom%jec
  do i = geom%isc,geom%iec
     jj = jj + 1
     field_loc(jj) = field_in(i,j)
  enddo
enddo

! Gather receive count from each processor
allocate(rcvcnt(geom%f_comm%size()))
call geom%f_comm%allgather(jj, rcvcnt)

! Diplacement for each processor
allocate(displs(geom%f_comm%size()))
displs(1) = 0
do j = 2,geom%f_comm%size()
   displs(j) = displs(j-1) + rcvcnt(j-1)
enddo

! Gather the lat/lons to all
call geom%f_comm%allgather(field_loc,field_all,jj,rcvcnt,displs)

! Deallocate
deallocate(rcvcnt,displs)
deallocate(field_loc)


! Set Field out using indices from the kdtree
! -------------------------------------------

do ob = 1,nobs
   field_out(ob,1) = field_all(interp_i(ob,1))
   field_out(ob,2) = field_all(interp_i(ob,2))
   field_out(ob,3) = field_all(interp_i(ob,3))
   field_out(ob,4) = field_all(interp_i(ob,4))
enddo

end subroutine barycentric_getneighbours

! ------------------------------------------------------------------------------

end module fv3jedi_interpolation_mod
