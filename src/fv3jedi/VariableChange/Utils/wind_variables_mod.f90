! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wind_vt_mod

! libs
use mpi
use netcdf

use fv3jedi_constants_mod, only: constant
use fv3jedi_geom_mod,  only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_communication_mod, only: gather_field, scatter_field
use fv3jedi_netcdf_utils_mod, only: nccheck

use fv_mp_mod, only: fill_corners
use fv_mp_adm_mod, only: mpp_update_domains_adm
use mpp_domains_mod, only: mpp_update_domains, dgrid_ne
use mpp_domains_mod, only: mpp_get_boundary, mpp_get_boundary_ad
use mpp_parameter_mod, only: CGRID_NE

use femps_fv3_mod, only: fv3field_to_ufield, ufield_to_fv3field
use femps_grid_mod, only: fempsgrid
use femps_operators_mod, only: fempsoprs
use femps_solve_mod, only: laplace, inverselaplace

implicit none
private

public sfc_10m_winds

! Stream function and velocity potential to wind (derivative)
public psichi_to_udvd
public psichi_to_udvd_adm

! Vorticity and divergence to stream function and velocity potential (inverse Laplacian)
public vortdivg_to_psichi

! Steam function and velocity potentail to voriticty and divergence (Laplacian)
public psichi_to_vortdivg

! Wind to vorticity and divergence
public udvd_to_vortdivg

! A to D grid winds (A is lonlat oriented)
public a_to_d
public a_to_d_ad

! D to A grid winds (A is lonlat oriented)
public d_to_a
public d_to_a_ad

! D to C grid winds (A is cubed sphere)
public d_to_a_to_c

! C to D grid winds (A is cubed sphere)
public c_to_a_to_d

contains

! --------------------------------------------------------------------------------------------------

subroutine sfc_10m_winds(geom,usrf,vsrf,f10r,spd10m,dir10m)

 !Arguments
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) ::   usrf(geom%isc:geom%iec,geom%jsc:geom%jec) !Lowest model level u m/s
 real(kind=kind_real), intent(in ) ::   vsrf(geom%isc:geom%iec,geom%jsc:geom%jec) !Lowest model level v m/s
 real(kind=kind_real), intent(in ) ::   f10r(geom%isc:geom%iec,geom%jsc:geom%jec) !Ratio of lowest level to 10m
 real(kind=kind_real), intent(out) :: spd10m(geom%isc:geom%iec,geom%jsc:geom%jec) !10m wind speed u m/s
 real(kind=kind_real), intent(out) :: dir10m(geom%isc:geom%iec,geom%jsc:geom%jec) !10m model wind direction

 !Locals
 integer :: isc,iec,jsc,jec,i,j
 integer :: iquad
 real(kind=kind_real) :: windangle, windratio
 real(kind=kind_real), parameter :: windscale = 999999.0_kind_real
 real(kind=kind_real), parameter :: windlimit = 0.0001_kind_real
 real(kind=kind_real), parameter :: quadcof(4,2) = reshape((/ 0.0_kind_real,  1.0_kind_real, &
                                                              1.0_kind_real,  2.0_kind_real, &
                                                              1.0_kind_real, -1.0_kind_real, &
                                                              1.0_kind_real, -1.0_kind_real /), (/4, 2/))
 real(kind=kind_real) :: pi, rad2deg

 !In GSI these calculations are done after interpolation to obs location

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

! Constants
 pi = constant('pi')
 rad2deg = constant('rad2deg')

 !10m wind speed
 spd10m(isc:iec,jsc:jec) = f10r(isc:iec,jsc:jec)*sqrt( usrf(isc:iec,jsc:jec)*usrf(isc:iec,jsc:jec) + &
                                                      vsrf(isc:iec,jsc:jec)*vsrf(isc:iec,jsc:jec) )

 !10m wind direction
 do j = jsc,jec
   do i = isc,iec

     if (usrf(i,j)*f10r(i,j) >= 0.0_kind_real .and. vsrf(i,j)*f10r(i,j) >= 0.0_kind_real) iquad = 1
     if (usrf(i,j)*f10r(i,j) >= 0.0_kind_real .and. vsrf(i,j)*f10r(i,j) <  0.0_kind_real) iquad = 2
     if (usrf(i,j)*f10r(i,j) <  0.0_kind_real .and. vsrf(i,j)*f10r(i,j) >= 0.0_kind_real) iquad = 3
     if (usrf(i,j)*f10r(i,j) <  0.0_kind_real .and. vsrf(i,j)*f10r(i,j) <  0.0_kind_real) iquad = 4

     if (abs(vsrf(i,j)*f10r(i,j)) >= windlimit) then
        windratio = (usrf(i,j)*f10r(i,j)) / (vsrf(i,j)*f10r(i,j))
     else
        windratio = 0.0_kind_real
        if (abs(usrf(i,j)*f10r(i,j)) > windlimit) then
          windratio = windscale * usrf(i,j)*f10r(i,j)
        endif
     endif

     windangle = atan(abs(windratio))

     dir10m(i,j) = rad2deg*(quadcof(iquad,1) * pi + windangle * quadcof(iquad, 2))

   enddo
 enddo

end subroutine sfc_10m_winds

! --------------------------------------------------------------------------------------------------

subroutine psichi_to_udvd(geom, psi_in, chi_in, ud, vd)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  :: psi_in(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(in)  :: chi_in(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(out) ::     ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(out) ::     vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)

integer :: i,j,k

real(kind=kind_real) :: ud_div(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: vd_div(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: uc_div(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: vc_div(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: ua_div(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: va_div(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)

real(kind=kind_real) ::   psi(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) ::   chi(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)


! Fill halos on psi and chi
! -------------------------
psi = 0.0_kind_real
chi = 0.0_kind_real
psi(geom%isc:geom%iec,geom%jsc:geom%jec,:) = psi_in(geom%isc:geom%iec,geom%jsc:geom%jec,:)
chi(geom%isc:geom%iec,geom%jsc:geom%jec,:) = chi_in(geom%isc:geom%iec,geom%jsc:geom%jec,:)

call mpp_update_domains(psi, geom%domain, complete=.true.)
call mpp_update_domains(chi, geom%domain, complete=.true.)


! Compute C-grid divergent winds
! ------------------------------
do k=1,geom%npz
  do j=geom%jsc,geom%jec
    do i=geom%isc,geom%iec
      uc_div(i,j,k) = (chi(i,j,k) - chi(i-1,j,k))/geom%dxc(i,j)
    enddo
  enddo

  do j=geom%jsc,geom%jec
    do i=geom%isc,geom%iec
      vc_div(i,j,k) = (chi(i,j,k) - chi(i,j-1,k))/geom%dyc(i,j)
    enddo
  enddo
enddo


! Interpolate C-grid divergent winds to D-grid (via cubed-sphere A-Grid)
! ----------------------------------------------------------------------
call fill_cgrid_winds(geom, uc_div, vc_div, fillhalo=.true.)

do k = 1,geom%npz
  call c_to_acs_domain_level(geom, uc_div(:,:,k), vc_div(:,:,k), ua_div(:,:,k), va_div(:,:,k))
enddo

! Refill halos of A-Grid winds
call mpp_update_domains(ua_div, geom%domain, complete=.true.)
call mpp_update_domains(va_div, geom%domain, complete=.true.)

! A -> D
do k = 1,geom%npz
  call acs_to_d_domain_level(geom, ua_div(:,:,k), va_div(:,:,k), ud_div(:,:,k), vd_div(:,:,k))
enddo


! Sum divergent and rotational winds
! ----------------------------------
ud = 0.0_kind_real
do k=1,geom%npz
  do j=geom%jsc,geom%jec
    do i=geom%isc,geom%iec
      ud(i,j,k) = -(psi(i,j,k) - psi(i,j-1,k))/geom%dyc(i,j) + ud_div(i,j,k)
    enddo
  enddo
enddo

vd = 0.0_kind_real
do k=1,geom%npz
  do j=geom%jsc,geom%jec
    do i=geom%isc,geom%iec
      vd(i,j,k) = (psi(i,j,k) - psi(i-1,j,k))/geom%dxc(i,j) + vd_div(i,j,k)
    enddo
  enddo
enddo

end subroutine psichi_to_udvd

! --------------------------------------------------------------------------------------------------

subroutine psichi_to_udvd_adm(geom, psi_in_ad, chi_in_ad, ud_ad, vd_ad)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(out)   :: psi_in_ad(geom%isc:geom%iec  , geom%jsc:geom%jec  , geom%npz)
real(kind=kind_real), intent(out)   :: chi_in_ad(geom%isc:geom%iec  , geom%jsc:geom%jec  , geom%npz)
real(kind=kind_real), intent(inout) ::     ud_ad(geom%isc:geom%iec  , geom%jsc:geom%jec+1, geom%npz)
real(kind=kind_real), intent(inout) ::     vd_ad(geom%isc:geom%iec+1, geom%jsc:geom%jec  , geom%npz)

integer :: i, j, k
real(kind=kind_real) :: ud_div_ad(geom%isd:geom%ied  , geom%jsd:geom%jed+1, geom%npz)
real(kind=kind_real) :: vd_div_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) :: uc_div_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) :: vc_div_ad(geom%isd:geom%ied  , geom%jsd:geom%jed+1, geom%npz)
real(kind=kind_real) :: ua_div_ad(geom%isd:geom%ied  , geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) :: va_div_ad(geom%isd:geom%ied  , geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) ::    psi_ad(geom%isd:geom%ied  , geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) ::    chi_ad(geom%isd:geom%ied  , geom%jsd:geom%jed  , geom%npz)
real(kind=kind_real) :: temp_ad

psi_ad = 0.0_kind_real
vd_div_ad = 0.0_kind_real
do k=geom%npz,1,-1
  do j=geom%jec,geom%jsc,-1
    do i=geom%iec,geom%isc,-1
      temp_ad = vd_ad(i,j,k)/geom%dxc(i,j)
      vd_div_ad(i,j,k) = vd_div_ad(i,j,k) + vd_ad(i,j,k)
      vd_ad(i,j,k) = 0.0_kind_real
      psi_ad(i,j,k) = psi_ad(i,j,k) + temp_ad
      psi_ad(i-1,j,k) = psi_ad(i-1,j,k) - temp_ad
    end do
  end do
end do

ud_div_ad = 0.0_kind_real
do k=geom%npz,1,-1
  do j=geom%jec,geom%jsc,-1
    do i=geom%iec,geom%isc,-1
      ud_div_ad(i,j,k) = ud_div_ad(i,j,k) + ud_ad(i,j,k)
      temp_ad = -(ud_ad(i,j,k)/geom%dyc(i,j))
      ud_ad(i,j,k) = 0.0_kind_real
      psi_ad(i,j,k) = psi_ad(i,j,k) + temp_ad
      psi_ad(i,j-1,k) = psi_ad(i,j-1,k) - temp_ad
    end do
  end do
end do

ua_div_ad = 0.0_kind_real
va_div_ad = 0.0_kind_real
do k=geom%npz,1,-1
  call acs_to_d_domain_level_adm(geom, ua_div_ad(:,:,k), va_div_ad(:,:,k), ud_div_ad(:,:,k), &
                                  vd_div_ad(:,:,k))
end do

call mpp_update_domains_adm(va_div_ad, va_div_ad, geom%domain, complete=.true.)
call mpp_update_domains_adm(va_div_ad, ua_div_ad, geom%domain, complete=.true.)

vc_div_ad = 0.0_kind_real
uc_div_ad = 0.0_kind_real
do k=geom%npz,1,-1
  call c_to_acs_domain_level_adm(geom, uc_div_ad(:,:,k), vc_div_ad(:,:,k), ua_div_ad(:,:,k), &
                                  va_div_ad(:,:,k))
end do
call fill_cgrid_winds_adm(geom, uc_div_ad, vc_div_ad, fillhalo=.true.)

chi_ad = 0.0_kind_real
do k=geom%npz,1,-1
  do j=geom%jec,geom%jsc,-1
    do i=geom%iec,geom%isc,-1
      temp_ad = vc_div_ad(i,j,k)/geom%dyc(i,j)
      vc_div_ad(i,j,k) = 0.0_kind_real
      chi_ad(i,j,k) = chi_ad(i,j,k) + temp_ad
      chi_ad(i,j-1,k) = chi_ad(i,j-1,k) - temp_ad
    end do
  end do
  do j=geom%jec,geom%jsc,-1
    do i=geom%iec,geom%isc,-1
      temp_ad = uc_div_ad(i,j,k)/geom%dxc(i,j)
      uc_div_ad(i,j,k) = 0.0_kind_real
      chi_ad(i,j,k) = chi_ad(i,j,k) + temp_ad
      chi_ad(i-1,j,k) = chi_ad(i-1,j,k) - temp_ad
    end do
  end do
end do

call mpp_update_domains_adm(chi_ad, chi_ad, geom%domain, complete=.true.)
call mpp_update_domains_adm(psi_ad, psi_ad, geom%domain, complete=.true.)

chi_in_ad = 0.0_kind_real
chi_in_ad(geom%isc:geom%iec,geom%jsc:geom%jec,:) = chi_ad(geom%isc:geom%iec,geom%jsc:geom%jec,:)
psi_in_ad = 0.0_kind_real
psi_in_ad(geom%isc:geom%iec,geom%jsc:geom%jec,:) = psi_ad(geom%isc:geom%iec,geom%jsc:geom%jec,:)
ud_ad = 0.0_kind_real
vd_ad = 0.0_kind_real

end subroutine psichi_to_udvd_adm

! --------------------------------------------------------------------------------------------------

subroutine vortdivg_to_psichi(geom, grid, oprs, vort, divg, psi, chi, lsize, lev_start, lev_final)

 type(fv3jedi_geom),   intent(in)  :: geom
 type(fempsgrid),      intent(in)  :: grid
 type(fempsoprs),      intent(in)  :: oprs
 real(kind=kind_real), intent(in)  :: vort(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(in)  :: divg(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(out) ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 real(kind=kind_real), intent(out) ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 integer,              intent(in)  :: lsize
 integer,              intent(in)  :: lev_start(lsize),lev_final(lsize)

 integer :: i,j,k
 real(kind=kind_real), allocatable, dimension(:,:,:) :: vorgcomm, divgcomm
 real(kind=kind_real), allocatable, dimension(:,:,:) :: psigcomm, chigcomm
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: vortg, divgg !Global level of vor and div
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: psig, chig !Global level of psi and chi

 real(kind=kind_real), allocatable, dimension(:) :: voru, divu     !Unstructured vor and div
 real(kind=kind_real), allocatable, dimension(:) :: psiu, chiu     !Unstructured psi and chi

 integer :: ranki, n
 integer, allocatable :: lev_proc(:)

 ! Convergence check
 real(kind=kind_real), allocatable, dimension(:,:) :: psi_convergence, chi_convergence
 integer :: comm, color, ncid, nlid, niid, vid(2), ierr
 integer :: starts(2), counts(2)


 ! Gather voricity and divergence to one processor and compute psi and chi
 ! -----------------------------------------------------------------------

 ranki = geom%f_comm%rank() + 1

 allocate(vorgcomm(1:geom%npx-1,1:geom%npy-1,6))
 allocate(divgcomm(1:geom%npx-1,1:geom%npy-1,6))

 if (ranki <= lsize) then
   allocate(vortg(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
   allocate(divgg(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
   allocate(psig(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
   allocate(chig(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))

   allocate(voru(grid%nface(grid%ngrids)))
   allocate(divu(grid%nface(grid%ngrids)))
   allocate(psiu(grid%nface(grid%ngrids)))
   allocate(chiu(grid%nface(grid%ngrids)))
 endif

 ! Processor that each level is associated with
 allocate(lev_proc(geom%npz))
 do k = 1,geom%npz
   do n = 1,lsize
     if (k >= lev_start(n) .and. k <= lev_final(n)) then
       lev_proc(k) = n
     endif
   enddo
 enddo

 ! Gather field to respective processor
 do k=1,geom%npz
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,vort(:,:,k),vorgcomm)
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,divg(:,:,k),divgcomm)
   if (ranki == lev_proc(k)) then
     vortg(:,:,:,k) = vorgcomm
     divgg(:,:,:,k) = divgcomm
   endif
 enddo

 deallocate(vorgcomm,divgcomm)

 ! Split comm
 if (grid%check_convergence) then
   color = 1
   if (ranki <= lsize) color = 0
   call mpi_comm_split(geom%f_comm%communicator(), color, geom%f_comm%rank(), comm, ierr)
 endif

 ! Loop over level and compute psi/chi
 if (ranki <= lsize) then ! Only processors with a level

   ! Arrays to hold convergence information
   if (grid%check_convergence) then
     allocate(psi_convergence(lev_start(ranki):lev_final(ranki),grid%niter))
     allocate(chi_convergence(lev_start(ranki):lev_final(ranki),grid%niter))
   endif

   do k = lev_start(ranki),lev_final(ranki)

     if (geom%f_comm%rank()==0) &
       print*, "Doing femps level ", k, ". Proc range:",lev_start(ranki),"-",lev_final(ranki)

     call fv3field_to_ufield(grid,geom%npx-1,vortg(:,:,:,k),voru)
     call fv3field_to_ufield(grid,geom%npx-1,divgg(:,:,:,k),divu)

     ! Convert to area integrals, required by femps
     voru = voru*grid%farea(:,grid%ngrids)
     divu = divu*grid%farea(:,grid%ngrids)

     ! Solve poisson equation (\psi=\nabla^{-2}\zeta, \chi=\nabla^{-2}D)
     if (grid%check_convergence) then
       call inverselaplace(grid,oprs,grid%ngrids,voru,psiu,field_conv=psi_convergence(k,:))
       call inverselaplace(grid,oprs,grid%ngrids,divu,chiu,field_conv=chi_convergence(k,:))
     else
       call inverselaplace(grid,oprs,grid%ngrids,voru,psiu)
       call inverselaplace(grid,oprs,grid%ngrids,divu,chiu)
     endif

     ! Convert from area integrals
     psiu = psiu/grid%farea(:,grid%ngrids)
     chiu = chiu/grid%farea(:,grid%ngrids)

     call ufield_to_fv3field(grid,geom%npx-1,psiu,psig(:,:,:,k))
     call ufield_to_fv3field(grid,geom%npx-1,chiu,chig(:,:,:,k))

   enddo

   ! Write files containing convergence information
   ! ----------------------------------------------

   if (grid%check_convergence) then

     ! Create file
     call nccheck( nf90_create( "femps_convergence.nc4", ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                                comm = comm, info = MPI_INFO_NULL), &
                   "nf90_create" )

     ! Define dimensions
     call nccheck( nf90_def_dim(ncid, "nlayers",     geom%npz,   nlid), "nf90_def_dim ngrids"  )
     call nccheck( nf90_def_dim(ncid, "niterations", grid%niter, niid), "nf90_def_dim nfacex"  )

     ! Define variables
     call nccheck( nf90_def_var(ncid, "psi", NF90_DOUBLE, (/ nlid, niid /), vid(1)), "nf90_def_var " )
     call nccheck( nf90_put_att(ncid, vid(1), "long_name", "convergence_psi" ), "nf90_put_att" )
     call nccheck( nf90_def_var(ncid, "chi", NF90_DOUBLE, (/ nlid, niid /), vid(2)), "nf90_def_var " )
     call nccheck( nf90_put_att(ncid, vid(2), "long_name", "convergence_chi" ), "nf90_put_att" )

     ! End define
     call nccheck( nf90_enddef(ncid), "nf90_enddef" )

     ! Write variables
     starts(1) = lev_start(ranki)
     starts(2) = 1
     counts(1) = lev_final(ranki) - lev_start(ranki) + 1
     counts(2) = grid%niter

     call nccheck( nf90_put_var( ncid, vid(1), psi_convergence, start = starts, count = counts), &
                   "nf90_put_var psi_convergence" )
     call nccheck( nf90_put_var( ncid, vid(2), chi_convergence, start = starts, count = counts), &
                   "nf90_put_var chi_convergence" )

     ! Close file
     call nccheck ( nf90_close(ncid), "nf90_close" )

   endif

 endif

 if (grid%check_convergence) call MPI_Comm_free(comm, ierr)

 ! Scatter field from respective processor
 allocate(psigcomm(1:geom%npx-1,1:geom%npy-1,6))
 allocate(chigcomm(1:geom%npx-1,1:geom%npy-1,6))

 do k=1,geom%npz
   if (ranki == lev_proc(k)) then
     psigcomm = psig(:,:,:,k)
     chigcomm = chig(:,:,:,k)
   endif
   call scatter_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,psigcomm,psi(:,:,k))
   call scatter_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,chigcomm,chi(:,:,k))
 enddo

end subroutine vortdivg_to_psichi

! --------------------------------------------------------------------------------------------------

subroutine psichi_to_vortdivg(geom, grid, oprs, psi, chi ,vor ,div, lsize, lev_start, lev_final)

 type(fv3jedi_geom),             intent(in)  :: geom
 type(fempsgrid),                intent(in)  :: grid
 type(fempsoprs),                intent(in)  :: oprs
 real(kind=kind_real),           intent(in)  :: psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real),           intent(in)  :: chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 real(kind=kind_real), optional, intent(out) :: vor(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Vorticity
 real(kind=kind_real), optional, intent(out) :: div(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Divergence
 integer,                        intent(in)  :: lsize
 integer,                        intent(in)  :: lev_start(lsize),lev_final(lsize)

 integer :: k
 real(kind=kind_real), allocatable, dimension(:,:,:) :: vorgcomm, divgcomm
 real(kind=kind_real), allocatable, dimension(:,:,:) :: psigcomm, chigcomm
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: vorg, divg !Global level of vor and div
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: psig, chig !Global level of psi and chi
 real(kind=kind_real), allocatable, dimension(:) :: voru, divu     !Unstructured vor and div
 real(kind=kind_real), allocatable, dimension(:) :: psiu, chiu     !Unstructured psi and chi

 integer :: ranki, n
 integer, allocatable :: lev_proc(:)

 ! Gather voricity and divergence to one processor and compute psi and chi
 ! -----------------------------------------------------------------------

 ranki = geom%f_comm%rank() + 1

  if (ranki <= lsize) then
    allocate(vorg(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
    allocate(divg(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
    allocate(psig(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
    allocate(chig(1:geom%npx-1,1:geom%npy-1,6,lev_start(ranki):lev_final(ranki)))
    allocate(voru(grid%nface(grid%ngrids)))
    allocate(divu(grid%nface(grid%ngrids)))
    allocate(psiu(grid%nface(grid%ngrids)))
    allocate(chiu(grid%nface(grid%ngrids)))
  endif

 ! Processor that each level is associated with
 allocate(lev_proc(geom%npz))
 do k = 1,geom%npz
   do n = 1,lsize
     if (k >= lev_start(n) .and. k <= lev_final(n)) then
       lev_proc(k) = n
     endif
   enddo
 enddo

 allocate(psigcomm(1:geom%npx-1,1:geom%npy-1,6))
 allocate(chigcomm(1:geom%npx-1,1:geom%npy-1,6))
 do k=1,geom%npz
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,psi(:,:,k),psigcomm)
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,chi(:,:,k),chigcomm)
   if (ranki == lev_proc(k)) then
     psig(:,:,:,k) = psigcomm
     chig(:,:,:,k) = chigcomm
   endif
 enddo
 deallocate(psigcomm,chigcomm)

 ! Loop over level and compute vort/divg
 if (ranki <= lsize) then ! Only processors with a level
   do k = lev_start(ranki),lev_final(ranki)

     if (geom%f_comm%rank()==0) &
       print*, "Doing femps level ", k, ". Proc range:",lev_start(ranki),"-",lev_final(ranki)

     call fv3field_to_ufield(grid,geom%npx-1,psig(:,:,:,k),psiu)
     call fv3field_to_ufield(grid,geom%npx-1,chig(:,:,:,k),chiu)

     ! Convert to area integrals, required by femps
     psiu = psiu*grid%farea(:,grid%ngrids)
     chiu = chiu*grid%farea(:,grid%ngrids)

     ! Solve poisson equation (\psi=\nabla^{-2}\zeta, \chi=\nabla^{-2}D)
     call laplace(grid,oprs,grid%ngrids,psiu,voru)
     call laplace(grid,oprs,grid%ngrids,chiu,divu)

     ! Convert from area integrals
     voru = voru/grid%farea(:,grid%ngrids)
     divu = divu/grid%farea(:,grid%ngrids)

     call ufield_to_fv3field(grid,geom%npx-1,voru,vorg(:,:,:,k))
     call ufield_to_fv3field(grid,geom%npx-1,divu,divg(:,:,:,k))

   enddo
 endif

 ! Scatter field from root to all
 allocate(vorgcomm(1:geom%npx-1,1:geom%npy-1,6))
 allocate(divgcomm(1:geom%npx-1,1:geom%npy-1,6))
 do k=1,geom%npz
   if (ranki == lev_proc(k)) then
     vorgcomm = vorg(:,:,:,k)
     divgcomm = divg(:,:,:,k)
   endif
   call scatter_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,vorgcomm,vor(:,:,k))
   call scatter_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,divgcomm,div(:,:,k))
 enddo
 deallocate(vorgcomm,divgcomm)

end subroutine psichi_to_vortdivg

! --------------------------------------------------------------------------------------------------

subroutine udvd_to_vortdivg(geom, ud_in, vd_in, vort, divg)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  :: ud_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(in)  :: vd_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(out) ::  vort(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(out) ::  divg(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)

! D-Grid
integer :: i, j, k
real(kind=kind_real) :: ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)

! C-Grid
real(kind=kind_real) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: ut(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: vt(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)


! Fill D-grid winds halo
! ----------------------
ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) = ud_in
vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) = vd_in
call fill_dgrid_winds(geom, ud, vd, fillhalo=.true.)


! D -> C grid
! -----------
do k = 1,geom%npz
  call d_to_a_to_c_domain_level(geom, ud(:,:,k), vd(:,:,k), uc(:,:,k), vc(:,:,k))
enddo


 ! Compute ut,vt
 ! -------------
 call fill_cgrid_winds(geom, uc, vc, fillhalo=.true.)

 do k = 1,geom%npz

   call c_to_t_domain_level(geom, uc(:,:,k), vc(:,:,k), ut(:,:,k), vt(:,:,k), 1.0_kind_real)

   do j=geom%jsc,geom%jec
      do i=geom%isc,geom%iec+1
         if ( ut(i,j,k) > 0. ) then
              ut(i,j,k) = geom%dy(i,j)*ut(i,j,k)*geom%sin_sg(i-1,j,3)
         else
              ut(i,j,k) = geom%dy(i,j)*ut(i,j,k)*geom%sin_sg(i,j,1)
        endif
      enddo
   enddo
   do j=geom%jsc,geom%jec+1
      do i=geom%isc,geom%iec
         if ( vt(i,j,k) > 0. ) then
              vt(i,j,k) = geom%dx(i,j)*vt(i,j,k)*geom%sin_sg(i,j-1,4)
         else
              vt(i,j,k) = geom%dx(i,j)*vt(i,j,k)*geom%sin_sg(i,j,2)
         endif
      enddo
   enddo

 enddo

 !D/C grid u and v to A grid vorticity and divergence
 !---------------------------------------------------
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       vort(i,j,k) = geom%rarea(i,j)*( ud(i,j,k)*geom%dx(i,j)-ud(i,j+1,k)*geom%dx(i,j+1) - &
                                       vd(i,j,k)*geom%dy(i,j)+vd(i+1,j,k)*geom%dy(i+1,j))
       divg(i,j,k) = geom%rarea(i,j)*( ut(i+1,j,k)-ut(i,j,k) + vt(i,j+1,k)-vt(i,j,k) )
     enddo
   enddo
 enddo

end subroutine udvd_to_vortdivg

! --------------------------------------------------------------------------------------------------

subroutine d_to_a_to_c(geom, ud_in, vd_in, uc_out, vc_out, ua_out, va_out, ut_out, vt_out)

type(fv3jedi_geom), target,     intent(in)  :: geom
real(kind=kind_real),           intent(in)  ::  ud_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Dgrid winds
real(kind=kind_real),           intent(in)  ::  vd_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Dgrid winds
real(kind=kind_real),           intent(out) :: uc_out(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Cgrid winds
real(kind=kind_real),           intent(out) :: vc_out(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Cgrid winds
real(kind=kind_real), optional, intent(out) :: ua_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Agrid winds
real(kind=kind_real), optional, intent(out) :: va_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Agrid winds
real(kind=kind_real), optional, intent(out) :: ut_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), optional, intent(out) :: vt_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)

! Locals
integer :: k
real(kind=kind_real) :: ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
real(kind=kind_real) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real) :: ut(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: vt(geom%isd:geom%ied  ,geom%jsd:geom%jed  )


! Fill halo for input (outside loop for communication efficiency)
! ---------------------------------------------------------------
ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) = ud_in
vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) = vd_in
call fill_dgrid_winds(geom, ud, vd, fillhalo=.true.)


! Call domain based calculation
! -----------------------------
do k = 1, geom%npz
  call d_to_a_to_c_domain_level(geom, ud(:,:,k), vd(:,:,k), uc, vc, ua, va, ut, vt)

  ! Fill outputs
  uc_out(:,:,k) = uc(geom%isc:geom%iec+1,geom%jsc:geom%jec  )
  vc_out(:,:,k) = vc(geom%isc:geom%iec  ,geom%jsc:geom%jec+1)
  if (present(ua_out)) ua_out(:,:,k) = ua(geom%isc:geom%iec,geom%jsc:geom%jec)
  if (present(va_out)) va_out(:,:,k) = va(geom%isc:geom%iec,geom%jsc:geom%jec)
  if (present(ut_out)) ut_out(:,:,k) = ut(geom%isc:geom%iec,geom%jsc:geom%jec)
  if (present(vt_out)) vt_out(:,:,k) = vt(geom%isc:geom%iec,geom%jsc:geom%jec)
enddo

end subroutine d_to_a_to_c

! --------------------------------------------------------------------------------------------------

subroutine d_to_a_to_c_domain_level(geom, u, v, uc, vc, ua_out, va_out, ut_out, vt_out)

type(fv3jedi_geom), target,     intent(in)  :: geom
real(kind=kind_real),           intent(in)  ::      u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real),           intent(in)  ::      v(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
real(kind=kind_real),           intent(out) ::     uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
real(kind=kind_real),           intent(out) ::     vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real), optional, intent(out) :: ua_out(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real), optional, intent(out) :: va_out(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real), optional, intent(out) :: ut_out(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real), optional, intent(out) :: vt_out(geom%isd:geom%ied  ,geom%jsd:geom%jed  )

! Local
real(kind=kind_real), dimension(geom%isd:geom%ied,geom%jsd:geom%jed):: utmp, vtmp
integer npt, i, j, ifirst, ilast, id, npx, npy
integer :: is,  ie,  js,  je
integer :: isd, ied, jsd, jed

real(kind=kind_real) :: ua(geom%isd:geom%ied,geom%jsd:geom%jed)
real(kind=kind_real) :: va(geom%isd:geom%ied,geom%jsd:geom%jed)
real(kind=kind_real) :: ut(geom%isd:geom%ied,geom%jsd:geom%jed)
real(kind=kind_real) :: vt(geom%isd:geom%ied,geom%jsd:geom%jed)

real(kind=kind_real), pointer, dimension(:,:,:) :: sin_sg
real(kind=kind_real), pointer, dimension(:,:)   :: cosa_u, cosa_v, cosa_s
real(kind=kind_real), pointer, dimension(:,:)   :: rsin_u, rsin_v, rsin2
real(kind=kind_real), pointer, dimension(:,:)   :: dxa,dya

real(kind=kind_real), parameter:: big_number = 1.E30_kind_real
real(kind=kind_real), parameter:: a1 =  0.5625_kind_real
real(kind=kind_real), parameter:: a2 = -0.0625_kind_real
real(kind=kind_real), parameter:: c1 = -2._kind_real/14._kind_real
real(kind=kind_real), parameter:: c2 = 11._kind_real/14._kind_real
real(kind=kind_real), parameter:: c3 =  5._kind_real/14._kind_real

is  = geom%isc
ie  = geom%iec
js  = geom%jsc
je  = geom%jec
isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed
npx = geom%npx
npy = geom%npy

sin_sg    => geom%sin_sg
cosa_u    => geom%cosa_u
cosa_v    => geom%cosa_v
cosa_s    => geom%cosa_s
rsin_u    => geom%rsin_u
rsin_v    => geom%rsin_v
rsin2     => geom%rsin2
dxa       => geom%dxa
dya       => geom%dya

 if ( geom%dord4 ) then
      id = 1
 else
      id = 0
 endif

 if (geom%grid_type < 3 .and. .not. geom%bounded_domain) then
    npt = 4
 else
    npt = -2
 endif

! Initialize the non-existing corner regions
 utmp(:,:) = big_number
 vtmp(:,:) = big_number

if ( geom%bounded_domain ) then

    do j=jsd+1,jed-1
       do i=isd,ied
          utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
       enddo
    enddo
    do i=isd,ied
       j = jsd
       utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
       j = jed
       utmp(i,j) = 0.5*(u(i,j)+u(i,j+1))
    end do

    do j=jsd,jed
       do i=isd+1,ied-1
          vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
       enddo
       i = isd
       vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
       i = ied
       vtmp(i,j) = 0.5*(v(i,j)+v(i+1,j))
    enddo

    do j=jsd,jed
       do i=isd,ied
          ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
          va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
       enddo
    enddo

else
    !----------
    ! Interior:
    !----------
    do j=max(npt,js-1),min(npy-npt,je+1)
       do i=max(npt,isd),min(npx-npt,ied)
          utmp(i,j) = a2*(u(i,j-1)+u(i,j+2)) + a1*(u(i,j)+u(i,j+1))
       enddo
    enddo
    do j=max(npt,jsd),min(npy-npt,jed)
       do i=max(npt,is-1),min(npx-npt,ie+1)
          vtmp(i,j) = a2*(v(i-1,j)+v(i+2,j)) + a1*(v(i,j)+v(i+1,j))
       enddo
    enddo

    !----------
    ! edges:
    !----------
    if (geom%grid_type < 3) then

       if ( js==1 .or. jsd<npt) then
          do j=jsd,npt-1
             do i=isd,ied
                utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
             enddo
          enddo
       endif
       if ( (je+1)==npy .or. jed>=(npy-npt)) then
          do j=npy-npt+1,jed
             do i=isd,ied
                utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
             enddo
          enddo
       endif

       if ( is==1 .or. isd<npt ) then
          do j=max(npt,jsd),min(npy-npt,jed)
             do i=isd,npt-1
                utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
             enddo
          enddo
       endif
       if ( (ie+1)==npx .or. ied>=(npx-npt)) then
          do j=max(npt,jsd),min(npy-npt,jed)
             do i=npx-npt+1,ied
                utmp(i,j) = 0.5*(u(i,j) + u(i,j+1))
                vtmp(i,j) = 0.5*(v(i,j) + v(i+1,j))
             enddo
          enddo
       endif

    endif

! Contra-variant components at cell center:
    do j=js-1-id,je+1+id
       do i=is-1-id,ie+1+id
          ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
          va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
       enddo
    enddo

end if

! A -> C
!--------------
! Fix the edges
!--------------
! Xdir:
    if( geom%sw_corner ) then
        do i=-2,0
           utmp(i,0) = -vtmp(0,1-i)
        enddo
    endif
    if( geom%se_corner ) then
        do i=0,2
           utmp(npx+i,0) = vtmp(npx,i+1)
        enddo
    endif
    if( geom%ne_corner ) then
        do i=0,2
           utmp(npx+i,npy) = -vtmp(npx,je-i)
        enddo
    endif
    if( geom%nw_corner ) then
        do i=-2,0
           utmp(i,npy) = vtmp(0,je+i)
        enddo
    endif

 if (geom%grid_type < 3 .and. .not. geom%bounded_domain) then
    ifirst = max(3,    is-1)
    ilast  = min(npx-2,ie+2)
 else
    ifirst = is-1
    ilast  = ie+2
 endif
!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
    do j=js-1,je+1
       do i=ifirst,ilast
          uc(i,j) = a2*(utmp(i-2,j)+utmp(i+1,j)) + a1*(utmp(i-1,j)+utmp(i,j))
          ut(i,j) = (uc(i,j) - v(i,j)*cosa_u(i,j))*rsin_u(i,j)
       enddo
    enddo

if (geom%grid_type < 3) then
! Xdir:
    if( geom%sw_corner ) then
        ua(-1,0) = -va(0,2)
        ua( 0,0) = -va(0,1)
    endif
    if( geom%se_corner ) then
        ua(npx,  0) = va(npx,1)
        ua(npx+1,0) = va(npx,2)
    endif
    if( geom%ne_corner ) then
        ua(npx,  npy) = -va(npx,npy-1)
        ua(npx+1,npy) = -va(npx,npy-2)
    endif
    if( geom%nw_corner ) then
        ua(-1,npy) = va(0,npy-2)
        ua( 0,npy) = va(0,npy-1)
    endif

    if( is==1 .and. .not. geom%bounded_domain  ) then
       do j=js-1,je+1
          uc(0,j) = c1*utmp(-2,j) + c2*utmp(-1,j) + c3*utmp(0,j)
          ut(1,j) = edge_interpolate4(ua(-1:2,j), dxa(-1:2,j))
          !Want to use the UPSTREAM value
          if (ut(1,j) > 0.) then
             uc(1,j) = ut(1,j)*sin_sg(0,j,3)
          else
             uc(1,j) = ut(1,j)*sin_sg(1,j,1)
          end if
          uc(2,j) = c1*utmp(3,j) + c2*utmp(2,j) + c3*utmp(1,j)
          ut(0,j) = (uc(0,j) - v(0,j)*cosa_u(0,j))*rsin_u(0,j)
          ut(2,j) = (uc(2,j) - v(2,j)*cosa_u(2,j))*rsin_u(2,j)
       enddo
    endif

    if( (ie+1)==npx  .and. .not. geom%bounded_domain ) then
       do j=js-1,je+1
          uc(npx-1,j) = c1*utmp(npx-3,j)+c2*utmp(npx-2,j)+c3*utmp(npx-1,j)
          ut(npx,  j) = edge_interpolate4(ua(npx-2:npx+1,j), dxa(npx-2:npx+1,j))
          if (ut(npx,j) > 0.) then
              uc(npx,j) = ut(npx,j)*sin_sg(npx-1,j,3)
          else
              uc(npx,j) = ut(npx,j)*sin_sg(npx,j,1)
          end if
          uc(npx+1,j) = c3*utmp(npx,j) + c2*utmp(npx+1,j) + c1*utmp(npx+2,j)
          ut(npx-1,j) = (uc(npx-1,j)-v(npx-1,j)*cosa_u(npx-1,j))*rsin_u(npx-1,j)
          ut(npx+1,j) = (uc(npx+1,j)-v(npx+1,j)*cosa_u(npx+1,j))*rsin_u(npx+1,j)
       enddo
    endif

endif

!------
! Ydir:
!------
    if( geom%sw_corner ) then
        do j=-2,0
           vtmp(0,j) = -utmp(1-j,0)
        enddo
    endif
    if( geom%nw_corner ) then
        do j=0,2
           vtmp(0,npy+j) = utmp(j+1,npy)
        enddo
    endif
    if( geom%se_corner ) then
        do j=-2,0
           vtmp(npx,j) = utmp(ie+j,0)
        enddo
    endif
    if( geom%ne_corner ) then
        do j=0,2
           vtmp(npx,npy+j) = -utmp(ie-j,npy)
        enddo
    endif
    if( geom%sw_corner ) then
        va(0,-1) = -ua(2,0)
        va(0, 0) = -ua(1,0)
    endif
    if( geom%se_corner ) then
        va(npx, 0) = ua(npx-1,0)
        va(npx,-1) = ua(npx-2,0)
    endif
    if( geom%ne_corner ) then
        va(npx,npy  ) = -ua(npx-1,npy)
        va(npx,npy+1) = -ua(npx-2,npy)
    endif
    if( geom%nw_corner ) then
        va(0,npy)   = ua(1,npy)
        va(0,npy+1) = ua(2,npy)
    endif

if (geom%grid_type < 3) then

    do j=js-1,je+2
     if ( j==1 .and. .not. geom%bounded_domain  ) then
       do i=is-1,ie+1
          vt(i,j) = edge_interpolate4(va(i,-1:2), dya(i,-1:2))
          if (vt(i,j) > 0.) then
             vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
          else
             vc(i,j) = vt(i,j)*sin_sg(i,j,2)
          end if
       enddo
     elseif ( j==0 .or. j==(npy-1) .and. .not. geom%bounded_domain  ) then
       do i=is-1,ie+1
          vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
          vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
       enddo
     elseif ( j==2 .or. j==(npy+1)  .and. .not. geom%bounded_domain ) then
       do i=is-1,ie+1
          vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
          vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
       enddo
     elseif ( j==npy .and. .not. geom%bounded_domain  ) then
       do i=is-1,ie+1
          vt(i,j) = edge_interpolate4(va(i,j-2:j+1), dya(i,j-2:j+1))
          if (vt(i,j) > 0.) then
             vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
          else
             vc(i,j) = vt(i,j)*sin_sg(i,j,2)
          end if
       enddo
     else
! 4th order interpolation for interior points:
       do i=is-1,ie+1
          vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1)) + a1*(vtmp(i,j-1)+vtmp(i,j))
          vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
       enddo
     endif
    enddo
else
! 4th order interpolation:
      do j=js-1,je+2
         do i=is-1,ie+1
            vc(i,j) = a2*(vtmp(i,j-2)+vtmp(i,j+1)) + a1*(vtmp(i,j-1)+vtmp(i,j))
            vt(i,j) = vc(i,j)
         enddo
      enddo
endif

! Fill outputs
if (present(ua_out)) ua_out = ua
if (present(va_out)) va_out = va
if (present(ut_out)) ut_out = ut
if (present(vt_out)) vt_out = vt

end subroutine d_to_a_to_c_domain_level

! --------------------------------------------------------------------------------------------------

subroutine c_to_a_to_d(geom, uc, vc, ud, vd, ua, va)

type(fv3jedi_geom), target,     intent(in)  :: geom
real(kind=kind_real),           intent(in)  :: uc(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Cgrid winds
real(kind=kind_real),           intent(in)  :: vc(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Cgrid winds
real(kind=kind_real),           intent(out) :: ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Dgrid winds
real(kind=kind_real),           intent(out) :: vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Dgrid winds
real(kind=kind_real), optional, intent(out) :: ua(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Agrid winds
real(kind=kind_real), optional, intent(out) :: va(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Agrid winds

integer :: k
real(kind=kind_real) :: ua_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds
real(kind=kind_real) :: va_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds

! C to A grid (cubed-sphere)
call c_to_acs(geom, uc, vc, ua_tmp, va_tmp)

! A (cubed-sphere) to D grid
call acs_to_d(geom, ua_tmp, va_tmp, ud, vd )

! Fill intermediates if requested.
if (present(ua)) ua = ua_tmp
if (present(va)) va = va_tmp

end subroutine c_to_a_to_d

! --------------------------------------------------------------------------------------------------

subroutine a_to_d(geom, ua, va, ud, vd)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  :: ua(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(in)  :: va(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(out) :: ud(geom%isc:geom%iec,  geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(out) :: vd(geom%isc:geom%iec+1,geom%jsc:geom%jec,  geom%npz)

integer :: is ,ie , js ,je
integer :: npx, npy, npz
integer :: i,j,k, im2,jm2

real(kind=kind_real) :: uatemp(geom%isd:geom%ied,geom%jsd:geom%jed,geom%npz)
real(kind=kind_real) :: vatemp(geom%isd:geom%ied,geom%jsd:geom%jed,geom%npz)

real(kind=kind_real) :: v3(geom%isc-1:geom%iec+1,geom%jsc-1:geom%jec+1,3)
real(kind=kind_real) :: ue(geom%isc-1:geom%iec+1,geom%jsc  :geom%jec+1,3)    ! 3D winds at edges
real(kind=kind_real) :: ve(geom%isc  :geom%iec+1,geom%jsc-1:geom%jec+1,3)    ! 3D winds at edges
real(kind=kind_real), dimension(geom%isc:geom%iec):: ut1, ut2, ut3
real(kind=kind_real), dimension(geom%jsc:geom%jec):: vt1, vt2, vt3
npx = geom%npx
npy = geom%npy
npz = geom%npz
is  = geom%isc
ie  = geom%iec
js  = geom%jsc
je  = geom%jec

im2 = (npx-1)/2
jm2 = (npy-1)/2

uatemp(:,:,:) = 0.0
vatemp(:,:,:) = 0.0

uatemp(is:ie,js:je,:) = ua
vatemp(is:ie,js:je,:) = va

call mpp_update_domains(uatemp, geom%domain, complete=.true.)
call mpp_update_domains(vatemp, geom%domain, complete=.true.)

do k=1, npz

  do j=js-1,je+1
    do i=is-1,ie+1
      v3(i,j,1) = uatemp(i,j,k)*geom%vlon(i,j,1) + vatemp(i,j,k)*geom%vlat(i,j,1)
      v3(i,j,2) = uatemp(i,j,k)*geom%vlon(i,j,2) + vatemp(i,j,k)*geom%vlat(i,j,2)
      v3(i,j,3) = uatemp(i,j,k)*geom%vlon(i,j,3) + vatemp(i,j,k)*geom%vlat(i,j,3)
    enddo
  enddo

  do j=js,je+1
    do i=is-1,ie+1
      ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
      ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
      ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
    enddo
  enddo

  do j=js-1,je+1
    do i=is,ie+1
      ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
      ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
      ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
    enddo
  enddo
 if (.not. geom%bounded_domain) then
  if ( is==1 ) then
    i = 1
    do j=js,je
      if ( j>jm2 ) then
        vt1(j) = geom%edge_vect_w(j)*ve(i,j-1,1)+(1.-geom%edge_vect_w(j))*ve(i,j,1)
        vt2(j) = geom%edge_vect_w(j)*ve(i,j-1,2)+(1.-geom%edge_vect_w(j))*ve(i,j,2)
        vt3(j) = geom%edge_vect_w(j)*ve(i,j-1,3)+(1.-geom%edge_vect_w(j))*ve(i,j,3)
      else
        vt1(j) = geom%edge_vect_w(j)*ve(i,j+1,1)+(1.-geom%edge_vect_w(j))*ve(i,j,1)
        vt2(j) = geom%edge_vect_w(j)*ve(i,j+1,2)+(1.-geom%edge_vect_w(j))*ve(i,j,2)
        vt3(j) = geom%edge_vect_w(j)*ve(i,j+1,3)+(1.-geom%edge_vect_w(j))*ve(i,j,3)
      endif
    enddo
    do j=js,je
      ve(i,j,1) = vt1(j)
      ve(i,j,2) = vt2(j)
      ve(i,j,3) = vt3(j)
    enddo
  endif

  if ( (ie+1)==npx ) then
    i = npx
    do j=js,je
      if ( j>jm2 ) then
        vt1(j) = geom%edge_vect_e(j)*ve(i,j-1,1)+(1.-geom%edge_vect_e(j))*ve(i,j,1)
        vt2(j) = geom%edge_vect_e(j)*ve(i,j-1,2)+(1.-geom%edge_vect_e(j))*ve(i,j,2)
        vt3(j) = geom%edge_vect_e(j)*ve(i,j-1,3)+(1.-geom%edge_vect_e(j))*ve(i,j,3)
      else
        vt1(j) = geom%edge_vect_e(j)*ve(i,j+1,1)+(1.-geom%edge_vect_e(j))*ve(i,j,1)
        vt2(j) = geom%edge_vect_e(j)*ve(i,j+1,2)+(1.-geom%edge_vect_e(j))*ve(i,j,2)
        vt3(j) = geom%edge_vect_e(j)*ve(i,j+1,3)+(1.-geom%edge_vect_e(j))*ve(i,j,3)
      endif
    enddo
    do j=js,je
      ve(i,j,1) = vt1(j)
      ve(i,j,2) = vt2(j)
      ve(i,j,3) = vt3(j)
    enddo
  endif

  if ( js==1 ) then
    j = 1
    do i=is,ie
      if ( i>im2 ) then
        ut1(i) = geom%edge_vect_s(i)*ue(i-1,j,1)+(1.-geom%edge_vect_s(i))*ue(i,j,1)
        ut2(i) = geom%edge_vect_s(i)*ue(i-1,j,2)+(1.-geom%edge_vect_s(i))*ue(i,j,2)
        ut3(i) = geom%edge_vect_s(i)*ue(i-1,j,3)+(1.-geom%edge_vect_s(i))*ue(i,j,3)
      else
        ut1(i) = geom%edge_vect_s(i)*ue(i+1,j,1)+(1.-geom%edge_vect_s(i))*ue(i,j,1)
        ut2(i) = geom%edge_vect_s(i)*ue(i+1,j,2)+(1.-geom%edge_vect_s(i))*ue(i,j,2)
        ut3(i) = geom%edge_vect_s(i)*ue(i+1,j,3)+(1.-geom%edge_vect_s(i))*ue(i,j,3)
      endif
    enddo
    do i=is,ie
      ue(i,j,1) = ut1(i)
      ue(i,j,2) = ut2(i)
      ue(i,j,3) = ut3(i)
    enddo
  endif

  if ( (je+1)==npy ) then
    j = npy
    do i=is,ie
      if ( i>im2 ) then
        ut1(i) = geom%edge_vect_n(i)*ue(i-1,j,1)+(1.-geom%edge_vect_n(i))*ue(i,j,1)
        ut2(i) = geom%edge_vect_n(i)*ue(i-1,j,2)+(1.-geom%edge_vect_n(i))*ue(i,j,2)
        ut3(i) = geom%edge_vect_n(i)*ue(i-1,j,3)+(1.-geom%edge_vect_n(i))*ue(i,j,3)
      else
        ut1(i) = geom%edge_vect_n(i)*ue(i+1,j,1)+(1.-geom%edge_vect_n(i))*ue(i,j,1)
        ut2(i) = geom%edge_vect_n(i)*ue(i+1,j,2)+(1.-geom%edge_vect_n(i))*ue(i,j,2)
        ut3(i) = geom%edge_vect_n(i)*ue(i+1,j,3)+(1.-geom%edge_vect_n(i))*ue(i,j,3)
      endif
    enddo
    do i=is,ie
      ue(i,j,1) = ut1(i)
      ue(i,j,2) = ut2(i)
      ue(i,j,3) = ut3(i)
    enddo
  endif

 endif ! .not. bounded_domain

  do j=js,je+1
    do i=is,ie
      ud(i,j,k) = 0.5*( ue(i,j,1)*geom%es(1,i,j,1) +  &
                        ue(i,j,2)*geom%es(2,i,j,1) +  &
                        ue(i,j,3)*geom%es(3,i,j,1) )
    enddo
  enddo

  do j=js,je
    do i=is,ie+1
      vd(i,j,k) = 0.5*( ve(i,j,1)*geom%ew(1,i,j,2) +  &
                        ve(i,j,2)*geom%ew(2,i,j,2) +  &
                        ve(i,j,3)*geom%ew(3,i,j,2) )
    enddo
  enddo

enddo

end subroutine a_to_d

! --------------------------------------------------------------------------------------------------

subroutine a_to_d_ad(geom, ua_ad, va_ad, ud_ad, vd_ad)

type(fv3jedi_geom), intent(in)      :: geom
real(kind=kind_real), intent(inout) :: ua_ad(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(inout) :: va_ad(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(out)   :: ud_ad(geom%isc:geom%iec,  geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(out)   :: vd_ad(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz)

integer :: is, ie, js, je
integer :: npx, npy, npz
integer :: i, j, k, im2, jm2

real(kind=kind_real) :: uatemp(geom%isd:geom%ied, geom%jsd:geom%jed, geom%npz)
real(kind=kind_real) :: vatemp(geom%isd:geom%ied, geom%jsd:geom%jed, geom%npz)

real(kind=kind_real) :: uatemp_ad(geom%isd:geom%ied, geom%jsd:geom%jed, geom%npz)
real(kind=kind_real) :: vatemp_ad(geom%isd:geom%ied, geom%jsd:geom%jed, geom%npz)
real(kind=kind_real) :: v3_ad(geom%isc-1:geom%iec+1, geom%jsc-1:geom%jec+1, 3)
real(kind=kind_real) :: ue_ad(geom%isc-1:geom%iec+1, geom%jsc:geom%jec+1, 3)
real(kind=kind_real) :: ve_ad(geom%isc:geom%iec+1, geom%jsc-1:geom%jec+1, 3)
real(kind=kind_real), dimension(geom%isc:geom%iec) :: ut1_ad, ut2_ad, ut3_ad
real(kind=kind_real), dimension(geom%jsc:geom%jec) :: vt1_ad, vt2_ad, vt3_ad
real(kind=kind_real) :: temp_ad
real(kind=kind_real) :: temp_ad0

npx = geom%npx
npy = geom%npy
npz = geom%npz
is = geom%isc
ie = geom%iec
js = geom%jsc
je = geom%jec
im2 = (npx-1)/2
jm2 = (npy-1)/2

uatemp(:, :, :) = 0.0
vatemp(:, :, :) = 0.0

v3_ad = 0.0_kind_real
uatemp_ad = 0.0_kind_real
ue_ad = 0.0_kind_real
vt1_ad = 0.0_kind_real
vt2_ad = 0.0_kind_real
vt3_ad = 0.0_kind_real
ve_ad = 0.0_kind_real
vatemp_ad = 0.0_kind_real
ut1_ad = 0.0_kind_real
ut2_ad = 0.0_kind_real
ut3_ad = 0.0_kind_real

do k=npz,1,-1

  do j=je,js,-1
    do i=ie+1,is,-1
      temp_ad0 = 0.5*vd_ad(i, j, k)
      ve_ad(i, j, 1) = ve_ad(i, j, 1) + geom%ew(1, i, j, 2)*temp_ad0
      ve_ad(i, j, 2) = ve_ad(i, j, 2) + geom%ew(2, i, j, 2)*temp_ad0
      ve_ad(i, j, 3) = ve_ad(i, j, 3) + geom%ew(3, i, j, 2)*temp_ad0
      vd_ad(i, j, k) = 0.0_kind_real
    end do
  end do

  do j=je+1,js,-1
    do i=ie,is,-1
      temp_ad = 0.5*ud_ad(i, j, k)
      ue_ad(i, j, 1) = ue_ad(i, j, 1) + geom%es(1, i, j, 1)*temp_ad
      ue_ad(i, j, 2) = ue_ad(i, j, 2) + geom%es(2, i, j, 1)*temp_ad
      ue_ad(i, j, 3) = ue_ad(i, j, 3) + geom%es(3, i, j, 1)*temp_ad
      ud_ad(i, j, k) = 0.0_kind_real
    end do
  end do
 if (.not. geom%bounded_domain) then
  if (je + 1 .eq. npy) then
    do i=ie,is,-1
      ut3_ad(i) = ut3_ad(i) + ue_ad(i, npy, 3)
      ue_ad(i, npy, 3) = 0.0_kind_real
      ut2_ad(i) = ut2_ad(i) + ue_ad(i, npy, 2)
      ue_ad(i, npy, 2) = 0.0_kind_real
      ut1_ad(i) = ut1_ad(i) + ue_ad(i, npy, 1)
      ue_ad(i, npy, 1) = 0.0_kind_real
    end do
    do i=ie,is,-1
      if (i .le. im2) then
        ue_ad(i+1, npy, 3) = ue_ad(i+1, npy, 3) + geom%edge_vect_n(i)*ut3_ad(i)
        ue_ad(i, npy, 3) = ue_ad(i, npy, 3) + (1.-geom%edge_vect_n(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_kind_real
        ue_ad(i+1, npy, 2) = ue_ad(i+1, npy, 2) + geom%edge_vect_n(i)*ut2_ad(i)
        ue_ad(i, npy, 2) = ue_ad(i, npy, 2) + (1.-geom%edge_vect_n(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_kind_real
        ue_ad(i+1, npy, 1) = ue_ad(i+1, npy, 1) + geom%edge_vect_n(i)*ut1_ad(i)
        ue_ad(i, npy, 1) = ue_ad(i, npy, 1) + (1.-geom%edge_vect_n(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_kind_real
      else
        ue_ad(i-1, npy, 3) = ue_ad(i-1, npy, 3) + geom%edge_vect_n(i)*ut3_ad(i)
        ue_ad(i, npy, 3) = ue_ad(i, npy, 3) + (1.-geom%edge_vect_n(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_kind_real
        ue_ad(i-1, npy, 2) = ue_ad(i-1, npy, 2) + geom%edge_vect_n(i)*ut2_ad(i)
        ue_ad(i, npy, 2) = ue_ad(i, npy, 2) + (1.-geom%edge_vect_n(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_kind_real
        ue_ad(i-1, npy, 1) = ue_ad(i-1, npy, 1) + geom%edge_vect_n(i)*ut1_ad(i)
        ue_ad(i, npy, 1) = ue_ad(i, npy, 1) + (1.-geom%edge_vect_n(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_kind_real
      end if
    end do
  end if

  if (js .eq. 1) then
    do i=ie,is,-1
      ut3_ad(i) = ut3_ad(i) + ue_ad(i, 1, 3)
      ue_ad(i, 1, 3) = 0.0_kind_real
      ut2_ad(i) = ut2_ad(i) + ue_ad(i, 1, 2)
      ue_ad(i, 1, 2) = 0.0_kind_real
      ut1_ad(i) = ut1_ad(i) + ue_ad(i, 1, 1)
      ue_ad(i, 1, 1) = 0.0_kind_real
    end do
    do i=ie,is,-1
      if (i .le. im2) then
        ue_ad(i+1, 1, 3) = ue_ad(i+1, 1, 3) + geom%edge_vect_s(i)*ut3_ad(i)
        ue_ad(i, 1, 3) = ue_ad(i, 1, 3) + (1.-geom%edge_vect_s(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_kind_real
        ue_ad(i+1, 1, 2) = ue_ad(i+1, 1, 2) + geom%edge_vect_s(i)*ut2_ad(i)
        ue_ad(i, 1, 2) = ue_ad(i, 1, 2) + (1.-geom%edge_vect_s(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_kind_real
        ue_ad(i+1, 1, 1) = ue_ad(i+1, 1, 1) + geom%edge_vect_s(i)*ut1_ad(i)
        ue_ad(i, 1, 1) = ue_ad(i, 1, 1) + (1.-geom%edge_vect_s(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_kind_real
      else
        ue_ad(i-1, 1, 3) = ue_ad(i-1, 1, 3) + geom%edge_vect_s(i)*ut3_ad(i)
        ue_ad(i, 1, 3) = ue_ad(i, 1, 3) + (1.-geom%edge_vect_s(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_kind_real
        ue_ad(i-1, 1, 2) = ue_ad(i-1, 1, 2) + geom%edge_vect_s(i)*ut2_ad(i)
        ue_ad(i, 1, 2) = ue_ad(i, 1, 2) + (1.-geom%edge_vect_s(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_kind_real
        ue_ad(i-1, 1, 1) = ue_ad(i-1, 1, 1) + geom%edge_vect_s(i)*ut1_ad(i)
        ue_ad(i, 1, 1) = ue_ad(i, 1, 1) + (1.-geom%edge_vect_s(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_kind_real
      end if
    end do
  end if

  if (ie + 1 .eq. npx) then
    do j=je,js,-1
      vt3_ad(j) = vt3_ad(j) + ve_ad(npx, j, 3)
      ve_ad(npx, j, 3) = 0.0_kind_real
      vt2_ad(j) = vt2_ad(j) + ve_ad(npx, j, 2)
      ve_ad(npx, j, 2) = 0.0_kind_real
      vt1_ad(j) = vt1_ad(j) + ve_ad(npx, j, 1)
      ve_ad(npx, j, 1) = 0.0_kind_real
    end do
    do j=je,js,-1
      if (j .le. jm2) then
        ve_ad(npx, j+1, 3) = ve_ad(npx, j+1, 3) + geom%edge_vect_e(j)*vt3_ad(j)
        ve_ad(npx, j, 3) = ve_ad(npx, j, 3) + (1.-geom%edge_vect_e(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_kind_real
        ve_ad(npx, j+1, 2) = ve_ad(npx, j+1, 2) + geom%edge_vect_e(j)*vt2_ad(j)
        ve_ad(npx, j, 2) = ve_ad(npx, j, 2) + (1.-geom%edge_vect_e(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_kind_real
        ve_ad(npx, j+1, 1) = ve_ad(npx, j+1, 1) + geom%edge_vect_e(j)*vt1_ad(j)
        ve_ad(npx, j, 1) = ve_ad(npx, j, 1) + (1.-geom%edge_vect_e(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_kind_real
      else
        ve_ad(npx, j-1, 3) = ve_ad(npx, j-1, 3) + geom%edge_vect_e(j)*vt3_ad(j)
        ve_ad(npx, j, 3) = ve_ad(npx, j, 3) + (1.-geom%edge_vect_e(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_kind_real
        ve_ad(npx, j-1, 2) = ve_ad(npx, j-1, 2) + geom%edge_vect_e(j)*vt2_ad(j)
        ve_ad(npx, j, 2) = ve_ad(npx, j, 2) + (1.-geom%edge_vect_e(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_kind_real
        ve_ad(npx, j-1, 1) = ve_ad(npx, j-1, 1) + geom%edge_vect_e(j)*vt1_ad(j)
        ve_ad(npx, j, 1) = ve_ad(npx, j, 1) + (1.-geom%edge_vect_e(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_kind_real
      end if
    end do
  end if

  if (is .eq. 1) then
    do j=je,js,-1
      vt3_ad(j) = vt3_ad(j) + ve_ad(1, j, 3)
      ve_ad(1, j, 3) = 0.0_kind_real
      vt2_ad(j) = vt2_ad(j) + ve_ad(1, j, 2)
      ve_ad(1, j, 2) = 0.0_kind_real
      vt1_ad(j) = vt1_ad(j) + ve_ad(1, j, 1)
      ve_ad(1, j, 1) = 0.0_kind_real
    end do
    do j=je,js,-1
      if (j .le. jm2) then
        ve_ad(1, j+1, 3) = ve_ad(1, j+1, 3) + geom%edge_vect_w(j)*vt3_ad(j)
        ve_ad(1, j, 3) = ve_ad(1, j, 3) + (1.-geom%edge_vect_w(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_kind_real
        ve_ad(1, j+1, 2) = ve_ad(1, j+1, 2) + geom%edge_vect_w(j)*vt2_ad(j)
        ve_ad(1, j, 2) = ve_ad(1, j, 2) + (1.-geom%edge_vect_w(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_kind_real
        ve_ad(1, j+1, 1) = ve_ad(1, j+1, 1) + geom%edge_vect_w(j)*vt1_ad(j)
        ve_ad(1, j, 1) = ve_ad(1, j, 1) + (1.-geom%edge_vect_w(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_kind_real
      else
        ve_ad(1, j-1, 3) = ve_ad(1, j-1, 3) + geom%edge_vect_w(j)*vt3_ad(j)
        ve_ad(1, j, 3) = ve_ad(1, j, 3) + (1.-geom%edge_vect_w(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_kind_real
        ve_ad(1, j-1, 2) = ve_ad(1, j-1, 2) + geom%edge_vect_w(j)*vt2_ad(j)
        ve_ad(1, j, 2) = ve_ad(1, j, 2) + (1.-geom%edge_vect_w(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_kind_real
        ve_ad(1, j-1, 1) = ve_ad(1, j-1, 1) + geom%edge_vect_w(j)*vt1_ad(j)
        ve_ad(1, j, 1) = ve_ad(1, j, 1) + (1.-geom%edge_vect_w(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_kind_real
      end if
    end do
  end if
 end if  !clt bounded=.false.

  do j=je+1,js-1,-1
    do i=ie+1,is,-1
      v3_ad(i-1, j, 3) = v3_ad(i-1, j, 3) + ve_ad(i, j, 3)
      v3_ad(i, j, 3) = v3_ad(i, j, 3) + ve_ad(i, j, 3)
      ve_ad(i, j, 3) = 0.0_kind_real
      v3_ad(i-1, j, 2) = v3_ad(i-1, j, 2) + ve_ad(i, j, 2)
      v3_ad(i, j, 2) = v3_ad(i, j, 2) + ve_ad(i, j, 2)
      ve_ad(i, j, 2) = 0.0_kind_real
      v3_ad(i-1, j, 1) = v3_ad(i-1, j, 1) + ve_ad(i, j, 1)
      v3_ad(i, j, 1) = v3_ad(i, j, 1) + ve_ad(i, j, 1)
      ve_ad(i, j, 1) = 0.0_kind_real
    end do
  end do

  do j=je+1,js,-1
    do i=ie+1,is-1,-1
      v3_ad(i, j-1, 3) = v3_ad(i, j-1, 3) + ue_ad(i, j, 3)
      v3_ad(i, j, 3) = v3_ad(i, j, 3) + ue_ad(i, j, 3)
      ue_ad(i, j, 3) = 0.0_kind_real
      v3_ad(i, j-1, 2) = v3_ad(i, j-1, 2) + ue_ad(i, j, 2)
      v3_ad(i, j, 2) = v3_ad(i, j, 2) + ue_ad(i, j, 2)
      ue_ad(i, j, 2) = 0.0_kind_real
      v3_ad(i, j-1, 1) = v3_ad(i, j-1, 1) + ue_ad(i, j, 1)
      v3_ad(i, j, 1) = v3_ad(i, j, 1) + ue_ad(i, j, 1)
      ue_ad(i, j, 1) = 0.0_kind_real
    end do
  end do

  do j=je+1,js-1,-1
    do i=ie+1,is-1,-1
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 3)*v3_ad(i, j, 3)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 3)*v3_ad(i, j, 3)
      v3_ad(i, j, 3) = 0.0_kind_real
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 2)*v3_ad(i, j, 2)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 2)*v3_ad(i, j, 2)
      v3_ad(i, j, 2) = 0.0_kind_real
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 1)*v3_ad(i, j, 1)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 1)*v3_ad(i, j, 1)
      v3_ad(i, j, 1) = 0.0_kind_real
    end do
  end do
end do

call mpp_update_domains_adm(vatemp, vatemp_ad, geom%domain, complete=.true.)
call mpp_update_domains_adm(uatemp, uatemp_ad, geom%domain, complete=.true.)

va_ad = va_ad + vatemp_ad(is:ie, js:je, :)
ua_ad = ua_ad + uatemp_ad(is:ie, js:je, :)

end subroutine a_to_d_ad

! --------------------------------------------------------------------------------------------------

subroutine d_to_a(geom, ud_in, vd_in, ua_out, va_out)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  ::  ud_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(in)  ::  vd_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz)
real(kind=kind_real), intent(out) :: ua_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)
real(kind=kind_real), intent(out) :: va_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)

integer :: k
real(kind=kind_real) :: ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,geom%npz)
real(kind=kind_real) :: vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,geom%npz)
real(kind=kind_real) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )

!Fill compute part from input
ud = 0.0_kind_real
vd = 0.0_kind_real
ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,:) = ud_in
vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,:) = vd_in

call fill_dgrid_winds(geom, ud, vd, .true.)

do k = 1,geom%npz

  call d_to_a_domain_level(geom, ud(:,:,k), vd(:,:,k), ua, va)

  ua_out(:,:,k) = ua(geom%isc:geom%iec,geom%jsc:geom%jec)
  va_out(:,:,k) = va(geom%isc:geom%iec,geom%jsc:geom%jec)

enddo

end subroutine d_to_a

! --------------------------------------------------------------------------------------------------

subroutine d_to_a_domain_level(geom, u, v, ua, va)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(in)    ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real), intent(in)    ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
real(kind=kind_real), intent(out)   :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real), intent(out)   :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )

integer i, j, k
integer :: is,  ie,  js,  je, npx, npy, npz
real(kind=kind_real) :: c1 =  1.125
real(kind=kind_real) :: c2 = -0.125
real(kind=kind_real) :: utmp(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
real(kind=kind_real) :: vtmp(geom%isc:geom%iec+1,geom%jsc:geom%jec)
real(kind=kind_real) :: wu(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
real(kind=kind_real) :: wv(geom%isc:geom%iec+1,geom%jsc:geom%jec)

!Shorten notation
is  = geom%isc
ie  = geom%iec
js  = geom%jsc
je  = geom%jec
npx = geom%npx
npy = geom%npy
npz = geom%npz

if (geom%bounded_domain) then

  do j=max(1,js),min(npy-1,je)
    do i=max(1,is),min(npx-1,ie)
      utmp(i,j) = c2*(u(i,j-1)+u(i,j+2)) + c1*(u(i,j)+u(i,j+1))
      vtmp(i,j) = c2*(v(i-1,j)+v(i+2,j)) + c1*(v(i,j)+v(i+1,j))
    enddo
  enddo

else

  do j=max(2,js),min(npy-2,je)
    do i=max(2,is),min(npx-2,ie)
      utmp(i,j) = c2*(u(i,j-1)+u(i,j+2)) + c1*(u(i,j)+u(i,j+1))
      vtmp(i,j) = c2*(v(i-1,j)+v(i+2,j)) + c1*(v(i,j)+v(i+1,j))
    enddo
  enddo

  if ( js==1  ) then
    do i=is,ie+1
      wv(i,1) = v(i,1)*geom%dy(i,1)
    enddo
    do i=is,ie
      vtmp(i,1) = 2.*(wv(i,1) + wv(i+1,1)) / (geom%dy(i,1)+geom%dy(i+1,1))
      utmp(i,1) = 2.*(u(i,1)*geom%dx(i,1) + u(i,2)*geom%dx(i,2)) / (geom%dx(i,1) + geom%dx(i,2))
    enddo
  endif

  if ( (je+1)==npy   ) then
    j = npy-1
    do i=is,ie+1
      wv(i,j) = v(i,j)*geom%dy(i,j)
    enddo
    do i=is,ie
      vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j)) / (geom%dy(i,j)+geom%dy(i+1,j))
      utmp(i,j) = 2.*(u(i,j)*geom%dx(i,j) + u(i,j+1)*geom%dx(i,j+1)) / (geom%dx(i,j) + geom%dx(i,j+1))
    enddo
  endif

  if ( is==1 ) then
    i = 1
    do j=js,je
      wv(1,j) = v(1,j)*geom%dy(1,j)
      wv(2,j) = v(2,j)*geom%dy(2,j)
    enddo
    do j=js,je+1
      wu(i,j) = u(i,j)*geom%dx(i,j)
    enddo
    do j=js,je
      utmp(i,j) = 2.*(wu(i,j) + wu(i,j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
      vtmp(i,j) = 2.*(wv(1,j) + wv(2,j  ))/(geom%dy(1,j)+geom%dy(2,j))
    enddo
  endif

  if ( (ie+1)==npx) then
    i = npx-1
    do j=js,je
      wv(i,  j) = v(i,  j)*geom%dy(i,  j)
      wv(i+1,j) = v(i+1,j)*geom%dy(i+1,j)
    enddo
    do j=js,je+1
      wu(i,j) = u(i,j)*geom%dx(i,j)
    enddo
    do j=js,je
      utmp(i,j) = 2.*(wu(i,j) + wu(i,  j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
      vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j  ))/(geom%dy(i,j)+geom%dy(i+1,j))
    enddo
  endif

endif

!Transform local a-grid winds into latitude-longitude coordinates
do j=js,je
  do i=is,ie
    ua(i,j) = geom%a11(i,j)*utmp(i,j) + geom%a12(i,j)*vtmp(i,j)
    va(i,j) = geom%a21(i,j)*utmp(i,j) + geom%a22(i,j)*vtmp(i,j)
  enddo
enddo

end subroutine d_to_a_domain_level

! --------------------------------------------------------------------------------------------------

subroutine d_to_a_ad(geom, u_ad_comp, v_ad_comp, ua_ad_comp, va_ad_comp)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(out) ::  u_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(out) ::  v_ad_comp(geom%isc:geom%iec+1,geom%jsc:geom%jec,  1:geom%npz)
real(kind=kind_real), intent(in)  :: ua_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec,  1:geom%npz)
real(kind=kind_real), intent(in)  :: va_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec,  1:geom%npz)

real(kind=kind_real) ::  u(geom%isd:geom%ied,  geom%jsd:geom%jed+1,geom%npz)
real(kind=kind_real) ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed,  geom%npz)
real(kind=kind_real) :: ua(geom%isd:geom%ied,  geom%jsd:geom%jed,  geom%npz)
real(kind=kind_real) :: va(geom%isd:geom%ied,  geom%jsd:geom%jed,  geom%npz)

real(kind=kind_real), save :: c1=1.125
real(kind=kind_real), save :: c2=-0.125
real(kind=kind_real) :: utmp(geom%isc:geom%iec, geom%jsc:geom%jec+1)
real(kind=kind_real) :: utmp_ad(geom%isc:geom%iec, geom%jsc:geom%jec+1)
real(kind=kind_real) :: vtmp(geom%isc:geom%iec+1, geom%jsc:geom%jec)
real(kind=kind_real) :: vtmp_ad(geom%isc:geom%iec+1, geom%jsc:geom%jec)
real(kind=kind_real) :: wu(geom%isc:geom%iec, geom%jsc:geom%jec+1)
real(kind=kind_real) :: wu_ad(geom%isc:geom%iec, geom%jsc:geom%jec+1)
real(kind=kind_real) :: wv(geom%isc:geom%iec+1, geom%jsc:geom%jec)
real(kind=kind_real) :: wv_ad(geom%isc:geom%iec+1, geom%jsc:geom%jec)
integer :: i, j, k
integer :: is, ie, js, je, npx, npy, npz
intrinsic max
intrinsic min
integer :: max1
integer :: max2
integer :: min1
integer :: min2
real(kind=kind_real) :: temp_ad
real(kind=kind_real) :: temp_ad0
real(kind=kind_real) :: temp_ad1
real(kind=kind_real) :: temp_ad2
real(kind=kind_real) :: temp_ad3
real(kind=kind_real) :: temp_ad4
real(kind=kind_real) :: temp_ad5
real(kind=kind_real) :: temp_ad6
integer :: ad_from
integer :: ad_to
integer :: ad_from0
integer :: ad_to0
integer :: branch

real(kind=kind_real) ::  u_ad(geom%isd:geom%ied,  geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) ::  v_ad(geom%isd:geom%ied+1,geom%jsd:geom%jed,  1:geom%npz)
real(kind=kind_real) :: ua_ad(geom%isd:geom%ied,  geom%jsd:geom%jed,  1:geom%npz)
real(kind=kind_real) :: va_ad(geom%isd:geom%ied,  geom%jsd:geom%jed,  1:geom%npz)

real(kind=kind_real) :: ebuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: nbufferx(geom%isc:geom%iec,1:geom%npz)
real(kind=kind_real) :: wbuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: sbufferx(geom%isc:geom%iec,1:geom%npz)

is = geom%isc
ie = geom%iec
js = geom%jsc
je = geom%jec
npx = geom%npx
npy = geom%npy
npz = geom%npz

ua_ad = 0.0_kind_real
va_ad = 0.0_kind_real
u_ad = 0.0_kind_real
v_ad = 0.0_kind_real

ua_ad(is:ie,js:je,:) = ua_ad_comp
va_ad(is:ie,js:je,:) = va_ad_comp

! TODO: OpenMP Re-enable this parallel section by fixing variable shared/private assignment
! !$omp parallel do default(none) shared(is,ie,js,je,npz,npx,npy,c2,c1, &
! !$omp                                  u,v,ua,va)         &
! !$omp                          private(utmp, vtmp, wu, wv)
do k=1,npz
  if (2 .lt. js) then
    max1 = js
  else
    max1 = 2
  end if
  if (npy - 2 .gt. je) then
    min1 = je
  else
    min1 = npy - 2
  end if
  ad_from0 = max1
  do j=ad_from0,min1
    if (2 .lt. is) then
      max2 = is
    else
      max2 = 2
    end if
    if (npx - 2 .gt. ie) then
      min2 = ie
    else
      min2 = npx - 2
    end if
    ad_from = max2
    call pushinteger4(i)
    i = min2 + 1
    call pushinteger4(i - 1)
    call pushinteger4(ad_from)
  end do
  call pushinteger4(j - 1)
  call pushinteger4(ad_from0)
  if (js .eq. 1) then
    call pushinteger4(i)
    call pushcontrol1b(0)
  else
    call pushcontrol1b(1)
  end if
  if (je + 1 .eq. npy) then
    j = npy - 1
    call pushinteger4(i)
    call pushcontrol1b(0)
  else
    call pushcontrol1b(1)
  end if
  if (is .eq. 1) then
    call pushinteger4(i)
    i = 1
    call pushinteger4(j)
    call pushcontrol1b(0)
  else
    call pushcontrol1b(1)
  end if
  if (ie + 1 .eq. npx) then
    call pushinteger4(i)
    i = npx - 1
    call pushinteger4(j)
    call pushcontrol1b(1)
  else
    call pushcontrol1b(0)
  end if
  call pushinteger4(j)
  do j=js,je
    call pushinteger4(i)
  end do
end do
vtmp_ad = 0.0_kind_real
wu_ad = 0.0_kind_real
wv_ad = 0.0_kind_real
utmp_ad = 0.0_kind_real
do k=npz,1,-1
  do j=je,js,-1
    do i=ie,is,-1
      utmp_ad(i, j) = utmp_ad(i, j) + geom%a11(i, j)*ua_ad(i, j, k) + geom%a21(i, j)*va_ad(i, j, k)
      vtmp_ad(i, j) = vtmp_ad(i, j) + geom%a12(i, j)*ua_ad(i, j, k) + geom%a22(i, j)*va_ad(i, j, k)
      va_ad(i, j, k) = 0.0_kind_real
      ua_ad(i, j, k) = 0.0_kind_real
    end do
    call popinteger4(i)
  end do
  call popinteger4(j)
  call popcontrol1b(branch)
  if (branch .ne. 0) then
    do j=je,js,-1
      temp_ad5 = 2.*vtmp_ad(i, j)/(geom%dy(i, j)+geom%dy(i+1, j))
      wv_ad(i, j) = wv_ad(i, j) + temp_ad5
      wv_ad(i+1, j) = wv_ad(i+1, j) + temp_ad5
      vtmp_ad(i, j) = 0.0_kind_real
      temp_ad6 = 2.*utmp_ad(i, j)/(geom%dx(i, j)+geom%dx(i, j+1))
      wu_ad(i, j) = wu_ad(i, j) + temp_ad6
      wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad6
      utmp_ad(i, j) = 0.0_kind_real
    end do
    do j=je+1,js,-1
      u_ad(i, j, k) = u_ad(i, j, k) + geom%dx(i, j)*wu_ad(i, j)
      wu_ad(i, j) = 0.0_kind_real
    end do
    do j=je,js,-1
      v_ad(i+1, j, k) = v_ad(i+1, j, k) + geom%dy(i+1, j)*wv_ad(i+1, j)
      wv_ad(i+1, j) = 0.0_kind_real
      v_ad(i, j, k) = v_ad(i, j, k) + geom%dy(i, j)*wv_ad(i, j)
      wv_ad(i, j) = 0.0_kind_real
    end do
    call popinteger4(j)
    call popinteger4(i)
  end if
  call popcontrol1b(branch)
  if (branch .eq. 0) then
    do j=je,js,-1
      temp_ad3 = 2.*vtmp_ad(i, j)/(geom%dy(1, j)+geom%dy(2, j))
      wv_ad(1, j) = wv_ad(1, j) + temp_ad3
      wv_ad(2, j) = wv_ad(2, j) + temp_ad3
      vtmp_ad(i, j) = 0.0_kind_real
      temp_ad4 = 2.*utmp_ad(i, j)/(geom%dx(i, j)+geom%dx(i, j+1))
      wu_ad(i, j) = wu_ad(i, j) + temp_ad4
      wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad4
      utmp_ad(i, j) = 0.0_kind_real
    end do
    do j=je+1,js,-1
      u_ad(i, j, k) = u_ad(i, j, k) + geom%dx(i, j)*wu_ad(i, j)
      wu_ad(i, j) = 0.0_kind_real
    end do
    do j=je,js,-1
      v_ad(2, j, k) = v_ad(2, j, k) + geom%dy(2, j)*wv_ad(2, j)
      wv_ad(2, j) = 0.0_kind_real
      v_ad(1, j, k) = v_ad(1, j, k) + geom%dy(1, j)*wv_ad(1, j)
      wv_ad(1, j) = 0.0_kind_real
    end do
    call popinteger4(j)
    call popinteger4(i)
  end if
  call popcontrol1b(branch)
  if (branch .eq. 0) then
    do i=ie,is,-1
      temp_ad1 = 2.*utmp_ad(i, j)/(geom%dx(i, j)+geom%dx(i, j+1))
      u_ad(i, j, k) = u_ad(i, j, k) + geom%dx(i, j)*temp_ad1
      u_ad(i, j+1, k) = u_ad(i, j+1, k) + geom%dx(i, j+1)*temp_ad1
      utmp_ad(i, j) = 0.0_kind_real
      temp_ad2 = 2.*vtmp_ad(i, j)/(geom%dy(i, j)+geom%dy(i+1, j))
      wv_ad(i, j) = wv_ad(i, j) + temp_ad2
      wv_ad(i+1, j) = wv_ad(i+1, j) + temp_ad2
      vtmp_ad(i, j) = 0.0_kind_real
    end do
    do i=ie+1,is,-1
      v_ad(i, j, k) = v_ad(i, j, k) + geom%dy(i, j)*wv_ad(i, j)
      wv_ad(i, j) = 0.0_kind_real
    end do
    call popinteger4(i)
  end if
  call popcontrol1b(branch)
  if (branch .eq. 0) then
    do i=ie,is,-1
      temp_ad = 2.*utmp_ad(i, 1)/(geom%dx(i, 1)+geom%dx(i, 2))
      u_ad(i, 1, k) = u_ad(i, 1, k) + geom%dx(i, 1)*temp_ad
      u_ad(i, 2, k) = u_ad(i, 2, k) + geom%dx(i, 2)*temp_ad
      utmp_ad(i, 1) = 0.0_kind_real
      temp_ad0 = 2.*vtmp_ad(i, 1)/(geom%dy(i, 1)+geom%dy(i+1, 1))
      wv_ad(i, 1) = wv_ad(i, 1) + temp_ad0
      wv_ad(i+1, 1) = wv_ad(i+1, 1) + temp_ad0
      vtmp_ad(i, 1) = 0.0_kind_real
    end do
    do i=ie+1,is,-1
      v_ad(i, 1, k) = v_ad(i, 1, k) + geom%dy(i, 1)*wv_ad(i, 1)
      wv_ad(i, 1) = 0.0_kind_real
    end do
    call popinteger4(i)
  end if
  call popinteger4(ad_from0)
  call popinteger4(ad_to0)
  do j=ad_to0,ad_from0,-1
    call popinteger4(ad_from)
    call popinteger4(ad_to)
    do i=ad_to,ad_from,-1
      v_ad(i-1, j, k) = v_ad(i-1, j, k) + c2*vtmp_ad(i, j)
      v_ad(i+2, j, k) = v_ad(i+2, j, k) + c2*vtmp_ad(i, j)
      v_ad(i, j, k) = v_ad(i, j, k) + c1*vtmp_ad(i, j)
      v_ad(i+1, j, k) = v_ad(i+1, j, k) + c1*vtmp_ad(i, j)
      vtmp_ad(i, j) = 0.0_kind_real
      u_ad(i, j-1, k) = u_ad(i, j-1, k) + c2*utmp_ad(i, j)
      u_ad(i, j+2, k) = u_ad(i, j+2, k) + c2*utmp_ad(i, j)
      u_ad(i, j, k) = u_ad(i, j, k) + c1*utmp_ad(i, j)
      u_ad(i, j+1, k) = u_ad(i, j+1, k) + c1*utmp_ad(i, j)
      utmp_ad(i, j) = 0.0_kind_real
    end do
    call popinteger4(i)
  end do
end do

!Adjoint fill halos
call mpp_update_domains_adm(u, u_ad, v, v_ad, geom%domain, gridtype=dgrid_ne)

!Adjoint fill edges
if (.not. geom%bounded_domain) then

  ebuffery = 0.0_kind_real
  nbufferx = 0.0_kind_real
  wbuffery = 0.0_kind_real
  sbufferx = 0.0_kind_real

  do k=1,npz
    do i=is,ie
      nbufferx(i,k) = u_ad(i,je+1,k)
    enddo
  enddo
  do k=1,npz
    do j=js,je
      ebuffery(j,k) = v_ad(ie+1,j,k)
    enddo
  enddo

  call mpp_get_boundary_ad( u_ad, v_ad, geom%domain, &
                            wbuffery=wbuffery, ebuffery=ebuffery, sbufferx=sbufferx, nbufferx=nbufferx, &
                            gridtype=DGRID_NE, complete=.true. )

endif

!Output compute part
u_ad_comp = u_ad(is:ie  ,js:je+1,:)
v_ad_comp = v_ad(is:ie+1,js:je  ,:)

end subroutine d_to_a_ad

! --------------------------------------------------------------------------------------------------

subroutine c_to_t_domain_level(geom, uc, vc, ut, vt, dt)

 type(fv3jedi_geom),   intent(in)  :: geom
 real(kind=kind_real), intent(in)  ::  uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
 real(kind=kind_real), intent(in)  ::  vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
 real(kind=kind_real), intent(out) ::  ut(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
 real(kind=kind_real), intent(out) ::  vt(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
 real(kind=kind_real), intent(in)  :: dt

! Local vars
  integer :: is,js,ie,je,isd,jsd,ied,jed,npx,npy
  real(kind=kind_real) :: damp
  integer i,j

  is  = geom%isc
  js  = geom%jsc
  ie  = geom%iec
  je  = geom%jec
  isd = geom%isd
  jsd = geom%jsd
  ied = geom%ied
  jed = geom%jed
  npx = geom%npx
  npy = geom%npy

  if ( geom%grid_type < 3 ) then

! Center of Cube Faces
        do j=jsd,jed
           if(j/=0 .and. j/=1 .and. j/=(npy-1) .and. j/=npy) then
             do i=is-1,ie+2
                ut(i,j) = ( uc(i,j) - 0.25 * geom%cosa_u(i,j) *     &
                    (vc(i-1,j)+vc(i,j)+vc(i-1,j+1)+vc(i,j+1)))*geom%rsin_u(i,j)
             enddo
           endif
        enddo
        do j=js-1,je+2
           if( j/=1 .and. j/=npy ) then
              do i=isd,ied
                 vt(i,j) = ( vc(i,j) - 0.25 * geom%cosa_v(i,j) *     &
                    (uc(i,j-1)+uc(i+1,j-1)+uc(i,j)+uc(i+1,j)))*geom%rsin_v(i,j)
              enddo
           endif
        enddo

! West edge:
       if ( is==1 ) then
          do j=jsd,jed
             if ( uc(1,j)*dt > 0. ) then
                ut(1,j) = uc(1,j) / geom%sin_sg(0,j,3)
             else
                ut(1,j) = uc(1,j) / geom%sin_sg(1,j,1)
             endif
          enddo
          do j=max(3,js), min(npy-2,je+1)
             vt(0,j) = vc(0,j) - 0.25*geom%cosa_v(0,j)*   &
                  (ut(0,j-1)+ut(1,j-1)+ut(0,j)+ut(1,j))
             vt(1,j) = vc(1,j) - 0.25*geom%cosa_v(1,j)*   &
                  (ut(1,j-1)+ut(2,j-1)+ut(1,j)+ut(2,j))
          enddo
       endif   ! West face

! East edge:
       if ( (ie+1)==npx ) then
          do j=jsd,jed
             if ( uc(npx,j)*dt > 0. ) then
                ut(npx,j) = uc(npx,j) / geom%sin_sg(npx-1,j,3)
             else
                ut(npx,j) = uc(npx,j) / geom%sin_sg(npx,j,1)
             endif
          enddo

           do j=max(3,js), min(npy-2,je+1)
              vt(npx-1,j) = vc(npx-1,j) - 0.25*geom%cosa_v(npx-1,j)*   &
                           (ut(npx-1,j-1)+ut(npx,j-1)+ut(npx-1,j)+ut(npx,j))
              vt(npx,j) = vc(npx,j) - 0.25*geom%cosa_v(npx,j)*   &
                         (ut(npx,j-1)+ut(npx+1,j-1)+ut(npx,j)+ut(npx+1,j))
           enddo
       endif

! South (Bottom) edge:
       if ( js==1 ) then

           do i=isd,ied
              if ( vc(i,1)*dt > 0. ) then
                   vt(i,1) = vc(i,1) / geom%sin_sg(i,0,4)
              else
                   vt(i,1) = vc(i,1) / geom%sin_sg(i,1,2)
              endif
           enddo

           do i=max(3,is),min(npx-2,ie+1)
              ut(i,0) = uc(i,0) - 0.25*geom%cosa_u(i,0)*   &
                       (vt(i-1,0)+vt(i,0)+vt(i-1,1)+vt(i,1))
              ut(i,1) = uc(i,1) - 0.25*geom%cosa_u(i,1)*   &
                       (vt(i-1,1)+vt(i,1)+vt(i-1,2)+vt(i,2))
           enddo
       endif

! North edge:
       if ( (je+1)==npy ) then
           do i=isd,ied
              if ( vc(i,npy)*dt > 0. ) then
                   vt(i,npy) = vc(i,npy) / geom%sin_sg(i,npy-1,4)
              else
                   vt(i,npy) = vc(i,npy) / geom%sin_sg(i,npy,2)
              endif
           enddo
           do i=max(3,is),min(npx-2,ie+1)
              ut(i,npy-1) = uc(i,npy-1) - 0.25*geom%cosa_u(i,npy-1)*   &
                           (vt(i-1,npy-1)+vt(i,npy-1)+vt(i-1,npy)+vt(i,npy))
              ut(i,npy) = uc(i,npy) - 0.25*geom%cosa_u(i,npy)*   &
                         (vt(i-1,npy)+vt(i,npy)+vt(i-1,npy+1)+vt(i,npy+1))
           enddo
       endif

!The following code solves a 2x2 system to get the interior parallel-to-edge
! uc,vc values near the corners (ex: for the sw corner ut(2,1) and vt(1,2) are solved for simultaneously).
! It then computes the halo uc, vc values so as to be consistent with the computations on the facing panel.

       !The system solved is:
       !  ut(2,1) = uc(2,1) - avg(vt)*geom%cosa_u(2,1)
       !  vt(1,2) = vc(1,2) - avg(ut)*geom%cosa_v(1,2)
       ! in which avg(vt) includes vt(1,2) and avg(ut) includes ut(2,1)

        if( geom%sw_corner ) then
            damp = 1. / (1.-0.0625*geom%cosa_u(2,0)*geom%cosa_v(1,0))
            ut(2,0) = (uc(2,0)-0.25*geom%cosa_u(2,0)*(vt(1,1)+vt(2,1)+vt(2,0) +vc(1,0) -   &
                      0.25*geom%cosa_v(1,0)*(ut(1,0)+ut(1,-1)+ut(2,-1))) ) * damp
            damp = 1. / (1.-0.0625*geom%cosa_u(0,1)*geom%cosa_v(0,2))
            vt(0,2) = (vc(0,2)-0.25*geom%cosa_v(0,2)*(ut(1,1)+ut(1,2)+ut(0,2)+uc(0,1) -   &
                      0.25*geom%cosa_u(0,1)*(vt(0,1)+vt(-1,1)+vt(-1,2))) ) * damp

            damp = 1. / (1.-0.0625*geom%cosa_u(2,1)*geom%cosa_v(1,2))
            ut(2,1) = (uc(2,1)-0.25*geom%cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2)+vc(1,2) -   &
                      0.25*geom%cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2))) ) * damp

            vt(1,2) = (vc(1,2)-0.25*geom%cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2)+uc(2,1) -   &
                      0.25*geom%cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2))) ) * damp
        endif

        if( geom%se_corner ) then
            damp = 1. / (1. - 0.0625*geom%cosa_u(npx-1,0)*geom%cosa_v(npx-1,0))
            ut(npx-1,0) = ( uc(npx-1,0)-0.25*geom%cosa_u(npx-1,0)*(   &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,0)+vc(npx-1,0) -   &
                      0.25*geom%cosa_v(npx-1,0)*(ut(npx,0)+ut(npx,-1)+ut(npx-1,-1))) ) * damp
            damp = 1. / (1. - 0.0625*geom%cosa_u(npx+1,1)*geom%cosa_v(npx,2))
            vt(npx,  2) = ( vc(npx,2)-0.25*geom%cosa_v(npx,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx+1,2)+uc(npx+1,1) -   &
                      0.25*geom%cosa_u(npx+1,1)*(vt(npx,1)+vt(npx+1,1)+vt(npx+1,2))) ) * damp

            damp = 1. / (1. - 0.0625*geom%cosa_u(npx-1,1)*geom%cosa_v(npx-1,2))
            ut(npx-1,1) = ( uc(npx-1,1)-0.25*geom%cosa_u(npx-1,1)*(  &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2)+vc(npx-1,2) -   &
                      0.25*geom%cosa_v(npx-1,2)*(ut(npx,1)+ut(npx,2)+ut(npx-1,2))) ) * damp
            vt(npx-1,2) = ( vc(npx-1,2)-0.25*geom%cosa_v(npx-1,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx-1,2)+uc(npx-1,1) -   &
                      0.25*geom%cosa_u(npx-1,1)*(vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2))) ) * damp
        endif

        if( geom%ne_corner ) then
            damp = 1. / (1. - 0.0625*geom%cosa_u(npx-1,npy)*geom%cosa_v(npx-1,npy+1))
            ut(npx-1,npy) = ( uc(npx-1,npy)-0.25*geom%cosa_u(npx-1,npy)*(   &
                              vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy+1)+vc(npx-1,npy+1) -   &
                0.25*geom%cosa_v(npx-1,npy+1)*(ut(npx,npy)+ut(npx,npy+1)+ut(npx-1,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*geom%cosa_u(npx+1,npy-1)*geom%cosa_v(npx,npy-1))
            vt(npx,  npy-1) = ( vc(npx,npy-1)-0.25*geom%cosa_v(npx,npy-1)*(   &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx+1,npy-2)+uc(npx+1,npy-1) -   &
                0.25*geom%cosa_u(npx+1,npy-1)*(vt(npx,npy)+vt(npx+1,npy)+vt(npx+1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*geom%cosa_u(npx-1,npy-1)*geom%cosa_v(npx-1,npy-1))
            ut(npx-1,npy-1) = ( uc(npx-1,npy-1)-0.25*geom%cosa_u(npx-1,npy-1)*(  &
                                vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1)+vc(npx-1,npy-1) -  &
                0.25*geom%cosa_v(npx-1,npy-1)*(ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2))) ) * damp
            vt(npx-1,npy-1) = ( vc(npx-1,npy-1)-0.25*geom%cosa_v(npx-1,npy-1)*(  &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2)+uc(npx-1,npy-1) -  &
                0.25*geom%cosa_u(npx-1,npy-1)*(vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1))) ) * damp
        endif

        if( geom%nw_corner ) then
            damp = 1. / (1. - 0.0625*geom%cosa_u(2,npy)*geom%cosa_v(1,npy+1))
            ut(2,npy) = ( uc(2,npy)-0.25*geom%cosa_u(2,npy)*(   &
                          vt(1,npy)+vt(2,npy)+vt(2,npy+1)+vc(1,npy+1) -   &
                      0.25*geom%cosa_v(1,npy+1)*(ut(1,npy)+ut(1,npy+1)+ut(2,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*geom%cosa_u(0,npy-1)*geom%cosa_v(0,npy-1))
            vt(0,npy-1) = ( vc(0,npy-1)-0.25*geom%cosa_v(0,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(0,npy-2)+uc(0,npy-1) -   &
                      0.25*geom%cosa_u(0,npy-1)*(vt(0,npy)+vt(-1,npy)+vt(-1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*geom%cosa_u(2,npy-1)*geom%cosa_v(1,npy-1))
            ut(2,npy-1) = ( uc(2,npy-1)-0.25*geom%cosa_u(2,npy-1)*(  &
                            vt(1,npy)+vt(2,npy)+vt(2,npy-1)+vc(1,npy-1) -   &
                      0.25*geom%cosa_v(1,npy-1)*(ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2))) ) * damp

            vt(1,npy-1) = ( vc(1,npy-1)-0.25*geom%cosa_v(1,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2)+uc(2,npy-1) -   &
                      0.25*geom%cosa_u(2,npy-1)*(vt(1,npy)+vt(2,npy)+vt(2,npy-1))) ) * damp
        endif

 else
! grid_type >= 3

        do j=jsd,jed
           do i=is,ie+1
              ut(i,j) =  uc(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              vt(i,j) = vc(i,j)
           enddo
        enddo
 endif

end subroutine c_to_t_domain_level

! --------------------------------------------------------------------------------------------------

subroutine c_to_acs(geom, uc_in, vc_in, ua_out, va_out, interp_order_in)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  ::  uc_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(in)  ::  vc_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(out) :: ua_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(out) :: va_out(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
integer, optional,    intent(in)  :: interp_order_in

integer :: k
real(kind=kind_real) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )

! Fill halo for input (outside loop for communication efficiency)
! ---------------------------------------------------------------
uc(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) = uc_in
vc(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) = vc_in
call fill_cgrid_winds(geom, uc, vc, fillhalo=.true.)


! Call domain based calculation
! -----------------------------
do k = 1, geom%npz
  call c_to_acs_domain_level(geom, uc(:,:,k), vc(:,:,k), ua, va, interp_order_in)

  ! Fill outputs
  ua_out(:,:,k) = ua(geom%isc:geom%iec,geom%jsc:geom%jec)
  va_out(:,:,k) = va(geom%isc:geom%iec,geom%jsc:geom%jec)
enddo

end subroutine c_to_acs

! --------------------------------------------------------------------------------------------------

subroutine c_to_acs_domain_level(geom, uc, vc, ua, va, interp_order_in)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
real(kind=kind_real), intent(in)  :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real), intent(out) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ) ! A-grid (cubed-sphere)
real(kind=kind_real), intent(out) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ) ! A-grid (cubed-sphere)
integer, optional,    intent(in)  :: interp_order_in

integer :: interp_order
integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
real(kind=kind_real) :: tmp1i(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp2i(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp3i(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp1j(geom%jsd:geom%jed+1)
real(kind=kind_real) :: tmp2j(geom%jsd:geom%jed+1)
real(kind=kind_real) :: tmp3j(geom%jsd:geom%jed+1)

interp_order = 1
if (present(interp_order_in)) interp_order = interp_order_in

is  = geom%isc
ie  = geom%iec
js  = geom%jsc
je  = geom%jec
isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed

! va
do i=isd,ied
  tmp1j(:) = 0.0
  tmp2j(:) = vc(i,:)*geom%dx(i,:)
  tmp3j(:) = geom%dyc(i,:)
  call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed+1, interp_order)
  va(i,jsd:jed) = tmp1j(jsd+1:jed+1)/geom%dxa(i,jsd:jed)
enddo

! ua
do j=jsd,jed
  tmp1i(:) = 0.0
  tmp2i(:) = uc(:,j)*geom%dy(:,j)
  tmp3i(:) = geom%dxc(:,j)
  call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied+1, interp_order)
  ua(isd:ied,j) = tmp1i(isd+1:ied+1)/geom%dya(isd:ied,j)
enddo

end subroutine c_to_acs_domain_level

! --------------------------------------------------------------------------------------------------

subroutine c_to_acs_domain_level_adm(geom, uc_ad, vc_ad, ua_ad, va_ad)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: uc_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed)
real(kind=kind_real), intent(inout) :: vc_ad(geom%isd:geom%ied, geom%jsd:geom%jed+1)
real(kind=kind_real), intent(inout) :: ua_ad(geom%isd:geom%ied, geom%jsd:geom%jed)
real(kind=kind_real), intent(inout) :: va_ad(geom%isd:geom%ied, geom%jsd:geom%jed)

! Locals
integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, iedp1, jedp1
real(kind=kind_real) :: tmp1i_ad(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp2i_ad(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp1j_ad(geom%jsd:geom%jed+1)
real(kind=kind_real) :: tmp2j_ad(geom%jsd:geom%jed+1)

isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed
jedp1 = jed + 1
iedp1 = ied + 1

do j=jed,jsd,-1
  tmp1i_ad = 0.0_kind_real
  tmp1i_ad(isd+1:ied+1) = tmp1i_ad(isd+1:ied+1) + ua_ad(isd:ied, j)/geom%dya(isd:ied, j)
  ua_ad(isd:ied, j) = 0.0_kind_real
  call interp_left_edge_1d_adm(isd, iedp1, iedp1, tmp1i_ad, tmp2i_ad)
  uc_ad(:, j) = uc_ad(:, j) + geom%dy(:, j)*tmp2i_ad
end do

do i=ied,isd,-1
  tmp1j_ad = 0.0_kind_real
  tmp1j_ad(jsd+1:jed+1) = tmp1j_ad(jsd+1:jed+1) + va_ad(i, jsd:jed)/geom%dxa(i, jsd:jed)
  va_ad(i, jsd:jed) = 0.0_kind_real
  call interp_left_edge_1d_adm(jsd, jedp1, jedp1, tmp1j_ad, tmp2j_ad)
  vc_ad(i, :) = vc_ad(i, :) + geom%dx(i, :)*tmp2j_ad
end do

end subroutine c_to_acs_domain_level_adm

! --------------------------------------------------------------------------------------------------

subroutine acs_to_d(geom, ua_in, va_in, ud_out, vd_out, interp_order_in)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  ::  ua_in(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(in)  ::  va_in(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz)
real(kind=kind_real), intent(out) :: ud_out(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(out) :: vd_out(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)
integer, optional,    intent(in)  :: interp_order_in

integer :: k
real(kind=kind_real) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real) :: ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real) :: vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  )

! Fill halo for input (outside loop for communication efficiency)
! ---------------------------------------------------------------
ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = ua_in
va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = va_in
call fill_agrid_winds(geom, ua, va)


! Call domain based calculation
! -----------------------------
do k = 1, geom%npz
  call acs_to_d_domain_level(geom, ua(:,:,k), va(:,:,k), ud, vd, interp_order_in)

  ! Fill outputs
  ud_out(:,:,k) = ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1)
  vd_out(:,:,k) = vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  )
enddo

end subroutine acs_to_d

! --------------------------------------------------------------------------------------------------

subroutine acs_to_d_domain_level(geom, ua, va, ud, vd, interp_order_in)

type(fv3jedi_geom),   intent(in)  :: geom
real(kind=kind_real), intent(in)  :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ) ! A-Grid (cubed-sphere)
real(kind=kind_real), intent(in)  :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ) ! A-Grid (cubed-sphere)
real(kind=kind_real), intent(out) :: ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
real(kind=kind_real), intent(out) :: vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
integer, optional,    intent(in)  :: interp_order_in

! Locals
integer :: interp_order
integer :: i, j, isd, ied, jsd, jed
real(kind=kind_real) :: tmp1i(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp2i(geom%isd:geom%ied)
real(kind=kind_real) :: tmp3i(geom%isd:geom%ied)
real(kind=kind_real) :: tmp1j(geom%jsd:geom%jed+1)
real(kind=kind_real) :: tmp2j(geom%jsd:geom%jed)
real(kind=kind_real) :: tmp3j(geom%jsd:geom%jed)

interp_order = 1
if (present(interp_order_in)) interp_order = interp_order_in

isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed

do j=jsd+1,jed
  tmp1i(:) = 0.0
  tmp2i(:) = va(:,j)*geom%dxa(:,j)
  tmp3i(:) = geom%dxa(:,j)
  call interp_left_edge_1d(tmp1i, tmp2i, tmp3i, isd, ied, interp_order)
  vd(:,j) = tmp1i(:)/geom%dxc(:,j)
enddo

do i=isd+1,ied
  tmp1j(:) = 0.0
  tmp2j(:) = ua(i,:)*geom%dya(i,:)
  tmp3j(:) = geom%dya(i,:)
  call interp_left_edge_1d(tmp1j, tmp2j, tmp3j, jsd, jed, interp_order)
  ud(i,:) = tmp1j(:)/geom%dyc(i,:)
enddo

call mpp_update_domains( ud, vd, geom%domain, gridtype=DGRID_NE, complete=.true.)
!if (.not. geom%bounded_domain) call fill_corners(ud, vd, geom%npx, geom%npy, VECTOR=.true., DGRID=.true.)

end subroutine acs_to_d_domain_level

! --------------------------------------------------------------------------------------------------

subroutine acs_to_d_domain_level_adm(geom, ua_ad, va_ad, ud_ad, vd_ad)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: ua_ad(geom%isd:geom%ied, geom%jsd:geom%jed)
real(kind=kind_real), intent(inout) :: va_ad(geom%isd:geom%ied, geom%jsd:geom%jed)
real(kind=kind_real), intent(inout) :: ud_ad(geom%isd:geom%ied, geom%jsd:geom%jed+1)
real(kind=kind_real), intent(inout) :: vd_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed)

! Locals
integer :: i, j, isd, ied, jsd, jed, iedp1, jedp1
real(kind=kind_real) :: tmp1i_ad(geom%isd:geom%ied+1)
real(kind=kind_real) :: tmp2i_ad(geom%isd:geom%ied)
real(kind=kind_real) :: tmp1j_ad(geom%jsd:geom%jed+1)
real(kind=kind_real) :: tmp2j_ad(geom%jsd:geom%jed)

isd = geom%isd
ied = geom%ied
jsd = geom%jsd
jed = geom%jed
jedp1 = jed+1
iedp1 = ied+1

call mpp_update_domains_adm(ud_ad, ud_ad, vd_ad, vd_ad, geom%domain, gridtype=dgrid_ne, complete=.true.)

do i=ied,isd+1,-1
  tmp1j_ad = 0.0_kind_real
  tmp1j_ad(:) = ud_ad(i, :)/geom%dyc(i, :)
  ud_ad(i, :) = 0.0_kind_real
  call interp_left_edge_1d_adm(jsd, jed, jedp1, tmp1j_ad, tmp2j_ad)
  ua_ad(i, :) = ua_ad(i, :) + geom%dya(i, :)*tmp2j_ad
end do

do j=jed,jsd+1,-1
  tmp1i_ad = 0.0_kind_real
  tmp1i_ad(:) = vd_ad(:, j)/geom%dxc(:, j)
  vd_ad(:, j) = 0.0_kind_real
  call interp_left_edge_1d_adm(isd, ied, iedp1, tmp1i_ad, tmp2i_ad)
  va_ad(:, j) = va_ad(:, j) + geom%dxa(:, j)*tmp2i_ad
end do

end subroutine acs_to_d_domain_level_adm

! --------------------------------------------------------------------------------------------------

subroutine fill_agrid_winds(geom, ua, va)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: ua(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)
real(kind=kind_real), intent(inout) :: va(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)

call mpp_update_domains(ua, geom%domain, complete=.true.)
call mpp_update_domains(va, geom%domain, complete=.true.)

end subroutine fill_agrid_winds

! --------------------------------------------------------------------------------------------------

subroutine fill_dgrid_winds(geom,u,v,fillhalo)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
real(kind=kind_real), intent(inout) :: v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
logical, optional,    intent(in)    :: fillhalo

integer :: isc, iec, jsc, jec, npz, i, j, k
real(kind=kind_real) :: ebuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: nbufferx(geom%isc:geom%iec,1:geom%npz)
real(kind=kind_real) :: wbuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: sbufferx(geom%isc:geom%iec,1:geom%npz)

! ---------------------------------------- !
! Fill edge and then halo of  D-grid winds !
! ---------------------------------------- !

! Shortcuts
! ---------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

! Fill north/east edges
! ---------------------
if (.not. geom%bounded_domain) then

  ebuffery = 0.0_kind_real
  nbufferx = 0.0_kind_real
  wbuffery = 0.0_kind_real
  sbufferx = 0.0_kind_real
  call mpp_get_boundary( u, v, geom%domain, &
                         wbuffery=wbuffery, ebuffery=ebuffery, &
                         sbufferx=sbufferx, nbufferx=nbufferx, &
                         gridtype=DGRID_NE, complete=.true. )
  do k=1,npz
     do i=isc,iec
        u(i,jec+1,k) = nbufferx(i,k)
     enddo
  enddo
  do k=1,npz
     do j=jsc,jec
        v(iec+1,j,k) = ebuffery(j,k)
     enddo
  enddo

endif

! Fill halos
! ----------
if (present(fillhalo)) then
  if (fillhalo) then
    call mpp_update_domains(u, v, geom%domain, gridtype=DGRID_NE)
  endif
endif

end subroutine fill_dgrid_winds

! --------------------------------------------------------------------------------------------------

subroutine fill_dgrid_winds_adm(geom, u, u_ad, v, v_ad, fillhalo)

  type(fv3jedi_geom), intent(in)      :: geom
  real(kind=kind_real), intent(inout) ::    u(geom%isd:geom%ied  , geom%jsd:geom%jed+1, geom%npz)
  real(kind=kind_real), intent(inout) :: u_ad(geom%isd:geom%ied  , geom%jsd:geom%jed+1, geom%npz)
  real(kind=kind_real), intent(inout) ::    v(geom%isd:geom%ied+1, geom%jsd:geom%jed  , geom%npz)
  real(kind=kind_real), intent(inout) :: v_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed  , geom%npz)
  logical, optional, intent(in) :: fillhalo

  integer :: isc, iec, jsc, jec, npz, i, j, k
  real(kind=kind_real) :: ebuffery(geom%jsc:geom%jec, geom%npz)
  real(kind=kind_real) :: nbufferx(geom%isc:geom%iec, geom%npz)
  real(kind=kind_real) :: wbuffery(geom%jsc:geom%jec, geom%npz)
  real(kind=kind_real) :: sbufferx(geom%isc:geom%iec, geom%npz)

  isc = geom%isc
  iec = geom%iec
  jsc = geom%jsc
  jec = geom%jec
  npz = geom%npz

  if (present(fillhalo)) then
    if (fillhalo) call mpp_update_domains_adm(u, u_ad, v, v_ad, geom%domain, gridtype=dgrid_ne)
  end if

  if (.not. geom%bounded_domain) then

    wbuffery = 0.0_kind_real
    sbufferx = 0.0_kind_real
    ebuffery = 0.0_kind_real
    nbufferx = 0.0_kind_real

    do k=npz,1,-1
      do j=jec,jsc,-1
        ebuffery(j, k) = v_ad(iec+1, j, k)
      end do
    end do
    do k=npz,1,-1
      do i=iec,isc,-1
        nbufferx(i, k) = u_ad(i, jec+1, k)
      end do
    end do

    call mpp_get_boundary_ad(u_ad, v_ad, geom%domain, wbuffery=wbuffery, ebuffery=ebuffery, &
                             sbufferx=sbufferx, nbufferx=nbufferx, gridtype=dgrid_ne, complete=.true.)

  endif

end subroutine fill_dgrid_winds_adm

! --------------------------------------------------------------------------------------------------

subroutine fill_cgrid_winds(geom, uc, vc, fillhalo)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
real(kind=kind_real), intent(inout) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
logical, optional,    intent(in)    :: fillhalo

integer :: isc, iec, jsc, jec, npz, i, j, k
real(kind=kind_real) :: wbufferx(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: sbuffery(geom%isc:geom%iec,1:geom%npz)
real(kind=kind_real) :: ebufferx(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: nbuffery(geom%isc:geom%iec,1:geom%npz)

! ---------------------------------------- !
! Fill edge and then halo of  C-grid winds !
! ---------------------------------------- !

! Shortcuts
! ---------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

! Fill north/east edges
! ---------------------
if (.not. geom%bounded_domain) then

  ebufferx = 0.0_kind_real
  nbuffery = 0.0_kind_real
  wbufferx = 0.0_kind_real
  sbuffery = 0.0_kind_real

  call mpp_get_boundary(uc, vc, geom%domain, &
                        wbufferx=wbufferx, ebufferx=ebufferx, &
                        sbuffery=sbuffery, nbuffery=nbuffery, &
                        gridtype=CGRID_NE, complete=.true.)
  do k=1,npz
     do j=jsc,jec
        uc(iec+1,j,k) = ebufferx(j,k)
     enddo
     do i=isc,iec
        vc(i,jec+1,k) = nbuffery(i,k)
     enddo
  enddo

endif

! Fill halos
! ----------
if (present(fillhalo)) then
  if (fillhalo) then
    call mpp_update_domains(uc, vc, geom%domain, gridtype=CGRID_NE)
  endif
endif

end subroutine fill_cgrid_winds

! --------------------------------------------------------------------------------------------------

subroutine fill_cgrid_winds_adm(geom, uc_ad, vc_ad, fillhalo)

type(fv3jedi_geom),   intent(in)    :: geom
real(kind=kind_real), intent(inout) :: uc_ad(geom%isd:geom%ied+1, geom%jsd:geom%jed, geom%npz)
real(kind=kind_real), intent(inout) :: vc_ad(geom%isd:geom%ied, geom%jsd:geom%jed+1, geom%npz)
logical, optional,    intent(in)    :: fillhalo

integer :: isc, iec, jsc, jec, npz, i, j, k
real(kind=kind_real) :: wbufferx(geom%jsc:geom%jec, geom%npz)
real(kind=kind_real) :: sbuffery(geom%isc:geom%iec, geom%npz)
real(kind=kind_real) :: ebufferx(geom%jsc:geom%jec, geom%npz)
real(kind=kind_real) :: nbuffery(geom%isc:geom%iec, geom%npz)

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

if (present(fillhalo)) then
  if (fillhalo) call mpp_update_domains_adm(uc_ad, uc_ad, vc_ad, vc_ad, geom%domain, gridtype=cgrid_ne)
end if

if (.not. geom%bounded_domain) then

  sbuffery = 0.0_kind_real
  wbufferx = 0.0_kind_real

  do k=npz,1,-1
    do i=iec,isc,-1
      nbuffery(i, k) = vc_ad(i,jec+1,k)
    end do
    do j=jec,jsc,-1
      ebufferx(j, k) = uc_ad(iec+1,j,k)
    end do
  end do

  call mpp_get_boundary_ad(uc_ad, vc_ad, geom%domain, sbuffery=sbuffery, nbuffery=nbuffery, &
                           wbufferx=wbufferx, ebufferx=ebufferx, gridtype=cgrid_ne, complete=.true.)

endif

end subroutine fill_cgrid_winds_adm

! --------------------------------------------------------------------------------------------------

subroutine interp_left_edge_1d(qout, qin, dx, ifirst, ilast, order)

real(kind=kind_real), intent(out) :: qout(ifirst:)
real(kind=kind_real), intent(in)  ::  qin(ifirst:)
real(kind=kind_real), intent(in)  ::   dx(ifirst:)
integer,              intent(in)  :: ifirst
integer,              intent(in)  :: ilast
integer,              intent(in)  :: order

integer :: i
real(kind=kind_real) :: dm(ifirst:ilast),qmax,qmin
real(kind=kind_real) :: r3, da1, da2, a6da, a6, al, ar
real(kind=kind_real) :: qLa, qLb1, qLb2
real(kind=kind_real) :: x

r3 = 1.0_kind_real/3.0_kind_real

qout(:) = 0.0_kind_real

if (order==1) then

  ! 1st order Uniform linear averaging
  do i=ifirst+1,ilast
    qout(i) = 0.5_kind_real * (qin(i-1) + qin(i))
  enddo

elseif (order==2) then

  ! Non-Uniform 1st order average
  do i=ifirst+1,ilast
    qout(i) = (dx(i-1)*qin(i-1) + dx(i)*qin(i))/(dx(i-1)+dx(i))
  enddo

elseif (order==3) then

  ! PPM - Uniform
   do i=ifirst+1,ilast-1
      dm(i) = 0.25_kind_real*(qin(i+1) - qin(i-1))
   enddo

  ! Applies monotonic slope constraint
  do i=ifirst+1,ilast-1
    qmax = max(qin(i-1),qin(i),qin(i+1)) - qin(i)
    qmin = qin(i) - min(qin(i-1),qin(i),qin(i+1))
    dm(i) = sign(min(abs(dm(i)),qmin,qmax),dm(i))
  enddo

  do i=ifirst+1,ilast-1
    qout(i) = 0.5_kind_real*(qin(i-1)+qin(i)) + r3*(dm(i-1) - dm(i))
  enddo

  ! First order average to fill in end points
  qout(ifirst+1) = 0.5_kind_real * (qin(ifirst) + qin(ifirst+1))
  qout(ilast) = 0.5_kind_real * (qin(ilast-1) + qin(ilast))

elseif (order==4) then

  ! Non-Uniform PPM
  do i=ifirst+1,ilast-1
    dm(i) = ( (2.0_kind_real*dx(i-1) + dx(i) ) /                         &
            (   dx(i+1) + dx(i) )  )  * ( qin(i+1) - qin(i) ) + &
            ( (dx(i)   + 2.0_kind_real*dx(i+1)) /                        &
            (dx(i-1) +    dx(i)  )  ) * ( qin(i) - qin(i-1) )
    dm(i) = ( dx(i) / ( dx(i-1) + dx(i) + dx(i+1) ) ) * dm(i)
    if ( (qin(i+1)-qin(i))*(qin(i)-qin(i-1)) > 0.0_kind_real) then
      dm(i) = SIGN( MIN( ABS(dm(i)), 2.0_kind_real*ABS(qin(i)-qin(i-1)), 2.0_kind_real*ABS(qin(i+1)-qin(i)) ) , dm(i) )
    else
      dm(i) = 0.
    endif
  enddo

  do i=ifirst+2,ilast-1
    qLa = ( (dx(i-2) + dx(i-1)) / (2.0_kind_real*dx(i-1) +  dx(i)) ) - &
          ( (dx(i+1) + dx(i)) / (2.0_kind_real*dx(i) +  dx(i-1)) )
    qLa = ( (2.0_kind_real*dx(i) * dx(i-1))  / (dx(i-1) + dx(i)) ) * qLa * &
          (qin(i) - qin(i-1))
    qLb1 = dx(i-1) * ( (dx(i-2) + dx(i-1)) / (2.0_kind_real*dx(i-1) + dx(i)) ) * &
           dm(i)
    qLb2 = dx(i) * ( (dx(i) + dx(i+1)) / (dx(i-1) + 2.0_kind_real*dx(i)) ) * &
           dm(i-1)

    qout(i) = 1. / ( dx(i-2) + dx(i-1) + dx(i) + dx(i+1) )
    qout(i) = qout(i) * ( qLa - qLb1 + qLb2 )
    qout(i) = qin(i-1) + ( dx(i-1) / ( dx(i-1) + dx(i) ) ) * (qin(i) - qin(i-1)) + qout(i)
  enddo

endif

end subroutine interp_left_edge_1d

! --------------------------------------------------------------------------------------------------

subroutine interp_left_edge_1d_adm(ifirst, ilast, ilastout, qout_ad, qin_ad)
  implicit none
  integer,              intent(in)    :: ifirst
  integer,              intent(in)    :: ilast
  integer,              intent(in)    :: ilastout
  real(kind=kind_real), intent(inout) :: qout_ad(ifirst:ilastout)
  real(kind=kind_real), intent(out)   :: qin_ad(ifirst:ilast)

  integer :: i

  qin_ad = 0.0_kind_real
  do i=ilast,ifirst+1,-1
    qin_ad(i-1) = qin_ad(i-1) + 0.5_kind_real*qout_ad(i)
    qin_ad(i) = qin_ad(i) + 0.5_kind_real*qout_ad(i)
    qout_ad(i) = 0.0_kind_real
  end do

end subroutine interp_left_edge_1d_adm

! --------------------------------------------------------------------------------------------------

real(kind=kind_real) function edge_interpolate4(ua, dxa)

real(kind=kind_real), intent(in) :: ua(4)
real(kind=kind_real), intent(in) :: dxa(4)
real(kind=kind_real):: t1, t2

t1 = dxa(1) + dxa(2)
t2 = dxa(3) + dxa(4)

edge_interpolate4 = 0.5*( ((t1+dxa(2))*ua(2)-dxa(2)*ua(1)) / t1 + &
                          ((t2+dxa(3))*ua(3)-dxa(3)*ua(4)) / t2 )

end function edge_interpolate4

! --------------------------------------------------------------------------------------------------

end module wind_vt_mod
