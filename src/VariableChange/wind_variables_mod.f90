! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wind_vt_mod

use fv3jedi_constants_mod, only: pi, rad2deg
use fv3jedi_geom_mod,  only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

use mpi
use fv_mp_adm_mod, only: mpp_update_domains_adm
use mpp_domains_mod, only: mpp_update_domains, dgrid_ne
use mpp_domains_mod, only: mpp_get_boundary, mpp_get_boundary_ad

implicit none
private

public sfc_10m_winds
public uv_to_vortdivg
public vortdivg_to_psichi
public gauss_seidel
public psichi_to_uava
public psichi_to_uava_adm
public psichi_to_udvd
public d_to_ac
public d2a2c_vect
public divergence_corner
public a2b_ord4
public PSICHI_TO_UDVD_ADM
public A2B_ORD4_ADM
public EXTRAP_CORNER_ADM
public a2d
public a2d_ad
public d2a
public D2A_AD

contains

!----------------------------------------------------------------------------
! Lowest model level winds to 10m speed and direction -----------------------
!----------------------------------------------------------------------------

subroutine sfc_10m_winds(geom,usrf,vsrf,f10r,spd10m,dir10m)

 implicit none

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
 real(kind=kind_real), parameter :: quadcof(4,2) = reshape((/ 0.0_kind_real,  1.0_kind_real, 1.0_kind_real,  2.0_kind_real,  &
                                                              1.0_kind_real, -1.0_kind_real, 1.0_kind_real, -1.0_kind_real /), (/4, 2/))

 !In GSI these calculations are done after interpolation to obs location

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

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

!----------------------------------------------------------------------------

subroutine uv_to_vortdivg(geom,u,v,ua,va,vort,divg)

 implicit none
 type(fv3jedi_geom),   intent(inout)  :: geom
 real(kind=kind_real), intent(inout) ::    u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::    v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real), intent(inout) :: vort(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Vorticity
 real(kind=kind_real), intent(inout) :: divg(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Divergence

 real(kind=kind_real), intent(inout) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real), intent(inout) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
 real(kind=kind_real) :: ut(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
 real(kind=kind_real) :: vt(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)

 real(kind=kind_real) :: wbuffer(geom%jsc:geom%jec,geom%npz)
 real(kind=kind_real) :: sbuffer(geom%isc:geom%iec,geom%npz)
 real(kind=kind_real) :: ebuffer(geom%jsc:geom%jec,geom%npz)
 real(kind=kind_real) :: nbuffer(geom%isc:geom%iec,geom%npz)

 integer :: i,j,k,isc,iec,jsc,jec,npz
 real(kind=kind_real) :: dt

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !  vd,uc|        x        |
 !       |   vort,divg     |
 !       |                 |
 !       |                 |
 !       x-----------------x
 !              ud,vc

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec
 npz = geom%npz

 ua = 0.0_kind_real
 va = 0.0_kind_real


 !Fill the halos of Dgrid winds and get A/C grid winds
 !----------------------------------------------------
 call d_to_ac(geom,u,v,ua,va,uc,vc)


 !Calculate vorticity (Agrid)
 !---------------------------
! vort = 0.0_kind_real
! do k=1,npz
!   do j=jsc,jec
!     do i=isc,iec
!       vort(i,j,k) = geom%rarea(i,j)*(u(i,j,k)*geom%dx(i,j)-u(i,j+1,k)*geom%dx(i,j+1)  - &
!                                      v(i,j,k)*geom%dy(i,j)+v(i+1,j,k)*geom%dy(i+1,j))
!     enddo
!   enddo
! enddo

 vort = 0.0_kind_real
 do k=1,npz
   do j=jsc,jec
     do i=isc,iec
       vort(i,j,k) = (va(i+1,j,k) - va(i-1,j,k))/(geom%dxc(i,j)+geom%dxc(i-1,j)) - (ua(i,j+1,k) - ua(i,j-1,k))/(geom%dyc(i,j)+geom%dyc(i,j-1))
     enddo
   enddo
 enddo

 !Calculate divergence (Agrid)
 !----------------------------
 !Using C grid winds to give divergence at the grid center
! divg = 0.0_kind_real
! do k=1,npz
!   do j=jsc,jec
!     do i=isc,iec
!       divg(i,j,k) = geom%rarea(i,j)*(uc(i,j,k)*geom%dy(i,j)-uc(i,j+1,k)*geom%dy(i,j+1)  + &
!                                      vc(i,j,k)*geom%dx(i,j)+vc(i+1,j,k)*geom%dx(i+1,j))
!     enddo
!   enddo
! enddo

 divg = 0.0_kind_real
 do k=1,npz
   do j=jsc,jec
     do i=isc,iec
       divg(i,j,k) = (ua(i+1,j,k) - ua(i-1,j,k))/(geom%dxc(i,j)+geom%dxc(i-1,j)) + (va(i,j+1,k) - va(i,j-1,k))/(geom%dyc(i,j)+geom%dyc(i,j-1))
     enddo
   enddo
 enddo

endsubroutine uv_to_vortdivg

!----------------------------------------------------------------------------

subroutine vortdivg_to_psichi(geom,vort,divg,psi,chi)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(in)    :: vort(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Vorticity
 real(kind=kind_real), intent(in)    :: divg(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Divergence
 real(kind=kind_real), intent(out)   :: psi (geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Stream function
 real(kind=kind_real), intent(out)   :: chi (geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Velocity potential

 real(kind=kind_real) :: psinew(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Stream function

 integer :: i,j,k,isc,iec,jsc,jec,npz
 integer :: maxiter, iter, ierr
 real(kind=kind_real) :: tolerance, converged, convergedg

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !       |        x        |
 !       |    vort,psi     |
 !       |    divg,chi     |
 !       |                 |
 !       x-----------------x

 call gauss_seidel(geom,psi,-vort)
 call gauss_seidel(geom,chi,-divg)

endsubroutine vortdivg_to_psichi

!----------------------------------------------------------------------------

subroutine gauss_seidel(geom,x,b)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: x(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)
 real(kind=kind_real), intent(in   ) :: b(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)

 integer :: i,j,k,isc,iec,jsc,jec,npz
 integer :: maxiter, iter, ierr
 real(kind=kind_real) :: tolerance, converged, convergedg
 real(kind=kind_real) :: xnew(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec
 npz = geom%npz

 x = 0.0_kind_real

 converged = 1.0_kind_real
 tolerance = 10e-4_kind_real

 maxiter = 10000
 iter = 0

 do while (iter < maxiter .and. converged > tolerance)

   iter = iter + 1

   xnew(isc:iec,jsc:jec,:) = x(isc:iec,jsc:jec,:)

   call mpp_update_domains(xnew, geom%domain, complete=.true.)

   !Rearranged Laplacian (not the corret dx's)
   do k=1,npz
     do j=jsc,jec
       do i=isc,iec
         xnew(i,j,k) = 0.25_kind_real*(xnew(i+1,j,k)+xnew(i,j+1,k)+xnew(i-1,j,k)+xnew(i,j-1,k)-geom%dxa(i,j)**2*b(i,j,k))
       enddo
     enddo
   enddo

   !Check for convergence
   converged = maxval(abs((xnew-x)))
   call mpi_allreduce(converged,convergedg,1,mpi_real8,mpi_max,mpi_comm_world,ierr)

   !Print info
   !print*, 'Gauss-Seidel iteration number and max difference, ', iter, convergedg

   !Update guess
   x = xnew

 enddo

endsubroutine gauss_seidel

!----------------------------------------------------------------------------

subroutine psichi_to_uava(geom,psi,chi,ua,va)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: psi(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) :: chi(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(out)   ::  ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (u)
 real(kind=kind_real), intent(out)   ::  va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (v)

 integer :: i,j,k

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !       |        x        |
 !       |     psi,chi     |
 !       |      ua,va      |
 !       |                 |
 !       |                 |
 !       x-----------------x

 call mpp_update_domains(psi, geom%domain, complete=.true.)
 call mpp_update_domains(chi, geom%domain, complete=.true.)

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec

        ua(i,j,k) =  (psi(i,j+1,k) - psi(i,j-1,k))/(geom%dyc(i,j) + geom%dyc(i,j+1)) + &
                     (chi(i+1,j,k) - chi(i-1,j,k))/(geom%dxc(i,j) + geom%dxc(i+1,j))
        va(i,j,k) = -(psi(i+1,j,k) - psi(i-1,j,k))/(geom%dxc(i,j) + geom%dxc(i+1,j)) + &
                     (chi(i,j+1,k) - chi(i,j-1,k))/(geom%dyc(i,j) + geom%dyc(i,j+1))

     enddo
   enddo
 enddo

end subroutine psichi_to_uava

!----------------------------------------------------------------------------

subroutine psichi_to_uava_adm(geom,psi_ad,chi_ad,ua_ad,va_ad)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: psi_ad(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) :: chi_ad(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(inout)   ::  ua_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (u)
 real(kind=kind_real), intent(inout)   ::  va_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Agrid winds (v)

 integer :: i,j,k
 real(kind=kind_real) :: temp1
 real(kind=kind_real) :: temp2
 real(kind=kind_real) :: temp3
 real(kind=kind_real) :: temp4

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !       |        x        |
 !       |     psi,chi     |
 !       |      ua,va      |
 !       |                 |
 !       |                 |
 !       x-----------------x

 do k=geom%npz,1,-1
   do j=geom%jec,geom%jsc,-1
     do i=geom%iec,geom%isc,-1

       temp1 =   va_ad(i,j,k)/(geom%dyc(i,j)+geom%dyc(i,j+1))
       temp2 = -(va_ad(i,j,k)/(geom%dxc(i,j)+geom%dxc(i+1,j)))
       chi_ad(i,j+1,k) = chi_ad(i,j+1,k) + temp1
       chi_ad(i,j-1,k) = chi_ad(i,j-1,k) - temp1
       psi_ad(i+1,j,k) = psi_ad(i+1,j,k) + temp2
       psi_ad(i-1,j,k) = psi_ad(i-1,j,k) - temp2
       va_ad(i,j,k) = 0.0_8

       temp3 = ua_ad(i,j,k)/(geom%dyc(i,j)+geom%dyc(i,j+1))
       temp4 = ua_ad(i,j,k)/(geom%dxc(i,j)+geom%dxc(i+1,j))
       psi_ad(i,j+1,k) = psi_ad(i,j+1,k) + temp3
       psi_ad(i,j-1,k) = psi_ad(i,j-1,k) - temp3
       chi_ad(i+1,j,k) = chi_ad(i+1,j,k) + temp4
       chi_ad(i-1,j,k) = chi_ad(i-1,j,k) - temp4
       ua_ad(i,j,k) = 0.0_8

     end do
   end do
 end do

 call mpp_update_domains_adm(ua_ad, chi_ad, geom%domain, complete=.true.)
 call mpp_update_domains_adm(va_ad, psi_ad, geom%domain, complete=.true.)

end subroutine psichi_to_uava_adm

!----------------------------------------------------------------------------

subroutine psichi_to_udvd(geom,psi,chi,u,v)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: psi(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Stream function
 real(kind=kind_real), intent(inout) :: chi(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(out)   ::   u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(out)   ::   v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz) !Dgrid winds (v)

 integer :: i,j,k
 real(kind=kind_real) :: chib(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz)

 !       x-----------------x
 !       |                 |
 !       |                 |
 !       |                 |
 !     vd|        x        |
 !       |     psi,chi     |
 !       |                 |
 !       |                 |
 !       |                 |
 !       x-----------------x
 !               ud

 !Fill halos of psi and chi
 call mpp_update_domains(psi, geom%domain, complete=.true.)
 call mpp_update_domains(chi, geom%domain, complete=.true.)

 !Interpolate chi to the B grid
 call a2b_ord4(chi, chib, geom, geom%npx, geom%npy, geom%isc, geom%iec, geom%jsc, geom%jec, geom%halo)

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec

        u(i,j,k) =  (psi (i,j+1,k) - psi (i,j,k))/(geom%dyc(i,j)) + &
                    (chib(i+1,j,k) - chib(i,j,k))/(geom%dx (i,j))
        v(i,j,k) = -(psi (i+1,j,k) - psi (i,j,k))/(geom%dxc(i,j)) + &
                    (chib(i,j+1,k) - chib(i,j,k))/(geom%dy (i,j))

     enddo
   enddo
 enddo

endsubroutine psichi_to_udvd

!----------------------------------------------------------------------------


!----------------------------------------------------------------------------

subroutine d_to_ac(geom, u, v, ua, va, uc, vc)

 !Convert D-grid winds to A-grid and C-grid
 !Also fills the halo of u and v on the D-grid


 implicit none

 type(fv3jedi_geom), intent(inout) :: geom
 real(kind=kind_real),           intent(inout) ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)
 real(kind=kind_real),           intent(inout) ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real),           intent(inout) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real),           intent(inout) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real), optional, intent(inout) :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real), optional, intent(inout) :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)

 real(kind=kind_real) :: wbuffer(geom%jsc:geom%jec,geom%npz)
 real(kind=kind_real) :: sbuffer(geom%isc:geom%iec,geom%npz)
 real(kind=kind_real) :: ebuffer(geom%jsc:geom%jec,geom%npz)
 real(kind=kind_real) :: nbuffer(geom%isc:geom%iec,geom%npz)

 integer isc,iec,jsc,jec
 integer isd,ied,jsd,jed
 integer npz
 integer i,j,k

 real(kind=kind_real) :: ut(geom%isd:geom%ied, geom%jsd:geom%jed)
 real(kind=kind_real) :: vt(geom%isd:geom%ied, geom%jsd:geom%jed)

 real(kind=kind_real) :: uctemp(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz)
 real(kind=kind_real) :: vctemp(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz)

 isc=geom%isc
 iec=geom%iec
 jsc=geom%jsc
 jec=geom%jec
 isd=geom%isd
 ied=geom%ied
 jsd=geom%jsd
 jed=geom%jed
 npz = geom%npz

 uctemp = 0.0_kind_real
 vctemp = 0.0_kind_real

 !Fill edge of patch winds
 call mpp_get_boundary(u, v, geom%domain, &
                       wbuffery=wbuffer, ebuffery=ebuffer, &
                       sbufferx=sbuffer, nbufferx=nbuffer, &
                       gridtype=DGRID_NE, complete=.true. )
 do k=1,npz
    do i=isc,iec
       u(i,jec+1,k) = nbuffer(i,k)
    enddo
    do j=jsc,jec
       v(iec+1,j,k) = ebuffer(j,k)
    enddo
 enddo

 !Fill halo on winds
 call mpp_update_domains(u, v, geom%domain, gridtype=DGRID_NE, complete=.true.)

 !Convert to C and A grids, winds perpendicular to grid cells
 do k=1,npz
   call d2a2c_vect(geom, u(:,:,k), v(:,:,k), &
                   ua(:,:,k), va(:,:,k), uctemp(:,:,k), vctemp(:,:,k), &
                   ut, vt, .true.)
 enddo

 if (present(uc)) uc = uctemp
 if (present(vc)) vc = vctemp

end subroutine d_to_ac

!----------------------------------------------------------------------------

subroutine d2a2c_vect(geom, u, v, ua, va, uc, vc, ut, vt, dord4)

  implicit none
  type(fv3jedi_geom), intent(IN), target :: geom
  logical, intent(in):: dord4
  real(kind=kind_real), intent(in) ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
  real(kind=kind_real), intent(in) ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed)
  real(kind=kind_real), intent(out), dimension(geom%isd:geom%ied+1,geom%jsd:geom%jed  ):: uc
  real(kind=kind_real), intent(out), dimension(geom%isd:geom%ied  ,geom%jsd:geom%jed+1):: vc
  real(kind=kind_real), intent(out), dimension(geom%isd:geom%ied  ,geom%jsd:geom%jed  ):: ua, va
  real(kind=kind_real), intent(out), dimension(geom%isd:geom%ied  ,geom%jsd:geom%jed  ):: ut, vt

 ! Local
  real(kind=kind_real), dimension(geom%isd:geom%ied,geom%jsd:geom%jed):: utmp, vtmp
  integer npt, i, j, ifirst, ilast, id
  integer :: npx, npy
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  real(kind=kind_real), pointer, dimension(:,:,:) :: sin_sg
  real(kind=kind_real), pointer, dimension(:,:)   :: cosa_u, cosa_v, cosa_s
  real(kind=kind_real), pointer, dimension(:,:)   :: rsin_u, rsin_v, rsin2
  real(kind=kind_real), pointer, dimension(:,:)   :: dxa,dya

  real(kind=kind_real), parameter :: a1 =  0.5625
  real(kind=kind_real), parameter :: a2 = -0.0625
  real(kind=kind_real), parameter :: c1 = -2./14.
  real(kind=kind_real), parameter :: c2 = 11./14.
  real(kind=kind_real), parameter :: c3 =  5./14.

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

  if ( dord4 ) then
       id = 1
  else
       id = 0
  endif

     npt = 4

! Initialize the non-existing corner regions
  utmp(:,:) = 1.0e8_kind_real
  vtmp(:,:) = 1.0e8_kind_real

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

! Contra-variant components at cell center:
     do j=js-1-id,je+1+id
        do i=is-1-id,ie+1+id
           ua(i,j) = (utmp(i,j)-vtmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
           va(i,j) = (vtmp(i,j)-utmp(i,j)*cosa_s(i,j)) * rsin2(i,j)
        enddo
     enddo

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

     ifirst = max(3,    is-1)
     ilast  = min(npx-2,ie+2)

!---------------------------------------------
! 4th order interpolation for interior points:
!---------------------------------------------
     do j=js-1,je+1
        do i=ifirst,ilast
           uc(i,j) = a2*(utmp(i-2,j)+utmp(i+1,j)) + a1*(utmp(i-1,j)+utmp(i,j))
           ut(i,j) = (uc(i,j) - v(i,j)*cosa_u(i,j))*rsin_u(i,j)
        enddo
     enddo

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

     if( is==1 ) then
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

     if( (ie+1)==npx ) then
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


     do j=js-1,je+2
      if ( j==1) then
        do i=is-1,ie+1
           vt(i,j) = edge_interpolate4(va(i,-1:2), dya(i,-1:2))
           if (vt(i,j) > 0.) then
              vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
           else
              vc(i,j) = vt(i,j)*sin_sg(i,j,2)
           end if
        enddo
      elseif ( j==0 .or. j==(npy-1)) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==2 .or. j==(npy+1)) then
        do i=is-1,ie+1
           vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
           vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
        enddo
      elseif ( j==npy ) then
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

 end subroutine d2a2c_vect

!----------------------------------------------------------------------------

real(kind=kind_real) function edge_interpolate4(ua, dxa)

 implicit none

 real(kind=kind_real), intent(in) :: ua(4)
 real(kind=kind_real), intent(in) :: dxa(4)

 real(kind=kind_real) :: t1, t2

 t1 = dxa(1) + dxa(2)
 t2 = dxa(3) + dxa(4)
 edge_interpolate4 = 0.5*( ((t1+dxa(2))*ua(2)-dxa(2)*ua(1)) / t1 + &
                           ((t2+dxa(3))*ua(3)-dxa(3)*ua(4)) / t2 )

end function edge_interpolate4

!----------------------------------------------------------------------------

 subroutine divergence_corner(geom, u, v, ua, va, divg_d)

 !> Refactored from sw_core.F90

 !> Divergence on the B-Grid from D and A grid winds

 implicit none
 type(fv3jedi_geom),   intent(in) :: geom
 real(kind=kind_real), intent(in),  dimension(geom%isd:geom%ied,  geom%jsd:geom%jed+1):: u
 real(kind=kind_real), intent(in),  dimension(geom%isd:geom%ied+1,geom%jsd:geom%jed  ):: v
 real(kind=kind_real), intent(in),  dimension(geom%isd:geom%ied  ,geom%jsd:geom%jed  ):: ua
 real(kind=kind_real), intent(in),  dimension(geom%isd:geom%ied  ,geom%jsd:geom%jed  ):: va
 real(kind=kind_real), intent(out), dimension(geom%isd:geom%ied+1,geom%jsd:geom%jed+1):: divg_d

 !Local
 integer :: is,  ie,  js,  je
 integer :: npx, npy
 integer i,j,is2,ie1
 real(kind=kind_real) :: uf(geom%isc-2:geom%iec+2,geom%jsc-1:geom%jec+2)
 real(kind=kind_real) :: vf(geom%isc-1:geom%iec+2,geom%jsc-2:geom%jec+2)

 is  = geom%isc
 ie  = geom%iec
 js  = geom%jsc
 je  = geom%jec

 npx = geom%npx
 npy = geom%npy

 is2 = max(2,is)
 ie1 = min(npx-1,ie+1)

 !uf
 do j=js,je+1
    if ( j==1 .or. j==npy ) then
      do i=is-1,ie+1
         uf(i,j) = u(i,j)*geom%dyc(i,j)*0.5*(geom%sin_sg(i,j-1,4)+geom%sin_sg(i,j,2))
      enddo
    else
      do i=is-1,ie+1
         uf(i,j) = (u(i,j)-0.25*(va(i,j-1)+va(i,j))*(geom%cos_sg(i,j-1,4)+geom%cos_sg(i,j,2)))   &
                                * geom%dyc(i,j)*0.5*(geom%sin_sg(i,j-1,4)+geom%sin_sg(i,j,2))
      enddo
    endif
 enddo

 !vf
 do j=js-1,je+1
    do i=is2,ie1
       vf(i,j) = (v(i,j) - 0.25*(ua(i-1,j)+ua(i,j))*(geom%cos_sg(i-1,j,3)+geom%cos_sg(i,j,1)))  &
                                      *geom%dxc(i,j)*0.5*(geom%sin_sg(i-1,j,3)+geom%sin_sg(i,j,1))
    enddo
    if (  is   ==  1 ) vf(1,  j) = v(1,  j)*geom%dxc(1,  j)*0.5*(geom%sin_sg(0,j,3)+geom%sin_sg(1,j,1))
    if ( (ie+1)==npx ) vf(npx,j) = v(npx,j)*geom%dxc(npx,j)*0.5*(geom%sin_sg(npx-1,j,3)+geom%sin_sg(npx,j,1))
 enddo

 !divg_d 1
 do j=js,je+1
    do i=is,ie+1
       divg_d(i,j) = (vf(i,j-1) - vf(i,j)) + (uf(i-1,j) - uf(i,j))
    enddo
 enddo

 !Remove the extra term at the corners:
 if (geom%sw_corner) divg_d(1,    1) = divg_d(1,    1) - vf(1,    0)
 if (geom%se_corner) divg_d(npx,  1) = divg_d(npx,  1) - vf(npx,  0)
 if (geom%ne_corner) divg_d(npx,npy) = divg_d(npx,npy) + vf(npx,npy)
 if (geom%nw_corner) divg_d(1,  npy) = divg_d(1,  npy) + vf(1,  npy)

 !Mult by grid area
 do j=js,je+1
    do i=is,ie+1
       divg_d(i,j) = geom%rarea_c(i,j)*divg_d(i,j)
    enddo
 enddo

end subroutine divergence_corner

!----------------------------------------------------------------------------

  subroutine a2b_ord4(qin, qout, geom, npx, npy, is, ie, js, je, ng)

  implicit none

  integer, intent(IN):: npx, npy, is, ie, js, je, ng
  real(kind=kind_real), intent(IN)::  qin(is-ng:ie+ng,js-ng:je+ng)   ! A-grid field
  real(kind=kind_real), intent(OUT):: qout(is-ng:ie+ng,js-ng:je+ng)   ! Output  B-grid field
  type(fv3jedi_geom), intent(IN), target :: geom

  real(kind=kind_real) qx(is:ie+1,js-ng:je+ng)
  real(kind=kind_real) qy(is-ng:ie+ng,js:je+1)
  real(kind=kind_real) qxx(is-ng:ie+ng,js-ng:je+ng)
  real(kind=kind_real) qyy(is-ng:ie+ng,js-ng:je+ng)
  real(kind=kind_real) g_in, g_ou
  real(kind=kind_real):: p0(2)
  real(kind=kind_real):: q1(is-1:ie+1), q2(js-1:je+1)
  integer:: i, j, is1, js1, is2, js2, ie1, je1

  real(kind=kind_real), parameter :: c1 =  2./3.
  real(kind=kind_real), parameter :: c2 = -1./6.
  real(kind=kind_real), parameter :: r3 = 1./3.
  real(kind=kind_real), parameter :: a1 =  0.5625  !  9/16
  real(kind=kind_real), parameter :: a2 = -0.0625  ! -1/16
  real(kind=kind_real), parameter :: b1 =  7./12.     ! 0.58333333
  real(kind=kind_real), parameter :: b2 = -1./12.

    is1 = max(1,is-1)
    js1 = max(1,js-1)
    is2 = max(2,is)
    js2 = max(2,js)

    ie1 = min(npx-1,ie+1)
    je1 = min(npy-1,je+1)

! Corners:
! 3-way extrapolation

    if ( geom%sw_corner ) then
          p0(1:2) = geom%grid(1,1,1:2)
        qout(1,1) = (extrap_corner(p0, geom%agrid(1,1,1:2), geom%agrid( 2, 2,1:2), qin(1,1), qin( 2, 2)) + &
                     extrap_corner(p0, geom%agrid(0,1,1:2), geom%agrid(-1, 2,1:2), qin(0,1), qin(-1, 2)) + &
                     extrap_corner(p0, geom%agrid(1,0,1:2), geom%agrid( 2,-1,1:2), qin(1,0), qin( 2,-1)))*r3

    endif
    if ( geom%se_corner ) then
            p0(1:2) = geom%grid(npx,1,1:2)
        qout(npx,1) = (extrap_corner(p0, geom%agrid(npx-1,1,1:2), geom%agrid(npx-2, 2,1:2), qin(npx-1,1), qin(npx-2, 2)) + &
                       extrap_corner(p0, geom%agrid(npx-1,0,1:2), geom%agrid(npx-2,-1,1:2), qin(npx-1,0), qin(npx-2,-1)) + &
                       extrap_corner(p0, geom%agrid(npx  ,1,1:2), geom%agrid(npx+1, 2,1:2), qin(npx  ,1), qin(npx+1, 2)))*r3
    endif
    if ( geom%ne_corner ) then
              p0(1:2) = geom%grid(npx,npy,1:2)
        qout(npx,npy) = (extrap_corner(p0, geom%agrid(npx-1,npy-1,1:2), geom%agrid(npx-2,npy-2,1:2), qin(npx-1,npy-1), qin(npx-2,npy-2)) + &
                         extrap_corner(p0, geom%agrid(npx  ,npy-1,1:2), geom%agrid(npx+1,npy-2,1:2), qin(npx  ,npy-1), qin(npx+1,npy-2)) + &
                         extrap_corner(p0, geom%agrid(npx-1,npy  ,1:2), geom%agrid(npx-2,npy+1,1:2), qin(npx-1,npy  ), qin(npx-2,npy+1)))*r3
    endif
    if ( geom%nw_corner ) then
            p0(1:2) = geom%grid(1,npy,1:2)
        qout(1,npy) = (extrap_corner(p0, geom%agrid(1,npy-1,1:2), geom%agrid( 2,npy-2,1:2), qin(1,npy-1), qin( 2,npy-2)) + &
                       extrap_corner(p0, geom%agrid(0,npy-1,1:2), geom%agrid(-1,npy-2,1:2), qin(0,npy-1), qin(-1,npy-2)) + &
                       extrap_corner(p0, geom%agrid(1,npy,  1:2), geom%agrid( 2,npy+1,1:2), qin(1,npy  ), qin( 2,npy+1)))*r3
    endif

!------------
! X-Interior:
!------------
    do j=max(1,js-2),min(npy-1,je+2)
       do i=max(3,is), min(npx-2,ie+1)
          qx(i,j) = b2*(qin(i-2,j)+qin(i+1,j)) + b1*(qin(i-1,j)+qin(i,j))
       enddo
    enddo

    ! *** West Edges:
    if ( is==1 ) then
       do j=js1, je1
          q2(j) = (qin(0,j)*geom%dxa(1,j) + qin(1,j)*geom%dxa(0,j))/(geom%dxa(0,j) + geom%dxa(1,j))
       enddo
       do j=js2, je1
          qout(1,j) = geom%edge_w(j)*q2(j-1) + (1.-geom%edge_w(j))*q2(j)
       enddo
!
       do j=max(1,js-2),min(npy-1,je+2)
             g_in = geom%dxa(2,j) / geom%dxa(1,j)
             g_ou = geom%dxa(-1,j) / geom%dxa(0,j)
          qx(1,j) = 0.5*( ((2.+g_in)*qin(1,j)-qin( 2,j))/(1.+g_in) +          &
                          ((2.+g_ou)*qin(0,j)-qin(-1,j))/(1.+g_ou) )
          qx(2,j) = ( 3.*(g_in*qin(1,j)+qin(2,j))-(g_in*qx(1,j)+qx(3,j)) ) / (2.+2.*g_in)
       enddo
    endif

    ! East Edges:
    if ( (ie+1)==npx ) then
       do j=js1, je1
          q2(j) = (qin(npx-1,j)*geom%dxa(npx,j) + qin(npx,j)*geom%dxa(npx-1,j))/(geom%dxa(npx-1,j) + geom%dxa(npx,j))
       enddo
       do j=js2, je1
          qout(npx,j) = geom%edge_e(j)*q2(j-1) + (1.-geom%edge_e(j))*q2(j)
       enddo
!
       do j=max(1,js-2),min(npy-1,je+2)
              g_in = geom%dxa(npx-2,j) / geom%dxa(npx-1,j)
              g_ou = geom%dxa(npx+1,j) / geom%dxa(npx,j)
          qx(npx,j) = 0.5*( ((2.+g_in)*qin(npx-1,j)-qin(npx-2,j))/(1.+g_in) +          &
                            ((2.+g_ou)*qin(npx,  j)-qin(npx+1,j))/(1.+g_ou) )
          qx(npx-1,j) = (3.*(qin(npx-2,j)+g_in*qin(npx-1,j)) - (g_in*qx(npx,j)+qx(npx-2,j)))/(2.+2.*g_in)
       enddo
    endif

!------------
! Y-Interior:
!------------

    do j=max(3,js),min(npy-2,je+1)
       do i=max(1,is-2), min(npx-1,ie+2)
          qy(i,j) = b2*(qin(i,j-2)+qin(i,j+1)) + b1*(qin(i,j-1) + qin(i,j))
       enddo
    enddo

    ! South Edges:
    if ( js==1 ) then
       do i=is1, ie1
          q1(i) = (qin(i,0)*geom%dya(i,1) + qin(i,1)*geom%dya(i,0))/(geom%dya(i,0) + geom%dya(i,1))
       enddo
       do i=is2, ie1
          qout(i,1) = geom%edge_s(i)*q1(i-1) + (1.-geom%edge_s(i))*q1(i)
       enddo
!
       do i=max(1,is-2),min(npx-1,ie+2)
             g_in = geom%dya(i,2) / geom%dya(i,1)
             g_ou = geom%dya(i,-1) / geom%dya(i,0)
          qy(i,1) = 0.5*( ((2.+g_in)*qin(i,1)-qin(i,2))/(1.+g_in) +          &
                          ((2.+g_ou)*qin(i,0)-qin(i,-1))/(1.+g_ou) )
          qy(i,2) = (3.*(g_in*qin(i,1)+qin(i,2)) - (g_in*qy(i,1)+qy(i,3)))/(2.+2.*g_in)
       enddo
    endif

    ! North Edges:
    if ( (je+1)==npy ) then
       do i=is1, ie1
          q1(i) = (qin(i,npy-1)*geom%dya(i,npy) + qin(i,npy)*geom%dya(i,npy-1))/(geom%dya(i,npy-1)+geom%dya(i,npy))
       enddo
       do i=is2, ie1
          qout(i,npy) = geom%edge_n(i)*q1(i-1) + (1.-geom%edge_n(i))*q1(i)
       enddo
!
       do i=max(1,is-2),min(npx-1,ie+2)
              g_in = geom%dya(i,npy-2) / geom%dya(i,npy-1)
              g_ou = geom%dya(i,npy+1) / geom%dya(i,npy)
          qy(i,npy) = 0.5*( ((2.+g_in)*qin(i,npy-1)-qin(i,npy-2))/(1.+g_in) +          &
                            ((2.+g_ou)*qin(i,npy  )-qin(i,npy+1))/(1.+g_ou) )
          qy(i,npy-1) = (3.*(qin(i,npy-2)+g_in*qin(i,npy-1)) - (g_in*qy(i,npy)+qy(i,npy-2)))/(2.+2.*g_in)
       enddo
    endif

!--------------------------------------

    do j=max(3,js),min(npy-2,je+1)
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,j) = a2*(qx(i,j-2)+qx(i,j+1)) + a1*(qx(i,j-1)+qx(i,j))
       enddo
    enddo

    if ( js==1 ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,2) = c1*(qx(i,1)+qx(i,2))+c2*(qout(i,1)+qxx(i,3))
       enddo
    endif
    if ( (je+1)==npy ) then
       do i=max(2,is),min(npx-1,ie+1)
          qxx(i,npy-1) = c1*(qx(i,npy-2)+qx(i,npy-1))+c2*(qout(i,npy)+qxx(i,npy-2))
       enddo
    endif


    do j=max(2,js),min(npy-1,je+1)
       do i=max(3,is),min(npx-2,ie+1)
          qyy(i,j) = a2*(qy(i-2,j)+qy(i+1,j)) + a1*(qy(i-1,j)+qy(i,j))
       enddo
       if ( is==1 ) qyy(2,j) = c1*(qy(1,j)+qy(2,j))+c2*(qout(1,j)+qyy(3,j))
       if((ie+1)==npx) qyy(npx-1,j) = c1*(qy(npx-2,j)+qy(npx-1,j))+c2*(qout(npx,j)+qyy(npx-2,j))

       do i=max(2,is),min(npx-1,ie+1)
          qout(i,j) = 0.5*(qxx(i,j) + qyy(i,j))   ! averaging
       enddo
    enddo

  end subroutine a2b_ord4

!----------------------------------------------------------------------------

real(kind=kind_real) function extrap_corner ( p0, p1, p2, q1, q2 )

    implicit none
    real(kind=kind_real), intent(in ), dimension(2) :: p0, p1, p2
    real(kind=kind_real), intent(in ) :: q1, q2
    real(kind=kind_real) :: x1, x2

    x1 = great_circle_dist( p1, p0 )
    x2 = great_circle_dist( p2, p0 )

    extrap_corner = q1 + x1/(x2-x1) * (q1-q2)

  end function extrap_corner

!----------------------------------------------------------------------------

  real(kind=kind_real) function great_circle_dist( q1, q2 )

       implicit none
       real(kind=kind_real), intent(IN)           :: q1(2), q2(2)

       real (kind=kind_real):: p1(2), p2(2)
       integer n

       do n=1,2
          p1(n) = q1(n)
          p2(n) = q2(n)
       enddo

       great_circle_dist = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                          sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

   end function great_circle_dist

!----------------------------------------------------------------------------

  SUBROUTINE PSICHI_TO_UDVD_ADM(geom, psi_ad, chi_ad, u_ad, v_ad)


    IMPLICIT NONE
    TYPE(FV3JEDI_GEOM), INTENT(INOUT) :: geom
!Stream function
    REAL(kind=kind_real) :: psi(geom%isd:geom%ied, &
&   geom%jsd:geom%jed, geom%npz)
    REAL(kind=kind_real), INTENT(INOUT) :: psi_ad(geom%isd:geom%&
&   ied, geom%jsd:geom%jed, geom%npz)
!Velocity potential
    REAL(kind=kind_real) :: chi(geom%isd:geom%ied, &
&   geom%jsd:geom%jed, geom%npz)
    REAL(kind=kind_real), INTENT(INOUT) :: chi_ad(geom%isd:geom%&
&   ied, geom%jsd:geom%jed, geom%npz)
!Dgrid winds (u)
    REAL(kind=kind_real) :: u(geom%isd:geom%ied, geom%jsd:geom%&
&   jed+1, geom%npz)
    REAL(kind=kind_real) :: u_ad(geom%isd:geom%ied, geom%jsd:&
&   geom%jed+1, geom%npz)
!Dgrid winds (v)
    REAL(kind=kind_real) :: v(geom%isd:geom%ied+1, geom%jsd:&
&   geom%jed, geom%npz)
    REAL(kind=kind_real) :: v_ad(geom%isd:geom%ied+1, geom%jsd:&
&   geom%jed, geom%npz)
    INTEGER :: i, j, k
    REAL(kind=kind_real) :: chib(geom%isd:geom%ied, geom%jsd:&
&   geom%jed, geom%npz)
    REAL(kind=kind_real) :: chib_ad(geom%isd:geom%ied, geom%jsd&
&   :geom%jed, geom%npz)
!       x-----------------x
!       |                 |
!       |                 |
!       |                 |
!     vd|        x        |
!       |     psi,chi     |
!       |                 |
!       |                 |
!       |                 |
!       x-----------------x
!               ud
!Fill halos of psi and chi
!Interpolate chi to the B grid
    REAL(kind=kind_real) :: temp_ad
    REAL(kind=kind_real) :: temp_ad0
    REAL(kind=kind_real) :: temp_ad1
    REAL(kind=kind_real) :: temp_ad2

    u = 0.0_kind_real
    v = 0.0_kind_real
    psi = 0.0_kind_real
    chi = 0.0_kind_real

    chib_ad = 0.0_8
    DO k=geom%npz,1,-1
      DO j=geom%jec,geom%jsc,-1
        DO i=geom%iec,geom%isc,-1
          temp_ad = v_ad(i, j, k)/geom%dy(i, j)
          temp_ad0 = -(v_ad(i, j, k)/geom%dxc(i, j))
          chib_ad(i, j+1, k) = chib_ad(i, j+1, k) + temp_ad
          chib_ad(i, j, k) = chib_ad(i, j, k) - temp_ad
          psi_ad(i+1, j, k) = psi_ad(i+1, j, k) + temp_ad0
          psi_ad(i, j, k) = psi_ad(i, j, k) - temp_ad0
          v_ad(i, j, k) = 0.0_8
          temp_ad1 = u_ad(i, j, k)/geom%dyc(i, j)
          temp_ad2 = u_ad(i, j, k)/geom%dx(i, j)
          psi_ad(i, j+1, k) = psi_ad(i, j+1, k) + temp_ad1
          psi_ad(i, j, k) = psi_ad(i, j, k) - temp_ad1
          chib_ad(i+1, j, k) = chib_ad(i+1, j, k) + temp_ad2
          chib_ad(i, j, k) = chib_ad(i, j, k) - temp_ad2
          u_ad(i, j, k) = 0.0_8
        END DO
      END DO
    END DO
    CALL A2B_ORD4_ADM(chi, chi_ad, chib, chib_ad, geom, geom%npx, geom%&
&               npy, geom%isc, geom%iec, geom%jsc, geom%jec&
&               , geom%halo)
    CALL MPP_UPDATE_DOMAINS_ADM(chi, chi_ad, geom%domain, complete=&
&                         .true.)
    CALL MPP_UPDATE_DOMAINS_ADM(psi, psi_ad, geom%domain, complete=&
&                         .true.)
  END SUBROUTINE PSICHI_TO_UDVD_ADM

!----------------------------------------------------------------------------

  SUBROUTINE A2B_ORD4_ADM(qin, qin_ad, qout, qout_ad, geom, npx, npy, is&
&   , ie, js, je, ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npx, npy, is, ie, js, je, ng
! A-grid field
    REAL(kind=kind_real), INTENT(IN) :: qin(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: qin_ad(is-ng:ie+ng, js-ng:je+ng)
! Output  B-grid field
    REAL(kind=kind_real) :: qout(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: qout_ad(is-ng:ie+ng, js-ng:je+ng)
    TYPE(FV3JEDI_GEOM), INTENT(IN), TARGET :: geom
    REAL(kind=kind_real) :: qx(is:ie+1, js-ng:je+ng)
    REAL(kind=kind_real) :: qx_ad(is:ie+1, js-ng:je+ng)
    REAL(kind=kind_real) :: qy(is-ng:ie+ng, js:je+1)
    REAL(kind=kind_real) :: qy_ad(is-ng:ie+ng, js:je+1)
    REAL(kind=kind_real) :: qxx(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: qxx_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: qyy(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: qyy_ad(is-ng:ie+ng, js-ng:je+ng)
    REAL(kind=kind_real) :: g_in, g_ou
    REAL(kind=kind_real) :: p0(2)
    REAL(kind=kind_real) :: q1(is-1:ie+1), q2(js-1:je+1)
    REAL(kind=kind_real) :: q1_ad(is-1:ie+1), q2_ad(js-1:je+1)
    INTEGER :: i, j, is1, js1, is2, js2, ie1, je1
    REAL(kind=kind_real), PARAMETER :: c1=2./3.
    REAL(kind=kind_real), PARAMETER :: c2=-(1./6.)
    REAL(kind=kind_real), PARAMETER :: r3=1./3.
!  9/16
    REAL(kind=kind_real), PARAMETER :: a1=0.5625
! -1/16
    REAL(kind=kind_real), PARAMETER :: a2=-0.0625
! 0.58333333
    REAL(kind=kind_real), PARAMETER :: b1=7./12.
    REAL(kind=kind_real), PARAMETER :: b2=-(1./12.)
    INTRINSIC MAX
    INTRINSIC MIN
    INTEGER :: max1
    INTEGER :: max2
    INTEGER :: max3
    INTEGER :: max4
    INTEGER :: max5
    INTEGER :: max6
    INTEGER :: max7
    INTEGER :: max8
    INTEGER :: max9
    INTEGER :: max10
    INTEGER :: max11
    INTEGER :: max12
    INTEGER :: max13
    INTEGER :: max14
    INTEGER :: max15
    INTEGER :: min1
    INTEGER :: min2
    INTEGER :: min3
    INTEGER :: min4
    INTEGER :: min5
    INTEGER :: min6
    INTEGER :: min7
    INTEGER :: min8
    INTEGER :: min9
    INTEGER :: min10
    INTEGER :: min11
    INTEGER :: min12
    INTEGER :: min13
    INTEGER :: min14
    INTEGER :: min15
    REAL(kind=kind_real) :: result1
    REAL(kind=kind_real) :: result1_ad
    REAL(kind=kind_real) :: result2
    REAL(kind=kind_real) :: result2_ad
    REAL(kind=kind_real) :: result3
    REAL(kind=kind_real) :: result3_ad
    REAL(kind=kind_real) :: temp_ad
    REAL(kind=kind_real) :: temp_ad0
    REAL(kind=kind_real) :: temp_ad1
    REAL(kind=kind_real) :: temp_ad2
    REAL(kind=kind_real) :: temp_ad3
    real(kind=kind_real) :: temp_ad4
    real(kind=kind_real) :: temp_ad5
    real(kind=kind_real) :: temp_ad6
    real(kind=kind_real) :: temp_ad7
    REAL(kind=kind_real) :: temp_ad8
    real(kind=kind_real) :: temp_ad9
    real(kind=kind_real) :: temp_ad10
    real(kind=kind_real) :: temp_ad11
    real(kind=kind_real) :: temp_ad12
    REAL(kind=kind_real) :: temp_ad13
    real(kind=kind_real) :: temp_ad14
    real(kind=kind_real) :: temp_ad15
    real(kind=kind_real) :: temp_ad16
    real(kind=kind_real) :: temp_ad17
    REAL(kind=kind_real) :: temp_ad18
    real(kind=kind_real) :: temp_ad19
    real(kind=kind_real) :: temp_ad20
    real(kind=kind_real) :: temp_ad21
    real(kind=kind_real) :: temp_ad22
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_from0
    INTEGER :: ad_to0
    INTEGER :: ad_from1
    INTEGER :: ad_to1
    INTEGER :: ad_from2
    INTEGER :: ad_to2
    INTEGER :: ad_from3
    INTEGER :: ad_to3
    INTEGER :: branch
    IF (1 .LT. is - 1) THEN
      is1 = is - 1
    ELSE
      is1 = 1
    END IF
    IF (1 .LT. js - 1) THEN
      js1 = js - 1
    ELSE
      js1 = 1
    END IF
    IF (2 .LT. is) THEN
      is2 = is
    ELSE
      is2 = 2
    END IF
    IF (2 .LT. js) THEN
      js2 = js
    ELSE
      js2 = 2
    END IF
    IF (npx - 1 .GT. ie + 1) THEN
      ie1 = ie + 1
    ELSE
      ie1 = npx - 1
    END IF
    IF (npy - 1 .GT. je + 1) THEN
      je1 = je + 1
    ELSE
      je1 = npy - 1
    END IF
! Corners:
! 3-way extrapolation
    IF (geom%sw_corner) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (geom%se_corner) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (geom%ne_corner) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (geom%nw_corner) THEN
      p0(1:2) = geom%grid(1, npy, 1:2)
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (1 .LT. js - 2) THEN
      max1 = js - 2
    ELSE
      max1 = 1
    END IF
    IF (npy - 1 .GT. je + 2) THEN
      min1 = je + 2
    ELSE
      min1 = npy - 1
    END IF
!------------
! X-Interior:
!------------
    DO j=max1,min1
      IF (3 .LT. is) THEN
        max2 = is
      ELSE
        max2 = 3
      END IF
      IF (npx - 2 .GT. ie + 1) THEN
        min2 = ie + 1
      ELSE
        min2 = npx - 2
      END IF
      ad_from = max2
      i = min2 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from)
    END DO
! *** West Edges:
    IF (is .EQ. 1) THEN
      IF (1 .LT. js - 2) THEN
        max3 = js - 2
      ELSE
        max3 = 1
      END IF
      IF (npy - 1 .GT. je + 2) THEN
        min3 = je + 2
      ELSE
        min3 = npy - 1
      END IF
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
! East Edges:
    IF (ie + 1 .EQ. npx) THEN
      IF (1 .LT. js - 2) THEN
        max4 = js - 2
      ELSE
        max4 = 1
      END IF
      IF (npy - 1 .GT. je + 2) THEN
        min4 = je + 2
      ELSE
        min4 = npy - 1
      END IF
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (3 .LT. js) THEN
      max5 = js
    ELSE
      max5 = 3
    END IF
    IF (npy - 2 .GT. je + 1) THEN
      min5 = je + 1
    ELSE
      min5 = npy - 2
    END IF
!------------
! Y-Interior:
!------------
    DO j=max5,min5
      IF (1 .LT. is - 2) THEN
        max6 = is - 2
      ELSE
        max6 = 1
      END IF
      IF (npx - 1 .GT. ie + 2) THEN
        min6 = ie + 2
      ELSE
        min6 = npx - 1
      END IF
      ad_from0 = max6
      i = min6 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from0)
    END DO
! South Edges:
    IF (js .EQ. 1) THEN
      IF (1 .LT. is - 2) THEN
        max7 = is - 2
      ELSE
        max7 = 1
      END IF
      IF (npx - 1 .GT. ie + 2) THEN
        min7 = ie + 2
      ELSE
        min7 = npx - 1
      END IF
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
! North Edges:
    IF (je + 1 .EQ. npy) THEN
      IF (1 .LT. is - 2) THEN
        max8 = is - 2
      ELSE
        max8 = 1
      END IF
      IF (npx - 1 .GT. ie + 2) THEN
        min8 = ie + 2
      ELSE
        min8 = npx - 1
      END IF
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (3 .LT. js) THEN
      max9 = js
    ELSE
      max9 = 3
    END IF
    IF (npy - 2 .GT. je + 1) THEN
      min9 = je + 1
    ELSE
      min9 = npy - 2
    END IF
!--------------------------------------
    DO j=max9,min9
      IF (2 .LT. is) THEN
        max10 = is
      ELSE
        max10 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        min10 = ie + 1
      ELSE
        min10 = npx - 1
      END IF
      ad_from1 = max10
      i = min10 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from1)
    END DO
    IF (js .EQ. 1) THEN
      IF (2 .LT. is) THEN
        max11 = is
      ELSE
        max11 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        min11 = ie + 1
      ELSE
        min11 = npx - 1
      END IF
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (je + 1 .EQ. npy) THEN
      IF (2 .LT. is) THEN
        max12 = is
      ELSE
        max12 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        min12 = ie + 1
      ELSE
        min12 = npx - 1
      END IF
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
    IF (2 .LT. js) THEN
      max13 = js
    ELSE
      max13 = 2
    END IF
    IF (npy - 1 .GT. je + 1) THEN
      min13 = je + 1
    ELSE
      min13 = npy - 1
    END IF
    DO j=max13,min13
      IF (3 .LT. is) THEN
        max14 = is
      ELSE
        max14 = 3
      END IF
      IF (npx - 2 .GT. ie + 1) THEN
        min14 = ie + 1
      ELSE
        min14 = npx - 2
      END IF
      ad_from2 = max14
      i = min14 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from2)
      IF (is .EQ. 1) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (ie + 1 .EQ. npx) THEN
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
      IF (2 .LT. is) THEN
        max15 = is
      ELSE
        max15 = 2
      END IF
      IF (npx - 1 .GT. ie + 1) THEN
        min15 = ie + 1
      ELSE
        min15 = npx - 1
      END IF
      ad_from3 = max15
      i = min15 + 1
      CALL PUSHINTEGER4(i - 1)
      CALL PUSHINTEGER4(ad_from3)
    END DO
    qy_ad = 0.0_8
    qxx_ad = 0.0_8
    qyy_ad = 0.0_8
    DO j=min13,max13,-1
      CALL POPINTEGER4(ad_from3)
      CALL POPINTEGER4(ad_to3)
      DO i=ad_to3,ad_from3,-1
        qxx_ad(i, j) = qxx_ad(i, j) + 0.5*qout_ad(i, j)
        qyy_ad(i, j) = qyy_ad(i, j) + 0.5*qout_ad(i, j)
        qout_ad(i, j) = 0.0_8
      END DO
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        qy_ad(npx-2, j) = qy_ad(npx-2, j) + c1*qyy_ad(npx-1, j)
        qy_ad(npx-1, j) = qy_ad(npx-1, j) + c1*qyy_ad(npx-1, j)
        qout_ad(npx, j) = qout_ad(npx, j) + c2*qyy_ad(npx-1, j)
        qyy_ad(npx-2, j) = qyy_ad(npx-2, j) + c2*qyy_ad(npx-1, j)
        qyy_ad(npx-1, j) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        qy_ad(1, j) = qy_ad(1, j) + c1*qyy_ad(2, j)
        qy_ad(2, j) = qy_ad(2, j) + c1*qyy_ad(2, j)
        qout_ad(1, j) = qout_ad(1, j) + c2*qyy_ad(2, j)
        qyy_ad(3, j) = qyy_ad(3, j) + c2*qyy_ad(2, j)
        qyy_ad(2, j) = 0.0_8
      END IF
      CALL POPINTEGER4(ad_from2)
      CALL POPINTEGER4(ad_to2)
      DO i=ad_to2,ad_from2,-1
        qy_ad(i-2, j) = qy_ad(i-2, j) + a2*qyy_ad(i, j)
        qy_ad(i+1, j) = qy_ad(i+1, j) + a2*qyy_ad(i, j)
        qy_ad(i-1, j) = qy_ad(i-1, j) + a1*qyy_ad(i, j)
        qy_ad(i, j) = qy_ad(i, j) + a1*qyy_ad(i, j)
        qyy_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      qx_ad = 0.0_8
    ELSE
      qx_ad = 0.0_8
      DO i=min12,max12,-1
        qx_ad(i, npy-2) = qx_ad(i, npy-2) + c1*qxx_ad(i, npy-1)
        qx_ad(i, npy-1) = qx_ad(i, npy-1) + c1*qxx_ad(i, npy-1)
        qout_ad(i, npy) = qout_ad(i, npy) + c2*qxx_ad(i, npy-1)
        qxx_ad(i, npy-2) = qxx_ad(i, npy-2) + c2*qxx_ad(i, npy-1)
        qxx_ad(i, npy-1) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=min11,max11,-1
        qx_ad(i, 1) = qx_ad(i, 1) + c1*qxx_ad(i, 2)
        qx_ad(i, 2) = qx_ad(i, 2) + c1*qxx_ad(i, 2)
        qout_ad(i, 1) = qout_ad(i, 1) + c2*qxx_ad(i, 2)
        qxx_ad(i, 3) = qxx_ad(i, 3) + c2*qxx_ad(i, 2)
        qxx_ad(i, 2) = 0.0_8
      END DO
    END IF
    DO j=min9,max9,-1
      CALL POPINTEGER4(ad_from1)
      CALL POPINTEGER4(ad_to1)
      DO i=ad_to1,ad_from1,-1
        qx_ad(i, j-2) = qx_ad(i, j-2) + a2*qxx_ad(i, j)
        qx_ad(i, j+1) = qx_ad(i, j+1) + a2*qxx_ad(i, j)
        qx_ad(i, j-1) = qx_ad(i, j-1) + a1*qxx_ad(i, j)
        qx_ad(i, j) = qx_ad(i, j) + a1*qxx_ad(i, j)
        qxx_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      q1_ad = 0.0_8
    ELSE
      DO i=min8,max8,-1
        g_in = geom%dya(i, npy-2)/geom%dya(i, npy-1)
        temp_ad19 = qy_ad(i, npy-1)/(g_in*2.+2.)
        qin_ad(i, npy-2) = qin_ad(i, npy-2) + 3.*temp_ad19
        qy_ad(i, npy) = qy_ad(i, npy) - g_in*temp_ad19
        qy_ad(i, npy-2) = qy_ad(i, npy-2) - temp_ad19
        qy_ad(i, npy-1) = 0.0_8
        g_ou = geom%dya(i, npy+1)/geom%dya(i, npy)
        temp_ad21 = 0.5*qy_ad(i, npy)
        temp_ad20 = temp_ad21/(g_in+1.)
        qin_ad(i, npy-1) = qin_ad(i, npy-1) + (g_in+2.)*temp_ad20 + 3.*&
&         g_in*temp_ad19
        temp_ad22 = temp_ad21/(g_ou+1.)
        qin_ad(i, npy-2) = qin_ad(i, npy-2) - temp_ad20
        qin_ad(i, npy) = qin_ad(i, npy) + (g_ou+2.)*temp_ad22
        qin_ad(i, npy+1) = qin_ad(i, npy+1) - temp_ad22
        qy_ad(i, npy) = 0.0_8
      END DO
      q1_ad = 0.0_8
      DO i=ie1,is2,-1
        q1_ad(i-1) = q1_ad(i-1) + geom%edge_n(i)*qout_ad(i, npy)
        q1_ad(i) = q1_ad(i) + (1.-geom%edge_n(i))*qout_ad(i, npy)
        qout_ad(i, npy) = 0.0_8
      END DO
      DO i=ie1,is1,-1
        temp_ad18 = q1_ad(i)/(geom%dya(i, npy-1)+geom%dya(i, npy))
        qin_ad(i, npy-1) = qin_ad(i, npy-1) + geom%dya(i, npy)*temp_ad18
        qin_ad(i, npy) = qin_ad(i, npy) + geom%dya(i, npy-1)*temp_ad18
        q1_ad(i) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO i=min7,max7,-1
        g_in = geom%dya(i, 2)/geom%dya(i, 1)
        temp_ad14 = qy_ad(i, 2)/(g_in*2.+2.)
        qin_ad(i, 1) = qin_ad(i, 1) + 3.*g_in*temp_ad14
        qin_ad(i, 2) = qin_ad(i, 2) + 3.*temp_ad14
        qy_ad(i, 1) = qy_ad(i, 1) - g_in*temp_ad14
        qy_ad(i, 3) = qy_ad(i, 3) - temp_ad14
        qy_ad(i, 2) = 0.0_8
        g_ou = geom%dya(i, -1)/geom%dya(i, 0)
        temp_ad15 = 0.5*qy_ad(i, 1)
        temp_ad16 = temp_ad15/(g_in+1.)
        temp_ad17 = temp_ad15/(g_ou+1.)
        qin_ad(i, 1) = qin_ad(i, 1) + (g_in+2.)*temp_ad16
        qin_ad(i, 2) = qin_ad(i, 2) - temp_ad16
        qin_ad(i, 0) = qin_ad(i, 0) + (g_ou+2.)*temp_ad17
        qin_ad(i, -1) = qin_ad(i, -1) - temp_ad17
        qy_ad(i, 1) = 0.0_8
      END DO
      DO i=ie1,is2,-1
        q1_ad(i-1) = q1_ad(i-1) + geom%edge_s(i)*qout_ad(i, 1)
        q1_ad(i) = q1_ad(i) + (1.-geom%edge_s(i))*qout_ad(i, 1)
        qout_ad(i, 1) = 0.0_8
      END DO
      DO i=ie1,is1,-1
        temp_ad13 = q1_ad(i)/(geom%dya(i, 0)+geom%dya(i, 1))
        qin_ad(i, 0) = qin_ad(i, 0) + geom%dya(i, 1)*temp_ad13
        qin_ad(i, 1) = qin_ad(i, 1) + geom%dya(i, 0)*temp_ad13
        q1_ad(i) = 0.0_8
      END DO
    END IF
    DO j=min5,max5,-1
      CALL POPINTEGER4(ad_from0)
      CALL POPINTEGER4(ad_to0)
      DO i=ad_to0,ad_from0,-1
        qin_ad(i, j-2) = qin_ad(i, j-2) + b2*qy_ad(i, j)
        qin_ad(i, j+1) = qin_ad(i, j+1) + b2*qy_ad(i, j)
        qin_ad(i, j-1) = qin_ad(i, j-1) + b1*qy_ad(i, j)
        qin_ad(i, j) = qin_ad(i, j) + b1*qy_ad(i, j)
        qy_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      q2_ad = 0.0_8
    ELSE
      DO j=min4,max4,-1
        g_in = geom%dxa(npx-2, j)/geom%dxa(npx-1, j)
        temp_ad9 = qx_ad(npx-1, j)/(g_in*2.+2.)
        qin_ad(npx-2, j) = qin_ad(npx-2, j) + 3.*temp_ad9
        qx_ad(npx, j) = qx_ad(npx, j) - g_in*temp_ad9
        qx_ad(npx-2, j) = qx_ad(npx-2, j) - temp_ad9
        qx_ad(npx-1, j) = 0.0_8
        g_ou = geom%dxa(npx+1, j)/geom%dxa(npx, j)
        temp_ad11 = 0.5*qx_ad(npx, j)
        temp_ad10 = temp_ad11/(g_in+1.)
        qin_ad(npx-1, j) = qin_ad(npx-1, j) + (g_in+2.)*temp_ad10 + 3.*&
&         g_in*temp_ad9
        temp_ad12 = temp_ad11/(g_ou+1.)
        qin_ad(npx-2, j) = qin_ad(npx-2, j) - temp_ad10
        qin_ad(npx, j) = qin_ad(npx, j) + (g_ou+2.)*temp_ad12
        qin_ad(npx+1, j) = qin_ad(npx+1, j) - temp_ad12
        qx_ad(npx, j) = 0.0_8
      END DO
      q2_ad = 0.0_8
      DO j=je1,js2,-1
        q2_ad(j-1) = q2_ad(j-1) + geom%edge_e(j)*qout_ad(npx, j)
        q2_ad(j) = q2_ad(j) + (1.-geom%edge_e(j))*qout_ad(npx, j)
        qout_ad(npx, j) = 0.0_8
      END DO
      DO j=je1,js1,-1
        temp_ad8 = q2_ad(j)/(geom%dxa(npx-1, j)+geom%dxa(npx, j))
        qin_ad(npx-1, j) = qin_ad(npx-1, j) + geom%dxa(npx, j)*temp_ad8
        qin_ad(npx, j) = qin_ad(npx, j) + geom%dxa(npx-1, j)*temp_ad8
        q2_ad(j) = 0.0_8
      END DO
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO j=min3,max3,-1
        g_in = geom%dxa(2, j)/geom%dxa(1, j)
        temp_ad4 = qx_ad(2, j)/(g_in*2.+2.)
        qin_ad(1, j) = qin_ad(1, j) + 3.*g_in*temp_ad4
        qin_ad(2, j) = qin_ad(2, j) + 3.*temp_ad4
        qx_ad(1, j) = qx_ad(1, j) - g_in*temp_ad4
        qx_ad(3, j) = qx_ad(3, j) - temp_ad4
        qx_ad(2, j) = 0.0_8
        g_ou = geom%dxa(-1, j)/geom%dxa(0, j)
        temp_ad5 = 0.5*qx_ad(1, j)
        temp_ad6 = temp_ad5/(g_in+1.)
        temp_ad7 = temp_ad5/(g_ou+1.)
        qin_ad(1, j) = qin_ad(1, j) + (g_in+2.)*temp_ad6
        qin_ad(2, j) = qin_ad(2, j) - temp_ad6
        qin_ad(0, j) = qin_ad(0, j) + (g_ou+2.)*temp_ad7
        qin_ad(-1, j) = qin_ad(-1, j) - temp_ad7
        qx_ad(1, j) = 0.0_8
      END DO
      DO j=je1,js2,-1
        q2_ad(j-1) = q2_ad(j-1) + geom%edge_w(j)*qout_ad(1, j)
        q2_ad(j) = q2_ad(j) + (1.-geom%edge_w(j))*qout_ad(1, j)
        qout_ad(1, j) = 0.0_8
      END DO
      DO j=je1,js1,-1
        temp_ad3 = q2_ad(j)/(geom%dxa(0, j)+geom%dxa(1, j))
        qin_ad(0, j) = qin_ad(0, j) + geom%dxa(1, j)*temp_ad3
        qin_ad(1, j) = qin_ad(1, j) + geom%dxa(0, j)*temp_ad3
        q2_ad(j) = 0.0_8
      END DO
    END IF
    DO j=min1,max1,-1
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO i=ad_to,ad_from,-1
        qin_ad(i-2, j) = qin_ad(i-2, j) + b2*qx_ad(i, j)
        qin_ad(i+1, j) = qin_ad(i+1, j) + b2*qx_ad(i, j)
        qin_ad(i-1, j) = qin_ad(i-1, j) + b1*qx_ad(i, j)
        qin_ad(i, j) = qin_ad(i, j) + b1*qx_ad(i, j)
        qx_ad(i, j) = 0.0_8
      END DO
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      temp_ad2 = r3*qout_ad(1, npy)
      result1_ad = temp_ad2
      result2_ad = temp_ad2
      result3_ad = temp_ad2
      qout_ad(1, npy) = 0.0_8
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(1, npy, 1:2), geom%agrid(2, &
&                      npy+1, 1:2), qin(1, npy), qin_ad(1, npy), qin(2, &
&                      npy+1), qin_ad(2, npy+1), result3_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(0, npy-1, 1:2), geom%agrid(-&
&                      1, npy-2, 1:2), qin(0, npy-1), qin_ad(0, npy-1), &
&                      qin(-1, npy-2), qin_ad(-1, npy-2), result2_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(1, npy-1, 1:2), geom%agrid(2&
&                      , npy-2, 1:2), qin(1, npy-1), qin_ad(1, npy-1), &
&                      qin(2, npy-2), qin_ad(2, npy-2), result1_ad)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad1 = r3*qout_ad(npx, npy)
      result1_ad = temp_ad1
      result2_ad = temp_ad1
      result3_ad = temp_ad1
      qout_ad(npx, npy) = 0.0_8
      p0(1:2) = geom%grid(npx, npy, 1:2)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx-1, npy, 1:2), geom%agrid&
&                      (npx-2, npy+1, 1:2), qin(npx-1, npy), qin_ad(npx-&
&                      1, npy), qin(npx-2, npy+1), qin_ad(npx-2, npy+1)&
&                      , result3_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx, npy-1, 1:2), geom%agrid&
&                      (npx+1, npy-2, 1:2), qin(npx, npy-1), qin_ad(npx&
&                      , npy-1), qin(npx+1, npy-2), qin_ad(npx+1, npy-2)&
&                      , result2_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx-1, npy-1, 1:2), geom%&
&                      agrid(npx-2, npy-2, 1:2), qin(npx-1, npy-1), &
&                      qin_ad(npx-1, npy-1), qin(npx-2, npy-2), qin_ad(&
&                      npx-2, npy-2), result1_ad)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad0 = r3*qout_ad(npx, 1)
      result1_ad = temp_ad0
      result2_ad = temp_ad0
      result3_ad = temp_ad0
      qout_ad(npx, 1) = 0.0_8
      p0(1:2) = geom%grid(npx, 1, 1:2)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx, 1, 1:2), geom%agrid(npx&
&                      +1, 2, 1:2), qin(npx, 1), qin_ad(npx, 1), qin(npx&
&                      +1, 2), qin_ad(npx+1, 2), result3_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx-1, 0, 1:2), geom%agrid(&
&                      npx-2, -1, 1:2), qin(npx-1, 0), qin_ad(npx-1, 0)&
&                      , qin(npx-2, -1), qin_ad(npx-2, -1), result2_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(npx-1, 1, 1:2), geom%agrid(&
&                      npx-2, 2, 1:2), qin(npx-1, 1), qin_ad(npx-1, 1), &
&                      qin(npx-2, 2), qin_ad(npx-2, 2), result1_ad)
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      temp_ad = r3*qout_ad(1, 1)
      result1_ad = temp_ad
      result2_ad = temp_ad
      result3_ad = temp_ad
      p0(1:2) = geom%grid(1, 1, 1:2)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(1, 0, 1:2), geom%agrid(2, -1&
&                      , 1:2), qin(1, 0), qin_ad(1, 0), qin(2, -1), &
&                      qin_ad(2, -1), result3_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(0, 1, 1:2), geom%agrid(-1, 2&
&                      , 1:2), qin(0, 1), qin_ad(0, 1), qin(-1, 2), &
&                      qin_ad(-1, 2), result2_ad)
      CALL EXTRAP_CORNER_ADM(p0, geom%agrid(1, 1, 1:2), geom%agrid(2, 2&
&                      , 1:2), qin(1, 1), qin_ad(1, 1), qin(2, 2), &
&                      qin_ad(2, 2), result1_ad)
    END IF
  END SUBROUTINE A2B_ORD4_ADM

!----------------------------------------------------------------------------

  SUBROUTINE EXTRAP_CORNER_ADM(p0, p1, p2, q1, q1_ad, q2, q2_ad, &
&   extrap_corner_ad)
    IMPLICIT NONE
    REAL(kind=kind_real), DIMENSION(2), INTENT(IN) :: p0, p1, p2
    REAL(kind=kind_real), INTENT(IN) :: q1, q2
    REAL(kind=kind_real) :: q1_ad, q2_ad
    REAL(kind=kind_real) :: x1, x2
    REAL(kind=kind_real) :: temp_ad
    REAL(kind=kind_real) :: extrap_corner_ad
    REAL(kind=kind_real) :: extrap_corner
    x1 = GREAT_CIRCLE_DIST(p1, p0)
    x2 = GREAT_CIRCLE_DIST(p2, p0)
    temp_ad = x1*extrap_corner_ad/(x2-x1)
    q1_ad = q1_ad + temp_ad + extrap_corner_ad
    q2_ad = q2_ad - temp_ad
  END SUBROUTINE EXTRAP_CORNER_ADM

! ------------------------------------------------------------------------------

subroutine a2d(geom, ua, va, ud, vd)

implicit none

type(fv3jedi_geom),   intent(inout) :: geom
real(kind=kind_real), intent(in)    :: ua(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(in)    :: va(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(inout) :: ud(geom%isc:geom%iec,  geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(inout) :: vd(geom%isc:geom%iec+1,geom%jsc:geom%jec,  geom%npz)

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

end subroutine a2d

! ------------------------------------------------------------------------------

subroutine a2d_ad(geom, ua_ad, va_ad, ud_ad, vd_ad)

implicit none
type(fv3jedi_geom), intent(inout) :: geom
real(kind=kind_real), intent(inout) :: ua_ad(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(inout) :: va_ad(geom%isc:geom%iec,  geom%jsc:geom%jec,  geom%npz)
real(kind=kind_real), intent(inout) :: ud_ad(geom%isc:geom%iec,  geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(inout) :: vd_ad(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz)

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

v3_ad = 0.0_8
uatemp_ad = 0.0_8
ue_ad = 0.0_8
vt1_ad = 0.0_8
vt2_ad = 0.0_8
vt3_ad = 0.0_8
ve_ad = 0.0_8
vatemp_ad = 0.0_8
ut1_ad = 0.0_8
ut2_ad = 0.0_8
ut3_ad = 0.0_8

do k=npz,1,-1

  do j=je,js,-1
    do i=ie+1,is,-1
      temp_ad0 = 0.5*vd_ad(i, j, k)
      ve_ad(i, j, 1) = ve_ad(i, j, 1) + geom%ew(1, i, j, 2)*temp_ad0
      ve_ad(i, j, 2) = ve_ad(i, j, 2) + geom%ew(2, i, j, 2)*temp_ad0
      ve_ad(i, j, 3) = ve_ad(i, j, 3) + geom%ew(3, i, j, 2)*temp_ad0
      vd_ad(i, j, k) = 0.0_8
    end do
  end do

  do j=je+1,js,-1
    do i=ie,is,-1
      temp_ad = 0.5*ud_ad(i, j, k)
      ue_ad(i, j, 1) = ue_ad(i, j, 1) + geom%es(1, i, j, 1)*temp_ad
      ue_ad(i, j, 2) = ue_ad(i, j, 2) + geom%es(2, i, j, 1)*temp_ad
      ue_ad(i, j, 3) = ue_ad(i, j, 3) + geom%es(3, i, j, 1)*temp_ad
      ud_ad(i, j, k) = 0.0_8
    end do
  end do

  if (je + 1 .eq. npy) then
    do i=ie,is,-1
      ut3_ad(i) = ut3_ad(i) + ue_ad(i, npy, 3)
      ue_ad(i, npy, 3) = 0.0_8
      ut2_ad(i) = ut2_ad(i) + ue_ad(i, npy, 2)
      ue_ad(i, npy, 2) = 0.0_8
      ut1_ad(i) = ut1_ad(i) + ue_ad(i, npy, 1)
      ue_ad(i, npy, 1) = 0.0_8
    end do
    do i=ie,is,-1
      if (i .le. im2) then
        ue_ad(i+1, npy, 3) = ue_ad(i+1, npy, 3) + geom%edge_vect_n(i)*ut3_ad(i)
        ue_ad(i, npy, 3) = ue_ad(i, npy, 3) + (1.-geom%edge_vect_n(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_8
        ue_ad(i+1, npy, 2) = ue_ad(i+1, npy, 2) + geom%edge_vect_n(i)*ut2_ad(i)
        ue_ad(i, npy, 2) = ue_ad(i, npy, 2) + (1.-geom%edge_vect_n(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_8
        ue_ad(i+1, npy, 1) = ue_ad(i+1, npy, 1) + geom%edge_vect_n(i)*ut1_ad(i)
        ue_ad(i, npy, 1) = ue_ad(i, npy, 1) + (1.-geom%edge_vect_n(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_8
      else
        ue_ad(i-1, npy, 3) = ue_ad(i-1, npy, 3) + geom%edge_vect_n(i)*ut3_ad(i)
        ue_ad(i, npy, 3) = ue_ad(i, npy, 3) + (1.-geom%edge_vect_n(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_8
        ue_ad(i-1, npy, 2) = ue_ad(i-1, npy, 2) + geom%edge_vect_n(i)*ut2_ad(i)
        ue_ad(i, npy, 2) = ue_ad(i, npy, 2) + (1.-geom%edge_vect_n(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_8
        ue_ad(i-1, npy, 1) = ue_ad(i-1, npy, 1) + geom%edge_vect_n(i)*ut1_ad(i)
        ue_ad(i, npy, 1) = ue_ad(i, npy, 1) + (1.-geom%edge_vect_n(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_8
      end if
    end do
  end if

  if (js .eq. 1) then
    do i=ie,is,-1
      ut3_ad(i) = ut3_ad(i) + ue_ad(i, 1, 3)
      ue_ad(i, 1, 3) = 0.0_8
      ut2_ad(i) = ut2_ad(i) + ue_ad(i, 1, 2)
      ue_ad(i, 1, 2) = 0.0_8
      ut1_ad(i) = ut1_ad(i) + ue_ad(i, 1, 1)
      ue_ad(i, 1, 1) = 0.0_8
    end do
    do i=ie,is,-1
      if (i .le. im2) then
        ue_ad(i+1, 1, 3) = ue_ad(i+1, 1, 3) + geom%edge_vect_s(i)*ut3_ad(i)
        ue_ad(i, 1, 3) = ue_ad(i, 1, 3) + (1.-geom%edge_vect_s(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_8
        ue_ad(i+1, 1, 2) = ue_ad(i+1, 1, 2) + geom%edge_vect_s(i)*ut2_ad(i)
        ue_ad(i, 1, 2) = ue_ad(i, 1, 2) + (1.-geom%edge_vect_s(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_8
        ue_ad(i+1, 1, 1) = ue_ad(i+1, 1, 1) + geom%edge_vect_s(i)*ut1_ad(i)
        ue_ad(i, 1, 1) = ue_ad(i, 1, 1) + (1.-geom%edge_vect_s(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_8
      else
        ue_ad(i-1, 1, 3) = ue_ad(i-1, 1, 3) + geom%edge_vect_s(i)*ut3_ad(i)
        ue_ad(i, 1, 3) = ue_ad(i, 1, 3) + (1.-geom%edge_vect_s(i))*ut3_ad(i)
        ut3_ad(i) = 0.0_8
        ue_ad(i-1, 1, 2) = ue_ad(i-1, 1, 2) + geom%edge_vect_s(i)*ut2_ad(i)
        ue_ad(i, 1, 2) = ue_ad(i, 1, 2) + (1.-geom%edge_vect_s(i))*ut2_ad(i)
        ut2_ad(i) = 0.0_8
        ue_ad(i-1, 1, 1) = ue_ad(i-1, 1, 1) + geom%edge_vect_s(i)*ut1_ad(i)
        ue_ad(i, 1, 1) = ue_ad(i, 1, 1) + (1.-geom%edge_vect_s(i))*ut1_ad(i)
        ut1_ad(i) = 0.0_8
      end if
    end do
  end if

  if (ie + 1 .eq. npx) then
    do j=je,js,-1
      vt3_ad(j) = vt3_ad(j) + ve_ad(npx, j, 3)
      ve_ad(npx, j, 3) = 0.0_8
      vt2_ad(j) = vt2_ad(j) + ve_ad(npx, j, 2)
      ve_ad(npx, j, 2) = 0.0_8
      vt1_ad(j) = vt1_ad(j) + ve_ad(npx, j, 1)
      ve_ad(npx, j, 1) = 0.0_8
    end do
    do j=je,js,-1
      if (j .le. jm2) then
        ve_ad(npx, j+1, 3) = ve_ad(npx, j+1, 3) + geom%edge_vect_e(j)*vt3_ad(j)
        ve_ad(npx, j, 3) = ve_ad(npx, j, 3) + (1.-geom%edge_vect_e(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_8
        ve_ad(npx, j+1, 2) = ve_ad(npx, j+1, 2) + geom%edge_vect_e(j)*vt2_ad(j)
        ve_ad(npx, j, 2) = ve_ad(npx, j, 2) + (1.-geom%edge_vect_e(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_8
        ve_ad(npx, j+1, 1) = ve_ad(npx, j+1, 1) + geom%edge_vect_e(j)*vt1_ad(j)
        ve_ad(npx, j, 1) = ve_ad(npx, j, 1) + (1.-geom%edge_vect_e(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_8
      else
        ve_ad(npx, j-1, 3) = ve_ad(npx, j-1, 3) + geom%edge_vect_e(j)*vt3_ad(j)
        ve_ad(npx, j, 3) = ve_ad(npx, j, 3) + (1.-geom%edge_vect_e(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_8
        ve_ad(npx, j-1, 2) = ve_ad(npx, j-1, 2) + geom%edge_vect_e(j)*vt2_ad(j)
        ve_ad(npx, j, 2) = ve_ad(npx, j, 2) + (1.-geom%edge_vect_e(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_8
        ve_ad(npx, j-1, 1) = ve_ad(npx, j-1, 1) + geom%edge_vect_e(j)*vt1_ad(j)
        ve_ad(npx, j, 1) = ve_ad(npx, j, 1) + (1.-geom%edge_vect_e(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_8
      end if
    end do
  end if

  if (is .eq. 1) then
    do j=je,js,-1
      vt3_ad(j) = vt3_ad(j) + ve_ad(1, j, 3)
      ve_ad(1, j, 3) = 0.0_8
      vt2_ad(j) = vt2_ad(j) + ve_ad(1, j, 2)
      ve_ad(1, j, 2) = 0.0_8
      vt1_ad(j) = vt1_ad(j) + ve_ad(1, j, 1)
      ve_ad(1, j, 1) = 0.0_8
    end do
    do j=je,js,-1
      if (j .le. jm2) then
        ve_ad(1, j+1, 3) = ve_ad(1, j+1, 3) + geom%edge_vect_w(j)*vt3_ad(j)
        ve_ad(1, j, 3) = ve_ad(1, j, 3) + (1.-geom%edge_vect_w(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_8
        ve_ad(1, j+1, 2) = ve_ad(1, j+1, 2) + geom%edge_vect_w(j)*vt2_ad(j)
        ve_ad(1, j, 2) = ve_ad(1, j, 2) + (1.-geom%edge_vect_w(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_8
        ve_ad(1, j+1, 1) = ve_ad(1, j+1, 1) + geom%edge_vect_w(j)*vt1_ad(j)
        ve_ad(1, j, 1) = ve_ad(1, j, 1) + (1.-geom%edge_vect_w(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_8
      else
        ve_ad(1, j-1, 3) = ve_ad(1, j-1, 3) + geom%edge_vect_w(j)*vt3_ad(j)
        ve_ad(1, j, 3) = ve_ad(1, j, 3) + (1.-geom%edge_vect_w(j))*vt3_ad(j)
        vt3_ad(j) = 0.0_8
        ve_ad(1, j-1, 2) = ve_ad(1, j-1, 2) + geom%edge_vect_w(j)*vt2_ad(j)
        ve_ad(1, j, 2) = ve_ad(1, j, 2) + (1.-geom%edge_vect_w(j))*vt2_ad(j)
        vt2_ad(j) = 0.0_8
        ve_ad(1, j-1, 1) = ve_ad(1, j-1, 1) + geom%edge_vect_w(j)*vt1_ad(j)
        ve_ad(1, j, 1) = ve_ad(1, j, 1) + (1.-geom%edge_vect_w(j))*vt1_ad(j)
        vt1_ad(j) = 0.0_8
      end if
    end do
  end if

  do j=je+1,js-1,-1
    do i=ie+1,is,-1
      v3_ad(i-1, j, 3) = v3_ad(i-1, j, 3) + ve_ad(i, j, 3)
      v3_ad(i, j, 3) = v3_ad(i, j, 3) + ve_ad(i, j, 3)
      ve_ad(i, j, 3) = 0.0_8
      v3_ad(i-1, j, 2) = v3_ad(i-1, j, 2) + ve_ad(i, j, 2)
      v3_ad(i, j, 2) = v3_ad(i, j, 2) + ve_ad(i, j, 2)
      ve_ad(i, j, 2) = 0.0_8
      v3_ad(i-1, j, 1) = v3_ad(i-1, j, 1) + ve_ad(i, j, 1)
      v3_ad(i, j, 1) = v3_ad(i, j, 1) + ve_ad(i, j, 1)
      ve_ad(i, j, 1) = 0.0_8
    end do
  end do

  do j=je+1,js,-1
    do i=ie+1,is-1,-1
      v3_ad(i, j-1, 3) = v3_ad(i, j-1, 3) + ue_ad(i, j, 3)
      v3_ad(i, j, 3) = v3_ad(i, j, 3) + ue_ad(i, j, 3)
      ue_ad(i, j, 3) = 0.0_8
      v3_ad(i, j-1, 2) = v3_ad(i, j-1, 2) + ue_ad(i, j, 2)
      v3_ad(i, j, 2) = v3_ad(i, j, 2) + ue_ad(i, j, 2)
      ue_ad(i, j, 2) = 0.0_8
      v3_ad(i, j-1, 1) = v3_ad(i, j-1, 1) + ue_ad(i, j, 1)
      v3_ad(i, j, 1) = v3_ad(i, j, 1) + ue_ad(i, j, 1)
      ue_ad(i, j, 1) = 0.0_8
    end do
  end do

  do j=je+1,js-1,-1
    do i=ie+1,is-1,-1
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 3)*v3_ad(i, j, 3)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 3)*v3_ad(i, j, 3)
      v3_ad(i, j, 3) = 0.0_8
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 2)*v3_ad(i, j, 2)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 2)*v3_ad(i, j, 2)
      v3_ad(i, j, 2) = 0.0_8
      uatemp_ad(i, j, k) = uatemp_ad(i, j, k) + geom%vlon(i, j, 1)*v3_ad(i, j, 1)
      vatemp_ad(i, j, k) = vatemp_ad(i, j, k) + geom%vlat(i, j, 1)*v3_ad(i, j, 1)
      v3_ad(i, j, 1) = 0.0_8
    end do
  end do
end do

call mpp_update_domains_adm(vatemp, vatemp_ad, geom%domain, complete=.true.)
call mpp_update_domains_adm(uatemp, uatemp_ad, geom%domain, complete=.true.)

va_ad = va_ad + vatemp_ad(is:ie, js:je, :)
ua_ad = ua_ad + uatemp_ad(is:ie, js:je, :)

end subroutine a2d_ad

! ------------------------------------------------------------------------------

subroutine d2a(geom, u_comp, v_comp, ua_comp, va_comp)

implicit none

type(fv3jedi_geom),   intent(inout) :: geom
real(kind=kind_real), intent(inout) ::  u_comp(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,geom%npz)
real(kind=kind_real), intent(inout) ::  v_comp(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,geom%npz)
real(kind=kind_real), intent(inout) :: ua_comp(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)
real(kind=kind_real), intent(inout) :: va_comp(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,geom%npz)

real(kind=kind_real) :: c1 =  1.125
real(kind=kind_real) :: c2 = -0.125
real(kind=kind_real) :: utmp(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
real(kind=kind_real) :: vtmp(geom%isc:geom%iec+1,geom%jsc:geom%jec)
real(kind=kind_real) :: wu(geom%isc:geom%iec,  geom%jsc:geom%jec+1)
real(kind=kind_real) :: wv(geom%isc:geom%iec+1,geom%jsc:geom%jec)

real(kind=kind_real) ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,geom%npz)
real(kind=kind_real) ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,geom%npz)
real(kind=kind_real) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,geom%npz)
real(kind=kind_real) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,geom%npz)

real(kind=kind_real) :: ebuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: nbufferx(geom%isc:geom%iec,1:geom%npz)
real(kind=kind_real) :: wbuffery(geom%jsc:geom%jec,1:geom%npz)
real(kind=kind_real) :: sbufferx(geom%isc:geom%iec,1:geom%npz)

integer i, j, k
integer :: is,  ie,  js,  je, npx, npy, npz

!Shorten notation
is  = geom%isc
ie  = geom%iec
js  = geom%jsc
je  = geom%jec
npx = geom%npx
npy = geom%npy
npz = geom%npz

!Fill compute part from input
u = 0.0_kind_real
v = 0.0_kind_real
u(is:ie  ,js:je+1,:) = u_comp
v(is:ie+1,js:je  ,:) = v_comp

!Fill edges
ebuffery = 0.0_kind_real
nbufferx = 0.0_kind_real
wbuffery = 0.0_kind_real
sbufferx = 0.0_kind_real
call mpp_get_boundary( u, v, geom%domain, &
                       wbuffery=wbuffery, ebuffery=ebuffery, &
                       sbufferx=sbufferx, nbufferx=nbufferx, &
                       gridtype=DGRID_NE, complete=.true. )
do k=1,npz
   do i=is,ie
      u(i,je+1,k) = nbufferx(i,k)
   enddo
enddo
do k=1,npz
   do j=js,je
      v(ie+1,j,k) = ebuffery(j,k)
   enddo
enddo

!Fill halos
call mpp_update_domains(u, v, geom%domain, gridtype=DGRID_NE)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,npx,npy,c2,c1, &
!$OMP                                  u,v,ua,va)         &
!$OMP                          private(utmp, vtmp, wu, wv)
do k=1,npz

  do j=max(2,js),min(npy-2,je)
    do i=max(2,is),min(npx-2,ie)
      utmp(i,j) = c2*(u(i,j-1,k)+u(i,j+2,k)) + c1*(u(i,j,k)+u(i,j+1,k))
      vtmp(i,j) = c2*(v(i-1,j,k)+v(i+2,j,k)) + c1*(v(i,j,k)+v(i+1,j,k))
    enddo
  enddo

  if ( js==1  ) then
    do i=is,ie+1
      wv(i,1) = v(i,1,k)*geom%dy(i,1)
    enddo
    do i=is,ie
      vtmp(i,1) = 2.*(wv(i,1) + wv(i+1,1)) / (geom%dy(i,1)+geom%dy(i+1,1))
      utmp(i,1) = 2.*(u(i,1,k)*geom%dx(i,1) + u(i,2,k)*geom%dx(i,2))   &
                   / (         geom%dx(i,1) +          geom%dx(i,2))
    enddo
  endif

  if ( (je+1)==npy   ) then
    j = npy-1
    do i=is,ie+1
      wv(i,j) = v(i,j,k)*geom%dy(i,j)
    enddo
    do i=is,ie
      vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j)) / (geom%dy(i,j)+geom%dy(i+1,j))
      utmp(i,j) = 2.*(u(i,j,k)*geom%dx(i,j) + u(i,j+1,k)*geom%dx(i,j+1))   &
                  / (         geom%dx(i,j) +            geom%dx(i,j+1))
    enddo
  endif

  if ( is==1 ) then
    i = 1
    do j=js,je
      wv(1,j) = v(1,j,k)*geom%dy(1,j)
      wv(2,j) = v(2,j,k)*geom%dy(2,j)
    enddo
    do j=js,je+1
      wu(i,j) = u(i,j,k)*geom%dx(i,j)
    enddo
    do j=js,je
      utmp(i,j) = 2.*(wu(i,j) + wu(i,j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
      vtmp(i,j) = 2.*(wv(1,j) + wv(2,j  ))/(geom%dy(1,j)+geom%dy(2,j))
    enddo
  endif

  if ( (ie+1)==npx) then
    i = npx-1
    do j=js,je
      wv(i,  j) = v(i,  j,k)*geom%dy(i,  j)
      wv(i+1,j) = v(i+1,j,k)*geom%dy(i+1,j)
    enddo
    do j=js,je+1
      wu(i,j) = u(i,j,k)*geom%dx(i,j)
    enddo
    do j=js,je
      utmp(i,j) = 2.*(wu(i,j) + wu(i,  j+1))/(geom%dx(i,j)+geom%dx(i,j+1))
      vtmp(i,j) = 2.*(wv(i,j) + wv(i+1,j  ))/(geom%dy(i,j)+geom%dy(i+1,j))
    enddo
  endif

  !Transform local a-grid winds into latitude-longitude coordinates
  do j=js,je
    do i=is,ie
      ua(i,j,k) = geom%a11(i,j)*utmp(i,j) + geom%a12(i,j)*vtmp(i,j)
      va(i,j,k) = geom%a21(i,j)*utmp(i,j) + geom%a22(i,j)*vtmp(i,j)
    enddo
  enddo

enddo

ua_comp = ua(is:ie,js:je,:)
va_comp = va(is:ie,js:je,:)

end subroutine d2a

! ------------------------------------------------------------------------------

subroutine d2a_ad(geom, u_ad_comp, v_ad_comp, ua_ad_comp, va_ad_comp)

implicit none
type(fv3jedi_geom),   intent(inout) :: geom
real(kind=kind_real), intent(inout) ::  u_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec+1,1:geom%npz)
real(kind=kind_real), intent(inout) ::  v_ad_comp(geom%isc:geom%iec+1,geom%jsc:geom%jec,  1:geom%npz)
real(kind=kind_real), intent(inout) :: ua_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec,  1:geom%npz)
real(kind=kind_real), intent(inout) :: va_ad_comp(geom%isc:geom%iec,  geom%jsc:geom%jec,  1:geom%npz)

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

!$omp parallel do default(none) shared(is,ie,js,je,npz,npx,npy,c2,c1, &
!$omp                                  u,v,ua,va)         &
!$omp                          private(utmp, vtmp, wu, wv)
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
vtmp_ad = 0.0_8
wu_ad = 0.0_8
wv_ad = 0.0_8
utmp_ad = 0.0_8
do k=npz,1,-1
  do j=je,js,-1
    do i=ie,is,-1
      utmp_ad(i, j) = utmp_ad(i, j) + geom%a11(i, j)*ua_ad(i, j, k) + geom%a21(i, j)*va_ad(i, j, k)
      vtmp_ad(i, j) = vtmp_ad(i, j) + geom%a12(i, j)*ua_ad(i, j, k) + geom%a22(i, j)*va_ad(i, j, k)
      va_ad(i, j, k) = 0.0_8
      ua_ad(i, j, k) = 0.0_8
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
      vtmp_ad(i, j) = 0.0_8
      temp_ad6 = 2.*utmp_ad(i, j)/(geom%dx(i, j)+geom%dx(i, j+1))
      wu_ad(i, j) = wu_ad(i, j) + temp_ad6
      wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad6
      utmp_ad(i, j) = 0.0_8
    end do
    do j=je+1,js,-1
      u_ad(i, j, k) = u_ad(i, j, k) + geom%dx(i, j)*wu_ad(i, j)
      wu_ad(i, j) = 0.0_8
    end do
    do j=je,js,-1
      v_ad(i+1, j, k) = v_ad(i+1, j, k) + geom%dy(i+1, j)*wv_ad(i+1, j)
      wv_ad(i+1, j) = 0.0_8
      v_ad(i, j, k) = v_ad(i, j, k) + geom%dy(i, j)*wv_ad(i, j)
      wv_ad(i, j) = 0.0_8
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
      vtmp_ad(i, j) = 0.0_8
      temp_ad4 = 2.*utmp_ad(i, j)/(geom%dx(i, j)+geom%dx(i, j+1))
      wu_ad(i, j) = wu_ad(i, j) + temp_ad4
      wu_ad(i, j+1) = wu_ad(i, j+1) + temp_ad4
      utmp_ad(i, j) = 0.0_8
    end do
    do j=je+1,js,-1
      u_ad(i, j, k) = u_ad(i, j, k) + geom%dx(i, j)*wu_ad(i, j)
      wu_ad(i, j) = 0.0_8
    end do
    do j=je,js,-1
      v_ad(2, j, k) = v_ad(2, j, k) + geom%dy(2, j)*wv_ad(2, j)
      wv_ad(2, j) = 0.0_8
      v_ad(1, j, k) = v_ad(1, j, k) + geom%dy(1, j)*wv_ad(1, j)
      wv_ad(1, j) = 0.0_8
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
      utmp_ad(i, j) = 0.0_8
      temp_ad2 = 2.*vtmp_ad(i, j)/(geom%dy(i, j)+geom%dy(i+1, j))
      wv_ad(i, j) = wv_ad(i, j) + temp_ad2
      wv_ad(i+1, j) = wv_ad(i+1, j) + temp_ad2
      vtmp_ad(i, j) = 0.0_8
    end do
    do i=ie+1,is,-1
      v_ad(i, j, k) = v_ad(i, j, k) + geom%dy(i, j)*wv_ad(i, j)
      wv_ad(i, j) = 0.0_8
    end do
    call popinteger4(i)
  end if
  call popcontrol1b(branch)
  if (branch .eq. 0) then
    do i=ie,is,-1
      temp_ad = 2.*utmp_ad(i, 1)/(geom%dx(i, 1)+geom%dx(i, 2))
      u_ad(i, 1, k) = u_ad(i, 1, k) + geom%dx(i, 1)*temp_ad
      u_ad(i, 2, k) = u_ad(i, 2, k) + geom%dx(i, 2)*temp_ad
      utmp_ad(i, 1) = 0.0_8
      temp_ad0 = 2.*vtmp_ad(i, 1)/(geom%dy(i, 1)+geom%dy(i+1, 1))
      wv_ad(i, 1) = wv_ad(i, 1) + temp_ad0
      wv_ad(i+1, 1) = wv_ad(i+1, 1) + temp_ad0
      vtmp_ad(i, 1) = 0.0_8
    end do
    do i=ie+1,is,-1
      v_ad(i, 1, k) = v_ad(i, 1, k) + geom%dy(i, 1)*wv_ad(i, 1)
      wv_ad(i, 1) = 0.0_8
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
      vtmp_ad(i, j) = 0.0_8
      u_ad(i, j-1, k) = u_ad(i, j-1, k) + c2*utmp_ad(i, j)
      u_ad(i, j+2, k) = u_ad(i, j+2, k) + c2*utmp_ad(i, j)
      u_ad(i, j, k) = u_ad(i, j, k) + c1*utmp_ad(i, j)
      u_ad(i, j+1, k) = u_ad(i, j+1, k) + c1*utmp_ad(i, j)
      utmp_ad(i, j) = 0.0_8
    end do
    call popinteger4(i)
  end do
end do

!Adjoint fill halos
call mpp_update_domains_adm(u, u_ad, v, v_ad, geom%domain, gridtype=dgrid_ne)

!Adjoint fill edges
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

!Output compute part
u_ad_comp = u_ad(is:ie  ,js:je+1,:)
v_ad_comp = v_ad(is:ie+1,js:je  ,:)

end subroutine d2a_ad

! ------------------------------------------------------------------------------

end module wind_vt_mod
