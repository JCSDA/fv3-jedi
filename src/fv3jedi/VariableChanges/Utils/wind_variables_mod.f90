! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wind_vt_mod

use fv3jedi_constants_mod, only: pi, rad2deg
use fv3jedi_geom_mod,  only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_communication_mod, only: gather_field, scatter_field

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

public psichi_to_uava
public psichi_to_uava_adm

public psichi_to_udvd
public udvd_to_psichi

public psichi_to_vortdivg

public udvd_to_vort
public uava_to_vort

public a2d
public a2d_ad
public d2a
public d2a_ad

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
 real(kind=kind_real), parameter :: quadcof(4,2) = reshape((/ 0.0_kind_real,  1.0_kind_real, &
                                                              1.0_kind_real,  2.0_kind_real, &
                                                              1.0_kind_real, -1.0_kind_real, &
                                                              1.0_kind_real, -1.0_kind_real /), (/4, 2/))

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

subroutine psichi_to_udvd(geom,psi,chi,u,v)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(in)    :: psi(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Stream function
 real(kind=kind_real), intent(in)    :: chi(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Velocity potential
 real(kind=kind_real), intent(out)   ::   u(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real), intent(out)   ::   v(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Dgrid winds (v)

 integer :: i,j,k
 real(kind=kind_real), allocatable, dimension(:,:,:) :: psid, chid
 real(kind=kind_real) :: u_rot(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz) !Rotational winds (u)
 real(kind=kind_real) :: v_rot(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz) !Rotational winds (v)
 real(kind=kind_real) :: u_div(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz) !Divergent winds (u)
 real(kind=kind_real) :: v_div(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz) !Divergent winds (v)

 real(kind=kind_real) :: uc_div(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz) !Divergent winds (u)
 real(kind=kind_real) :: vc_div(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz) !Divergent winds (v)
 real(kind=kind_real) :: ua_div(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Divergent winds (u)
 real(kind=kind_real) :: va_div(geom%isd:geom%ied  ,geom%jsd:geom%jed  ,1:geom%npz) !Divergent winds (v)

 ! Fill halos on psi and chi
 !--------------------------
 allocate(psid(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 allocate(chid(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%npz))
 psid = 0.0_kind_real
 chid = 0.0_kind_real

 psid(geom%isc:geom%iec,geom%jsc:geom%jec,:) = psi(geom%isc:geom%iec,geom%jsc:geom%jec,:)
 chid(geom%isc:geom%iec,geom%jsc:geom%jec,:) = chi(geom%isc:geom%iec,geom%jsc:geom%jec,:)

 call mpp_update_domains(psid, geom%domain, complete=.true.)
 call mpp_update_domains(chid, geom%domain, complete=.true.)


 ! Compute rotational and divergent winds
 ! --------------------------------------

 do k=1,geom%npz
   do j=geom%jsc,geom%jec+1
     do i=geom%isc,geom%iec
        u_rot(i,j,k) = -(psid(i,j,k) - psid(i,j-1,k))/geom%dyc(i,j)
     enddo
   enddo
 enddo

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec+1
        v_rot(i,j,k) = (psid(i,j,k) - psid(i-1,j,k))/geom%dxc(i,j)
     enddo
   enddo
 enddo

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec+1
        uc_div(i,j,k) = (chid(i,j,k) - chid(i-1,j,k))/geom%dxc(i,j)
     enddo
   enddo
 enddo

 do k=1,geom%npz
   do j=geom%jsc,geom%jec+1
     do i=geom%isc,geom%iec
        vc_div(i,j,k) = (chid(i,j,k) - chid(i,j-1,k))/geom%dyc(i,j)
     enddo
   enddo
 enddo

 ! Convert C-grid divergent winds to A-grid divergent winds
 ! --------------------------------------------------------
 call fill_cgrid_winds(geom,uc_div,vc_div,fillhalo=.true.)

 do k = 1,geom%npz
   call ctoa(geom, uc_div(:,:,k), vc_div(:,:,k), ua_div(:,:,k), va_div(:,:,k))
 enddo

 ! Convert A-grid divergent winds to D-grid divergent winds
 ! --------------------------------------------------------

 do k = 1,geom%npz
   call atod(geom, ua_div(:,:,k), va_div(:,:,k), u_div(:,:,k), v_div(:,:,k))
 enddo


 ! Sum rotational and divergent parts
 ! ----------------------------------

 call fill_dgrid_winds(geom,u_div,v_div,fillhalo=.true.)

 do k=1,geom%npz
   do j=geom%jsc,geom%jec+1
     do i=geom%isc,geom%iec
        u(i,j,k) = u_rot(i,j,k) + u_div(i,j,k)
     enddo
   enddo
 enddo

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec+1
        v(i,j,k) = v_rot(i,j,k) + v_div(i,j,k)
     enddo
   enddo
 enddo

end subroutine psichi_to_udvd

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

subroutine udvd_to_psichi(geom,grid,oprs,u_in,v_in,psi,chi,lsize,lev_start,lev_final,vor_out,div_out)

 implicit none
 type(fv3jedi_geom),             intent(inout) :: geom
 type(fempsgrid),                intent(in)    :: grid
 type(fempsoprs),                intent(in)    :: oprs
 real(kind=kind_real),           intent(inout) :: u_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) !Dgrid winds (u)
 real(kind=kind_real),           intent(inout) :: v_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) !Dgrid winds (v)
 real(kind=kind_real),           intent(out)   ::  psi(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Stream function
 real(kind=kind_real),           intent(out)   ::  chi(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) !Velocity potential
 integer,                        intent(in)    :: lsize
 integer,                        intent(in)    :: lev_start(lsize),lev_final(lsize)
 real(kind=kind_real), optional, intent(out)   :: vor_out(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real), optional, intent(out)   :: div_out(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential

 integer :: i,j,k
 real(kind=kind_real), allocatable, dimension(:,:,:) :: u, v
 real(kind=kind_real), allocatable, dimension(:,:,:) :: uc, vc
 real(kind=kind_real), allocatable, dimension(:,:,:) :: ut, vt

 real(kind=kind_real), allocatable, dimension(:,:,:) :: vor, div

 real(kind=kind_real), allocatable, dimension(:,:,:) :: vorgcomm, divgcomm
 real(kind=kind_real), allocatable, dimension(:,:,:) :: psigcomm, chigcomm
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: vorg, divg !Global level of vor and div
 real(kind=kind_real), allocatable, dimension(:,:,:,:) :: psig, chig !Global level of psi and chi

 real(kind=kind_real), allocatable, dimension(:) :: voru, divu     !Unstructured vor and div
 real(kind=kind_real), allocatable, dimension(:) :: psiu, chiu     !Unstructured psi and chi

 integer :: ranki, n
 integer, allocatable :: lev_proc(:)

 ! ------------------------------------------ !
 ! Convert D-grid winds to A-grid psi and chi !
 ! ------------------------------------------ !


 ! Fill edge of D grid winds
 ! -------------------------
 allocate(u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz))
 allocate(v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz))

 ! Copy internal part
 u(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = u_in(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 v(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = v_in(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 call fill_dgrid_winds(geom,u,v,fillhalo=.true.)


 ! Get C grid winds
 ! ----------------
 allocate(uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz))
 allocate(vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz))

 do k = 1,geom%npz
   call udvd_to_ucvc(geom,u(:,:,k),v(:,:,k),uc(:,:,k),vc(:,:,k))
 enddo

 call fill_cgrid_winds(geom,uc,vc,fillhalo=.true.)


 ! Compute ut,vt
 ! -------------
 allocate(ut(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz))
 allocate(vt(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz))

 do k = 1,geom%npz

   call compute_utvt(geom, uc(:,:,k), vc(:,:,k), ut(:,:,k), vt(:,:,k), 1.0_kind_real)

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
 allocate(vor(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
 allocate(div(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
 vor = 0.0_kind_real
 div = 0.0_kind_real

 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       vor(i,j,k) = geom%rarea(i,j)*( u(i,j,k)*geom%dx(i,j)-u(i,j+1,k)*geom%dx(i,j+1) - &
                                      v(i,j,k)*geom%dy(i,j)+v(i+1,j,k)*geom%dy(i+1,j))
       div(i,j,k) = geom%rarea(i,j)*( ut(i+1,j,k)-ut(i,j,k) + vt(i,j+1,k)-vt(i,j,k) )
                                  !geom%rarea(i,j)*( uc(i+1,j,k)*geom%dy(i+1,j)-uc(i,j,k)*geom%dy(i,j) + &
                                  !    vc(i,j+1,k)*geom%dx(i,j+1)-vc(i,j,k)*geom%dx(i,j+1) )
     enddo
   enddo
 enddo

 if (present(vor_out)) vor_out = vor
 if (present(div_out)) div_out = div

 ! Gather voricity and divergence to one processor and compute psi and chi
 ! -----------------------------------------------------------------------

 ranki = geom%f_comm%rank() + 1

 allocate(vorgcomm(1:geom%npx-1,1:geom%npy-1,6))
 allocate(divgcomm(1:geom%npx-1,1:geom%npy-1,6))

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

 ! Gather field to respective processor
 do k=1,geom%npz
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,vor(:,:,k),vorgcomm)
   call gather_field(geom,geom%f_comm%communicator(),lev_proc(k)-1,div(:,:,k),divgcomm)
   if (ranki == lev_proc(k)) then
     vorg(:,:,:,k) = vorgcomm
     divg(:,:,:,k) = divgcomm
   endif
 enddo

 deallocate(vorgcomm,divgcomm)

 ! Loop over level and compute psi/chi
 if (ranki <= lsize) then ! Only processors with a level
   do k = lev_start(ranki),lev_final(ranki)

     if (geom%f_comm%rank()==0) &
       print*, "Doing femps level ", k, ". Proc range:",lev_start(ranki),"-",lev_final(ranki)

     call fv3field_to_ufield(grid,geom%npx-1,vorg(:,:,:,k),voru)
     call fv3field_to_ufield(grid,geom%npx-1,divg(:,:,:,k),divu)

     ! Convert to area integrals, required by femps
     voru = voru*grid%farea(:,grid%ngrids)
     divu = divu*grid%farea(:,grid%ngrids)

     ! Solve poisson equation (\psi=\nabla^{-2}\zeta, \chi=\nabla^{-2}D)
     call inverselaplace(grid,oprs,grid%ngrids,voru,psiu,level=k)
     call inverselaplace(grid,oprs,grid%ngrids,divu,chiu,level=k)

     ! Convert from area integrals
     psiu = psiu/grid%farea(:,grid%ngrids)
     chiu = chiu/grid%farea(:,grid%ngrids)

     call ufield_to_fv3field(grid,geom%npx-1,psiu,psig(:,:,:,k))
     call ufield_to_fv3field(grid,geom%npx-1,chiu,chig(:,:,:,k))

   enddo
 endif

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

end subroutine udvd_to_psichi

! ------------------------------------------------------------------------------

subroutine psichi_to_vortdivg(geom,grid,oprs,psi,chi,lsize,lev_start,lev_final,vor,div)

 implicit none
 type(fv3jedi_geom),             intent(inout) :: geom
 type(fempsgrid),                intent(in)    :: grid
 type(fempsoprs),                intent(in)    :: oprs
 real(kind=kind_real),           intent(in)    ::  psi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Stream function
 real(kind=kind_real),           intent(in)    ::  chi(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Velocity potential
 integer,                        intent(in)    :: lsize
 integer,                        intent(in)    :: lev_start(lsize),lev_final(lsize)
 real(kind=kind_real), optional, intent(out)   ::  vor(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Vorticity
 real(kind=kind_real), optional, intent(out)   ::  div(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !Divergence

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

! ------------------------------------------------------------------------------

subroutine udvd_to_vort(geom, ud_in, vd_in, vort)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) :: ud_in(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz) ! Dgrid winds (u)
 real(kind=kind_real), intent(inout) :: vd_in(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz) ! Dgrid winds (v)
 real(kind=kind_real), intent(inout) ::  vort(geom%isc:geom%iec  ,geom%jsc:geom%jec  ,1:geom%npz) ! Vorticity

 real(kind=kind_real), allocatable :: ud(:,:,:)
 real(kind=kind_real), allocatable :: vd(:,:,:)
 integer :: i, j, k

 ! --------------------------------- !
 ! Convert D-grid winds to vorticity !
 ! --------------------------------- !

 ! Fill edge of D-grid winds
 ! -------------------------
 allocate(ud(geom%isd:geom%ied  ,geom%jsd:geom%jed+1,1:geom%npz))
 allocate(vd(geom%isd:geom%ied+1,geom%jsd:geom%jed  ,1:geom%npz))

 ! Copy internal part
 ud(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = ud_in(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)
 vd(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) = vd_in(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)

 call fill_dgrid_winds(geom, ud, vd, fillhalo=.true.)

 !D-grid u and v to A-grid vorticity
 !----------------------------------
 vort = 0.0_kind_real
 do k=1,geom%npz
   do j=geom%jsc,geom%jec
     do i=geom%isc,geom%iec
       vort(i,j,k) = geom%rarea(i,j)*( ud(i,j,k)*geom%dx(i,j)-ud(i,j+1,k)*geom%dx(i,j+1) - &
                                       vd(i,j,k)*geom%dy(i,j)+vd(i+1,j,k)*geom%dy(i+1,j))
     enddo
   enddo
 enddo

end subroutine udvd_to_vort

! ------------------------------------------------------------------------------

subroutine uava_to_vort(geom, ua, va, vort)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(inout) ::   ua(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) ! Dgrid winds (u)
 real(kind=kind_real), intent(inout) ::   va(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) ! Dgrid winds (v)
 real(kind=kind_real), intent(inout) :: vort(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) ! Vorticity

 real(kind=kind_real) :: ud(geom%isc:geom%iec  ,geom%jsc:geom%jec+1,1:geom%npz)
 real(kind=kind_real) :: vd(geom%isc:geom%iec+1,geom%jsc:geom%jec  ,1:geom%npz)
 integer :: i, j, k

 ! --------------------------------- !
 ! Convert A-grid winds to vorticity !
 ! --------------------------------- !

 ! A to D grid winds
 ! -----------------
 call a2d(geom, ua, va, ud, vd)

 ! Compute vorticity
 ! -----------------
 call udvd_to_vort(geom, ud, vd, vort)

end subroutine uava_to_vort

! ------------------------------------------------------------------------------

subroutine udvd_to_ucvc(geom,u,v,uc,vc,ua_out,va_out)

implicit none
type(fv3jedi_geom), target,     intent(inout) :: geom
real(kind=kind_real),           intent(in)    ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1) !Dgrid winds
real(kind=kind_real),           intent(in)    ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed  ) !Dgrid winds
real(kind=kind_real),           intent(out)   :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  ) !Cgrid winds
real(kind=kind_real),           intent(out)   :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1) !Cgrid winds
real(kind=kind_real), optional, intent(out)   :: ua_out(geom%isd:geom%ied  ,geom%jsd:geom%jed) !Agrid winds
real(kind=kind_real), optional, intent(out)   :: va_out(geom%isd:geom%ied  ,geom%jsd:geom%jed) !Agrid winds

! Normally input/output
real(kind=kind_real) :: ut  (geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: vt  (geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: utmp(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: vtmp(geom%isd:geom%ied  ,geom%jsd:geom%jed  )
real(kind=kind_real) :: ua  (geom%isd:geom%ied  ,geom%jsd:geom%jed )
real(kind=kind_real) :: va  (geom%isd:geom%ied  ,geom%jsd:geom%jed )
logical :: dord4
logical :: nested
integer :: grid_type
integer :: npx, npy

! Locals
real(kind=kind_real), parameter:: a1 =  0.5625_kind_real
real(kind=kind_real), parameter:: a2 = -0.0625_kind_real
real(kind=kind_real), parameter:: c1 = -2.0_kind_real/14.0_kind_real
real(kind=kind_real), parameter:: c2 = 11.0_kind_real/14.0_kind_real
real(kind=kind_real), parameter:: c3 =  5.0_kind_real/14.0_kind_real

integer npt, i, j, ifirst, ilast, id
integer :: is,  ie,  js,  je
integer :: isd, ied, jsd, jed
real(kind=kind_real), pointer, dimension(:,:,:) :: sin_sg
real(kind=kind_real), pointer, dimension(:,:)   :: cosa_u, cosa_v, cosa_s
real(kind=kind_real), pointer, dimension(:,:)   :: rsin_u, rsin_v, rsin2
real(kind=kind_real), pointer, dimension(:,:)   :: dxa,dya

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

dord4 = geom%dord4
nested = geom%nested
grid_type = geom%grid_type

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

if (grid_type < 3 .and. .not. nested) then
  npt = 4
else
  npt = -2
endif

! Initialize the non-existing corner regions
utmp(:,:) = 1.0e30_kind_real
vtmp(:,:) = 1.0e30_kind_real

if ( nested) then

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
  if (grid_type < 3) then

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

if (grid_type < 3 .and. .not. nested) then
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

if (grid_type < 3) then
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

  if( is==1 .and. .not. nested  ) then
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

  if( (ie+1)==npx  .and. .not. nested ) then
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

if (grid_type < 3) then

  do j=js-1,je+2
    if ( j==1 .and. .not. nested  ) then
      do i=is-1,ie+1
        vt(i,j) = edge_interpolate4(va(i,-1:2), dya(i,-1:2))
        if (vt(i,j) > 0.) then
          vc(i,j) = vt(i,j)*sin_sg(i,j-1,4)
        else
          vc(i,j) = vt(i,j)*sin_sg(i,j,2)
        end if
      enddo
    elseif ( j==0 .or. j==(npy-1) .and. .not. nested  ) then
      do i=is-1,ie+1
        vc(i,j) = c1*vtmp(i,j-2) + c2*vtmp(i,j-1) + c3*vtmp(i,j)
        vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
      enddo
    elseif ( j==2 .or. j==(npy+1)  .and. .not. nested ) then
      do i=is-1,ie+1
        vc(i,j) = c1*vtmp(i,j+1) + c2*vtmp(i,j) + c3*vtmp(i,j-1)
        vt(i,j) = (vc(i,j) - u(i,j)*cosa_v(i,j))*rsin_v(i,j)
      enddo
    elseif ( j==npy .and. .not. nested  ) then
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

if (present(ua_out) .and. present(va_out)) then
  ua_out(is:ie,js:je) = ua(is:ie,js:je)
  va_out(is:ie,js:je) = va(is:ie,js:je)
endif

end subroutine udvd_to_ucvc

! ------------------------------------------------------------------------------

real(kind=kind_real) function edge_interpolate4(ua, dxa)

implicit none
real(kind=kind_real), intent(in) :: ua(4)
real(kind=kind_real), intent(in) :: dxa(4)
real(kind=kind_real):: t1, t2

t1 = dxa(1) + dxa(2)
t2 = dxa(3) + dxa(4)

edge_interpolate4 = 0.5*( ((t1+dxa(2))*ua(2)-dxa(2)*ua(1)) / t1 + &
                          ((t2+dxa(3))*ua(3)-dxa(3)*ua(4)) / t2 )

end function edge_interpolate4

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
 if (.not. geom%bounded_domain) then
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
 end if  !clt bounded=.false.
 
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

! TODO: OpenMP Re-enable this parallel section by fixing variable shared/private assignment
! !$OMP parallel do default(none) shared(is,ie,js,je,npz,npx,npy,c2,c1, &
! !$OMP                                  u,v,ua,va)         &
! !$OMP                          private(utmp, vtmp, wu, wv)
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

subroutine fill_dgrid_winds(geom,u,v,fillhalo)

implicit none
type(fv3jedi_geom),   intent(inout) :: geom
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

! Fill halos
! ----------
if (present(fillhalo)) then
  if (fillhalo) then
    call mpp_update_domains(u, v, geom%domain, gridtype=DGRID_NE)
  endif
endif

end subroutine fill_dgrid_winds

! ------------------------------------------------------------------------------

subroutine fill_cgrid_winds(geom,uc,vc,fillhalo)

implicit none
type(fv3jedi_geom),   intent(inout) :: geom
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

! Fill halos
! ----------
if (present(fillhalo)) then
  if (fillhalo) then
    call mpp_update_domains(uc, vc, geom%domain, gridtype=CGRID_NE)
  endif
endif

end subroutine fill_cgrid_winds

! ------------------------------------------------------------------------------

subroutine compute_utvt(geom, uc, vc, ut, vt, dt)

 implicit none
 type(fv3jedi_geom),   intent(inout) :: geom
 real(kind=kind_real), intent(in   ) ::  uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
 real(kind=kind_real), intent(in   ) ::  vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
 real(kind=kind_real), intent(inout) ::  ut(geom%isd:geom%ied+1,geom%jsd:geom%jed  )
 real(kind=kind_real), intent(inout) ::  vt(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)
 real(kind=kind_real), intent(in   ) :: dt

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

end subroutine compute_utvt

! ------------------------------------------------------------------------------

subroutine ctoa(geom, uc, vc, ua, va)

 implicit none
 type(fv3jedi_geom),   intent(in)  :: geom
 real(kind=kind_real), intent(in)  :: uc(geom%isd:geom%ied+1,geom%jsd:geom%jed  )    ! c-grid u-wind field
 real(kind=kind_real), intent(in)  :: vc(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)    ! c-grid v-wind field
 real(kind=kind_real), intent(out) :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )    ! a-grid u-wind field
 real(kind=kind_real), intent(out) :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )    ! a-grid v-wind field

 integer :: i,j
 real(kind=kind_real) :: tmp1i(geom%isd:geom%ied+1)
 real(kind=kind_real) :: tmp2i(geom%isd:geom%ied+1)
 real(kind=kind_real) :: tmp1j(geom%jsd:geom%jed+1)
 real(kind=kind_real) :: tmp2j(geom%jsd:geom%jed+1)

 do i=geom%isd,geom%ied
    tmp1j(:) = 0.0_kind_real
    tmp2j(:) = vc(i,:)*geom%dx(i,:)
    call interp_left_edge_1d(geom%jsd, geom%jed+1, tmp2j, tmp1j)
    va(i,geom%jsd:geom%jed) = tmp1j(geom%jsd+1:geom%jed+1)/geom%dxa(i,geom%jsd:geom%jed)
 enddo

 do j=geom%jsd,geom%jed
    tmp1i(:) = 0.0_kind_real
    tmp2i(:) = uc(:,j)*geom%dy(:,j)
    call interp_left_edge_1d(geom%isd, geom%ied+1, tmp2i, tmp1i)
    ua(geom%isd:geom%ied,j) = tmp1i(geom%isd+1:geom%ied+1)/geom%dya(geom%isd:geom%ied,j)
 enddo

end subroutine ctoa

! ------------------------------------------------------------------------------

subroutine atod(geom, ua, va, u, v)

 implicit none
 type(fv3jedi_geom),   intent(in)  :: geom
 real(kind=kind_real), intent(in)  :: ua(geom%isd:geom%ied  ,geom%jsd:geom%jed  )    ! a-grid u-wind field
 real(kind=kind_real), intent(in)  :: va(geom%isd:geom%ied  ,geom%jsd:geom%jed  )    ! a-grid v-wind field
 real(kind=kind_real), intent(out) ::  u(geom%isd:geom%ied  ,geom%jsd:geom%jed+1)    ! d-grid u-wind field
 real(kind=kind_real), intent(out) ::  v(geom%isd:geom%ied+1,geom%jsd:geom%jed  )    ! d-grid v-wind field

 integer :: i,j
 real(kind=kind_real) :: tmp1i(geom%isd:geom%ied+1)
 real(kind=kind_real) :: tmp2i(geom%isd:geom%ied)
 real(kind=kind_real) :: tmp1j(geom%jsd:geom%jed+1)
 real(kind=kind_real) :: tmp2j(geom%jsd:geom%jed)

 do j=geom%jsd+1,geom%jed
    tmp1i(:) = 0.0
    tmp2i(:) = va(:,j)*geom%dxa(:,j)
    call interp_left_edge_1d(geom%isd, geom%ied, tmp1i, tmp2i)
    v(:,j) = tmp1i(:)/geom%dxc(:,j)
 enddo
 do i=geom%isd+1,geom%ied
    tmp1j(:) = 0.0
    tmp2j(:) = ua(i,:)*geom%dya(i,:)
    call interp_left_edge_1d(geom%jsd, geom%jed, tmp1j, tmp2j)
    u(i,:) = tmp1j(:)/geom%dyc(i,:)
 enddo

end subroutine atod

! ------------------------------------------------------------------------------

subroutine interp_left_edge_1d(ifirst, ilast, qin, qout)

 implicit none
 integer,              intent(in)  :: ifirst,ilast
 real(kind=kind_real), intent(in)  ::  qin(ifirst:ilast)
 real(kind=kind_real), intent(out) :: qout(ifirst:ilast)

 integer :: i

 qout(:) = 0.0_kind_real
 do i=ifirst+1,ilast
   qout(i) = 0.5_kind_real * (qin(i-1) + qin(i))
 enddo

end subroutine interp_left_edge_1d

! ------------------------------------------------------------------------------

end module wind_vt_mod
