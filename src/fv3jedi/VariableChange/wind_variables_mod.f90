! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wind_vt_mod

use fv3jedi_constants_mod, only: pi, rad2deg
use fv3jedi_geom_mod,  only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

use fv_mp_adm_mod, only: mpp_update_domains_adm
use mpp_domains_mod, only: mpp_update_domains, dgrid_ne
use mpp_domains_mod, only: mpp_get_boundary, mpp_get_boundary_ad

implicit none
private

public sfc_10m_winds

public psichi_to_uava
public psichi_to_uava_adm

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
