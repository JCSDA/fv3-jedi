! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module pressure_vt_mod

use fv3jedi_constants_mod, only: kappa
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

implicit none
private

public compute_fv3_pressures
public compute_fv3_pressures_tlm
public compute_fv3_pressures_bwd
public delp_to_pe_p_logp
public pe_to_delp
public delp_to_pe
public pe_to_pk
public ps_to_delp
public ps_to_pe
public ps_to_delp_tl
public ps_to_delp_ad

contains

!----------------------------------------------------------------------------
! Compute all pressures needed as input to fv3 ------------------------------
!----------------------------------------------------------------------------

subroutine compute_fv3_pressures( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                  kappa, ptop, delp, pe, pk, pkz, peln )

 implicit none
 integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
 real(kind=kind_real), intent(in) :: kappa, ptop
 real(kind=kind_real), intent(in) :: delp(isd:ied, jsd:jed, npz)
 real(kind=kind_real), intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
 real(kind=kind_real), intent(out) :: pk(is:ie, js:je, npz+1)
 real(kind=kind_real), intent(out) :: peln(is:ie, npz+1, js:je)
 real(kind=kind_real), intent(out) :: pkz(is:ie, js:je, npz)
 integer :: i, j, k

 pe(:, :, :) = 0.0
 pe(:, 1, :) = ptop
 do k=2,npz+1
   do j=js,je
     do i=is,ie
       pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
     end do
   end do
 end do

 do k=1,npz+1
   do j=js,je
     do i=is,ie
       peln(i, k, j) = log(pe(i, k, j))
     end do
   end do
 end do

 do k=1,npz+1
   do j=js,je
     do i=is,ie
       pk(i, j, k) = exp(kappa*peln(i, k, j))
     end do
   end do
 end do

 do k=1,npz
   do j=js,je
     do i=is,ie
       pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1, j)-peln(i, k, j)))
     end do
   end do
 end do

end subroutine compute_fv3_pressures

subroutine compute_fv3_pressures_tlm( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                      kappa, ptop, delp, delp_tl, &
                                      pe, pe_tl, pk, pk_tl, pkz, pkz_tl, peln, peln_tl )

 implicit none
 integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
 real(kind=kind_real), intent(in) :: kappa, ptop
 real(kind=kind_real), intent(in) :: delp(isd:ied, jsd:jed, npz)
 real(kind=kind_real), intent(in) :: delp_tl(isd:ied, jsd:jed, npz)
 real(kind=kind_real), intent(out) :: pe(is-1:ie+1, npz+1, js-1:je+1)
 real(kind=kind_real), intent(out) :: pe_tl(is-1:ie+1, npz+1, js-1:je+1)
 real(kind=kind_real), intent(out) :: pk(is:ie, js:je, npz+1)
 real(kind=kind_real), intent(out) :: pk_tl(is:ie, js:je, npz+1)
 real(kind=kind_real), intent(out) :: peln(is:ie, npz+1, js:je)
 real(kind=kind_real), intent(out) :: peln_tl(is:ie, npz+1, js:je)
 real(kind=kind_real), intent(out) :: pkz(is:ie, js:je, npz)
 real(kind=kind_real), intent(out) :: pkz_tl(is:ie, js:je, npz)
 integer :: i, j, k

 pe(:, :, :) = 0.0
 pe(:, 1, :) = ptop
 pe_tl = 0.0
 do k=2,npz+1
   do j=js,je
     do i=is,ie
       pe_tl(i, k, j) = pe_tl(i, k-1, j) + delp_tl(i, j, k-1)
       pe(i, k, j) = pe(i, k-1, j) + delp(i, j, k-1)
     end do
   end do
 end do

 peln_tl = 0.0
 do k=1,npz+1
   do j=js,je
     do i=is,ie
       peln_tl(i, k, j) = pe_tl(i, k, j)/pe(i, k, j)
       peln(i, k, j) = log(pe(i, k, j))
     end do
   end do
 end do

 pk_tl = 0.0
 do k=1,npz+1
   do j=js,je
     do i=is,ie
       pk_tl(i, j, k) = kappa*peln_tl(i, k, j)*exp(kappa*peln(i, k, j))
       pk(i, j, k) = exp(kappa*peln(i, k, j))
     end do
   end do
 end do

 pkz_tl = 0.0
 do k=1,npz
   do j=js,je
     do i=is,ie
       pkz_tl(i, j, k) = ((pk_tl(i, j, k+1)-pk_tl(i, j, k))*kappa*(peln   (i, k+1, j)-peln   (i, k, j)) &
                         -(pk   (i, j, k+1)-pk   (i, j, k))*kappa*(peln_tl(i, k+1, j)-peln_tl(i, k, j))) &
                         /(kappa*(peln(i, k+1, j)-peln(i, k, j)))**2
       pkz(i, j, k) = (pk(i, j, k+1)-pk(i, j, k))/(kappa*(peln(i, k+1, j)-peln(i, k, j)))
     end do
   end do
 end do

end subroutine compute_fv3_pressures_tlm

subroutine compute_fv3_pressures_bwd( is, ie, js, je, isd, ied, jsd, jed, npz, &
                                      kappa, ptop, delp, delp_ad, &
                                      pe, pe_ad, pk, pk_ad, pkz, pkz_ad, peln, peln_ad)

  implicit none
  integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed, npz
  real(kind=kind_real), intent(in) :: kappa, ptop
  real(kind=kind_real), intent(in)    :: delp(isd:ied, jsd:jed, npz)
  real(kind=kind_real), intent(inout) :: delp_ad(isd:ied, jsd:jed, npz)
  real(kind=kind_real), intent(in)    :: pe(is-1:ie+1, npz+1, js-1:je+1)
  real(kind=kind_real), intent(inout) :: pe_ad(is-1:ie+1, npz+1, js-1:je+1)
  real(kind=kind_real), intent(in)    :: pk(is:ie, js:je, npz+1)
  real(kind=kind_real), intent(inout) :: pk_ad(is:ie, js:je, npz+1)
  real(kind=kind_real), intent(in)    :: peln(is:ie, npz+1, js:je)
  real(kind=kind_real), intent(inout) :: peln_ad(is:ie, npz+1, js:je)
  real(kind=kind_real), intent(in)    :: pkz(is:ie, js:je, npz)
  real(kind=kind_real), intent(inout) :: pkz_ad(is:ie, js:je, npz)
  integer :: i, j, k
  real(kind=kind_real) :: temp_tj1, temp_ad1, temp_ad2

  do k=npz,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        temp_tj1 = kappa*(peln(i, k+1, j)-peln(i, k, j))
        temp_ad1 = pkz_ad(i, j, k)/temp_tj1
        temp_ad2 = -((pk(i, j, k+1)-pk(i, j, k))*kappa*temp_ad1/temp_tj1)
        pk_ad(i, j, k+1) = pk_ad(i, j, k+1) + temp_ad1
        pk_ad(i, j, k) = pk_ad(i, j, k) - temp_ad1
        peln_ad(i, k+1, j) = peln_ad(i, k+1, j) + temp_ad2
        peln_ad(i, k, j) = peln_ad(i, k, j) - temp_ad2
        pkz_ad(i, j, k) = 0.0
      end do
    end do
  end do

  do k=npz+1,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        peln_ad(i, k, j) = peln_ad(i, k, j) + exp(kappa*peln(i, k, j))*kappa*pk_ad(i, j, k)
        pk_ad(i, j, k) = 0.0
      end do
    end do
  end do

  do k=npz+1,1,-1
    do j=je,js,-1
      do i=ie,is,-1
        pe_ad(i, k, j) = pe_ad(i, k, j) + peln_ad(i, k, j)/pe(i, k, j)
        peln_ad(i, k, j) = 0.0
      end do
    end do
  end do

  do k=npz+1,2,-1
    do j=je,js,-1
      do i=ie,is,-1
        pe_ad(i, k-1, j) = pe_ad(i, k-1, j) + pe_ad(i, k, j)
        delp_ad(i, j, k-1) = delp_ad(i, j, k-1) + pe_ad(i, k, j)
        pe_ad(i, k, j) = 0.0
      end do
    end do
  end do

end subroutine compute_fv3_pressures_bwd

!----------------------------------------------------------------------------
! Pressure thickness to pressure (edge), pressure (mid) and log p (mid) -----
!----------------------------------------------------------------------------

subroutine delp_to_pe_p_logp(geom,delp,pe,p,logp)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness
 real(kind=kind_real), intent(out) ::   pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure edge/interface
 real(kind=kind_real), intent(out) ::    p(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure mid point
 real(kind=kind_real), optional, intent(out) :: logp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Log of pressure mid point

 !Locals
 integer :: isc,iec,jsc,jec,k

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

 !Pressure at layer edge
 pe(isc:iec,jsc:jec,1) = geom%ptop
 do k = 2,geom%npz+1
   pe(isc:iec,jsc:jec,k) = pe(isc:iec,jsc:jec,k-1) + delp(isc:iec,jsc:jec,k-1)
 enddo

 !Midpoint pressure
 p(isc:iec,jsc:jec,:) = 0.5*(pe(isc:iec,jsc:jec,2:geom%npz+1) + pe(isc:iec,jsc:jec,1:geom%npz))

 if (present(logp)) then
   !Log pressure
   logp(isc:iec,jsc:jec,:) = log(p(isc:iec,jsc:jec,:))
 endif

end subroutine delp_to_pe_p_logp

!----------------------------------------------------------------------------

subroutine pe_to_delp(geom,pe,delp)

 implicit none
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in ) ::   pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure edge/interface
 real(kind=kind_real), intent(out) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness

 !Locals
 integer :: isc,iec,jsc,jec,k

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

 !Pressure at layer edge
 do k = 1,geom%npz
   delp(isc:iec,jsc:jec,k) = pe(isc:iec,jsc:jec,k+1) - pe(isc:iec,jsc:jec,k)
 enddo

end subroutine pe_to_delp

!----------------------------------------------------------------------------

subroutine delp_to_pe( geom, delp, pe )

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
real(kind=kind_real), intent(in ) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness
real(kind=kind_real), intent(out) ::   pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure edge/interface

!Locals
integer :: isc,iec,jsc,jec,k

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec

!Pressure at layer edge
pe(isc:iec,jsc:jec,1) = geom%ptop
do k = 2,geom%npz+1
  pe(isc:iec,jsc:jec,k) = pe(isc:iec,jsc:jec,k-1) + delp(isc:iec,jsc:jec,k-1)
enddo

end subroutine delp_to_pe

!----------------------------------------------------------------------------

subroutine pe_to_pk( geom, pe, pk )

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
real(kind=kind_real), intent(in ) :: pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure edge/interface
real(kind=kind_real), intent(out) :: pk(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure to the kappa

!Locals
integer :: i, j, k
real(kind=kind_real) :: pel1, pel2
real(kind=kind_real) :: pek1, pek2

do k=1,geom%npz
  do j = geom%jsc,geom%jec
    do i = geom%isc,geom%iec

      pel1 = log(pe(i,j,k+1))
      pel2 = log(pe(i,j,k))

      pek1 = exp(kappa*pel1)
      pek2 = exp(kappa*pel2)

      pk(i,j,k) = (pek1-pek2)/(kappa*(pel1-pel2))

    end do
  end do
end do

end subroutine pe_to_pk

!----------------------------------------------------------------------------

subroutine ps_to_pe(geom,ps,pe)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in   ) :: ps(geom%isc:geom%iec,geom%jsc:geom%jec,1           )   !Surface pressure
 real(kind=kind_real), intent(inout) :: pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1)   !Pressure thickness

 integer :: k

 do k = 1,geom%npz+1
   pe(:,:,k) = geom%ak(k) + geom%bk(k) * ps(:,:,1)
 enddo

endsubroutine ps_to_pe

!----------------------------------------------------------------------------

subroutine ps_to_delp(geom,ps,delp)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in   ) ::   ps(geom%isc:geom%iec,geom%jsc:geom%jec           )   !Surface pressure
 real(kind=kind_real), intent(inout) :: delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness

 integer :: isc,iec,jsc,jec,i,j,k

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

 do k = 1,geom%npz
   do j = jsc,jec
     do i = isc,iec
       delp(i,j,k) = geom%ak(k+1) + geom%bk(k+1)*ps(i,j) - (geom%ak(k) + geom%bk(k)*ps(i,j))
     enddo
   enddo
 enddo

endsubroutine ps_to_delp

subroutine ps_to_delp_tl(geom,ps_tl,delp_tl)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in   ) ::   ps_tl(geom%isc:geom%iec,geom%jsc:geom%jec           )   !Surface pressure
 real(kind=kind_real), intent(inout) :: delp_tl(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness

 integer :: isc,iec,jsc,jec,i,j,k

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

 delp_tl = 0.0_kind_real
 do k = 1,geom%npz
   do j = jsc,jec
     do i = isc,iec
       delp_tl(i,j,k) = geom%bk(k+1)*ps_tl(i,j) - geom%bk(k)*ps_tl(i,j)
     enddo
   enddo
 enddo

endsubroutine ps_to_delp_tl

subroutine ps_to_delp_ad(geom,ps_ad,delp_ad)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom !Geometry for the model
 real(kind=kind_real), intent(inout) ::   ps_ad(geom%isc:geom%iec,geom%jsc:geom%jec           )   !Surface pressure
 real(kind=kind_real), intent(inout) :: delp_ad(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure thickness

 integer :: isc,iec,jsc,jec,i,j,k

 isc = geom%isc
 iec = geom%iec
 jsc = geom%jsc
 jec = geom%jec

 ps_ad = 0.0_kind_real

 do k = geom%npz,1,-1
   do j = jec,jsc,-1
     do i = iec,isc,-1
       ps_ad(i,j) = ps_ad(i,j) + (geom%bk(k+1) - geom%bk(k))*delp_ad(i,j,k)
     enddo
   enddo
 enddo

endsubroutine ps_to_delp_ad

end module pressure_vt_mod
