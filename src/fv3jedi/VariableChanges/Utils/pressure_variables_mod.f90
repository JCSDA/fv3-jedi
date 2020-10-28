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

public delp_to_pe_p_logp
public pe_to_pkz
public ps_to_pkz
public pe_to_delp
public delp_to_pe
public pe_to_pk
public ps_to_delp
public ps_to_pe
public ps_to_delp_tl
public ps_to_delp_ad

contains

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

subroutine pe_to_pkz(geom,pe,pkz)

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
real(kind=kind_real), intent(in ) ::  pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure edge/interface
real(kind=kind_real), intent(out) :: pkz(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure to the kappa

integer :: k
real(kind=kind_real) :: peln(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1)
real(kind=kind_real) ::   pk(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1)

peln = log(pe)
pk = exp(kappa*peln)

do k=1,geom%npz
  pkz(:,:,k) = (pk(:,:,k+1)-pk(:,:,k)) / (kappa*(peln(:,:,k+1)-peln(:,:,k)))
enddo

end subroutine pe_to_pkz

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

subroutine ps_to_pkz(geom,ps,pkz)

 implicit none
 type(fv3jedi_geom)  , intent(in   ) :: geom !Geometry for the model
 real(kind=kind_real), intent(in   ) ::  ps(geom%isc:geom%iec,geom%jsc:geom%jec,1         )   !Surface pressure
 real(kind=kind_real), intent(inout) :: pkz(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz)   !Pressure ^ kappa

 integer :: k
 real(kind=kind_real) :: pe(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1) !Pressure

 call ps_to_pe(geom, ps, pe)
 call pe_to_pkz(geom, pe, pkz)

endsubroutine ps_to_pkz

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
