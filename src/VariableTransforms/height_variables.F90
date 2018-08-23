module height_vt_mod

use fv3jedi_constants, only: grav, rvap, tice,rdry, zvir 
use fv3jedi_geom_mod,  only: fv3jedi_geom
use kinds,             only: kind_real

implicit none
public

contains
subroutine geop_height(geom,prs,prsi,T,q,phis,use_compress,gph)

implicit none
type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
real(kind_real), intent(in ) :: prs(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)    !mid layerpressure
real(kind_real), intent(in ) :: prsi(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz+1) ! interface pressure
real(kind_real), intent(in ) :: phis(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)              !Surface geopotential (grav*Z_sfc)
real(kind_real), intent(in ) :: T(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)            
real(kind_real), intent(in ) :: q(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)     ! specific humidity
real(kind_real), intent(out) :: gph(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)   !geopotential height (meters)

!locals
real(kind_real)       :: Tv(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz)
real(kind_real)       :: qmr(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,1:geom%npz) ! mixing ratio|kg/kg
logical               :: use_compress
integer               :: isc,iec,jsc,jec,npz,i,j,k
real(kind=kind_real)  :: Tkk,Tvk,Tc, qmk,Pak,dpk,dz
real(kind=kind_real)  :: prs_sv, prs_v           
real(kind=kind_real)  :: ehn_fct,x_v,cmpr

! Constants  from GSI
! Constants for compressibility factor (Davis et al 1992)
real(kind_real),parameter ::  cpf_a0 =  1.58123e-6_kind_real ! K/Pa
real(kind_real),parameter ::  cpf_a1 = -2.9331e-8_kind_real  ! 1/Pa
real(kind_real),parameter ::  cpf_a2 =  1.1043e-10_kind_real ! 1/K 1/Pa
real(kind_real),parameter ::  cpf_b0 =  5.707e-6_kind_real   ! K/Pa
real(kind_real),parameter ::  cpf_b1 = -2.051e-8_kind_real   ! 1/Pa
real(kind_real),parameter ::  cpf_c0 =  1.9898e-4_kind_real  ! K/Pa
real(kind_real),parameter ::  cpf_c1 = -2.376e-6_kind_real   ! 1/Pa
real(kind_real),parameter ::  cpf_d  =  1.83e-11_kind_real   ! K2/Pa2
real(kind_real),parameter ::  cpf_e  = -0.765e-8_kind_real   ! K2/Pa2

real(kind_real),parameter ::  psv_a =  1.2378847e-5_kind_real       !  (1/K2)
real(kind_real),parameter ::  psv_b = -1.9121316e-2_kind_real       !  (1/K)
real(kind_real),parameter ::  psv_c = 33.93711047_kind_real         !
real(kind_real),parameter ::  psv_d = -6.3431645e+3_kind_real       !  (K)
! Constants for enhancement factor to calculating the mole fraction of water vapor
real(kind_real),parameter ::  ef_alpha = 1.00062_kind_real           !
real(kind_real),parameter ::  ef_beta  = 3.14e-8_kind_real           !  (1/Pa)
real(kind_real),parameter ::  ef_gamma = 5.6e-7_kind_real            !  (1/K2)

isc = geom%bd%isc
iec = geom%bd%iec
jsc = geom%bd%jsc
jec = geom%bd%jec
npz = geom%npz


!get qmr--mixing ratio and virtual temeprature
qmr(isc:iec,jsc:jec,:) = q(isc:iec,jsc:jec,:)/(1-q(isc:iec,jsc:jec,:))
Tv(isc:iec,jsc:jec,:)  = T(isc:iec,jsc:jec,:)*(1.0 + zvir*qmr(isc:iec,jsc:jec,:))

if (use_compress) then
 
!  Compute compressibility factor (Picard et al 2008) and geopotential heights at midpoint 
  do j = jsc,jec
  do i = isc,iec

     do k = geom%npz, 1, -1
        if ( k == geom%npz) then
           Tkk  = T(i,j,k)
           Tvk  = Tv(i,j,k)
           Pak  = exp(0.5_kind_real*(log(prsi(i,j,k+1))+log(prs(i,j,k))))
           dpk  = prsi(i,j,k+1)/prs(i,j,k)
        else 
           Tkk  = 0.5_kind_real * ( T(i,j,k+1) +  T(i,j,k) )
           Tvk  = 0.5_kind_real * (Tv(i,j,k+1) + Tv(i,j,k) )
           Pak  = exp(0.5_kind_real*(log(prs(i,j,k+1))+log(prs(i,j,k)))) 
           dpk  = prs(i,j,k+1)/prs(i,j,k)
        end if

        Tc   = Tkk - tice
        qmk  = qmr(i,j,k)
        prs_sv  = exp(psv_a*Tkk**2 + psv_b*Tkk + psv_c + psv_d/Tkk ) ! Pvap sat, eq A1.1 (Pa)
        ehn_fct = ef_alpha + ef_beta*Pak + ef_gamma*Tc**2 ! enhancement factor (eq. A1.2)
        prs_v   = qmk* Pak/(1.0+qmk*rdry/rvap)            ! vapor pressure (Pa)
        x_v     = prs_v/prs_sv * ehn_fct * prs_sv/Pak     ! molar fraction of water vapor (eq. A1.3)
!       Compressibility factor (eq A1.4 from Picard et al 2008)
        cmpr    = 1.0_kind_real - (Pak/Tkk) * (cpf_a0 + cpf_a1*Tc + cpf_a2*Tc**2 &
                    + (cpf_b0 + cpf_b1*Tc)*x_v + (cpf_c0 + cpf_c1*Tc)*x_v**2 ) &
                    + (Pak**2/Tkk**2) * (cpf_d + cpf_e*x_v**2)
        dz      = rdry/grav * Tvk * cmpr * log(dpk)
        if ( k == geom%npz) then
          gph(i,j,k) = phis(i,j)/grav + dz
        else
          gph(i,j,k) = gph(i,j,k+1) + dz
        end if
      end do ! end k loop

    end do ! end i loop
    end do ! end j loop

else  ! not use compressivity

  do j = jsc,jec
  do i = isc,iec

     k  = geom%npz
     dz         = rdry/grav * Tv(i,j,k) * log(prsi(i,j,k+1)/prs(i,j,k))
     gph(i,j,k) = phis(i,j)/grav + dz
 
     do k = geom%npz-1, 1, -1
        dz         = rdry/grav * 0.5_kind_real * (Tv(i,j,k+1)+Tv(i,j,k)) * log(prs(i,j,k+1)/prs(i,j,k))
        gph(i,j,k) = gph(i,j,k+1) + dz
     end do

   end do
   end do

end if

end subroutine geop_height


end module height_vt_mod
