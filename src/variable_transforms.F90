module variable_transforms

use fv3jedi_constants, only: kappa, epsilon
use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real

implicit none

public

contains

subroutine dpppk_tj(geom, ptop, delp, p, pe, pk, pke, pkco)

  !Save the trajectory for when converting pressure variable

  implicit none

  type(fv3jedi_geom)  , intent(in)   :: geom !Geometry for the model

  real(kind=kind_real), intent(in ) :: ptop
  real(kind=kind_real), intent(in ) :: delp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )   !Pressure thickness
  real(kind=kind_real), intent(out) ::    p(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )   !Pressure (mid point)
  real(kind=kind_real), intent(out) ::   pe(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)   !Pressure (edges)
  real(kind=kind_real), intent(out) ::   pk(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )   !Pressure to the kappa (mid point)
  real(kind=kind_real), intent(out) ::  pke(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)   !Pressure to the kappa (edges)
  real(kind=kind_real), intent(out) :: pkco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,2)   !Saved coeffs for pkp

  integer :: k

  pe = 0.0
  pk = 0.0
  pke = 0.0
  pkco = 0.0

  !Pressure at layer edge
  pe(:,:,1) = ptop
  do k = 2,geom%nlevs+1
    pe(:,:,k) = pe(:,:,k-1) + delp(:,:,k-1)
  enddo

  !Pressure at layer mid point
  p = 0.5*(pe(:,:,2:geom%nlevs+1) + pe(:,:,1:geom%nlevs))

  !Pressure to the kappa at layer edge
  pke = pe**kappa

  !Multiplication coefficients used in tl and ad of pk computation
  pkco(:,:,:,1)  = 1.0/(kappa*(log(pe(:,:,2:geom%nlevs+1)) - log(pe(:,:,1:geom%nlevs))))
  pkco(:,:,:,2)  = pkco(:,:,:,1)*pk*kappa

  !Pressure to the kappa at layer mid point
  pk = pkco(:,:,:,1)*(pke(:,:,2:geom%nlevs+1) - pke(:,:,1:geom%nlevs))

end subroutine dpppk_tj

! ------------------------------------------------------------------------------

subroutine delp2lnp(geom,ptop,delp)

  !Tangent linear of conversion from delp to log(p)

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom
  real(kind=kind_real), intent(in)  :: ptop
  real(kind=kind_real), intent(inout)  :: delp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  p(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure (mid)
  real(kind=kind_real) :: pe(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Pressure at the edges
  pe(:,:,1) = ptop
  do k = 2,geom%nlevs+1
    pe(:,:,k) = pe(:,:,k-1) + delp(:,:,k-1)
  enddo

  !Pressure at mid points
  p = 0.5*(pe(:,:,2:geom%nlevs+1) + pe(:,:,1:geom%nlevs))

  !Overwrite with log pressure
  delp = log(p)

endsubroutine delp2lnp

! ------------------------------------------------------------------------------

subroutine delp2lnp_tl(geom,p,delpp)

  !Tangent linear of conversion from delp to log(p)

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Trajectory
  real(kind=kind_real), intent(in)  ::     p(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure (edges)
  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) ::   pp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure (edges)
  integer :: k

  !Linearized pressure at the edges
  pep(:,:,1) = 0.0
  do k = 2,geom%nlevs+1
    pep(:,:,k) = pep(:,:,k-1) + delpp(:,:,k-1)
  enddo

  !Linearized pressure at mid points
  pp = 0.5*(pep(:,:,2:geom%nlevs+1) + pep(:,:,1:geom%nlevs))

  !Linearized overwrite with log pressure
  delpp = pp / p

endsubroutine delp2lnp_tl

! ------------------------------------------------------------------------------

subroutine delp2lnp_ad(geom,p,delpp)

  !Adjoint of conversion from delp to log(p)

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Trajectory
  real(kind=kind_real), intent(in)  ::     p(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure (edges)
  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) ::   pp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure (edges)
  real(kind=kind_real) :: lnpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Log pressure
  integer :: k

  !Adjoint overwrite with log pressure
  pp = delpp / p

  !Adjoint pressure at mid points
  pep = 0.0_8
  pep(:, :, 2:geom%nlevs+1) = pep(:, :, 2:geom%nlevs+1) + 0.5*pp
  pep(:, :, 1:geom%nlevs  ) = pep(:, :, 1:geom%nlevs  ) + 0.5*pp

  !Adjoint pressure at the edges
  delpp = 0.0_8
  DO k=geom%nlevs+1,2,-1
    pep  (:, :, k-1) = pep  (:, :, k-1) + pep(:, :, k)
    delpp(:, :, k-1) =                    pep(:, :, k)
  END DO

endsubroutine delp2lnp_ad

! ------------------------------------------------------------------------------

subroutine delp2p(geom,ptop,delp)

  !Tangent linear of conversion from delp to p at layer mid points 

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom
  real(kind=kind_real), intent(in)  :: ptop

  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pe(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Linearized pressure at the edges
  pe(:,:,1) = ptop
  do k = 2,geom%nlevs+1
    pe(:,:,k) = pe(:,:,k-1) + delp(:,:,k-1)
  enddo

  !Linearized overwrite with pressure at mid points
  delp = 0.5*(pe(:,:,2:geom%nlevs+1) + pe(:,:,1:geom%nlevs))

endsubroutine delp2p

! ------------------------------------------------------------------------------

subroutine delp2p_tl(geom,delpp)

  !Tangent linear of conversion from delp to p at layer mid points 

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Linearized pressure at the edges
  pep(:,:,1) = 0.0
  do k = 2,geom%nlevs+1
    pep(:,:,k) = pep(:,:,k-1) + delpp(:,:,k-1)
  enddo

  !Linearized overwrite with pressure at mid points
  delpp = 0.5*(pep(:,:,2:geom%nlevs+1) + pep(:,:,1:geom%nlevs))

endsubroutine delp2p_tl

! ------------------------------------------------------------------------------

subroutine delp2p_ad(geom,delpp)

  !Adjoint of conversion from delp to p at layer mid points 

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Adjoint overwrite with pressure at mid points
  pep = 0.0_8
  pep(:, :, 2:geom%nlevs+1) = pep(:, :, 2:geom%nlevs+1) + 0.5*delpp
  pep(:, :, 1:geom%nlevs  ) = pep(:, :, 1:geom%nlevs  ) + 0.5*delpp

  !Adjoint pressure at the edges
  delpp = 0.0_8
  DO k=geom%nlevs+1,2,-1
    pep  (:, :, k-1) = pep  (:, :, k-1) + pep(:, :, k)
    delpp(:, :, k-1) =                    pep(:, :, k)
  END DO

endsubroutine delp2p_ad

! ------------------------------------------------------------------------------

subroutine delp2pe_tl(geom,delpp)

  !Tangent linear of conversion from delp to p at layer edges 
  !NOTE: only returns p from layer 2 to nlevs + 1. At top of domain pertubtaion is zero

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Linearized pressure at the edges
  pep(:,:,1) = 0.0
  do k = 2,geom%nlevs+1
    pep(:,:,k) = pep(:,:,k-1) + delpp(:,:,k-1)
  enddo

  !Linearized overwrite with pressure at mid points
  delpp = pep(:,:,2:geom%nlevs+1)

endsubroutine delp2pe_tl

! ------------------------------------------------------------------------------

subroutine delp2pe_ad(geom,delpp)

  !Adjoint of conversion from delp to p at layer edges 
  !NOTE: only returns p from layer 2 to nlevs + 1. At top of domain pertubtaion is zero

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in)  :: geom

  !Perturbations  
  real(kind=kind_real), intent(inout)  :: delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs  )  !Pressure thickness

  !Locals
  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  integer :: k

  !Adjoint overwrite with pressure at mid points
  pep = 0.0_8
  pep(:, :, 2:geom%nlevs+1) = delpp

  !Adjoint pressure at the edges
  delpp = 0.0_8
  DO k=geom%nlevs+1,2,-1
    pep  (:, :, k-1) = pep  (:, :, k-1) + pep(:, :, k)
    delpp(:, :, k-1) =                    pep(:, :, k)
  END DO

endsubroutine delp2pe_ad

! ------------------------------------------------------------------------------

subroutine pt2tv(geom, p0, pt, q, pk)

  !Nonlinear calculation of virtual temperature from potential temperature.

  implicit none

  !Convert potential temperature to virtual temperature
  !Tv = (1 + eps q) * T
  !   = (1 + eps q) * (P/P0)**kappa * pt

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in ) :: geom

  real(kind=kind_real), intent(in   ) :: p0
  real(kind=kind_real), intent(inout) :: pt(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(in   ) ::  q(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(in   ) :: pk(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)

  pt = (1.0 + epsilon * q) * p0**(-kappa) * pk * pt

end subroutine pt2tv

! ------------------------------------------------------------------------------

subroutine pt2tv_tj(geom, p0, pt, q, pk, pt2tvco)

  !Trajectory componenents for tl and ad of potential temperature to virtual temperature

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in ) :: geom

  !Trajecotory
  real(kind=kind_real), intent(in)  :: p0
  real(kind=kind_real), intent(in)  :: pt(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(in)  ::  q(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(in)  :: pk(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(out) :: pt2tvco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,3)

  real(kind=kind_real) :: p0mk

  p0mk = p0**(-kappa)
  pt2tvco(:,:,:,1) = p0mk*(1.0 + epsilon*q)*pk
  pt2tvco(:,:,:,2) = p0mk*epsilon*pk*pt
  pt2tvco(:,:,:,3) = p0mk*(1.0 + epsilon*q)*pt

end subroutine pt2tv_tj

! ------------------------------------------------------------------------------

subroutine pt2tv_tl(geom, pe, pke, pkco, pt2tvco, delpp, ptp, qp)

  !Tangent linear of potential temperature to virtual temperature

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in ) :: geom

  !Trajecotory
  real(kind=kind_real), intent(in) ::      pe(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real), intent(in) ::     pke(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure to the kappa (edges)
  real(kind=kind_real), intent(in) ::    pkco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,2)
  real(kind=kind_real), intent(in) :: pt2tvco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,3)
  !Perturbations
  real(kind=kind_real), intent(in   ) ::  delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(inout) ::    ptp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(in   ) ::     qp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)

  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) :: lpep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) :: pkep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) ::  pkp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real) ::  tvp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)

  integer :: k

  !Linearized pressure at the edges
  pep(:,:,1) = 0.0
  do k = 2,geom%nlevs+1
    pep(:,:,k) = pep(:,:,k-1) + delpp(:,:,k-1)
  enddo

  !Linearized log pressure (edges)
  lpep = pep / pe

  !Linearized p to the kappa at edges
  pkep = kappa * (pke/pe) * pep 

  !Linearized p to the kappa at midpoint
  pkp = pkco(:,:,:,1) * (pkep(:,:,2:geom%nlevs+1) - pkep(:,:,1:geom%nlevs)) - pkco(:,:,:,2) * (lpep(:,:,2:geom%nlevs+1) - lpep(:,:,1:geom%nlevs))

  !Linearized potential temperature to virtual temperature
  tvp = pt2tvco(:,:,:,1) * ptp + pt2tvco(:,:,:,2) * qp + pt2tvco(:,:,:,3) * pkp

  !Linearized overwrite
  ptp = tvp

end subroutine pt2tv_tl

! ------------------------------------------------------------------------------

subroutine pt2tv_ad(geom, pe, pke, pkco, pt2tvco, delpp, ptp, qp)

  !Adjoint of potential temperature to virtual temperature

  implicit none

  !Geometry for the model
  type(fv3jedi_geom)  , intent(in ) :: geom

  !Trajecotory
  real(kind=kind_real), intent(in) ::      pe(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real), intent(in) ::     pke(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure to the kappa (edges)
  real(kind=kind_real), intent(in) ::    pkco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,2)
  real(kind=kind_real), intent(in) :: pt2tvco(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs,3)
  !Perturbations
  real(kind=kind_real), intent(inout) ::  delpp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(inout) ::    ptp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real), intent(inout) ::     qp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)

  real(kind=kind_real) ::  pep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) :: lpep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) :: pkep(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs+1)  !Pressure (edges)
  real(kind=kind_real) ::  pkp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)
  real(kind=kind_real) ::  tvp(geom%bd%isc:geom%bd%iec,geom%bd%jsc:geom%bd%jec,1:geom%nlevs)

  integer :: k

  !Adjoint of overwrite
  tvp = ptp

  !Adjoint of potential temperature to virtual temperature
  ptp =      pt2tvco(:,:,:,1) * tvp
   qp = qp + pt2tvco(:,:,:,2) * tvp
  pkp =      pt2tvco(:,:,:,3) * tvp

  !Adjoint of p to the kappa at midpoint
  pkep = 0.0
  lpep = 0.0
  pkep(:,:,2:geom%nlevs+1) = pkep(:,:,2:geom%nlevs+1) + pkco(:,:,:,1) * pkp
  pkep(:,:,1:geom%nlevs  ) = pkep(:,:,1:geom%nlevs  ) - pkco(:,:,:,1) * pkp
  lpep(:,:,2:geom%nlevs+1) = lpep(:,:,2:geom%nlevs+1) - pkco(:,:,:,2) * pkp
  lpep(:,:,1:geom%nlevs  ) = lpep(:,:,1:geom%nlevs  ) + pkco(:,:,:,2) * pkp

  !Adjoint of p to the kappa at edges + adjoint of linearized log pressure
  pep = kappa * (pke/pe) * pkep + lpep/pe

  !Adjoint of pressure at edges
  do k = geom%nlevs+1,2,-1
       pep(:,:,k-1) =   pep(:,:,k-1) + pep(:,:,k)
     delpp(:,:,k-1) = delpp(:,:,k-1) + pep(:,:,k)
  enddo

end subroutine pt2tv_ad

! ------------------------------------------------------------------------------

! Convert d-grid winds to a-grid.
! -------------------------------

 subroutine d2a2c_vect(u, v, ua, va, dord4, gridstruct, &
                       bd, npx, npy, nested, grid_type)

  use fv3jedi_mod, only: fv_grid_bounds_type, fv_grid_type

!There is a limit to how far this routine can fill uc and vc in the
! halo, and so either mpp_update_domains or some sort of boundary
!  routine (extrapolation, outflow, interpolation from a nested grid)
!   is needed after c_sw is completed if these variables are needed
!    in the halo

 !Modified from sw_core.F90
 !dont need the c grid winds in the output

 !This routine takes in d-grid winds and requires that the halos are full
 !subroutine d2a2c_vect(u, v, ua, va, uc, vc, ut, vt, dord4, gridstruct, &
 !                      bd, npx, npy, nested, grid_type)

  type(fv_grid_bounds_type), intent(IN) :: bd
  logical, intent(in):: dord4
  real(kind=kind_real), intent(in) ::  u(bd%isd:bd%ied,bd%jsd:bd%jed+1)
  real(kind=kind_real), intent(in) ::  v(bd%isd:bd%ied+1,bd%jsd:bd%jed)
!dh  real(kind=kind_real), intent(out), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ):: uc
!dh  real(kind=kind_real), intent(out), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1):: vc
  real(kind=kind_real), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ):: uc
  real(kind=kind_real), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1):: vc
!dh  real(kind=kind_real), intent(out), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed  ):: ua, va, ut, vt
  real(kind=kind_real), intent(out), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed  ):: ua, va
  real(kind=kind_real), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed  ):: ut, vt
  integer, intent(IN) :: npx, npy, grid_type
  logical, intent(IN) :: nested
  type(fv_grid_type), intent(IN), target :: gridstruct
! Local 
  real(kind=kind_real), dimension(bd%isd:bd%ied,bd%jsd:bd%jed):: utmp, vtmp
  integer npt, i, j, ifirst, ilast, id
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed

  real(kind=kind_real), pointer, dimension(:,:,:) :: sin_sg
  real(kind=kind_real), pointer, dimension(:,:)   :: cosa_u, cosa_v, cosa_s
  real(kind=kind_real), pointer, dimension(:,:)   :: rsin_u, rsin_v, rsin2
  real(kind=kind_real), pointer, dimension(:,:)   :: dxa,dya

  real(kind=kind_real), dimension(4) :: etmp1, etmp2

!----------------------------
! 4-pt Lagrange interpolation
!----------------------------
  real(kind=kind_real), parameter:: a1 =  0.5625
  real(kind=kind_real), parameter:: a2 = -0.0625
!----------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
  real(kind=kind_real), parameter:: c1 = -2./14.
  real(kind=kind_real), parameter:: c2 = 11./14.
  real(kind=kind_real), parameter:: c3 =  5./14.


      is  = bd%isc
      ie  = bd%iec
      js  = bd%jsc
      je  = bd%jec
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      sin_sg    => gridstruct%sin_sg  
      cosa_u    => gridstruct%cosa_u  
      cosa_v    => gridstruct%cosa_v  
      cosa_s    => gridstruct%cosa_s  
      rsin_u    => gridstruct%rsin_u  
      rsin_v    => gridstruct%rsin_v  
      rsin2     => gridstruct%rsin2   
      dxa       => gridstruct%dxa     
      dya       => gridstruct%dya     

  if ( dord4 ) then
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
  utmp(:,:) = 0.0
  vtmp(:,:) = 0.0 

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
     if( gridstruct%sw_corner ) then
         do i=-2,0
            utmp(i,0) = -vtmp(0,1-i)
         enddo
     endif
     if( gridstruct%se_corner ) then
         do i=0,2
            utmp(npx+i,0) = vtmp(npx,i+1)
         enddo
     endif
     if( gridstruct%ne_corner ) then
         do i=0,2
            utmp(npx+i,npy) = -vtmp(npx,je-i)
         enddo
     endif
     if( gridstruct%nw_corner ) then
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
     if( gridstruct%sw_corner ) then
         ua(-1,0) = -va(0,2)
         ua( 0,0) = -va(0,1) 
     endif
     if( gridstruct%se_corner ) then
         ua(npx,  0) = va(npx,1)
         ua(npx+1,0) = va(npx,2) 
     endif
     if( gridstruct%ne_corner ) then
         ua(npx,  npy) = -va(npx,npy-1)
         ua(npx+1,npy) = -va(npx,npy-2) 
     endif
     if( gridstruct%nw_corner ) then
         ua(-1,npy) = va(0,npy-2)
         ua( 0,npy) = va(0,npy-1) 
     endif

     if( is==1 .and. .not. nested  ) then
        do j=js-1,je+1
           uc(0,j) = c1*utmp(-2,j) + c2*utmp(-1,j) + c3*utmp(0,j)
           etmp1 = ua(-1:2,j)
           etmp2 = dxa(-1:2,j)
           ut(1,j) = edge_interpolate4(etmp1,etmp2)
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
           etmp1 = ua(npx-2:npx+1,j)
           etmp2 = dxa(npx-2:npx+1,j)
           ut(npx,  j) = edge_interpolate4(etmp1,etmp2)
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
     if( gridstruct%sw_corner ) then
         do j=-2,0
            vtmp(0,j) = -utmp(1-j,0)
         enddo
     endif
     if( gridstruct%nw_corner ) then
         do j=0,2
            vtmp(0,npy+j) = utmp(j+1,npy)
         enddo
     endif
     if( gridstruct%se_corner ) then
         do j=-2,0
            vtmp(npx,j) = utmp(ie+j,0)
         enddo
     endif
     if( gridstruct%ne_corner ) then
         do j=0,2
            vtmp(npx,npy+j) = -utmp(ie-j,npy)
         enddo
     endif
     if( gridstruct%sw_corner ) then
         va(0,-1) = -ua(2,0)
         va(0, 0) = -ua(1,0)
     endif
     if( gridstruct%se_corner ) then
         va(npx, 0) = ua(npx-1,0)
         va(npx,-1) = ua(npx-2,0)
     endif
     if( gridstruct%ne_corner ) then
         va(npx,npy  ) = -ua(npx-1,npy)
         va(npx,npy+1) = -ua(npx-2,npy)
     endif
     if( gridstruct%nw_corner ) then
         va(0,npy)   = ua(1,npy)
         va(0,npy+1) = ua(2,npy)
     endif

 if (grid_type < 3) then

     do j=js-1,je+2
      if ( j==1 .and. .not. nested  ) then
        do i=is-1,ie+1
           etmp1 = va(i,-1:2)
           etmp2 = dya(i,-1:2)
           vt(i,j) = edge_interpolate4(etmp1,etmp2)
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
           etmp1 = va(i,j-2:j+1)
           etmp2 = dya(i,j-2:j+1)
           vt(i,j) = edge_interpolate4(etmp1,etmp2)
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

 end subroutine d2a2c_vect

 real function edge_interpolate4(ua, dxa)

   implicit none

   real(kind=kind_real), intent(in) :: ua(4)
   real(kind=kind_real), intent(in) :: dxa(4)
   real(kind=kind_real) :: t1, t2

   t1 = dxa(1) + dxa(2)
   t2 = dxa(3) + dxa(4)
   edge_interpolate4 = 0.5*( ((t1+dxa(2))*ua(2)-dxa(2)*ua(1)) / t1 + &
                             ((t2+dxa(3))*ua(3)-dxa(3)*ua(4)) / t2 )

 end function edge_interpolate4

end module variable_transforms
