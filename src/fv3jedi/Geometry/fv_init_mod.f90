! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv_init_mod

! fv3 uses
use fv_arrays_mod,     only: fv_atmos_type
use fv_control_mod,    only: fv_init_fv3 => fv_init, pelist_all

! fv3jedi uses
use fv3jedi_kinds_mod, only: kind_real
use fv_prec_mod,       only: kind_fv3

implicit none
private
public fv_init

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fv_init(Atm, dt_atmos_in, grids_on_this_pe, p_split, gtile, skip_nml_read)

type(fv_atmos_type), allocatable, intent(inout) :: Atm(:)
real(kind=kind_real),             intent(in)    :: dt_atmos_in
logical, allocatable,             intent(inout) :: grids_on_this_pe(:)
integer,                          intent(inout) :: p_split
integer, optional,                intent(out)   :: gtile
logical,                          intent(in)    :: skip_nml_read

real(kind_fv3) :: dt_atmos

dt_atmos = dt_atmos_in

call fv_init_fv3(Atm, dt_atmos, grids_on_this_pe, p_split)
deallocate(pelist_all)

if (present(gtile)) gtile = Atm(1)%tile

end subroutine fv_init

! --------------------------------------------------------------------------------------------------

end module fv_init_mod
