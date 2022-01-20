! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_vc_coldstartwinds_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! fv3
use fv_arrays_mod,        only: R_GRID
use fv_grid_utils_mod,    only: mid_pt_sphere, get_unit_vect2, get_latlon_vector, inner_prod

! fv3jedi
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_field_mod,     only: copy_subset, field_clen
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state

implicit none
private
public :: fv3jedi_vc_coldstartwinds

type :: fv3jedi_vc_coldstartwinds
 integer :: isc, iec, jsc, jec
 real(kind=kind_real), allocatable, dimension(:,:,:) :: grid
 contains
   procedure :: create
   procedure :: delete
   procedure :: changevar
end type fv3jedi_vc_coldstartwinds

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_vc_coldstartwinds), intent(inout) :: self
type(fv3jedi_geom),               intent(in)    :: geom
type(fckit_configuration),        intent(in)    :: conf

allocate(self%grid(geom%isd:geom%ied+1,geom%jsd:geom%jed+1,2))
self%grid = geom%grid

self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_vc_coldstartwinds), intent(inout) :: self

deallocate(self%grid)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self, xin, xout)

class(fv3jedi_vc_coldstartwinds), intent(inout) :: self
type(fv3jedi_state),              intent(in)    :: xin
type(fv3jedi_state),              intent(inout) :: xout

integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Winds
logical :: have_d_winds
integer :: i, j, k, levp
real(kind=kind_real), allocatable :: ud_cold (:,:,:)         !u component D-grid
real(kind=kind_real), allocatable :: vd_cold (:,:,:)         !v component D-grid
real(kind=kind_real), pointer     :: u_w_cold(:,:,:)
real(kind=kind_real), pointer     :: v_w_cold(:,:,:)
real(kind=kind_real), pointer     :: u_s_cold(:,:,:)
real(kind=kind_real), pointer     :: v_s_cold(:,:,:)

real(kind=R_GRID), dimension(2) :: p1, p2, p3
real(kind=R_GRID), dimension(3) :: e1, e2, ex, ey

! Identity part of the change of variables
! ----------------------------------------
call copy_subset(xin%fields, xout%fields, fields_to_do)


! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return


! D-Grid winds
! ------------
have_d_winds = .false.
if ( xin%has_field('u_w_cold') .and. xin%has_field('v_w_cold') .and. &
     xin%has_field('u_s_cold') .and. xin%has_field('v_s_cold') ) then

  call xin%get_field('u_w_cold', u_w_cold)
  call xin%get_field('v_w_cold', v_w_cold)
  call xin%get_field('u_s_cold', u_s_cold)
  call xin%get_field('v_s_cold', v_s_cold)

  levp = size(u_w_cold,3)

  allocate(ud_cold(self%isc:self%iec,  self%jsc:self%jec+1,1:levp))
  allocate(vd_cold(self%isc:self%iec+1,self%jsc:self%jec,  1:levp))

  do k = 1, levp
    do j = self%jsc, self%jec+1
      do i = self%isc, self%iec
        p1(:) = self%grid(i,  j,1:2)
        p2(:) = self%grid(i+1,j,1:2)
        call  mid_pt_sphere(p1, p2, p3)
        call get_unit_vect2(p1, p2, e1)
        call get_latlon_vector(p3, ex, ey)
        ud_cold(i,j,k) = u_s_cold(i,j,k)*inner_prod(e1, ex) + v_s_cold(i,j,k)*inner_prod(e1, ey)
      enddo
    enddo
    do j = self%jsc, self%jec
      do i = self%isc, self%iec+1
        p1(:) = self%grid(i,j  ,1:2)
        p2(:) = self%grid(i,j+1,1:2)
        call  mid_pt_sphere(p1, p2, p3)
        call get_unit_vect2(p1, p2, e2)
        call get_latlon_vector(p3, ex, ey)
        vd_cold(i,j,k) = u_w_cold(i,j,k)*inner_prod(e2, ex) + v_w_cold(i,j,k)*inner_prod(e2, ey)
      enddo
    enddo
  enddo

  have_d_winds = .true.

endif


! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call xout%get_field(trim(fields_to_do(f)), field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ud_cold")

    if (.not. have_d_winds) call field_fail(fields_to_do(f))
    field_ptr = ud_cold

  case ("vd_cold")

    if (.not. have_d_winds) call field_fail(fields_to_do(f))
    field_ptr = vd_cold

  case default

    call abor1_ftn("fv3jedi_vc_coldstartwinds_mod.changevar unknown field: "//trim(fields_to_do(f))&
                   //". Not in input field and no transform case specified.")

  end select

enddo

end subroutine changevar

! --------------------------------------------------------------------------------------------------

end module fv3jedi_vc_coldstartwinds_mod
