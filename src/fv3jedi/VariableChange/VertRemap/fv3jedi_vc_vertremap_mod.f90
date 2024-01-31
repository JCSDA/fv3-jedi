! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_vc_vertremap_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! fms
use constants_mod,           only: grav
use field_manager_mod,       only: MODEL_ATMOS
use mpp_domains_mod,         only: mpp_update_domains
use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, get_tracer_index, NO_TRACER, &
                                   set_tracer_profile

! fv3
use external_ic_mod,         only: remap_scalar, remap_dwinds
use fv_arrays_mod,           only: fv_atmos_type, deallocate_fv_atmos_type, R_GRID
use fv_grid_utils_mod,       only: mid_pt_sphere, get_unit_vect2, get_latlon_vector, inner_prod
use test_cases_mod,          only: checker_tracers

! fv3jedi
use fv_prec_mod,             only: kind_fv3
use fv_init_mod,             only: fv_init
use fv3jedi_fmsnamelist_mod, only: fv3jedi_fmsnamelist
use fv3jedi_geom_mod,        only: fv3jedi_geom
use fv3jedi_fieldfail_mod,   only: field_fail
use fv3jedi_field_mod,       only: copy_subset, field_clen, fv3jedi_field
use fv3jedi_kinds_mod,       only: kind_real
use fv3jedi_state_mod,       only: fv3jedi_state

implicit none
private
public :: fv3jedi_vc_vertremap

type :: fv3jedi_vc_vertremap
 type(fv_atmos_type), allocatable :: Atm(:)
 logical, allocatable :: grids_on_this_pe(:)
 logical :: from_cold_start, checker_tr
 integer :: nt_checker
 contains
   procedure :: create
   procedure :: delete
   procedure :: changevar
end type fv3jedi_vc_vertremap

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_vc_vertremap), intent(inout) :: self
type(fv3jedi_geom),          intent(in)    :: geom
type(fckit_configuration),   intent(in)    :: conf

integer :: gtile, p_split = 1, n
logical :: checks_passed
character(len=:), allocatable :: str
type(fv3jedi_fmsnamelist) :: fmsnamelist

! Prepare namelist
call fmsnamelist%replace_namelist(conf)

! Create Atm structure
call fv_init(self%Atm, 300.0_kind_real, self%grids_on_this_pe, p_split, gtile, .true.)

! Flag to use cold starts
if( .not. conf%get('input is cold starts', self%from_cold_start) ) self%from_cold_start = .true.

! Check tracers
if( .not. conf%get('check tracers',    self%checker_tr) ) self%checker_tr = .false.
if( .not. conf%get('check tracers nt', self%nt_checker) ) self%nt_checker = 0

! Set global source variable in external_ic
if (.not. conf%get("source of inputs", str)) then
  str = 'FV3GFS GAUSSIAN NETCDF FILE'
endif

! Remapping needs nggps_ic to be true
self%Atm(1)%flagstruct%nggps_ic = .true.

! ak/bk
self%Atm(1)%ak = real(geom%ak,kind_fv3)
self%Atm(1)%bk = real(geom%bk,kind_fv3)
self%Atm%ptop = real(geom%ak(1),kind_fv3)

! Sanity checks
checks_passed = .true.
if (checks_passed) checks_passed = geom%npx == self%Atm(1)%npx
if (checks_passed) checks_passed = geom%npy == self%Atm(1)%npy
if (checks_passed) checks_passed = geom%npz == self%Atm(1)%npz
if (checks_passed) checks_passed = geom%isd == self%Atm(1)%bd%isd
if (checks_passed) checks_passed = geom%ied == self%Atm(1)%bd%ied
if (checks_passed) checks_passed = geom%jsd == self%Atm(1)%bd%jsd
if (checks_passed) checks_passed = geom%jed == self%Atm(1)%bd%jed
if (.not.checks_passed) call abor1_ftn("fv3jedi_vc_vertremap_mod.field_fail: Geometry generated"// &
                                       " here does not match fv3-jedi geometry.")

if (.not.size(self%Atm)==1) call abor1_ftn("fv3jedi_vc_vertremap_mod.field_fail: Atm strucutre"// &
                                           " with size > 1 not supported.")

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_vc_vertremap), intent(inout) :: self

integer :: n

do n = 1,size(self%Atm)
  call deallocate_fv_atmos_type(self%Atm(n))
enddo
deallocate(self%Atm)
deallocate(self%grids_on_this_pe)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self, xin, xout)

class(fv3jedi_vc_vertremap), target, intent(inout) :: self
type(fv3jedi_state),                 intent(in)    :: xin
type(fv3jedi_state),                 intent(inout) :: xout

! Copy versus transform control
integer :: f
character(len=field_clen), allocatable :: fields_to_do(:)
type(fv3jedi_field), pointer :: field_ptr

! Local versions of fv3jedi state
real(kind=kind_real), allocatable :: orog_filt(:,:,:)
real(kind=kind_real), allocatable :: ps_cold(:,:,:)
real(kind=kind_real), allocatable :: zh_cold(:,:,:)
real(kind=kind_real), allocatable ::  w_cold(:,:,:)
real(kind=kind_real), allocatable ::  t_cold(:,:,:)
real(kind=kind_real), allocatable :: q_tmp(:,:,:)
real(kind=kind_real), allocatable :: ud_cold(:,:,:)
real(kind=kind_real), allocatable :: vd_cold(:,:,:)

logical :: have_remapped
type(fv_atmos_type), pointer :: Atm
integer:: i, j, k, nt, ntracers, ntprog, itoa, levp, isc, iec, jsc, jec, npz, nts
integer:: liq_wat, ice_wat, rainwat, snowwat, graupel, ntclamt
character(len=64) :: tracer_name
real(kind=kind_fv3), allocatable :: ak(:), bk(:)
real(kind=kind_fv3), allocatable ::  q_cold(:,:,:,:)
real(kind=kind_fv3) :: wt, qt, m_fac

! Field names
character(len=field_clen) :: ps_fname
character(len=field_clen) :: zh_fname
character(len=field_clen) :: w_fname
character(len=field_clen) :: t_fname
character(len=field_clen) :: ud_fname
character(len=field_clen) :: vd_fname


! Identity part of the change of variables
! ----------------------------------------
call copy_subset(xin%fields, xout%fields, fields_to_do)


! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return


! Field names
! -----------
ps_fname = 'ps'
zh_fname = 'zh'
w_fname = 'w'
t_fname = 't'
ud_fname = 'ud'
vd_fname = 'vd'
if (self%from_cold_start) then
  ps_fname = trim(ps_fname)//'_cold'
  zh_fname = trim(zh_fname)//'_cold'
  w_fname = trim(w_fname)//'_cold'
  t_fname = trim(t_fname)//'_cold'
  ud_fname = trim(ud_fname)//'_cold'
  vd_fname = trim(vd_fname)//'_cold'
endif


! Remap to new Lagrangian vertical coordinate
! -------------------------------------------
have_remapped = .false.
if ( xin%has_field(ps_fname) .and. xin%has_field(zh_fname) .and. &
     xin%has_field(w_fname) .and. &
     xin%has_field('orog_filt') .and. &
     xin%has_field(ud_fname) .and. xin%has_field(vd_fname) ) then

  ! Shortcuts
  Atm => self%Atm(1)
  isc = Atm%bd%is
  iec = Atm%bd%ie
  jsc = Atm%bd%js
  jec = Atm%bd%je
  npz = Atm%npz

  ! Orography
  call xin%get_field('orog_filt', orog_filt)
  Atm%phis(isc:iec,jsc:jec) = real(orog_filt(isc:iec,jsc:jec,1),kind_fv3)*grav  ! Convert to phis

  ! Dynamics fields
  call xin%get_field(ps_fname, ps_cold)
  call xin%get_field(zh_fname, zh_cold)
  call xin%get_field( w_fname,  w_cold)

  ! akbk from cold starts have different levels (for now extra levels are zero)
  levp = size(w_cold,3)
  itoa = levp - npz + 1

  allocate(ak(levp+1))
  allocate(bk(levp+1))
  ak = 0.0_kind_fv3
  bk = 0.0_kind_fv3
  ak(itoa:levp+1) = Atm%ak(1:npz+1)
  bk(itoa:levp+1) = Atm%bk(1:npz+1)
  ak(1) = max(1.e-9_kind_fv3, ak(1))

  ! Tracers
  ! -------
  ! Number of tracers in Atm
  call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)
  call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers, num_prog=ntprog)

  ! initialize all tracers to default values prior to being input
  do nt = 1, ntprog
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    ! set all tracers to an initial profile value
    call set_tracer_profile (MODEL_ATMOS, nt, Atm%q(:,:,:,nt)  )
  enddo
  do nt = ntprog+1, ntracers
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    ! set all tracers to an initial profile value
    call set_tracer_profile (MODEL_ATMOS, nt, Atm%qdiag(:,:,:,nt)  )
  enddo

  ! Copy from fv3-jedi to temporary cold q structure
  allocate (q_cold(isc:iec, jsc:jec, levp, ntracers))
  q_cold = 0.0_kind_fv3
  do nt = 1, ntracers
    call get_tracer_names(MODEL_ATMOS, nt, tracer_name)
    if (self%from_cold_start) tracer_name = trim(tracer_name)//"_cold"
    if (xin%has_field(trim(tracer_name))) then
      call xin%get_field(trim(tracer_name), q_tmp)
      q_cold(:,:,:,nt) = real(q_tmp, kind_fv3)
    endif
  enddo

  ! Call remapping non-wind variables
  ! ---------------------------------
  if (xin%has_field(t_fname)) then
    call xin%get_field( t_fname,  t_cold)
    call remap_scalar(Atm, levp, npz, ntracers, ak, bk, real(ps_cold(:,:,1),kind_fv3), q_cold, &
                      real(zh_cold,kind_fv3), real(w_cold,kind_fv3), real(t_cold,kind_fv3))
  else
    call remap_scalar(Atm, levp, npz, ntracers, ak, bk, real(ps_cold(:,:,1),kind_fv3), q_cold, &
                      real(zh_cold,kind_fv3), real(w_cold,kind_fv3))
  endif

  ! Call remapping wind variables
  ! -----------------------------
  call xin%get_field(ud_fname, ud_cold)
  call xin%get_field(vd_fname, vd_cold)

  call remap_dwinds(levp, npz, ak, bk, real(ps_cold(:,:,1),kind_fv3), real(ud_cold,kind_fv3), &
                    real(vd_cold,kind_fv3), Atm)

  ! Tracer weighting
  ! ----------------
  liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
  ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
  rainwat = get_tracer_index(MODEL_ATMOS, 'rainwat')
  snowwat = get_tracer_index(MODEL_ATMOS, 'snowwat')
  graupel = get_tracer_index(MODEL_ATMOS, 'graupel')
  ntclamt = get_tracer_index(MODEL_ATMOS, 'cld_amt')

  if (self%from_cold_start) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          wt = Atm%delp(i,j,k)
          if ( Atm%flagstruct%nwat == 6 ) then
            qt = wt*(1.0_kind_fv3 + Atm%q(i,j,k,liq_wat) + Atm%q(i,j,k,ice_wat) + &
                                    Atm%q(i,j,k,rainwat) + Atm%q(i,j,k,snowwat) + &
                                    Atm%q(i,j,k,graupel))
          else
            qt = wt*(1.0_kind_fv3 + sum(Atm%q(i,j,k,2:Atm%flagstruct%nwat)))
          endif
          Atm%delp(i,j,k) = qt
          if (ntclamt > 0) Atm%q(i,j,k,ntclamt) = 0.0
        enddo
      enddo
    enddo
  else
    ! TODO Question is do we need this if just remapping after adding the increment, say?
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          wt = Atm%delp(i,j,k)
          if ( Atm%flagstruct%nwat == 6 ) then
            qt = wt*(1.0_kind_fv3 + Atm%q(i,j,k,liq_wat) + Atm%q(i,j,k,ice_wat) + &
                                    Atm%q(i,j,k,rainwat) + Atm%q(i,j,k,snowwat) + &
                                    Atm%q(i,j,k,graupel))
          else
             qt = wt*(1.0_kind_fv3 + sum(Atm%q(i,j,k,2:Atm%flagstruct%nwat)))
          endif
          m_fac = wt / qt
          do nt=1,ntracers
            Atm%q(i,j,k,nt) = m_fac * Atm%q(i,j,k,nt)
          enddo
          Atm%delp(i,j,k) = qt
          if (ntclamt > 0) Atm%q(i,j,k,ntclamt) = 0.0
        enddo
      enddo
    enddo
  endif

  if (self%checker_tr) then
    nts = ntracers - self%nt_checker+1
    call checker_tracers(isc, iec, jsc, jec, Atm%bd%isd, Atm%bd%ied, Atm%bd%jsd, Atm%bd%jed, &
                         self%nt_checker, npz, Atm%q(:,:,:,nts:ntracers), &
                         Atm%gridstruct%agrid_64(isc:iec,jsc:jec,1),     &
                         Atm%gridstruct%agrid_64(isc:iec,jsc:jec,2), &
                         real(9.0_kind_real,kind_fv3), real(9.0_kind_real,kind_fv3))
  endif

  have_remapped = .true.

endif


! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  ! Remapping needs to have happened in order to copy anything
  if (.not. have_remapped) call field_fail(fields_to_do(f))

  call xout%get_field(trim(fields_to_do(f)), field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ud")

    field_ptr%array(isc:iec,jsc:jec+1,:) = Atm%u(isc:iec,jsc:jec+1,:)

  case ("vd")

    field_ptr%array(isc:iec+1,jsc:jec,:) = Atm%v(isc:iec+1,jsc:jec,:)

  case ("t")

    field_ptr%array(isc:iec,jsc:jec,:) = Atm%pt(isc:iec,jsc:jec,:)

  case ("delp")

    field_ptr%array(isc:iec,jsc:jec,:) = Atm%delp(isc:iec,jsc:jec,:)

  case ("ps")

    field_ptr%array(isc:iec,jsc:jec,1) = Atm%ps(isc:iec,jsc:jec)

  case ("delz")

    field_ptr%array(isc:iec,jsc:jec,:) = Atm%delz(isc:iec,jsc:jec,:)

  case ("w")

    field_ptr%array(isc:iec,jsc:jec,:) = Atm%w(isc:iec,jsc:jec,:)

  case ("phis")

    field_ptr%array(isc:iec,jsc:jec,1) = Atm%phis(isc:iec,jsc:jec)

  case ('sphum','liq_wat','ice_wat','o3mr','graupel','snowwat','rainwat')

    nt = get_tracer_index(MODEL_ATMOS, trim(fields_to_do(f)))
    if (nt == NO_TRACER) call field_fail(fields_to_do(f))
    field_ptr%array(isc:iec,jsc:jec,:) = Atm%q(isc:iec,jsc:jec,:,nt)

  case ('sgs_tke','cld_amt')

    field_ptr%array = 0.0_kind_real

  case default

    call abor1_ftn("fv3jedi_vc_coldstartwinds_mod.changevar unknown field: "//trim(fields_to_do(f))&
                   //". Not in input field and no transform case specified.")

  end select

enddo

end subroutine changevar

! --------------------------------------------------------------------------------------------------

end module fv3jedi_vc_vertremap_mod
