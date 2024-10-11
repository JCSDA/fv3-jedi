! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_lvc_model2geovals_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,           only: fckit_log

use datetime_mod

use fv3jedi_constants_mod, only: constant
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_field_mod,     only: copy_subset, field_clen
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state

use radii_vt_mod
use height_vt_mod
use moisture_vt_mod
use pressure_vt_mod
use surface_vt_mod
use temperature_vt_mod
use wind_vt_mod

implicit none

private
public :: fv3jedi_lvc_model2geovals

type :: fv3jedi_lvc_model2geovals
  integer :: isc, iec, jsc, jec, npz
  real(kind=kind_real), allocatable ::      t(:,:,:)
  real(kind=kind_real), allocatable ::      q(:,:,:)
  real(kind=kind_real), allocatable ::     o3(:,:,:)
  real(kind=kind_real), allocatable ::     ql(:,:,:)
  real(kind=kind_real), allocatable ::     qi(:,:,:)
  real(kind=kind_real), allocatable ::     qr(:,:,:)
  real(kind=kind_real), allocatable ::     qs(:,:,:)
  real(kind=kind_real), allocatable ::     qg(:,:,:)
  real(kind=kind_real), allocatable ::   delp(:,:,:)
  real(kind=kind_real), allocatable ::  slmsk(:,:,:)
  real(kind=kind_real), allocatable :: sheleg(:,:,:)
  real(kind=kind_real), allocatable :: frseaice(:,:,:)
  real(kind=kind_real), allocatable :: frland(:,:,:)
  real(kind=kind_real), allocatable :: frsnow(:,:,:)
  real(kind=kind_real), allocatable ::frocean(:,:,:)
  real(kind=kind_real), allocatable :: frlake(:,:,:)
  real(kind=kind_real), allocatable ::     ts(:,:,:)
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: multiply
    procedure, public :: multiplyadjoint
end type fv3jedi_lvc_model2geovals

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, dummyconf)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(in)    :: geom
type(fv3jedi_state),              intent(in)    :: bg
type(fv3jedi_state),              intent(in)    :: fg
type(fckit_configuration),        intent(in)    :: dummyconf

integer :: i,j
logical :: have_fractions,have_slmsk,have_ts
real(kind=kind_real), allocatable :: local_swe(:,:,:)

!Locals
real(kind=kind_real), parameter :: minswe = 1.0_kind_real / 10.0_kind_real

!!! DO NOT USE CONF !!!

! Grid convenience
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz

! Trajectory fields
if (bg%has_field('t'     )) call bg%get_field('t'     , self%t )
if (bg%has_field('sphum' )) call bg%get_field('sphum' , self%q )
if (bg%has_field('o3mr'  )) call bg%get_field('o3mr'  , self%o3)
if (bg%has_field('o3ppmv')) call bg%get_field('o3ppmv', self%o3)
if (bg%has_field('delp'  )) call bg%get_field('delp'  , self%delp)
if (bg%has_field('cloud_liquid_water')) call bg%get_field('cloud_liquid_water', self%ql)
if (bg%has_field('cloud_liquid_ice'))   call bg%get_field('cloud_liquid_ice'  , self%qi)
if (bg%has_field('rain_water'))         call bg%get_field('rain_water'        , self%qr)
if (bg%has_field('snow_water'))         call bg%get_field('snow_water'        , self%qs)
if (bg%has_field('graupel'))            call bg%get_field('graupel'           , self%qg)
if (bg%has_field('fraction_of_ice'))    call bg%get_field('fraction_of_ice'   , self%frseaice)
if (bg%has_field('fraction_of_snow'))   call bg%get_field('fraction_of_snow'  , self%frsnow)
if (bg%has_field('fraction_of_land'))   call bg%get_field('fraction_of_land'  , self%frland)
if (bg%has_field('fraction_of_lake'))   call bg%get_field('fraction_of_lake'  , self%frlake)
if (bg%has_field('fraction_of_ocean'))  call bg%get_field('fraction_of_ocean' , self%frocean)
if (bg%has_field('slmsk'))              call bg%get_field('slmsk'             , self%slmsk)
if (bg%has_field('sheleg'))             call bg%get_field('sheleg'            , self%sheleg)

have_ts=.false.
if (bg%has_field('ts')) then
   call bg%get_field('ts', self%ts)
   have_ts=.true.
else if (bg%has_field('tsea')) then
   call bg%get_field('tsea'              , self%ts)
   have_ts=.true.
endif

have_fractions=allocated(self%frlake) .and. allocated(self%frocean) .and. &
               allocated(self%frseaice)

! Land sea mask
! -------------
have_slmsk = .false.
if (allocated(self%slmsk)) then
  have_slmsk = .true.
elseif ( have_fractions .and. have_ts ) then

  allocate(self%slmsk(self%isc:self%iec,self%jsc:self%jec,1))
  self%slmsk = 1.0_kind_real !Land
  do j = self%jsc,self%jec
    do i = self%isc,self%iec
      if ( self%frocean(i,j,1) + self%frlake(i,j,1) >= 0.6_kind_real) then
        self%slmsk(i,j,1) = 0.0_kind_real ! Water
      endif
      if ( self%slmsk(i,j,1) == 0.0_kind_real .and. self%frseaice(i,j,1) > 0.5_kind_real) then
        self%slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
      if ( self%slmsk(i,j,1) == 0.0_kind_real .and. self%ts(i,j,1) < 271.4_kind_real ) then
        self%slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
    enddo
  enddo
  have_slmsk = .true.
endif

! Land sea mask
! -------------
if (have_slmsk) then

  if(.not.allocated(self%frocean)) allocate(self%frocean(self%isc:self%iec,self%jsc:self%jec,1))
  if(.not.allocated(self%frland)) allocate(self%frland(self%isc:self%iec,self%jsc:self%jec,1))
  if(.not.allocated(self%frseaice)) allocate(self%frseaice(self%isc:self%iec,self%jsc:self%jec,1))
  if(.not.allocated(self%frsnow)) allocate(self%frsnow(self%isc:self%iec,self%jsc:self%jec,1))
  allocate(local_swe(geom%isc:geom%iec,geom%jsc:geom%jec,1))

  ! Potential for missing values in snow water equivalent (if missing set to 0.0)
  local_swe = self%sheleg  ! SWE is named "sheleg" in backgrounds
  where (abs(local_swe) > 10.0e10_kind_real) local_swe = 0.0_kind_real

  ! Note: The GFS slmsk has values {0,1,2} denoting {sea,land,ice}.
  !       Locally within this function, we also use an additional value (3) to denote snow.
  self%slmsk = nint(self%slmsk)
  where (self%slmsk >= 1 .and. local_swe > minswe) self%slmsk = 3

  do j = self%jsc,self%jec
    do i = self%isc,self%iec
      if ( self%slmsk(i,j,1) == 0.0_kind_real) then
        self%frocean(i,j,1) = 1.0_kind_real
      elseif ( self%slmsk(i,j,1) == 1.0_kind_real) then
        self%frland(i,j,1) = 1.0_kind_real
      elseif ( self%slmsk(i,j,1) == 2.0_kind_real) then
        self%frseaice(i,j,1) = 1.0_kind_real
      elseif ( self%slmsk(i,j,1) == 3.0_kind_real) then
        self%frsnow(i,j,1) = 1.0_kind_real
      endif
    enddo
  enddo

  deallocate(local_swe)

endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self

if (allocated(self%slmsk)) deallocate(self%slmsk)
if (allocated(self%sheleg)) deallocate(self%sheleg)
if (allocated(self%frocean)) deallocate(self%frocean  )
if (allocated(self%frlake)) deallocate(self%frlake  )
if (allocated(self%frland)) deallocate(self%frland  )
if (allocated(self%frsnow)) deallocate(self%frsnow  )
if (allocated(self%frseaice)) deallocate(self%frseaice  )
if (allocated(self%qg  )) deallocate(self%qg  )
if (allocated(self%qs  )) deallocate(self%qs  )
if (allocated(self%qr  )) deallocate(self%qr  )
if (allocated(self%qi  )) deallocate(self%qi  )
if (allocated(self%ql  )) deallocate(self%ql  )
if (allocated(self%delp)) deallocate(self%delp)
if (allocated(self%o3  )) deallocate(self%o3  )
if (allocated(self%q   )) deallocate(self%q   )
if (allocated(self%t   )) deallocate(self%t   )

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, geom, dxm, dxg)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(inout) :: geom
type(fv3jedi_increment),          intent(in)    :: dxm
type(fv3jedi_increment),          intent(inout) :: dxg

integer :: f, i, j, k, nf2do
character(len=field_clen), allocatable :: fields_to_do_(:)
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

!Winds
logical :: have_winds
real(kind=kind_real), allocatable :: ua  (:,:,:)         !A-grid wind u component
real(kind=kind_real), allocatable :: va  (:,:,:)         !A-grid wind v component
real(kind=kind_real), pointer     :: ud  (:,:,:)         !D-grid wind u component
real(kind=kind_real), pointer     :: vd  (:,:,:)         !D-grid wind v component

!Virtual temperature
logical :: have_tv
real(kind=kind_real), pointer     :: q   (:,:,:)         !Specific humidity
real(kind=kind_real), pointer     :: t   (:,:,:)         !Temperature
real(kind=kind_real), allocatable :: tv  (:,:,:)         !Virtual temperature

!Humidity mixing ratio
logical :: have_qmr
real(kind=kind_real), allocatable :: qmr (:,:,:)         !Humidity mixing ratio

!Cloud liquid water mixing ratio
logical :: have_ql,have_qi,have_qr,have_qs,have_qg
real(kind=kind_real), allocatable :: clwpath (:,:,:)     !Cloud liquid  water path
real(kind=kind_real), allocatable :: ciwpath (:,:,:)     !Cloud ice     water path
real(kind=kind_real), allocatable :: crwpath (:,:,:)     !Cloud rain    water path
real(kind=kind_real), allocatable :: cswpath (:,:,:)     !Cloud snow    water path
real(kind=kind_real), allocatable :: cgwpath (:,:,:)     !Cloud graupel water path

real(kind=kind_real), allocatable :: cmxr (:,:,:)        !Cloud mixing ratio

!Ozone mixing ratio
logical :: have_o3
real(kind=kind_real), allocatable :: o3mr  (:,:,:)       !Ozone mixing ratio
real(kind=kind_real), allocatable :: o3ppmv(:,:,:)       !Ozone ppmv

!Surface pressure
logical :: have_ps
real(kind=kind_real), pointer     :: ps  (:,:,:)         !Surface pressure
real(kind=kind_real), pointer     :: delp(:,:,:)         !Pressure thickness

!Skin temperature
logical :: have_tskin
real(kind=kind_real), pointer     :: tskin(:,:,:)        !Skin temperature

! Identity part of the change of fields
! -------------------------------------
call copy_subset(dxm%fields, dxg%fields, fields_to_do_)


! Winds are always special case
! -----------------------------
nf2do = 0
if (allocated(fields_to_do_)) nf2do = size(fields_to_do_)

if (dxg%has_field('ua')) then
  allocate(fields_to_do(nf2do+2))
  if (allocated(fields_to_do_)) fields_to_do(1:nf2do) = fields_to_do_
  fields_to_do(nf2do+1) = 'ua'
  fields_to_do(nf2do+2) = 'va'
else
  if (allocated(fields_to_do_)) then
    allocate(fields_to_do(nf2do))
    fields_to_do(1:nf2do) = fields_to_do_
  endif
endif


! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return


! Assertion on D-Grid winds
! -------------------------
if (dxg%has_field('ud')) call abor1_ftn("GeoVaLs state should not have D-Grid winds")


! Winds
! -----
have_winds = .false.
if (dxm%has_field('ud')) then
  call dxm%get_field('ud', ud)
  call dxm%get_field('vd', vd)
  allocate(ua(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(va(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call d_to_a(geom, ud, vd, ua, va)
  have_winds = .true.
elseif (dxm%has_field('ua')) then
    call dxm%get_field('ua', ua)
    call dxm%get_field('va', va)
    have_winds = .true.
endif



! Virtual temperature needed but now done in VADER
! ------------------------------------------------


! Humidity mixing ratio
! ---------------------
have_qmr = .false.
if (allocated(self%q) .and. dxm%has_field('sphum')) then
  call dxm%get_field('sphum', q)
  allocate(qmr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call crtm_mixratio_tl(geom, self%q, q, qmr)
  have_qmr = .true.
endif

! Cloud liquid water
! ------------------
have_ql = .false.
if (allocated(self%ql).and.allocated(self%delp).and.dxm%has_field('liq_wat').and.&
   dxg%has_field('mass_content_of_cloud_liquid_water_in_atmosphere_layer')) then
  call dxm%get_field('liq_wat', cmxr)
  allocate(clwpath(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call hydro_mixr_to_wpath_tl (geom, self%delp, self%ql, cmxr, clwpath)
  have_ql = .true.
endif

! Cloud ice water
! ---------------
have_qi = .false.
if (allocated(self%qi).and.allocated(self%delp).and.dxm%has_field('ice_wat').and.&
   dxg%has_field('mass_content_of_cloud_ice_in_atmosphere_layer')) then
  call dxm%get_field('ice_wat', cmxr)
  allocate(ciwpath(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call hydro_mixr_to_wpath_tl (geom, self%delp, self%qi, cmxr, ciwpath)
  have_qi = .true.
endif

! Rain
! ----
have_qr = .false.
if (allocated(self%qr).and.allocated(self%delp).and.dxm%has_field('rainwat').and.&
   dxg%has_field('mass_content_of_rain_in_atmosphere_layer')) then
  call dxm%get_field('rainwat', cmxr)
  allocate(crwpath(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call hydro_mixr_to_wpath_tl (geom, self%delp, self%qr, cmxr, crwpath)
  have_qr = .true.
endif

! Snow
! ----
have_qs = .false.
if (allocated(self%qs).and.allocated(self%delp).and.dxm%has_field('snowwat').and.&
   dxg%has_field('mass_content_of_snow_in_atmosphere_layer')) then
  call dxm%get_field('snowwat', cmxr)
  allocate(cswpath(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call hydro_mixr_to_wpath_tl (geom, self%delp, self%qs, cmxr, cswpath)
  have_qs = .true.
endif

! Graupel
! -------
have_qg = .false.
if (allocated(self%qg).and.allocated(self%delp).and.dxm%has_field('graupel').and.&
   dxg%has_field('mass_content_of_graupel_in_atmosphere_layer')) then
  call dxm%get_field('graupel', cmxr)
  allocate(cgwpath(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call hydro_mixr_to_wpath_tl (geom, self%delp, self%qg, cmxr, cgwpath)
  have_qg = .true.
endif

! Ozone
! -----
have_o3   = .false.
if (dxm%has_field( 'o3mr')) then
  call dxm%get_field('o3mr', o3mr)
  allocate(o3ppmv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  o3ppmv = o3mr * constant('constoz')
  have_o3 = .true.
elseif (dxm%has_field('o3ppmv')) then
  call dxm%get_field('o3ppmv', o3ppmv)
  have_o3 = .true.
endif

if (have_o3) then
  if (.not.allocated(self%o3)) call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiply no ozone" // &
                                              "found in trajectory")
  do k = 1, self%npz
    do j = self%jsc, self%jec
      do i = self%isc, self%iec
        if (self%o3(i,j,k) < 0.0_kind_real ) then
          o3ppmv(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif

! Surface pressure
! ----------------
have_ps = .false.
if (dxm%has_field( 'ps')) then
  call dxm%get_field('ps', ps)
  have_ps = .true.
elseif (dxm%has_field('delp')) then
  call dxm%get_field('delp', delp)
  allocate(ps(self%isc:self%iec,self%jsc:self%jec,1))
  ps(:,:,1) = sum(delp,3)
  have_ps = .true.
endif

! Skin temperature
! ----------------
have_tskin = .false.
if (dxm%has_field( 'skin_temperature_at_surface')) then
  call dxm%get_field('skin_temperature_at_surface', tskin)
  have_tskin = .true.
endif

! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call dxg%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ua")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = ua

  case ("va")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = va

! Virtual temperature needed but now done in VADER
!   case ("tv")

  case ("ps")

    if (.not. have_ps) call field_fail(fields_to_do(f))
    field_ptr = ps

  case ("skin_temperature_at_surface")

    if (.not. have_tskin) call field_fail(fields_to_do(f))
    field_ptr = tskin

  case ("water_vapor_mixing_ratio_wrt_dry_air")

    if (.not. have_qmr) call field_fail(fields_to_do(f))
    field_ptr = qmr

  case ("o3ppmv", "mole_fraction_of_ozone_in_air")

    if (.not. have_o3) call field_fail(fields_to_do(f))
    field_ptr = o3ppmv

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

    if (have_ql) field_ptr = clwpath

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")
    if (have_qi) field_ptr = ciwpath

  case ("mass_content_of_rain_in_atmosphere_layer")
    if (have_qr) field_ptr = crwpath

  case ("mass_content_of_snow_in_atmosphere_layer")
    if (have_qs) field_ptr = cswpath

  case ("mass_content_of_graupel_in_atmosphere_layer")

    if (have_qg) field_ptr = cgwpath

  ! Simulated but not assimilated
  case ("skin_temperature_at_surface_where_sea")
  case ("skin_temperature_at_surface_where_land")
  case ("skin_temperature_at_surface_where_ice")
  case ("skin_temperature_at_surface_where_snow")
  case ("pe")
  case ("p")

  case default

    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiply unknown field: "//trim(fields_to_do(f)) &
                   //". Not in input field and no transform case specified.")

  end select

enddo

if(allocated(ua)) deallocate(ua)
if(allocated(va)) deallocate(va)
if(allocated(qmr)) deallocate(qmr)
if(allocated(o3ppmv)) deallocate(o3ppmv)
if(allocated(clwpath)) deallocate(clwpath)
if(allocated(ciwpath)) deallocate(ciwpath)
if(allocated(crwpath)) deallocate(crwpath)
if(allocated(cswpath)) deallocate(cswpath)
if(allocated(cgwpath)) deallocate(cgwpath)

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiplyadjoint(self, geom, dxg, dxm)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(inout) :: geom
type(fv3jedi_increment),          intent(in)    :: dxg
type(fv3jedi_increment),          intent(inout) :: dxm

integer :: fg, fm, i, j, k, dxg_index, num_not_copied
real(kind=kind_real), pointer :: field_ptr(:,:,:)
character(len=field_clen), allocatable :: fields_to_do(:)
character(len=field_clen) :: not_copied_(size(dxm%fields))
logical, allocatable :: field_passed(:)
integer :: noassim_index

!Winds
logical :: have_awinds, have_dwinds
integer :: ua_index, va_index
real(kind=kind_real), pointer     :: ua   (:,:,:)         !A-grid wind u component
real(kind=kind_real), pointer     :: va   (:,:,:)         !A-grid wind v component
real(kind=kind_real), allocatable :: ud   (:,:,:)         !D-grid wind u component
real(kind=kind_real), allocatable :: vd   (:,:,:)         !D-grid wind v component

!Virtual temperature
logical :: have_tv
integer :: tv_index
real(kind=kind_real), pointer     :: tv   (:,:,:)         !Virtual temperature
real(kind=kind_real), allocatable :: q_tv (:,:,:)         !Specific humidity
real(kind=kind_real), allocatable :: t_tv (:,:,:)         !Temperature
real(kind=kind_real), pointer     :: tptr (:,:,:)         !Temperature

!Humidity mixing ratio
logical :: have_qmr
integer :: qmr_index
real(kind=kind_real), pointer     :: qmr  (:,:,:)         !Virtual temperature
real(kind=kind_real), allocatable :: q_qmr(:,:,:)         !Specific humidity
real(kind=kind_real), pointer     :: qptr (:,:,:)         !Specific humidity

!Cloud liquid water mixing ratio
logical :: have_ql,have_qi,have_qr,have_qs,have_qg
real(kind=kind_real), pointer     :: wpath (:,:,:)        !Water path

real(kind=kind_real), allocatable :: dql (:,:,:)          !Cloud liq water mixing ratio ad
real(kind=kind_real), allocatable :: dqi (:,:,:)          !Cloud ice water mixing ratio ad
real(kind=kind_real), allocatable :: dqr (:,:,:)          !Rain water mixing ratio ad
real(kind=kind_real), allocatable :: dqs (:,:,:)          !Snow water mixing ratio ad
real(kind=kind_real), allocatable :: dqg (:,:,:)          !Graupel mixing ratio ad
integer :: ql_index
integer :: qi_index
integer :: qr_index
integer :: qs_index
integer :: qg_index

!Ozone mixing ratio
logical :: have_o3mr, have_o3ppmv
integer :: o3_index
real(kind=kind_real), allocatable :: o3mr  (:,:,:)        !Ozone mixing ratio
real(kind=kind_real), allocatable :: o3ppmv(:,:,:)        !Ozone ppmv

!Surface pressure
logical :: have_ps
integer :: ps_index
real(kind=kind_real), pointer     :: ps   (:,:,:)         !Surface pressure
real(kind=kind_real), allocatable :: delp (:,:,:)         !Pressure thickness

!Skin temperature
logical :: have_tsea,have,have_tland,have_tice,have_tsnow,have_tskin
integer :: tskin_index
real(kind=kind_real), pointer     :: dtsea   (:,:,:)       !Sea surface temperature
real(kind=kind_real), pointer     :: dtland  (:,:,:)       !Land surface temperature
real(kind=kind_real), pointer     :: dtice   (:,:,:)       !Ice surface temperature
real(kind=kind_real), pointer     :: dtsnow  (:,:,:)       !Snow surface temperature

! initialize pointers
ps_index=0
ql_index=0
qi_index=0
qr_index=0
qs_index=0
qg_index=0

! Print information
!if (geom%f_comm%rank()==0) then
!  do fg = 1, size(dxg%fields)
!    print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields IN: ", trim(dxg%fields(fg)%short_name), &
!            minval(dxg%fields(fg)%array), maxval(dxg%fields(fg)%array)
!  enddo
!  do fm = 1, size(dxm%fields)
!    print*, "Model2GeoVaLs.multiplyAD, Model fields IN:   ", trim(dxm%fields(fm)%short_name), &
!            minval(dxm%fields(fm)%array), maxval(dxm%fields(fm)%array)
!  enddo
!endif

! Keep track of input fields passed to output
allocate(field_passed(size(dxg%fields)))
field_passed = .false.

! Loop over model fields
num_not_copied = 0
do fm = 1, size(dxm%fields)
  ! Identity if found and not winds
  if (.not.trim(dxm%fields(fm)%short_name) == 'ua' .and. &
      .not.trim(dxm%fields(fm)%short_name) == 'va' .and. &
      dxg%has_field( dxm%fields(fm)%short_name, dxg_index)) then
    call dxg%get_field(dxm%fields(fm)%short_name, field_ptr)
    dxm%fields(fm)%array = dxm%fields(fm)%array + field_ptr
    field_passed(dxg_index) = .true.
  else
    num_not_copied = num_not_copied + 1
    not_copied_(num_not_copied) = dxm%fields(fm)%short_name
  endif
enddo

allocate(fields_to_do(num_not_copied))
fields_to_do(1:num_not_copied) = not_copied_(1:num_not_copied)


! Winds
! -----
have_awinds = .false.
have_dwinds = .false.
if (dxg%has_field( "ua", ua_index) .and. dxg%has_field( "va", va_index)) then
  call dxg%get_field('ua', ua)
  call dxg%get_field('va', va)
  if (dxm%has_field('ud')) then
    allocate(ud(self%isc:self%iec  ,self%jsc:self%jec+1,self%npz))
    allocate(vd(self%isc:self%iec+1,self%jsc:self%jec  ,self%npz))
    ud = 0.0_kind_real
    vd = 0.0_kind_real
    call d_to_a_ad(geom, ud, vd, ua, va)
    have_dwinds = .true.
  elseif (dxm%has_field('ua')) then
    have_awinds = .true.
  else
    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiplyadjoint: Winds found in GeoVaLs but"// &
                   " not in the model.")
  endif
endif


! Virtual temperature needed but now done in VADER
! ------------------------------------------------

! Humidity mixing ratio
! ---------------------
have_qmr = .false.
if (allocated(self%q) .and. dxg%has_field('water_vapor_mixing_ratio_wrt_dry_air', qmr_index)) then
  call dxg%get_field('water_vapor_mixing_ratio_wrt_dry_air', qmr)
  allocate(q_qmr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  q_qmr = 0.0_kind_real
  call crtm_mixratio_ad(geom, self%q, q_qmr, qmr)
  have_qmr = .true.
endif

! Pressure
! --------
have_ps = .false.
if (dxg%has_field( "ps", ps_index)) then
  call dxg%get_field('ps', ps)
  allocate(delp(self%isc:self%iec,self%jsc:self%jec,self%npz))
  delp = 0.0_kind_real
  do k = 1, self%npz
!   delp(:,:,k) = delp(:,:,k) + (geom%bk(k+1)-geom%bk(k))*ps(:,:,1)
    delp(:,:,k) = delp(:,:,k) + ps(:,:,1)
  enddo
  have_ps = .true.
endif

! Cloud liquid water
! ------------------
have_ql = .false.
if (allocated(self%ql).and.allocated(self%delp).and.dxm%has_field('liq_wat').and.&
    dxg%has_field('mass_content_of_cloud_liquid_water_in_atmosphere_layer',ql_index)) then
  call dxg%get_field('mass_content_of_cloud_liquid_water_in_atmosphere_layer', wpath)
  allocate(dql(self%isc:self%iec,self%jsc:self%jec,self%npz))
  dql=0.0_kind_real
  call hydro_mixr_to_wpath_ad (geom, self%delp, self%ql, dql, wpath)
  have_ql = .true.
endif

! Cloud ice water
! ---------------
have_qi = .false.
if (allocated(self%qi).and.allocated(self%delp).and.dxm%has_field('ice_wat').and.&
    dxg%has_field('mass_content_of_cloud_ice_in_atmosphere_layer',qi_index)) then
  call dxg%get_field('mass_content_of_cloud_ice_in_atmosphere_layer', wpath)
  allocate(dqi(self%isc:self%iec,self%jsc:self%jec,self%npz))
  dqi=0.0_kind_real
  call hydro_mixr_to_wpath_ad (geom, self%delp, self%qi, dqi, wpath)
  have_qi = .true.
endif

! Rain
! ----
have_qr = .false.
if (allocated(self%qr).and.allocated(self%delp).and.dxm%has_field('rainwat').and.&
  dxg%has_field('mass_content_of_rain_in_atmosphere_layer',qr_index)) then
  call dxg%get_field('mass_content_of_rain_in_atmosphere_layer', wpath)
  allocate(dqr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  dqr=0.0_kind_real
  call hydro_mixr_to_wpath_ad (geom, self%delp, self%qr, dqr, wpath)
  have_qr = .true.
endif

! Snow
! ----
have_qs = .false.
if (allocated(self%qs).and.allocated(self%delp).and.dxm%has_field('snowwat').and.&
  dxg%has_field('mass_content_of_snow_in_atmosphere_layer',qs_index)) then
  call dxg%get_field('mass_content_of_snow_in_atmosphere_layer', wpath)
  allocate(dqs(self%isc:self%iec,self%jsc:self%jec,self%npz))
  dqs=0.0_kind_real
  call hydro_mixr_to_wpath_ad (geom, self%delp, self%qs, dqs, wpath)
  have_qs = .true.
endif

! Graupel
! -------
have_qg = .false.
if (allocated(self%qg).and.allocated(self%delp).and.dxm%has_field('graupel').and.&
  dxg%has_field('mass_content_of_graupel_in_atmosphere_layer',qg_index)) then
  call dxg%get_field('mass_content_of_graupel_in_atmosphere_layer', wpath)
  allocate(dqg(self%isc:self%iec,self%jsc:self%jec,self%npz))
  dqg=0.0_kind_real
  call hydro_mixr_to_wpath_ad (geom, self%delp, self%qg, dqg, wpath)
  have_qg = .true.
endif

! Skin temperature
! ----------------
have_tland = .false.
have_tsea  = .false.
have_tice  = .false.
have_tsnow = .false.
have_tskin = .false.
if (dxm%has_field('skin_temperature_at_surface',tskin_index)) then
  have_tskin = .true.
else if (dxm%has_field('tsea',tskin_index)) then
  have_tskin = .true.
else if (dxm%has_field('ts',tskin_index)) then
  have_tskin = .true.
endif
if (dxm%has_field('skin_temperature_at_surface',tskin_index)) then
   if ( allocated(self%frland).and.dxg%has_field('skin_temperature_at_surface_where_land') ) then
     call dxg%get_field('skin_temperature_at_surface_where_land',dtland)
     have_tland = .true.
   endif
   if ( allocated(self%frocean).and.dxg%has_field('skin_temperature_at_surface_where_sea') ) then
     call dxg%get_field('skin_temperature_at_surface_where_sea',dtsea)
     have_tsea = .true.
   endif
   if ( allocated(self%frseaice).and.dxg%has_field('skin_temperature_at_surface_where_ice') ) then
     call dxg%get_field('skin_temperature_at_surface_where_ice',dtice)
     have_tice = .true.
   endif
   if ( allocated(self%frsnow).and.dxg%has_field('skin_temperature_at_surface_where_snow') ) then
     call dxg%get_field('skin_temperature_at_surface_where_snow',dtsnow)
     have_tsnow = .true.
   endif
endif

! Ozone
! -----
have_o3ppmv = .false.
if (dxg%has_field('o3ppmv', o3_index)) then
  call dxg%get_field('o3ppmv', o3ppmv)
  have_o3ppmv = .true.
elseif (dxg%has_field('mole_fraction_of_ozone_in_air', o3_index)) then
  call dxg%get_field('mole_fraction_of_ozone_in_air', o3ppmv)
  have_o3ppmv = .true.
endif

have_o3mr = .false.
if (have_o3ppmv) then
  if (.not.allocated(self%o3)) call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiply no ozone" // &
                                              "found in trajectory")
  do k = 1, self%npz
    do j = self%jsc, self%jec
      do i = self%isc, self%iec
        if (self%o3(i,j,k) < 0.0_kind_real ) then
          o3ppmv(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo

  ! Adjoint of ppmv to mixing ratio
  allocate(o3mr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  o3mr = o3ppmv * constant('constoz')
  have_o3mr = .true.

endif

! Simulated but not assimilated
if (dxg%has_field( "mass_content_of_cloud_liquid_water_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "mass_content_of_cloud_ice_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "mass_content_of_rain_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "mass_content_of_snow_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "mass_content_of_graupel_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "skin_temperature_at_surface_where_sea", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "skin_temperature_at_surface_where_land", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "skin_temperature_at_surface_where_ice", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "skin_temperature_at_surface_where_snow", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "p", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "pe", noassim_index)) &
  field_passed(noassim_index) = .true.


! Print information
!if (geom%f_comm%rank()==0) then
! do fg = 1, size(dxg%fields)
!   print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields: ", trim(dxg%fields(fg)%short_name), &
!           ", passed: ", field_passed(fg)
! enddo
!endif


! Loop over the fields not obtainable from the input
! --------------------------------------------------
do fm = 1, size(fields_to_do)

  call dxm%get_field(trim(fields_to_do(fm)), field_ptr)

  select case(trim(fields_to_do(fm)))

  case ("ud")

    if (have_dwinds) then
      field_passed(ua_index) = .true.
      field_ptr = field_ptr + ud
    endif

  case ("vd")

    if (have_dwinds) then
      field_passed(va_index) = .true.
      field_ptr = field_ptr + vd
    endif

  case ("ua")

    if (have_awinds) then
      field_passed(ua_index) = .true.
      field_ptr = field_ptr + ua
    endif

  case ("va")

    if (have_awinds) then
      field_passed(va_index) = .true.
      field_ptr = field_ptr + va
    endif

  case ("delp")

    if (have_ps) then
      field_passed(ps_index) = .true.
      field_ptr = field_ptr + delp
    endif

  case ("o3mr")

    if (have_o3mr) then
      field_passed(o3_index) = .true.
      field_ptr = field_ptr + o3mr
    endif

  case ("o3ppmv")

    if (have_o3ppmv) then
      field_passed(o3_index) = .true.
      field_ptr = field_ptr + o3ppmv
    endif

  case ("liq_wat")

    if (have_ql) then
      field_passed(ql_index) = .true.
      field_ptr = field_ptr + dql
    endif

  case ("ice_wat")

    if (have_qi) then
      field_passed(qi_index) = .true.
      field_ptr = field_ptr + dqi
    endif

  case ("rainwat")

    if (have_qr) then
      field_passed(qr_index) = .true.
      field_ptr = field_ptr + dqr
    endif

  case ("snowwat")

    if (have_qs) then
      field_passed(qs_index) = .true.
      field_ptr = field_ptr + dqs
    endif

  case ("graupel")

    if (have_qg) then
      field_passed(qg_index) = .true.
      field_ptr = field_ptr + dqg
    endif

  case ("ts","skin_temperature_at_surface", "tsea")

    if (have_tsea) then
      field_passed(tskin_index) = .true.
      where (self%frocean>0.0_kind_real)
        field_ptr = field_ptr + dtsea
      endwhere
    endif

    if (have_tland) then
      field_passed(tskin_index) = .true.
      where (self%frland>0.0_kind_real)
        field_ptr = field_ptr + dtland
      endwhere
    endif

    if (have_tice) then
      field_passed(tskin_index) = .true.
      where (self%frseaice>0.0_kind_real)
        field_ptr = field_ptr + dtice
      endwhere
    endif

    if (have_tsnow) then
      field_passed(tskin_index) = .true.
      where (self%frsnow>0.0_kind_real)
        field_ptr = field_ptr + dtsnow
      endwhere
    endif

  end select

enddo

! Print information
!if (geom%f_comm%rank()==0) then
! do fg = 1, size(dxg%fields)
!   print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields: ", trim(dxg%fields(fg)%short_name), &
!           "Passed: ", field_passed(fg)
! enddo
!endif


! Loop over fields in input not yet assigned an output field
do fg = 1, size(dxg%fields)

  if (.not. field_passed(fg)) then

    select case(trim(dxg%fields(fg)%short_name))

    case ("water_vapor_mixing_ratio_wrt_dry_air")

      if (.not. have_qmr) call field_fail(trim(dxg%fields(fg)%short_name))
      field_passed(qmr_index) = .true.
      call dxm%get_field("sphum", qptr)
      qptr = qptr + q_qmr

    case default

      call abor1_ftn("GeoVaLs field "//trim(dxg%fields(fg)%short_name)//" has no known link "// &
                      "to fields in model state")

    end select

  endif

enddo


! Check all fields have been linked to an output field
do fg = 1, size(dxg%fields)
  if (.not. field_passed(fg)) then
    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiplyadjoint failed to send all geoval "// &
                   "fields to a model field")
  endif
enddo

! ! Print information
! if (geom%f_comm%rank()==0) then
!   do fg = 1, size(dxg%fields)
!     print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields OUT: ", trim(dxg%fields(fg)%short_name), minval(dxg%fields(fg)%array), maxval(dxg%fields(fg)%array)
!   enddo
!   do fm = 1, size(dxm%fields)
!     print*, "Model2GeoVaLs.multiplyAD, Model fields OUT:   ", trim(dxm%fields(fm)%short_name), minval(dxm%fields(fg)%array), maxval(dxm%fields(fg)%array)
!   enddo
! endif

if(allocated(delp)) deallocate(delp)
if(allocated(o3mr)) deallocate(o3mr)
if(allocated(dql)) deallocate(dql)
if(allocated(dqi)) deallocate(dqi)
if(allocated(dqr)) deallocate(dqr)
if(allocated(dqs)) deallocate(dqs)
if(allocated(dqg)) deallocate(dqg)
if(allocated(q_qmr)) deallocate(q_qmr)

end subroutine multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module fv3jedi_lvc_model2geovals_mod
