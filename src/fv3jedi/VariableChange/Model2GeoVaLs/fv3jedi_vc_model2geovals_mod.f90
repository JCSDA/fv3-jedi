! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_vc_model2geovals_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,           only: fckit_log

use datetime_mod

use fv3jedi_constants_mod, only: constoz, grav
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_field_mod,     only: copy_subset, field_clen
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state

use height_vt_mod
use moisture_vt_mod
use pressure_vt_mod
use surface_vt_mod
use temperature_vt_mod
use wind_vt_mod

implicit none

private
public :: fv3jedi_vc_model2geovals

type :: fv3jedi_vc_model2geovals
  integer :: isc, iec, jsc, jec, npz
  character(len=10) :: tropprs_method
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: changevar
end type fv3jedi_vc_model2geovals

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_vc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),              intent(in)    :: geom
type(fckit_configuration),       intent(in)    :: conf

character(len=:), allocatable :: str

! Method to use for tropopause pressure ([gsi] or thompson)
if (.not. conf%get("tropopause pressure method", str)) str = 'gsi'
self%tropprs_method = trim(str)

! Grid convenience
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_vc_model2geovals), intent(inout) :: self

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self, geom, xm, xg)

class(fv3jedi_vc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),              intent(inout) :: geom
type(fv3jedi_state),             intent(in)    :: xm
type(fv3jedi_state),             intent(inout) :: xg

integer :: f, jlev, i, j, k
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

! Specific humidity
logical :: have_q
real(kind=kind_real), pointer     :: q     (:,:,:)         !Specific humidity

! Relative humidity
logical :: have_rh
real(kind=kind_real), allocatable :: rh    (:,:,:)         !Relative humidity
real(kind=kind_real), pointer     :: qsat   (:,:,:)        !Saturation specific humidity

! Pressure fields
logical :: have_pressures
real(kind=kind_real), allocatable :: ps    (:,:,:)         !Surface pressure
real(kind=kind_real), allocatable :: delp  (:,:,:)         !Pressure thickness
real(kind=kind_real), allocatable :: prsi  (:,:,:)         !Pressure, interfaces
real(kind=kind_real), allocatable :: prs   (:,:,:)         !Pressure, midpoint

! Temperature fields
logical :: have_t
real(kind=kind_real), allocatable :: t     (:,:,:)         !Temperature
real(kind=kind_real), pointer     :: pt    (:,:,:)         !Potential temperature
real(kind=kind_real), allocatable :: pkz   (:,:,:)         !Pressure to the kapaa

! Vitual temperature
logical :: have_tv
real(kind=kind_real), allocatable :: tv    (:,:,:)         !Virtual temperature

! Geopotential heights
logical :: have_geoph
real(kind=kind_real), allocatable :: geophi(:,:,:)         !Geopotential height, interfaces
real(kind=kind_real), allocatable :: geoph (:,:,:)         !Geopotential height
real(kind=kind_real), allocatable :: suralt(:,:,:)         !Surface altitude
real(kind=kind_real), pointer     :: phis  (:,:,:)         !Surface geopotential height times gravity
logical,  parameter :: use_compress = .true.

! Ozone
logical :: have_o3
real(kind=kind_real), allocatable :: o3mr  (:,:,:)         !Ozone mixing ratio
real(kind=kind_real), allocatable :: o3ppmv(:,:,:)         !Ozone ppmv

! Winds
logical :: have_winds
real(kind=kind_real), allocatable :: ua    (:,:,:)         !Eastward wind
real(kind=kind_real), allocatable :: va    (:,:,:)         !Northward wind
real(kind=kind_real), pointer     :: ud    (:,:,:)         !u component D-grid
real(kind=kind_real), pointer     :: vd    (:,:,:)         !v component D-grid

! Surface roughness length
logical :: have_zorl
real(kind=kind_real), pointer     :: zorl     (:,:,:)      !Incoming in centimeters
real(kind=kind_real), allocatable :: sfc_rough(:,:,:)      !Outgoing in meters

! Sea-land mask
logical :: have_slmsk
real(kind=kind_real), allocatable :: slmsk   (:,:,:)       !Land-sea mask
real(kind=kind_real), pointer     :: frocean (:,:,:)       !Fraction ocean
real(kind=kind_real), pointer     :: frlake  (:,:,:)       !Fraction lake
real(kind=kind_real), pointer     :: frseaice(:,:,:)       !Fraction seaice
real(kind=kind_real), pointer     :: tsea    (:,:,:)       !Surface temperature

!f10m
logical :: have_f10m
real(kind=kind_real), allocatable :: f10m    (:,:,:)       !Surface wind reduction factor
real(kind=kind_real), pointer     :: u_srf   (:,:,:)
real(kind=kind_real), pointer     :: v_srf   (:,:,:)
real(kind=kind_real) :: wspd

!qiql
logical :: have_qiql
real(kind=kind_real), allocatable :: qi      (:,:,:)
real(kind=kind_real), allocatable :: ql      (:,:,:)
real(kind=kind_real), pointer     :: qils    (:,:,:)
real(kind=kind_real), pointer     :: qicn    (:,:,:)
real(kind=kind_real), pointer     :: qlls    (:,:,:)
real(kind=kind_real), pointer     :: qlcn    (:,:,:)

!qrqs
logical :: have_qrqs
real(kind=kind_real), allocatable :: qr      (:,:,:)
real(kind=kind_real), allocatable :: qs      (:,:,:)
real(kind=kind_real), pointer     :: qrls    (:,:,:)
real(kind=kind_real), pointer     :: qrcn    (:,:,:)
real(kind=kind_real), pointer     :: qsls    (:,:,:)
real(kind=kind_real), pointer     :: qscn    (:,:,:)

!CRTM mixing ratio
logical :: have_qmr
real(kind=kind_real), allocatable :: qmr     (:,:,:)       !Land-sea mask

!CRTM moisture fields
logical :: have_crtm_cld
real(kind=kind_real), allocatable :: ql_ade  (:,:,:)
real(kind=kind_real), allocatable :: qi_ade  (:,:,:)
real(kind=kind_real), allocatable :: qr_ade  (:,:,:)
real(kind=kind_real), allocatable :: qs_ade  (:,:,:)
real(kind=kind_real), allocatable :: ql_efr  (:,:,:)
real(kind=kind_real), allocatable :: qi_efr  (:,:,:)
real(kind=kind_real), allocatable :: qr_efr  (:,:,:)
real(kind=kind_real), allocatable :: qs_efr  (:,:,:)
real(kind=kind_real), allocatable :: watercov(:,:)

!Salinity
logical :: have_sss
real(kind=kind_real), allocatable :: sss     (:,:,:)

!CRTM surface
logical :: have_crtm_surface, have_soilt, have_soilm
real(kind=kind_real), pointer     :: sheleg   (:,:,:)
real(kind=kind_real), pointer     :: vtype    (:,:,:)
real(kind=kind_real), pointer     :: stype    (:,:,:)
real(kind=kind_real), pointer     :: vfrac    (:,:,:)
real(kind=kind_real), allocatable :: soilt    (:,:,:)
real(kind=kind_real), allocatable :: soilm    (:,:,:)
real(kind=kind_real), pointer     :: soil_tmp (:,:,:)
real(kind=kind_real), allocatable :: land_type_index_npoess                    (:,:,:)
real(kind=kind_real), allocatable :: land_type_index_igbp                      (:,:,:)
real(kind=kind_real), allocatable :: vegetation_type_index                     (:,:,:)
real(kind=kind_real), allocatable :: soil_type                                 (:,:,:)
real(kind=kind_real), allocatable :: water_area_fraction                       (:,:,:)
real(kind=kind_real), allocatable :: land_area_fraction                        (:,:,:)
real(kind=kind_real), allocatable :: ice_area_fraction                         (:,:,:)
real(kind=kind_real), allocatable :: surface_snow_area_fraction                (:,:,:)
real(kind=kind_real), allocatable :: leaf_area_index                           (:,:,:)
real(kind=kind_real), allocatable :: surface_temperature_where_sea             (:,:,:)
real(kind=kind_real), allocatable :: surface_temperature_where_land            (:,:,:)
real(kind=kind_real), allocatable :: surface_temperature_where_ice             (:,:,:)
real(kind=kind_real), allocatable :: surface_temperature_where_snow            (:,:,:)
real(kind=kind_real), allocatable :: volume_fraction_of_condensed_water_in_soil(:,:,:)
real(kind=kind_real), allocatable :: vegetation_area_fraction                  (:,:,:)
real(kind=kind_real), allocatable :: soil_temperature                          (:,:,:)
real(kind=kind_real), allocatable :: surface_snow_thickness                    (:,:,:)
real(kind=kind_real), allocatable :: surface_wind_speed                        (:,:,:)
real(kind=kind_real), allocatable :: surface_wind_from_direction               (:,:,:)
real(kind=kind_real), allocatable :: sea_surface_salinity                      (:,:,:)

! Vorticity
logical :: have_vort
real(kind=kind_real), allocatable :: vort     (:,:,:)
real(kind=kind_real), allocatable :: divg     (:,:,:)

! Tropopause pressure
logical :: have_tropprs
real(kind=kind_real), allocatable :: tprs     (:,:,:)



! Identity part of the change of fields
! -------------------------------------
call copy_subset(xm%fields, xg%fields, fields_to_do)


! if (geom%f_comm%rank()==0) then
!   do f = 1, size(xm%fields)
!     print*, "Model2GeoVaLs.changeVar, Model fields:   ", trim(xm%fields(f)%fv3jedi_name)
!   enddo
!   do f = 1, size(xg%fields)
!     print*, "Model2GeoVaLs.changeVar, GeoVaLs fields: ", trim(xg%fields(f)%fv3jedi_name)
!   enddo
!   do f = 1, size(fields_to_do)
!     print*, "Model2GeoVaLs.changeVar, GeoVaLs needed by transform: ", trim(fields_to_do(f))
!   enddo
! endif


! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return


! Get pressures at edge, center & log center
! ------------------------------------------
have_pressures = .false.

if (xm%has_field('delp')) then
  call xm%get_field('delp', delp)
  allocate(ps(self%isc:self%iec, self%jsc:self%jec, 1))
  ps(:,:,1) = sum(delp,3)
  have_pressures = .true.
elseif (xm%has_field('ps')) then
  call xm%get_field('ps', ps)
  allocate(delp(self%isc:self%iec, self%jsc:self%jec, self%npz))
  do jlev = 1,self%npz
    delp(:,:,jlev) = (geom%ak(jlev+1)-geom%ak(jlev))+(geom%bk(jlev+1)-geom%bk(jlev))*ps(:,:,1)
  enddo
  have_pressures = .true.
elseif (xm%has_field('pe')) then
  call xm%get_field('pe', prsi)
  allocate(ps(self%isc:self%iec, self%jsc:self%jec, 1))
  ps(:,:,1) = prsi(:,:,self%npz+1)
  allocate(delp(self%isc:self%iec, self%jsc:self%jec, self%npz))
  do jlev = 1,self%npz
    delp(:,:,jlev) = prsi(:,:,jlev+1) - prsi(:,:,jlev)
  enddo
  have_pressures = .true.
endif

if (have_pressures) then
  if (.not.allocated(prsi)) allocate(prsi(self%isc:self%iec,self%jsc:self%jec,self%npz+1))
  if (.not.allocated(prs )) allocate(prs (self%isc:self%iec,self%jsc:self%jec,self%npz  ))
  call delp_to_pe_p_logp(geom, delp, prsi, prs)
endif


! Temperature
! -----------
have_t = .false.
if (xm%has_field( 't')) then
  call xm%get_field('t', t)
  have_t = .true.
elseif (xm%has_field( 'pt')) then
  if (.not. have_pressures) &
    call abor1_ftn("fv3jedi_vc_model2geovals_mod.changevar: a state with pt needs pressures")
  allocate(t  (self%isc:self%iec, self%jsc:self%jec, self%npz))
  allocate(pkz(self%isc:self%iec, self%jsc:self%jec, self%npz))
  call xm%get_field('pt',  pt)
  call pe_to_pkz(geom, prsi, pkz)
  call pt_to_t(geom, pkz, pt, t)
  have_t = .true.
endif


! Specific humidity
! -----------------
have_q = .false.
if (xm%has_field( 'sphum')) then
  call xm%get_field('sphum',  q)
  have_q = .true.
endif


! Relative humidity
! -----------------
have_rh = .false.
if (xm%has_field('rh')) then
  call xm%get_field('rh', rh)
  have_rh = .true.
elseif (have_t .and. have_pressures .and. have_q) then
  allocate(qsat(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(rh  (self%isc:self%iec,self%jsc:self%jec,self%npz))
  call get_qsat(geom,delp,t,q,qsat)
  call q_to_rh(geom,qsat,q,rh)
  deallocate(qsat)
  have_rh = .true.
endif


! Geopotential height
! -------------------
have_geoph = .false.
if (have_t .and. have_pressures .and. have_q .and. xm%has_field( 'phis')) then
  call xm%get_field('phis',  phis)
  if (.not.allocated(geophi)) allocate(geophi(self%isc:self%iec,self%jsc:self%jec,self%npz+1))
  if (.not.allocated(geoph )) allocate(geoph (self%isc:self%iec,self%jsc:self%jec,self%npz  ))
  if (.not.allocated(suralt)) allocate(suralt(self%isc:self%iec,self%jsc:self%jec,1))
  call geop_height(geom, prs, prsi, t, q, phis(:,:,1), use_compress, geoph)
  call geop_height_levels(geom, prs, prsi, t, q, phis(:,:,1), use_compress, geophi)
  suralt = phis / grav
  have_geoph = .true.
endif


! Virtual temperature
! -------------------
have_tv = .false.
if (have_t .and. have_q) then
  allocate(tv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call t_to_tv(geom, t, q, tv)
  have_tv = .true.
endif


! Ozone
! -----
have_o3   = .false.
if (xm%has_field( 'o3mr')) then
  call xm%get_field('o3mr', o3mr)
  allocate(o3ppmv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  o3ppmv = o3mr * constoz
  have_o3 = .true.
elseif (xm%has_field('o3ppmv')) then
  call xm%get_field('o3ppmv', o3ppmv)
  allocate(o3mr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  o3mr = o3ppmv / constoz
  have_o3 = .true.
endif

if (have_o3) then
  do k = 1, self%npz
    do j = self%jsc, self%jec
      do i = self%isc, self%iec
        if (o3mr(i,j,k) < 0.0_kind_real ) then
          o3mr(i,j,k)   = 0.0_kind_real
          o3ppmv(i,j,k) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif

! Wind transforms
! ---------------
have_winds = .false.
if (xm%has_field('ua')) then
  call xm%get_field('ua', ua)
  call xm%get_field('va', va)
  have_winds = .true.
elseif (xm%has_field('ud')) then
  call xm%get_field('ud', ud)
  call xm%get_field('vd', vd)
  allocate(ua(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(va(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call d_to_a(geom, ud, vd, ua, va)
  have_winds = .true.
  nullify(ud,vd)
endif

! Land sea mask
! -------------
have_slmsk = .false.
if (xm%has_field( 'slmsk')) then
  call xm%get_field('slmsk', slmsk)
  have_slmsk = .true.
elseif ( xm%has_field('frocean' ) .and. xm%has_field('frlake'  ) .and. &
         xm%has_field('frseaice') .and. xm%has_field('tsea'    ) ) then
  call xm%get_field('frocean' , frocean )
  call xm%get_field('frlake'  , frlake  )
  call xm%get_field('frseaice', frseaice)
  call xm%get_field('tsea'    , tsea    )

  allocate(slmsk(self%isc:self%iec,self%jsc:self%jec,1))
  slmsk = 1.0_kind_real !Land
  do j = self%jsc,self%jec
    do i = self%isc,self%iec
      if ( frocean(i,j,1) + frlake(i,j,1) >= 0.6_kind_real) then
        slmsk(i,j,1) = 0.0_kind_real ! Water
      endif
      if ( slmsk(i,j,1) == 0.0_kind_real .and. frseaice(i,j,1) > 0.5_kind_real) then
        slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
      if ( slmsk(i,j,1) == 0.0_kind_real .and. tsea(i,j,1) < 271.4_kind_real ) then
        slmsk(i,j,1) = 2.0_kind_real ! Ice
      endif
    enddo
  enddo
  have_slmsk = .true.
endif


! Transform surface roughness length (zorl) from cm to meters.
! -------------
have_zorl = .false.
allocate(sfc_rough(self%isc:self%iec,self%jsc:self%jec,1))
if (xm%has_field( 'zorl')) then
  call xm%get_field('zorl', zorl)
  have_zorl = .true.
  sfc_rough = zorl*0.01
else
  sfc_rough = 0.01_kind_real
endif


! f10m
! ----
have_f10m = .false.
if (xm%has_field('f10m')) then
  call xm%get_field('f10m', f10m)
  have_f10m = .true.
elseif ( xm%has_field( 'u_srf') .and. xm%has_field( 'v_srf') .and. have_winds ) then
  call xm%get_field('u_srf' , u_srf)
  call xm%get_field('v_srf' , v_srf)

  allocate(f10m(self%isc:self%iec,self%jsc:self%jec,1))
  f10m = sqrt(u_srf**2 + v_srf**2)

  do j = self%jsc,self%jec
    do i = self%isc,self%iec
      wspd = sqrt(ua(i,j,self%npz)**2 +  va(i,j,self%npz)**2)
      if (f10m(i,j,1) > 0.0_kind_real) then
        f10m(i,j,1) = f10m(i,j,1)/wspd
      else
        f10m(i,j,1) = 1.0_kind_real
      endif
    enddo
  enddo
  have_f10m = .true.
endif


! CRTM mixing ratio
! -----------------
have_qmr = .false.
if (have_q) then
  allocate(qmr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call crtm_mixratio(geom, q, qmr)
  have_qmr = .true.
endif


! Clouds
! ------
have_qiql = .false.
if (xm%has_field( 'ice_wat') .and. xm%has_field( 'liq_wat')) then
  call xm%get_field('ice_wat', qi)
  call xm%get_field('liq_wat', ql)
  have_qiql = .true.
elseif (xm%has_field( 'qils') .and. xm%has_field( 'qicn') .and. &
        xm%has_field( 'qlls') .and. xm%has_field( 'qlcn')) then
  call xm%get_field('qils', qils)
  call xm%get_field('qicn', qicn)
  call xm%get_field('qlls', qlls)
  call xm%get_field('qlcn', qlcn)
  allocate(qi(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(ql(self%isc:self%iec,self%jsc:self%jec,self%npz))
  qi = qils + qicn
  ql = qlls + qlcn
  have_qiql = .true.
endif
! ------
have_qrqs = .false.
if (xm%has_field( 'rainwat') .and. xm%has_field( 'snowwat')) then
  call xm%get_field('rainwat', qr)
  call xm%get_field('snowwat', qs)
  have_qrqs = .true.
elseif (xm%has_field( 'qrls') .and. xm%has_field( 'qrcn') .and. &
        xm%has_field( 'qsls') .and. xm%has_field( 'qscn')) then
  call xm%get_field('qrls', qrls)
  call xm%get_field('qrcn', qrcn)
  call xm%get_field('qsls', qsls)
  call xm%get_field('qscn', qscn)
  allocate(qr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(qs(self%isc:self%iec,self%jsc:self%jec,self%npz))
  qr = qrls + qrcn
  qs = qsls + qscn
  have_qrqs = .true.
endif


! Get CRTM moisture fields
! ------------------------
have_crtm_cld = .false.
if (have_slmsk .and. have_t .and. have_pressures .and. have_q .and. have_qiql ) then
  allocate(ql_ade(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(qi_ade(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(ql_efr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(qi_efr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(watercov(self%isc:self%iec,self%jsc:self%jec))
  ql_ade = 0.0_kind_real
  qi_ade = 0.0_kind_real
  ql_efr = 0.0_kind_real
  qi_efr = 0.0_kind_real

  !TODO Is it water_area_fraction or sea_coverage fed in here?
  watercov = 0.0_kind_real
  do j = self%jsc,self%jec
    do i = self%isc,self%iec
      if (slmsk(i,j,1) == 0) watercov(i,j) = 1.0_kind_real
    enddo
  enddo
  if (have_qrqs ) then
    allocate(qr_ade(self%isc:self%iec,self%jsc:self%jec,self%npz))
    allocate(qs_ade(self%isc:self%iec,self%jsc:self%jec,self%npz))
    allocate(qr_efr(self%isc:self%iec,self%jsc:self%jec,self%npz))
    allocate(qs_efr(self%isc:self%iec,self%jsc:self%jec,self%npz))
    qr_ade = 0.0_kind_real
    qs_ade = 0.0_kind_real
    qr_efr = 0.0_kind_real
    qs_efr = 0.0_kind_real
    call crtm_ade_efr( geom, prs, t, delp, watercov, q, ql, qi, &
                     ql_ade, qi_ade, ql_efr, qi_efr, &
                     qr, qs, qr_ade, qs_ade, qr_efr, qs_efr )
  else
    call crtm_ade_efr( geom, prs, t, delp, watercov, q, ql, qi, &
                     ql_ade, qi_ade, ql_efr, qi_efr )
  endif
  have_crtm_cld = .true.
endif


! CRTM surface fields
! -------------------

! Soil temperature
have_soilt = .false.
if (xm%has_field( 'soilt' )) then
  call xm%get_field('soilt' , soilt )
  have_soilt = .true.
elseif (xm%has_field( 'stc' )) then
  call xm%get_field('stc' , soil_tmp )
  allocate(soilt(self%isc:self%iec,self%jsc:self%jec,1))
  soilt(:,:,1) = soil_tmp(:,:,1) ! Which of the 4 levels should we use?
  have_soilt = .true.
endif

! Soil moisture
have_soilm = .false.
if (xm%has_field( 'soilm' )) then
  call xm%get_field('soilm' , soilm )
  have_soilm = .true.
elseif (xm%has_field( 'smc' )) then
  call xm%get_field('smc' , soil_tmp )
  allocate(soilm(self%isc:self%iec,self%jsc:self%jec,1))
  soilm(:,:,1) = soil_tmp(:,:,1) ! Which of the 4 levels should we use?
  have_soilm = .true.
endif

have_crtm_surface = .false.
have_sss = .false.
if ( have_slmsk .and. have_f10m .and. xm%has_field( 'sheleg') .and. &
     xm%has_field( 'tsea'  ) .and. xm%has_field( 'vtype' ) .and. &
     xm%has_field( 'stype' ) .and. xm%has_field( 'vfrac' ) .and. &
     have_soilt .and. have_soilm .and. &
     xm%has_field( 'u_srf' ) .and. xm%has_field( 'v_srf' ) ) then

  call xm%get_field('sheleg', sheleg)
  call xm%get_field('tsea'  , tsea  )
  call xm%get_field('vtype' , vtype )
  call xm%get_field('stype' , stype )
  call xm%get_field('vfrac' , vfrac )
  call xm%get_field('u_srf' , u_srf )
  call xm%get_field('v_srf' , v_srf )

  allocate(land_type_index_npoess                    (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(land_type_index_igbp                      (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(vegetation_type_index                     (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(soil_type                                 (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(water_area_fraction                       (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(land_area_fraction                        (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(ice_area_fraction                         (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_snow_area_fraction                (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(leaf_area_index                           (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_temperature_where_sea             (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_temperature_where_land            (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_temperature_where_ice             (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_temperature_where_snow            (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(volume_fraction_of_condensed_water_in_soil(self%isc:self%iec,self%jsc:self%jec,1))
  allocate(vegetation_area_fraction                  (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(soil_temperature                          (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_snow_thickness                    (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_wind_speed                        (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(surface_wind_from_direction               (self%isc:self%iec,self%jsc:self%jec,1))
  allocate(sea_surface_salinity                      (self%isc:self%iec,self%jsc:self%jec,1))

  allocate(sss(self%isc:self%iec,self%jsc:self%jec,1))
  sss = 0.0_kind_real
  if (xm%has_field( 'sss')) then
    call xm%get_field('sss', sss)
    have_sss = .true.
  endif

  ! Here we compute the surface data for CRTM. A particular note about the surface types:
  !
  ! We compute at once the CRTM surface types supported in fv3-jedi: land types in two possible
  ! classification schemes, soil type, and vegetation type (*). This simplifies the logic, and
  ! should only be a negligible additional expense. A minor optimization would be to pass to
  ! `crtm_surface` information about which surface type was requested from this function
  ! (changevar), so only the requested type would be computed. Typically, this would be land type
  ! (in just a single classification scheme) for vis/IR obs, or soil and vegetation types for
  ! microwave obs.
  !
  ! (*) Note: CRTM also supports the USGS land type classification -- to use this with fv3-jedi, the
  ! mapping from the fv3-jedi land type to the USGS land type must be added to `crtm_surface`. This
  ! would exactly follow the code currently in place for the NPOESS and IGBP classifications.
  call crtm_surface( geom, slmsk, sheleg, tsea, vtype, stype, vfrac, soilt, soilm, u_srf, v_srf, &
                      f10m, sss, land_type_index_npoess, land_type_index_igbp, &
                      vegetation_type_index, soil_type, &
                      water_area_fraction, land_area_fraction, ice_area_fraction, &
                      surface_snow_area_fraction, leaf_area_index, surface_temperature_where_sea, &
                      surface_temperature_where_land, surface_temperature_where_ice, &
                      surface_temperature_where_snow, volume_fraction_of_condensed_water_in_soil, &
                      vegetation_area_fraction, soil_temperature, surface_snow_thickness, &
                      surface_wind_speed, surface_wind_from_direction, sea_surface_salinity)

  have_crtm_surface = .true.

endif


! Vorticity
! ---------
have_vort = .false.
if (xm%has_field('ud') .and. xm%has_field('vd')) then
  call xm%get_field('ud', ud)
  call xm%get_field('vd', vd)
  allocate(vort(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(divg(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call udvd_to_vortdivg(geom, ud, vd, vort, divg)
  have_vort = .true.
  nullify(ud, vd)
elseif (have_winds) then
  allocate(vort(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(divg(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(ud(self%isc:self%iec  ,self%jsc:self%jec+1,self%npz))
  allocate(vd(self%isc:self%iec+1,self%jsc:self%jec  ,self%npz))
  call a_to_d(geom, ua, va, ud, vd)
  call udvd_to_vortdivg(geom, ud, vd, vort, divg)
  have_vort = .true.
endif


! Tropopause pressure
! -------------------
have_tropprs = .false.
if (trim(self%tropprs_method) == "gsi") then
  if (have_vort .and. have_tv .and. have_pressures .and. have_o3) then
    allocate(tprs(self%isc:self%iec,self%jsc:self%jec,1))
    call tropprs(geom, ps, prs, tv, o3ppmv, vort, tprs)
    have_tropprs = .true.
  endif
elseif (trim(self%tropprs_method) == "thompson") then
  if (have_pressures .and. have_geoph .and. have_t) then
    allocate(tprs(self%isc:self%iec,self%jsc:self%jec,1))
    call tropprs_th(geom, prs, geoph, t, tprs)
    have_tropprs = .true.
  endif
endif


! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call xg%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ua")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = ua

  case ("va")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = va

  case ("q")

    if (.not. have_q) call field_fail(fields_to_do(f))
    field_ptr = q

  case ("rh")

    if (.not. have_rh) call field_fail(fields_to_do(f))
    field_ptr = rh

  case ("p")

    if (.not. have_pressures) call field_fail(fields_to_do(f))
    field_ptr = prs

  case ("pe")

    if (.not. have_pressures) call field_fail(fields_to_do(f))
    field_ptr = prsi

  case ("delp")

    if (.not. have_pressures) call field_fail(fields_to_do(f))
    field_ptr = delp

  case ("ps")

    if (.not. have_pressures) call field_fail(fields_to_do(f))
    field_ptr = ps

  case ("t")

    if (.not. have_t) call field_fail(fields_to_do(f))
    field_ptr = t

  case ("tv")

    if (.not. have_tv) call field_fail(fields_to_do(f))
    field_ptr = tv

  case ("mole_fraction_of_ozone_in_air", "o3ppmv")

    if (.not. have_o3) call field_fail(fields_to_do(f))
    field_ptr = o3ppmv

  case ("geopotential_height", "height")

    if (.not. have_geoph) call field_fail(fields_to_do(f))
    field_ptr = geoph

  case ("geopotential_height_levels")

    if (.not. have_geoph) call field_fail(fields_to_do(f))
    field_ptr = geophi

  case ("surface_altitude", "surface_geopotential_height")

    if (.not. have_geoph) call field_fail(fields_to_do(f))
    field_ptr = suralt

  case ("mole_fraction_of_carbon_dioxide_in_air")

    field_ptr = 407.0_kind_real

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) call field_fail(fields_to_do(f))
    field_ptr = qmr

  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = ql_ade

  case ("mass_content_of_cloud_ice_in_atmosphere_layer")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qi_ade

  case ("mass_content_of_rain_in_atmosphere_layer")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qr_ade

  case ("mass_content_of_snow_in_atmosphere_layer")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qs_ade

  case ("effective_radius_of_cloud_liquid_water_particle")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = ql_efr

  case ("effective_radius_of_cloud_ice_particle")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qi_efr

  case ("effective_radius_of_rain_particle")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qr_efr

  case ("effective_radius_of_snow_particle")

    if (.not. have_crtm_cld) call field_fail(fields_to_do(f))
    field_ptr = qs_efr

  case ("water_area_fraction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = water_area_fraction

  case ("land_area_fraction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = land_area_fraction

  case ("ice_area_fraction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = ice_area_fraction

  case ("surface_snow_area_fraction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_snow_area_fraction

  case ("surface_temperature_where_sea")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_temperature_where_sea

  case ("sea_surface_salinity")

    if (.not. have_sss) call field_fail(fields_to_do(f))
    field_ptr = sea_surface_salinity

  case ("surface_temperature_where_land")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_temperature_where_land

  case ("surface_temperature_where_ice")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_temperature_where_ice

  case ("surface_temperature_where_snow")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_temperature_where_snow

  case ("surface_snow_thickness")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_snow_thickness

  case ("vegetation_area_fraction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = vegetation_area_fraction

  case ("surface_wind_speed")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_wind_speed

  case ("surface_wind_from_direction")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_wind_from_direction

  case ("wind_reduction_factor_at_10m")

    if (.not. have_f10m) call field_fail(fields_to_do(f))
    field_ptr = f10m

  case ("leaf_area_index")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = leaf_area_index

  case ("volume_fraction_of_condensed_water_in_soil")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = volume_fraction_of_condensed_water_in_soil

  case ("soil_temperature")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = soil_temperature

  case ("land_type_index_NPOESS")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = land_type_index_npoess

  case ("land_type_index_IGBP")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = land_type_index_igbp

  case ("land_type_index_USGS")

    call abor1_ftn("fv3jedi_vc_model2geovals_mod.changevar does not currently implement " &
                   //"the variable change to the USGS land type classification. " &
                   //"This variable change should be added to surface_variables_mod.")

  case ("surface_roughness_length")

    if (.not. have_zorl) call field_fail(fields_to_do(f))
    field_ptr = sfc_rough

  case ("vegetation_type_index")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = vegetation_type_index

  case ("soil_type")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = soil_type

  case ("vorticity", "vort")

    if (.not. have_vort) call field_fail(fields_to_do(f))
    field_ptr = vort

  case ("tropopause_pressure")

    if (.not. have_tropprs) call field_fail(fields_to_do(f))
    field_ptr = tprs

  case ("average_surface_temperature_within_field_of_view")

    if (.not. have_crtm_surface) call field_fail(fields_to_do(f))
    field_ptr = surface_temperature_where_sea

  case default

    call abor1_ftn("fv3jedi_vc_model2geovals_mod.changevar unknown field: "//trim(fields_to_do(f)) &
                   //". Not in input field and no transform case specified.")

  end select

enddo

if (associated(sheleg)) nullify(sheleg)
if (associated(vtype)) nullify(vtype)
if (associated(stype)) nullify(stype)
if (associated(vfrac)) nullify(vfrac)
if (allocated(soilt)) deallocate(soilt)
if (allocated(soilm)) deallocate(soilm)
if (associated(zorl)) nullify(zorl)
if (associated(field_ptr)) nullify(field_ptr)
if (associated(q)) nullify(q)
if (associated(qsat)) nullify(qsat)
if (associated(pt)) nullify(pt)
if (associated(phis)) nullify(phis)
if (associated(ud)) nullify(ud)
if (associated(vd)) nullify(vd)
if (associated(frocean)) nullify(frocean)
if (associated(frlake)) nullify(frlake)
if (associated(frseaice)) nullify(frseaice)
if (associated(tsea)) nullify(tsea)
if (associated(u_srf)) nullify(u_srf)
if (associated(v_srf)) nullify(v_srf)
if (associated(qils)) nullify(qils)
if (associated(qlls)) nullify(qlls)
if (associated(qrls)) nullify(qrls)
if (associated(qsls)) nullify(qsls)
if (associated(qicn)) nullify(qicn)
if (associated(qlcn)) nullify(qlcn)
if (associated(qrcn)) nullify(qrcn)
if (associated(qscn)) nullify(qscn)

if (allocated(fields_to_do)) deallocate(fields_to_do)
if (allocated(rh)) deallocate(rh)
if (allocated(t)) deallocate(t)
if (allocated(tv)) deallocate(tv)
if (allocated(ps)) deallocate(ps)
if (allocated(delp)) deallocate(delp)
if (allocated(prsi)) deallocate(prsi)
if (allocated(prs)) deallocate(prs)
if (allocated(pkz)) deallocate(pkz)
if (allocated(geophi)) deallocate(geophi)
if (allocated(geoph)) deallocate(geoph)
if (allocated(suralt)) deallocate(suralt)
if (allocated(o3mr)) deallocate(o3mr)
if (allocated(o3ppmv)) deallocate(o3ppmv)
if (allocated(ua)) deallocate(ua)
if (allocated(va)) deallocate(va)
if (allocated(slmsk)) deallocate(slmsk)
if (allocated(sfc_rough)) deallocate(sfc_rough)
if (allocated(f10m)) deallocate(f10m)
if (allocated(ql)) deallocate(ql)
if (allocated(qi)) deallocate(qi)
if (allocated(qr)) deallocate(qr)
if (allocated(qs)) deallocate(qs)
if (allocated(qmr)) deallocate(qmr)
if (allocated(ql_ade)) deallocate(ql_ade)
if (allocated(qi_ade)) deallocate(qi_ade)
if (allocated(qr_ade)) deallocate(qr_ade)
if (allocated(qs_ade)) deallocate(qs_ade)
if (allocated(ql_efr)) deallocate(ql_efr)
if (allocated(qi_efr)) deallocate(qi_efr)
if (allocated(qr_efr)) deallocate(qr_efr)
if (allocated(qs_efr)) deallocate(qs_efr)
if (allocated(watercov)) deallocate(watercov)
if (allocated(sss)) deallocate(sss)
if (allocated(land_type_index_npoess)) deallocate(land_type_index_npoess)
if (allocated(land_type_index_igbp)) deallocate(land_type_index_igbp)
if (allocated(vegetation_type_index)) deallocate(vegetation_type_index)
if (allocated(soil_type)) deallocate(soil_type)
if (allocated(water_area_fraction)) deallocate(water_area_fraction)
if (allocated(land_area_fraction)) deallocate(land_area_fraction)
if (allocated(ice_area_fraction)) deallocate(ice_area_fraction)
if (allocated(surface_snow_area_fraction)) deallocate(surface_snow_area_fraction)
if (allocated(leaf_area_index)) deallocate(leaf_area_index)
if (allocated(surface_temperature_where_sea)) deallocate(surface_temperature_where_sea)
if (allocated(surface_temperature_where_land)) deallocate(surface_temperature_where_land)
if (allocated(surface_temperature_where_ice)) deallocate(surface_temperature_where_ice)
if (allocated(surface_temperature_where_snow)) deallocate(surface_temperature_where_snow)
if (allocated(volume_fraction_of_condensed_water_in_soil)) &
                                              deallocate(volume_fraction_of_condensed_water_in_soil)
if (allocated(vegetation_area_fraction)) deallocate(vegetation_area_fraction)
if (allocated(soil_temperature)) deallocate(soil_temperature)
if (allocated(surface_snow_thickness)) deallocate(surface_snow_thickness)
if (allocated(surface_wind_speed)) deallocate(surface_wind_speed)
if (allocated(surface_wind_from_direction)) deallocate(surface_wind_from_direction)
if (allocated(sea_surface_salinity)) deallocate(sea_surface_salinity)
if (allocated(vort)) deallocate(vort)
if (allocated(tprs)) deallocate(tprs)

end subroutine changevar

! --------------------------------------------------------------------------------------------------

end module fv3jedi_vc_model2geovals_mod
