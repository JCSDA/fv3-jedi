! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module radii_vt_mod

      use fv3jedi_kinds_mod, only: kind_real
      use fv3jedi_geom_mod, only: fv3jedi_geom
      use fv3jedi_constants_mod, only: constant
      use fckit_log_module, only: fckit_log

      implicit none
      private

      public crtm_ade_efr
      public effectiveRadius
      public calculateNumber
      public thompson08_snow
      public wgamma

      character(len=256):: debug_msg
      real, parameter :: PI = ACOS(-1.)

      !..Snow, rain, and graupel are large so their contribution to Vis/IR is minimal,
      !.. but larger particles are much more influential to microwave and have a
      !.. characteristic size more closely related to mass-squared.
      logical :: assume_microwave = .true.

      !..Pre-computed ratio of gamma functions for calculating moments of distribution.
      real, dimension(0:15), parameter:: g_ratio=(/6,24,60,120,210,336, &
     &                 504,720,990,1320,1716,2184,2730,3360,4080,4896/)

      !..Any hydrometeor with less than min_qx is ignored as zero cloud/precip.
      real(kind=kind_real), parameter :: min_qx = 1.0E-8_kind_real

      !..A user can override the spherical water drops using
      !.. auxillary variables that are optionally passed into functions
      !.. where they are called a_mass, b_mass, and mu_x.
      real, parameter, private:: am_x = PI * 1000./6.   ! Spherical water drops
      real, parameter, private:: bm_x = 3.
      integer, parameter, private:: nu_x = 0            ! Inverse exponential
      double precision:: intercept, lambda

      !..A pretty good approximation of ice crystal size as a function of
      !.. temperature from -94 to 0C by Jon Egill Kristjansson and coauthors.
      !.. Only used when ice number is absent and not using GFDL microphys option.
      real retab(95)
      data retab /5.92779, 6.26422, 6.61973, 6.99539, 7.39234,          &
                  7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930, &
                  10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
                  15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955, &
                  20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125, &
                  27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
                  31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
                  34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078, &
                  38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635, &
                  42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221, &
                  50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898, &
                  65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833, &
                  93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
                  124.954, 130.630, 136.457, 142.446, 148.608, 154.956, &
                  161.503, 168.262, 175.248, 182.473, 189.952, 197.699, &
                  205.728, 214.055, 222.694, 231.661, 240.971, 250.639/

contains

!>----------------------------------------------------------------------------
!> Compute cloud area density and effective radius for the crtm --------------
!>----------------------------------------------------------------------------

subroutine crtm_ade_efr( geom,p,T,delp,sea_frac,q,ql,qi,qr,qs,qg,nc,ni,nr,ns,ng, &
                         ql_ade,qi_ade,qr_ade,qs_ade,qg_ade,                     &
                         ql_efr,qi_efr,qr_efr,qs_efr,qg_efr, method, use_mask)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom
real(kind=kind_real), intent(in)  :: p(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Pressure | Pa
real(kind=kind_real), intent(in)  :: t(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Temperature | K
real(kind=kind_real), intent(in)  :: delp(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)  !Layer thickness | Pa
real(kind=kind_real), intent(in)  :: sea_frac(geom%isc:geom%iec,geom%jsc:geom%jec)          !Sea fraction | 1
real(kind=kind_real), intent(in)  :: q(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)     !Specific humidity | kg/kg
real(kind=kind_real), intent(in)  :: ql(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of cloud liquid water | kg/kg
real(kind=kind_real), intent(in)  :: qi(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of cloud ice | kg/kg
real(kind=kind_real), intent(in), optional  :: qr(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of rain | kg/kg
real(kind=kind_real), intent(in), optional  :: qs(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of snow | kg/kg
real(kind=kind_real), intent(in), optional  :: qg(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Mixing ratio of graupel | kg/kg
real(kind=kind_real), intent(in), optional  :: nc(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Number of cloud droplets | /kg
real(kind=kind_real), intent(in), optional  :: ni(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Number of cloud ice | /kg
real(kind=kind_real), intent(in), optional  :: nr(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Number of rain drops | /kg
real(kind=kind_real), intent(in), optional  :: ns(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Number of snow particles | /kg
real(kind=kind_real), intent(in), optional  :: ng(geom%isc:geom%iec,geom%jsc:geom%jec, 1:geom%npz)    !Number of graupel particles | /kg

real(kind=kind_real), intent(out) :: ql_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for cloud liquid water | kg/m^2
real(kind=kind_real), intent(out) :: qi_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for cloud ice | kg/m^2
real(kind=kind_real), intent(out) :: ql_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for cloud liquid water | microns
real(kind=kind_real), intent(out) :: qi_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for cloud ice | microns
real(kind=kind_real), intent(out), optional :: qr_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for rain | kg/m^2
real(kind=kind_real), intent(out), optional :: qs_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for snow | kg/m^2
real(kind=kind_real), intent(out), optional :: qg_ade(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !area density for graupel | kg/m^2
real(kind=kind_real), intent(out), optional :: qr_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for rain | microns
real(kind=kind_real), intent(out), optional :: qs_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for snow | microns
real(kind=kind_real), intent(out), optional :: qg_efr(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz) !efr for graupel | microns

character(len=*):: method              ! Method for radii calculation (thompson, gfdl, or gsi)
character(len=*), optional:: use_mask  ! Set to land or sea to mask either one, otherwise no mask

!Locals
integer :: isc,iec,jsc,jec,npz,i,j,k
logical, allocatable :: seamask(:,:,:)
real(kind=kind_real) :: tem1, tem2, kgkg_to_kgm2, xqq
logical :: have_qr, have_qs, have_qg, have_nc, have_ni, have_nr, have_ns, have_ng
logical :: mask_land, mask_sea
real :: tempK, wcontent, nconc, answer, ygra1, zans1
integer :: mu, idx_rei
real(kind=kind_real), allocatable :: rho_air(:,:,:)
real, allocatable :: nnc(:,:,:), nnr(:,:,:), nng(:,:,:)
real(kind=kind_real) :: rdry, grav, tice, zvir
!+---+

if (geom%f_comm%rank() == 0 ) then
  debug_msg = 'Inside crtm_ade_efr routine to compute effective radii'
  call fckit_log%debug(debug_msg)
endif

! Constants
rdry = constant('rdry')
grav = constant('grav')
tice = constant('tice')
zvir = constant('zvir')

mask_land = .false.
mask_sea = .false.
have_qr = .false.
have_qs = .false.
have_qg = .false.
have_nc = .false.
have_ni = .false.
have_nr = .false.
have_ns = .false.
have_ng = .false.
if(present(qr)) have_qr = .true.
if(present(qs)) have_qs = .true.
if(present(qg)) have_qg = .true.
if(present(nc)) have_nc = .true.
if(present(ni)) have_ni = .true.
if(present(nr)) have_nr = .true.
if(present(ns)) have_ns = .true.
if(present(ng)) have_ng = .true.
if(present(use_mask)) then
  if (use_mask .eq. 'land') then
    mask_land = .true.
  elseif (use_mask .eq. 'sea') then
    mask_sea = .true.
  elseif (use_mask .eq. 'none') then
    if (geom%f_comm%rank() == 0 ) then
      debug_msg = 'DEBUG,  not masking either land or sea'
      call fckit_log%debug(debug_msg)
    endif
  else
    if (geom%f_comm%rank() == 0 ) then
      debug_msg = 'ERROR,  mask over must be set to either land or sea'
      call fckit_log%debug(debug_msg)
    endif
  endif
endif

! Grid convenience
! ----------------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

! Allocate some convenient arrays
! -------------------------------
allocate(seamask(isc:iec,jsc:jec,1:npz))  ! 2-D really, but 3-D is convenient for where statement
allocate(rho_air(isc:iec,jsc:jec,1:npz))


! Set outputs to zero
! -------------------
ql_ade = 0.0_kind_real
qi_ade = 0.0_kind_real
ql_efr = 0.0_kind_real
qi_efr = 0.0_kind_real
if(have_qr) then
  qr_ade = 0.0_kind_real
  qr_efr = 0.0_kind_real
endif
if(have_qs) then
  qs_ade = 0.0_kind_real
  qs_efr = 0.0_kind_real
endif
if(have_qg) then
  qg_ade = 0.0_kind_real
  qg_efr = 0.0_kind_real
endif

! Sea mask
! --------
seamask = .false.
do j = jsc,jec
  do i = isc,iec
     seamask(i,j,:) = min(max(0.0_kind_real,sea_frac(i,j)),1.0_kind_real)  >= 0.99_kind_real
  enddo
enddo

! Calculate air density
! ---------------------
do k = 1,npz
   do j = jsc,jec
     do i = isc,iec
        rho_air(i,j,k) = p(i,j,k)/(rdry*t(i,j,k)* (1.0_kind_real + zvir * max(q(i,j,k),0.0_kind_real)))
     enddo
   enddo
enddo

! Convert hydrometeor mixing ratio to liquid/ice water path (kg/kg to kg/m^2)
! ---------------------------------------------------------------------------
do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
         if (ql(i,j,k) .lt. min_qx) CYCLE
         kgkg_to_kgm2 = delp(i,j,k) / grav
         ql_ade(i,j,k) = ql(i,j,k) * kgkg_to_kgm2
    enddo
  enddo
enddo

do k = 1,npz
  do j = jsc,jec
    do i = isc,iec
         if (qi(i,j,k) .lt. min_qx) CYCLE
         kgkg_to_kgm2 = delp(i,j,k) / grav
         qi_ade(i,j,k) = qi(i,j,k) * kgkg_to_kgm2
    enddo
  enddo
enddo

if( have_qr ) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
         if (qr(i,j,k) .lt. min_qx) CYCLE
         kgkg_to_kgm2 = delp(i,j,k) / grav
         qr_ade(i,j,k) = qr(i,j,k) * kgkg_to_kgm2
      enddo
    enddo
  enddo
endif

if( have_qs ) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
         if (qs(i,j,k) .lt. min_qx) CYCLE
         kgkg_to_kgm2 = delp(i,j,k) / grav
         qs_ade(i,j,k) = qs(i,j,k) * kgkg_to_kgm2
      enddo
    enddo
  enddo
endif

if( have_qg ) then
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
         if (qg(i,j,k) .lt. min_qx) CYCLE
         kgkg_to_kgm2 = delp(i,j,k) / grav
         qg_ade(i,j,k) = qg(i,j,k) * kgkg_to_kgm2
      enddo
    enddo
  enddo
endif

! Compute effective radius using various methods.
! -----------------------------------------------

!..The Thompson method follows Thompson and Eidhammer (2014), Thompson et al (2008, 2004).
!.. This method is based on Slingo (1989) that relates the radiative effective radius
!.. as the 3rd momenent divided by 2nd moment of the distribution, which should be
!.. appropriate for shortwave and longwave IR wavelengths.  Rain and graupel are
!.. considered too big compared to those wavelengths, but radii are computed anyway.
!.. If the host model passes the number concentrations of cloud water, cloud ice, and
!.. rain, then the following code utilizes each properly.  When lacking number concentration
!.. variables, it is not true to the scheme, but we make a few assumptions to estimate
!.. particle number before computing size.

if (method .eq. 'thompson') then
  if (geom%f_comm%rank() == 0 ) then
    debug_msg = 'DEBUG,  using Thompson method for radii calculations'
    call fckit_log%debug(debug_msg)
  endif
  allocate(nnc(isc:iec,jsc:jec,1:npz))
  if (.not. have_nc) then
    where (seamask)
      nnc = 75.E6        ! 75 drops per cc over ocean
    elsewhere
      nnc = 250.E6       ! 250 drops per cc over land
    end where
  else
    nnc = nc
  endif
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        if (ql(i,j,k) .lt. min_qx) CYCLE
        wcontent = ql(i,j,k)*rho_air(i,j,k)
        nconc    = nnc(i,j,k)*rho_air(i,j,k)
        mu = MAX(2, MIN((NINT(1000.E6/nconc) + 2), 15))
        answer = effectiveRadius(rx=wcontent, nx=nconc, mu=mu)
        ql_efr(i,j,k) = max(2.0_kind_real, min(answer*1.0E6_kind_real, 25.0_kind_real))
      enddo
    enddo
  enddo
  deallocate(nnc)

  if (.not. have_ni) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qi(i,j,k) .lt. min_qx) CYCLE
          idx_rei = int(t(i,j,k)-179._kind_real)
          idx_rei = min(max(idx_rei,1),94)
          qi_efr(i,j,k) = retab(idx_rei)*1.0_kind_real
        enddo
      enddo
    enddo
  else
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qi(i,j,k) .lt. min_qx) CYCLE
          wcontent = qi(i,j,k)*rho_air(i,j,k)
          nconc    = ni(i,j,k)*rho_air(i,j,k)
          answer = effectiveRadius(rx=wcontent, nx=nconc)
          qi_efr(i,j,k) = max(2.5_kind_real, min(answer*1.0E6_kind_real, 250.0_kind_real))
        enddo
      enddo
    enddo
  endif

  if (have_qr) then
    allocate(nnr(isc:iec,jsc:jec,1:npz))
    if (.not. have_nr) then
      do k = 1,npz
        do j = jsc,jec
          do i = isc,iec
            if (qr(i,j,k) .lt. min_qx) CYCLE
            wcontent = qr(i,j,k)*rho_air(i,j,k)
            nnr(i,j,k) = calculateNumber(rx=wcontent, N0_exp=8.E6)/rho_air(i,j,k)
          enddo
        enddo
      enddo
    else
      nnr = nr
    endif
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qr(i,j,k) .lt. min_qx) CYCLE
          wcontent = qr(i,j,k)*rho_air(i,j,k)
          nconc    = nnr(i,j,k)*rho_air(i,j,k)
          answer = effectiveRadius(rx=wcontent, nx=nconc)
          qr_efr(i,j,k) = max(50.0_kind_real, min(answer*1.0E6_kind_real, 1000.0_kind_real))
        enddo
      enddo
    enddo
    if (assume_microwave) qr_efr = qr_efr*2.0_kind_real
    deallocate(nnr)
  endif

  if (have_qs) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qs(i,j,k) .lt. min_qx) CYCLE
          wcontent = qs(i,j,k)*rho_air(i,j,k)
          tempK = t(i,j,k)
          answer = thompson08_snow(rs=wcontent, temp=tempK)
          qs_efr(i,j,k) = max(5.0_kind_real, min(answer*1.0E6_kind_real, 5000.0_kind_real))
        enddo
      enddo
    enddo
    if (assume_microwave) qs_efr = qs_efr*2.0_kind_real
  endif

  if (have_qg) then
    allocate(nng(isc:iec,jsc:jec,1:npz))
    if (.not. have_ng) then
      do k = 1,npz
        do j = jsc,jec
          do i = isc,iec
            if (qg(i,j,k) .lt. min_qx) CYCLE
            wcontent = qg(i,j,k)*rho_air(i,j,k)
            ygra1 = alog10(max(1.E-9, wcontent))
            zans1 = (3.4 + 2./7. * (ygra1+8.))
            zans1 = MAX(2., MIN(zans1, 7.))
            nng(i,j,k) = calculateNumber(rx=wcontent, N0_exp=10.**(zans1))/rho_air(i,j,k)
          enddo
        enddo
      enddo
    else
      nng = ng
    endif
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qg(i,j,k) .lt. min_qx) CYCLE
          wcontent = qg(i,j,k)*rho_air(i,j,k)
          nconc    = nng(i,j,k)*rho_air(i,j,k)
          answer = effectiveRadius(rx=wcontent, nx=nconc)
          qg_efr(i,j,k) = max(150.0_kind_real, min(answer*1.0E6_kind_real, 5000.0_kind_real))
        enddo
      enddo
    enddo
    if (assume_microwave) qg_efr = qg_efr*2.0_kind_real
    deallocate(nng)
  endif

!..GFDL basically has one-moment rain, snow, and graupel with constant intercept
!.. parameters. Their final calculator of size involves strange factors that
!.. deviate from 3rd/2nd moment.  It results in larger sizes than the pure
!.. mathematics of the classical 3rd/2nd moments the following manner:
!.  a) rain: instead of (3+Mu) with mu=0, it is as if (3+0.9038) prefactor;
!.  b) snow: instead of (3+Mu) with mu=0, it is as if (3+0.6356) prefactor;
!.  c) graupel: instead of (3+Mu) with mu=0, it is as if (3+0.7583) prefactor;
!.  d) cloud water: instead of (3+Mu), it is as if mu is set to -2 in prefactor;
!.  e) cloud ice: completely lacking traditional 3rd/2nd moment calculator and
!.   is based on temperature and ice mixing ratio only.

elseif (method .eq. 'gfdl') then

  if (geom%f_comm%rank() == 0 ) then
    debug_msg = 'DEBUG,  using GFDL method for radii calculations'
    call fckit_log%debug(debug_msg)
  endif
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        if (ql(i,j,k) .lt. min_qx) CYCLE
        wcontent = ql(i,j,k)*rho_air(i,j,k)
        nconc    = 100.E6
        answer = exp(1./3. * log((3*wcontent)/(4.*PI*1000.*nconc)))
        ql_efr(i,j,k) = max(2.0_kind_real, min(answer*1.0E6_kind_real, 25.0_kind_real))
      enddo
    enddo
  enddo

  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        if (qi(i,j,k) .lt. min_qx) CYCLE
        wcontent = qi(i,j,k)*rho_air(i,j,k)
        if (t(i,j,k)-tice .lt. -50) then
          answer  = 1.22/9.917 * exp((1 - 0.891) * log(1.0e3*wcontent))*1.0e-3
        elseif (t(i,j,k)-tice .lt. -40) then
          answer  = 1.22/9.337 * exp((1 - 0.920) * log(1.0e3*wcontent))*1.0e-3
        elseif (t(i,j,k)-tice .lt. -30) then
          answer  = 1.22/9.208 * exp((1 - 0.945) * log(1.0e3*wcontent))*1.0e-3
        else
          answer  = 1.22/9.387 * exp((1 - 0.969) * log(1.0e3*wcontent))*1.0e-3
        endif
        qi_efr(i,j,k) = max(2.0_kind_real, min(answer*1.0E6_kind_real, 250.0_kind_real))
      enddo
    enddo
  enddo

  if (have_qr) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qr(i,j,k) .lt. min_qx) CYCLE
          wcontent = qr(i,j,k)*rho_air(i,j,k)
          nconc = calculateNumber(rx=wcontent, N0_exp=8.E6)
          answer = effectiveRadius(rx=wcontent, nx=nconc)
          answer = 3.9038/3.0 * answer    ! Adjust for exact match to GFDL scheme.
          qr_efr(i,j,k) = max(50.0_kind_real, min(answer*1.0E6_kind_real, 1000.0_kind_real))
        enddo
      enddo
    enddo
  endif

  if (have_qs) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qs(i,j,k) .lt. min_qx) CYCLE
          wcontent = qs(i,j,k)*rho_air(i,j,k)
          nconc =  calculateNumber(rx=wcontent, N0_exp=2.E6)
          answer = effectiveRadius(rx=wcontent, nx=nconc, a=PI*100./6., b=3.0)
          answer = 3.6356/3.0 * answer    ! Adjust for exact match to GFDL scheme.
          qs_efr(i,j,k) = max(10.0_kind_real, min(answer*1.0E6_kind_real, 2000.0_kind_real))
        enddo
      enddo
    enddo
  endif

  if (have_qg) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qg(i,j,k) .lt. min_qx) CYCLE
          wcontent = qg(i,j,k)*rho_air(i,j,k)
          nconc =  calculateNumber(rx=wcontent, N0_exp=4.E6)
          answer = effectiveRadius(rx=wcontent, nx=nconc, a=PI*500./6., b=3.0)
          answer = 3.7583/3.0 * answer    ! Adjust for exact match to GFDL scheme.
          qg_efr(i,j,k) = max(50.0_kind_real, min(answer*1.0E6_kind_real, 5000.0_kind_real))
        enddo
      enddo
    enddo
  endif

elseif (method .eq. 'gsi') then

  if (geom%f_comm%rank() == 0 ) then
    debug_msg = 'DEBUG,  using GSI method for radii calculations'
    call fckit_log%debug(debug_msg)
  endif
  ! Cloud liquid water effective radius
  ! -----------------------------------
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        if (ql(i,j,k) .lt. min_qx) CYCLE
        tem1 = (tice-t(i,j,k))*0.05_kind_real
        ql_efr(i,j,k) = 5.0_kind_real + 5.0_kind_real * min(1.0_kind_real, tem1)
        ql_efr(i,j,k) = max(1.0_kind_real,ql_efr(i,j,k))
      enddo
    enddo
  enddo

  ! Cloud ice water effective radius
  ! ---------------------------------
  do k = 1,npz
    do j = jsc,jec
      do i = isc,iec
        if (qi(i,j,k) .lt. min_qx) CYCLE
        tem2 = t(i,j,k) - tice
        wcontent = qi(i,j,k)*rho_air(i,j,k)
        if (tem2 < -50.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.917_kind_real)*wcontent**0.109_kind_real
        elseif (tem2 < -40.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.337_kind_real)*wcontent**0.08_kind_real
        elseif (tem2 < -30.0_kind_real ) then
           qi_efr(i,j,k) =  (1250._kind_real/9.208_kind_real)*wcontent**0.055_kind_real
        else
           qi_efr(i,j,k) =  (1250._kind_real/9.387_kind_real)*wcontent**0.031_kind_real
        endif
        qi_efr(i,j,k) = max(5.0_kind_real,qi_efr(i,j,k))
      enddo
    enddo
  enddo

  ! Rain water effective radius (Taken from set_crtm_cloudmod.f90 in GSI for GEOS qr & qs .)
  ! ---------------------------------
  if( have_qr ) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qr(i,j,k) .lt. min_qx) CYCLE
          wcontent = 1000.0_kind_real * qr(i,j,k)*rho_air(i,j,k)
          xqq = log10(wcontent)
          qr_efr(i,j,k) = 7.934_kind_real*xqq*xqq*xqq + 90.858_kind_real*xqq*xqq &
                        + 387.807_kind_real*xqq + 679.939_kind_real
          qr_efr(i,j,k) = max(qr_efr(i,j,k), 100.0_kind_real)
        enddo
      enddo
    enddo
  endif

  ! Snow water effective radius (Taken from set_crtm_cloudmod.f90 in GSI for GEOS qr & qs.)
  ! ---------------------------------
  if( have_qs ) then
    do k = 1,npz
      do j = jsc,jec
        do i = isc,iec
          if (qs(i,j,k) .lt. min_qx) CYCLE
          wcontent = 1000.0_kind_real * qs(i,j,k)*rho_air(i,j,k)
          xqq = log10(wcontent)
          qs_efr(i,j,k) = 9.33_kind_real*xqq*xqq*xqq +  84.779_kind_real*xqq*xqq &
                        + 351.1345_kind_real*xqq + 691.391_kind_real !Liu DDA_type5
          qs_efr(i,j,k) = max(qs_efr(i,j,k), 100.0_kind_real)
        enddo
      enddo
    enddo
  endif

else

  if (geom%f_comm%rank() == 0 ) then
    debug_msg = 'WARNING, hydrometeor effective radii method must be one of: gfdl, thompson, gsi'
    call fckit_log%debug(debug_msg)
  endif

endif

! If a land/sea mask is desired, just zero out water/ice path variables.
! ----------------------------------------------------------------------
if (mask_land) then
   where (.not. seamask)
      ql_ade = 0.0_kind_real
      qi_ade = 0.0_kind_real
   end where
   if (have_qr) then
      where (.not. seamask)
         qr_ade = 0.0_kind_real
      end where
   endif
   if (have_qs) then
      where (.not. seamask)
         qs_ade = 0.0_kind_real
      end where
   endif
   if (have_qg) then
      where (.not. seamask)
         qg_ade = 0.0_kind_real
      end where
   endif
elseif (mask_sea) then
   where (seamask)
      ql_ade = 0.0_kind_real
      qi_ade = 0.0_kind_real
   end where
   if (have_qr) then
      where (seamask)
         qr_ade = 0.0_kind_real
      end where
   endif
   if (have_qs) then
      where (seamask)
         qs_ade = 0.0_kind_real
      end where
   endif
   if (have_qg) then
      where (seamask)
         qg_ade = 0.0_kind_real
      end where
   endif
endif

deallocate(seamask)
deallocate(rho_air)

end subroutine crtm_ade_efr

!+---+-----------------------------------------------------------------+

      real function effectiveRadius(rx, nx, a, b, mu)

      real, intent(in):: rx, nx
      real, optional, intent(in):: a, b
      integer, optional, intent(in):: mu
      real:: obmx, xcontent, xnumber
      real:: a_mass, b_mass
      integer:: mu_x, n
      real, dimension(0:15):: g_rat

      if (present(a)) then
         a_mass = a
      else
         a_mass = am_x
      endif
      if (present(b)) then
         b_mass = b
      else
         b_mass = bm_x
      endif
      obmx = 1./b_mass
      if (present(mu)) then
         mu_x = mu
         mu_x = MAX(0, MIN(mu, 15))
      else
         mu_x = nu_x
      endif
      if ((b_mass.lt.2.9999) .or. (b_mass.gt.3.0001)) then
         do n = 0, 15
            g_rat(n) = WGAMMA(b_mass + mu_x + 1.)/WGAMMA(mu_x + 1.)
         enddo
      else
         g_rat = g_ratio
      endif

      xnumber = MAX(1.E-6, nx)
      xcontent = MAX(1.E-12, rx)

      lambda = (a_mass*g_rat(mu_x)*xnumber/xcontent)**obmx
      effectiveRadius = SNGL(0.5D0 * DBLE(3.+mu_x)/lambda)

   end function effectiveRadius

!+---+-----------------------------------------------------------------+

      real function calculateNumber(rx, N0_exp, a, b, mu)

      real, intent(in):: rx, N0_exp
      real, optional, intent(in):: a, b
      integer, optional, intent(in):: mu
      real:: obmx
      double precision:: lam_exp
      integer:: n
      real:: a_mass, b_mass
      integer:: mu_x
      !..Variables to hold exponents and gamma values.
      real, dimension(3,0:15):: ce, cg
      real, dimension(0:15)::  ocg1, ocg2, ocg3
      save ce, cg, ocg1, ocg2, ocg3

      if (present(a)) then
         a_mass = a
      else
         a_mass = am_x
      endif
      if (present(b)) then
         b_mass = b
      else
         b_mass = bm_x
      endif
      obmx = 1./b_mass
      if (present(mu)) then
         mu_x = mu
         mu_x = MAX(0, MIN(mu, 15))
      else
         mu_x = nu_x
      endif

      do n = 0, 15
         ce(1,n) = n + 1.
         ce(2,n) = b_mass + n + 1.
         ce(3,n) = b_mass + 1.
         cg(1,n) = WGAMMA(ce(1,n))
         cg(2,n) = WGAMMA(ce(2,n))
         cg(3,n) = WGAMMA(ce(3,n))
         ocg1(n) = 1./cg(1,n)
         ocg2(n) = 1./cg(2,n)
         ocg3(n) = 1./cg(3,n)
      enddo

      !..Always assuming the intercept is from exponential distribution,
      !.. then if user wants a gamma distrib with N0_exp reference, we
      !.. can convert to the proper lambda to get final number.

      lam_exp = (N0_exp*a_mass*cg(3,1)/rx)**(1./ce(3,1))
      lambda = lam_exp * (cg(2,mu_x)*ocg1(mu_x)*ocg3(mu_x))**obmx

      !..Double-check final intercept when incoming mu=0, this should match N0_exp
      !..intercept = N0_exp/(cg(1,mu_x)*lam_exp) * lambda**ce(1,mu_x)
      !..print*, '    double-CHECK intercept = ', intercept

      calculateNumber = cg(1,mu_x)*ocg2(mu_x)*rx*lambda**b_mass / a_mass

   end function calculateNumber

!+---+-----------------------------------------------------------------+

      real function thompson08_snow(rs, temp)

      real, intent(in):: rs, temp
      !..For snow moments conversions (from Field et al. 2005)
      REAL, DIMENSION(10), PARAMETER :: &
      sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285,      &
     &        0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
      REAL, DIMENSION(10), PARAMETER :: &
      sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141,      &
     &        0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)
      real, parameter :: am_s = 0.069
      real, parameter :: bm_s = 2.0
      real, dimension(1), parameter :: cse = (/ bm_s + 1. /)
      real :: tc0, smo2, smob, smoc, smoz, a_, b_, loga_, rs_eff

      !..Calculate bm_s+1 (th) moment, smoc.  Useful for diameter calcs.
      tc0 = MIN(-0.1, temp-273.15)
      smob = rs/am_s
      loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
            + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
            + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
            + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
            + sa(10)*cse(1)*cse(1)*cse(1)
      a_ = 10.0**loga_
      b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
         + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
         + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
         + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
      smoc = a_ * smob**b_

      rs_eff = MAX(2.51E-6, MIN(0.5*(smoc/smob), 1999.E-6))

      thompson08_snow = rs_eff

      end function thompson08_snow

!+---+-----------------------------------------------------------------+
      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL, INTENT(IN):: XX
      DOUBLE PRECISION, PARAMETER:: STP = 2.5066282746310005D0
      DOUBLE PRECISION, DIMENSION(6), PARAMETER:: &
               COF = (/76.18009172947146D0, -86.50532032941677D0, &
                       24.01409824083091D0, -1.231739572450155D0, &
                      .1208650973866179D-2, -.5395239384953D-5/)
      DOUBLE PRECISION:: SER,TMP,X,Y
      INTEGER:: J

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION WGAMMA(y)

      IMPLICIT NONE
      REAL, INTENT(IN):: y

      WGAMMA = EXP(GAMMLN(y))

      END FUNCTION WGAMMA
!+---+-----------------------------------------------------------------+

end module radii_vt_mod
