! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_constants_mod

! ieee
use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

! fv3jedi
use fv3jedi_kinds_mod, only: kind_real

implicit none

real(kind=kind_real), parameter :: rad2deg = 57.2957779186820_kind_real
real(kind=kind_real), parameter :: deg2rad = 0.01745329300562541_kind_real
real(kind=8)        , parameter :: pi_r8   = 3.14159265358979323846
real(kind=kind_real), parameter :: pi      = real(pi_r8,kind_real)
real(kind=kind_real), parameter :: grav    = 9.80665_kind_real
real(kind=kind_real), parameter :: radius  = 6371.0e3_kind_real
real(kind=kind_real), parameter :: omega   = 2.0_kind_real*pi/86164.0_kind_real
real(kind=kind_real), parameter :: stfbol  = 5.6734e-8_kind_real
real(kind=kind_real), parameter :: airmw   = 28.965_kind_real
real(kind=kind_real), parameter :: h2omw   = 18.015_kind_real
real(kind=kind_real), parameter :: o3mw    = 47.9982_kind_real
real(kind=kind_real), parameter :: runiv   = 8314.47_kind_real
real(kind=kind_real), parameter :: alhl    = 2.4665e6_kind_real
real(kind=kind_real), parameter :: alhf    = 3.3370e5_kind_real
real(kind=kind_real), parameter :: alhs    = alhl+alhf
real(kind=kind_real), parameter :: rdry    = runiv/airmw
real(kind=kind_real), parameter :: cpdry   = 3.5_kind_real*rdry
real(kind=kind_real), parameter :: cvdry   = cpdry-rdry
real(kind=kind_real), parameter :: rvap    = runiv/h2omw
real(kind=kind_real), parameter :: cpvap   = 4._kind_real*rvap
real(kind=kind_real), parameter :: cvvap   = cpvap-rvap
real(kind=kind_real), parameter :: kappa   = rdry/cpdry
real(kind=kind_real), parameter :: epsilon = h2omw/airmw
real(kind=kind_real), parameter :: deltap  = cpvap/cpdry
real(kind=kind_real), parameter :: deltav  = cvvap/cvdry
real(kind=kind_real), parameter :: gammad  = cpdry/cvdry
real(kind=kind_real), parameter :: rgas    = rdry
real(kind=kind_real), parameter :: cp      = rgas/kappa
real(kind=kind_real), parameter :: zvir    = rvap/rgas - 1._kind_real
real(kind=kind_real), parameter :: vireps  = 1.0_kind_real/epsilon-1.0_kind_real
real(kind=kind_real), parameter :: p00     = 100000.0_kind_real
real(kind=kind_real), parameter :: ps      = 101300.0_kind_real
real(kind=kind_real), parameter :: capice  = 2000._kind_real
real(kind=kind_real), parameter :: capwtr  = 4218._kind_real
real(kind=kind_real), parameter :: rhowtr  = 1000._kind_real
real(kind=kind_real), parameter :: nuair   = 1.533e-5_kind_real
real(kind=kind_real), parameter :: tice    = 273.16_kind_real
real(kind=kind_real), parameter :: srfprs  = 98470_kind_real
real(kind=kind_real), parameter :: karman  = 0.40_kind_real
real(kind=kind_real), parameter :: usmin   = 1.00_kind_real
real(kind=kind_real), parameter :: avogad  = 6.023e26_kind_real
real(kind=kind_real), parameter :: rho_seawater  = 1026.0_kind_real
real(kind=kind_real), parameter :: rho_seaice    = 917.0_kind_real
real(kind=kind_real), parameter :: rho_snow      = 330.0_kind_real
real(kind=kind_real), parameter :: f_coriolis_angle = 0.0_kind_real
!real(kind=kind_real), parameter :: kg_kg_to_ppmv_o3 = 1.0e6_kind_real*airmw/o3mw ! kg/kg -> ppmv
real(kind=kind_real), parameter :: constoz = 604229.0_kind_real ! kg/kg -> ppmv (what the GSI uses, but of unknown origin)
real(kind=kind_real), parameter :: kap1 = kappa + 1.0_kind_real
real(kind=kind_real), parameter :: kapr = 1.0_kind_real/kappa

contains

real(kind=kind_real) function nan()
  real(kind=kind_real) :: nan_in
  nan = ieee_value(nan_in, ieee_quiet_nan)
end function nan

end module fv3jedi_constants_mod
