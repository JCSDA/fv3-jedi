! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_constants_mod

use fv3jedi_kinds_mod, only: kind_real

implicit none
private

real(kind=kind_real), parameter, public :: rad2deg = 57.2957779186820_kind_real
real(kind=kind_real), parameter, public :: deg2rad = 0.01745329300562541_kind_real
real(kind=8)        , parameter, public :: pi_r8   = 3.14159265358979323846
real(kind=kind_real), parameter, public :: pi      = real(pi_r8,kind_real)
real(kind=kind_real), parameter, public :: grav    = 9.80665_kind_real
real(kind=kind_real), parameter, public :: radius  = 6371.0e3_kind_real
real(kind=kind_real), parameter, public :: omega   = 2.0_kind_real*pi/86164.0_kind_real
real(kind=kind_real), parameter, public :: stfbol  = 5.6734e-8_kind_real
real(kind=kind_real), parameter, public :: airmw   = 28.965_kind_real
real(kind=kind_real), parameter, public :: h2omw   = 18.015_kind_real
real(kind=kind_real), parameter, public :: o3mw    = 47.9982_kind_real
real(kind=kind_real), parameter, public :: runiv   = 8314.47_kind_real
real(kind=kind_real), parameter, public :: alhl    = 2.4665e6_kind_real
real(kind=kind_real), parameter, public :: alhf    = 3.3370e5_kind_real
real(kind=kind_real), parameter, public :: alhs    = alhl+alhf
real(kind=kind_real), parameter, public :: rdry    = runiv/airmw
real(kind=kind_real), parameter, public :: cpdry   = 3.5_kind_real*rdry
real(kind=kind_real), parameter, public :: cvdry   = cpdry-rdry
real(kind=kind_real), parameter, public :: rvap    = runiv/h2omw
real(kind=kind_real), parameter, public :: cpvap   = 4._kind_real*rvap
real(kind=kind_real), parameter, public :: cvvap   = cpvap-rvap
real(kind=kind_real), parameter, public :: kappa   = rdry/cpdry
real(kind=kind_real), parameter, public :: epsilon = h2omw/airmw
real(kind=kind_real), parameter, public :: deltap  = cpvap/cpdry
real(kind=kind_real), parameter, public :: deltav  = cvvap/cvdry
real(kind=kind_real), parameter, public :: gammad  = cpdry/cvdry
real(kind=kind_real), parameter, public :: rgas    = rdry
real(kind=kind_real), parameter, public :: cp      = rgas/kappa
real(kind=kind_real), parameter, public :: zvir    = rvap/rgas - 1._kind_real
real(kind=kind_real), parameter, public :: vireps  = 1.0_kind_real/epsilon-1.0_kind_real
real(kind=kind_real), parameter, public :: p00     = 100000.0_kind_real
real(kind=kind_real), parameter, public :: ps      = 101300.0_kind_real
real(kind=kind_real), parameter, public :: capice  = 2000._kind_real
real(kind=kind_real), parameter, public :: capwtr  = 4218._kind_real
real(kind=kind_real), parameter, public :: rhowtr  = 1000._kind_real
real(kind=kind_real), parameter, public :: nuair   = 1.533e-5_kind_real
real(kind=kind_real), parameter, public :: tice    = 273.16_kind_real
real(kind=kind_real), parameter, public :: srfprs  = 98470_kind_real
real(kind=kind_real), parameter, public :: karman  = 0.40_kind_real
real(kind=kind_real), parameter, public :: usmin   = 1.00_kind_real
real(kind=kind_real), parameter, public :: avogad  = 6.023e26_kind_real
real(kind=kind_real), parameter, public :: rho_seawater  = 1026.0_kind_real
real(kind=kind_real), parameter, public :: rho_seaice    = 917.0_kind_real
real(kind=kind_real), parameter, public :: rho_snow      = 330.0_kind_real
real(kind=kind_real), parameter, public :: f_coriolis_angle = 0.0_kind_real
!real(kind=kind_real), parameter, public :: kg_kg_to_ppmv_o3 = 1.0e6_kind_real*airmw/o3mw ! kg/kg -> ppmv 
real(kind=kind_real), parameter, public :: constoz = 604229.0_kind_real ! kg/kg -> ppmv (what the GSI uses, but of unknown origin) 
real(kind=kind_real), parameter, public :: kap1 = kappa + 1.0_kind_real
real(kind=kind_real), parameter, public :: kapr = 1.0_kind_real/kappa

end module fv3jedi_constants_mod
