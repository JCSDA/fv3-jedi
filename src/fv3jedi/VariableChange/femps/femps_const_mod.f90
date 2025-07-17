! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_const_mod

use femps_kinds_mod

implicit none
private

! Pi
real(kind=kind_real), parameter, public :: piby4 = atan(1.0_kind_real)
real(kind=kind_real), parameter, public :: pi = 4.0_kind_real*piby4
real(kind=kind_real), parameter, public :: piby2 = 0.5_kind_real*pi

! Earth's radius
real(kind=kind_real), parameter, public :: rearth = 6371220.0_kind_real

! Radians to degrees
real(kind=kind_real), parameter, public :: rad2deg = 57.2957779186820_kind_real
real(kind=kind_real), parameter, public :: deg2rad = 0.01745329300562541_kind_real

end module femps_const_mod
