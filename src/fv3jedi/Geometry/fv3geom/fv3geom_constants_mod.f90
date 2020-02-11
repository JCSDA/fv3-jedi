module fv3geom_constants_mod

use, intrinsic :: iso_c_binding

implicit none

integer,              public, parameter :: kind_real = c_double
real(kind=kind_real), public, parameter :: pi = 3.1415926535897931_kind_real
real(kind=kind_real), public, parameter :: torad = pi/180.0_kind_real
real(kind=kind_real), public, parameter :: radius = 6.3712e+6_kind_real

end module fv3geom_constants_mod
