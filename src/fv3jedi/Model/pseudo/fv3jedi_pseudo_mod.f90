! (C) Copyright 2019-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_pseudo_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use datetime_mod
use duration_mod

! fv3-jedi
use fv3jedi_kinds_mod,   only: kind_real
use fv3jedi_geom_mod,    only: fv3jedi_geom
use fv3jedi_state_mod,   only: fv3jedi_state
use fv3jedi_io_gfs_mod,  only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos

implicit none
private
public :: pseudo_model

! --------------------------------------------------------------------------------------------------

!> Fortran derived type to hold model definition
type :: pseudo_model
  character(len=255)    :: pseudo_type !geos of gfs
  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos
  contains
    procedure :: create
    procedure :: delete
    procedure :: initialize
    procedure :: step
    procedure :: finalize
end type pseudo_model

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
class(pseudo_model),       intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

character(len=:), allocatable :: str

! Get IO type to use
call conf%get_or_die("pseudo_type",str)
self%pseudo_type = str

! Setup IO from config
if (trim(self%pseudo_type) == "geos") then
  call self%geos%setup_conf(geom, conf)
elseif (trim(self%pseudo_type) == "gfs") then
  call self%gfs%setup_conf(conf)
else
  call abor1_ftn("fv3jedi_pseudo_mod: pseudo_type must be geos or gfs")
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(pseudo_model), intent(inout) :: self

if (trim(self%pseudo_type) == "geos") call self%geos%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine initialize(self, state)

implicit none
class(pseudo_model), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: state

end subroutine initialize

! --------------------------------------------------------------------------------------------------

subroutine step(self, state, geom, vdate)

implicit none
class(pseudo_model), intent(inout) :: self
type(fv3jedi_state), intent(inout) :: state
type(fv3jedi_geom),  intent(inout) :: geom
type(datetime),      intent(inout) :: vdate !< Valid datetime after step

if (trim(self%pseudo_type) == "gfs") then
  call self%gfs%setup_date(vdate)
  call self%gfs%read_fields(state%fields, geom%domain, geom%npz)
elseif (trim(self%pseudo_type) == "geos") then
  call self%geos%setup_date(vdate)
  call self%geos%read_fields(geom, state%fields)
else
  call abor1_ftn("fv3jedi_pseudo_mod: pseudo_model, model choice must be geos or gfs")
endif

end subroutine step

! --------------------------------------------------------------------------------------------------

subroutine finalize(self, state)

implicit none
class(pseudo_model) :: self
type(fv3jedi_state) :: state

end subroutine finalize

! --------------------------------------------------------------------------------------------------

end module fv3jedi_pseudo_mod
