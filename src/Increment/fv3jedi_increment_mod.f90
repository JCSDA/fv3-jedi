! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_increment_mod

use iso_c_binding
use config_mod
use datetime_mod
use variables_mod, only: oops_vars

use random_mod
use fckit_mpi_module
use unstructured_grid_mod

use fv3jedi_field_mod,           only: fv3jedi_field, get_field, fields_rms, fields_print, fields_gpnorm
use fv3jedi_constants_mod,       only: rad2deg, constoz, cp, alhl, rgas
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_interpolation_mod,   only: bilinear_bump_interp
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_state_utils_mod,     only: fv3jedi_state
use fv3jedi_getvalues_mod,       only: getvalues_tl, getvalues_ad

use wind_vt_mod, only: d2a

use mpp_domains_mod, only: mpp_global_sum, bitwise_efp_sum, center, east, north, center

implicit none
private
public :: fv3jedi_increment, create, delete, zeros, random, copy, &
          self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          dot_prod, diff_incr, &
          read_file, write_file, gpnorm, rms, &
          change_resol, getvalues_tl, getvalues_ad, &
          ug_coord, increment_to_ug, increment_from_ug, dirac, jnormgrad, &
          increment_print

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(oops_vars),         intent(in)    :: vars

integer :: var, vcount

! Total fields
! ------------
self%nf = vars%nv

! Allocate fields structure
! -------------------------
allocate(self%fields(self%nf))

! Loop through and allocate main increment fields
! -----------------------------------------------
vcount = 0
do var = 1, vars%nv
  select case (trim(vars%fldnames(var)))
    case("u","ud")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_d_grid_u_wind', &
           fv3jedi_name = 'ud', units = 'm s-1', staggerloc = north, arraypointer = self%ud)
    case("v","vd")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_d_grid_v_wind', &
           fv3jedi_name = 'vd', units = 'm s-1', staggerloc = east, arraypointer = self%vd)
    case("ua")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_eastward_wind', &
           fv3jedi_name = 'ua', units = 'm s-1', staggerloc = center, arraypointer = self%ua)
    case("va")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_northward_wind', &
           fv3jedi_name = 'va', units = 'm s-1', staggerloc = center, arraypointer = self%va)
    case("t","T")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_air_temperature', &
           fv3jedi_name = 't', units = 'K', staggerloc = center, arraypointer = self%t)
    case("ps")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
           short_name = vars%fldnames(var), long_name = 'increment_of_surface_pressure', &
           fv3jedi_name = 'ps', units = 'Pa', staggerloc = center, arraypointer = self%ps)
    case("q","sphum")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_specific_humidity', &
           fv3jedi_name = 'q', units = 'kg kg-1', staggerloc = center, arraypointer = self%q)
    case("qi","ice_wat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_cloud_liquid_ice', &
           fv3jedi_name = 'qi', units = 'kg kg-1', staggerloc = center, arraypointer = self%qi)
    case("ql","liq_wat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_cloud_liquid_ice_water', &
           fv3jedi_name = 'ql', units = 'kg kg-1', staggerloc = center, arraypointer = self%ql)
    case("qs","snowwat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_snow_water', &
           fv3jedi_name = 'qs', units = 'kg kg-1', staggerloc = center, arraypointer = self%qs)
    case("qr","rainwat")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_rain_water', &
           fv3jedi_name = 'qr', units = 'kg kg-1', staggerloc = center, arraypointer = self%qr)
    case("o3","o3mr")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_ozone_mass_mixing_ratio', &
           fv3jedi_name = 'o3', units = 'kg kg-1', staggerloc = center, arraypointer = self%o3)
    case("psi")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_stream_function', &
           fv3jedi_name = 'psi', units = 'm+2 s', staggerloc = center, arraypointer = self%psi)
    case("chi")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_velocity_potential', &
           fv3jedi_name = 'chi', units = 'm+2 s', staggerloc = center, arraypointer = self%chi)
    case("tv")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_virtual_temperature', &
           fv3jedi_name = 'tv', units = 'K', staggerloc = center, arraypointer = self%tv)
    case("rh")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_relative_humidity', &
           fv3jedi_name = 'rh', units = '1', staggerloc = center, arraypointer = self%rh)
    case("delp","DELP")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_pressure_thickness', &
           fv3jedi_name = 'delp', units = 'Pa', staggerloc = center, arraypointer = self%delp)
    case("w","W")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_vertical_wind', &
           fv3jedi_name = 'w', units = 'm s-1', staggerloc = center, arraypointer = self%w)
    case("delz","DZ")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_layer_thickness', &
           fv3jedi_name = 'delz', units = 'm', staggerloc = center, arraypointer = self%delz)
    !Aerosols
    case("du001","DU001","dust1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_dust_mixing_ratio_bin1', &
           fv3jedi_name = 'du001', units = 'kg kg-1', staggerloc = center, arraypointer = self%du001, &
           tracer = .true.)
    case("du002","DU002","dust2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_dust_mixing_ratio_bin2', &
           fv3jedi_name = 'du002', units = 'kg kg-1', staggerloc = center, arraypointer = self%du002, &
           tracer = .true.)
    case("du003","DU003","dust3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_dust_mixing_ratio_bin3', &
           fv3jedi_name = 'du003', units = 'kg kg-1', staggerloc = center, arraypointer = self%du003, &
           tracer = .true.)
    case("du004","DU004","dust4")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_dust_mixing_ratio_bin4', &
           fv3jedi_name = 'du004', units = 'kg kg-1', staggerloc = center, arraypointer = self%du004, &
           tracer = .true.)
    case("du005","DU005","dust5")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_dust_mixing_ratio_bin5', &
           fv3jedi_name = 'du005', units = 'kg kg-1', staggerloc = center, arraypointer = self%du005, &
           tracer = .true.)
    case("ss001","SS001","seas1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sea_salt_mixing_ratio_bin1', &
           fv3jedi_name = 'ss001', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss001)
    case("ss002", "SS002","seas2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sea_salt_mixing_ratio_bin2', &
           fv3jedi_name = 'ss002', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss002, &
           tracer = .true.)
    case("ss003","SS003","seas3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sea_salt_mixing_ratio_bin3', &
           fv3jedi_name = 'ss003', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss003, &
           tracer = .true.)
    case("ss004","SS004","seas4")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sea_salt_mixing_ratio_bin4', &
           fv3jedi_name = 'ss004', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss004, &
           tracer = .true.)
    case("ss005","SS005","seas5")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sea_salt_mixing_ratio_bin5', &
           fv3jedi_name = 'ss005', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss005, &
           tracer = .true.)
    case("bcphobic","BCPHOBIC","bc1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_hydrophobic_black_carbon', &
           fv3jedi_name = 'bcphobic', units = 'kg kg-1', staggerloc = center, arraypointer = self%bcphobic, &
           tracer = .true.)
    case("bcphilic","BCPHILIC","bc2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_hydrophilic_black_carbon', &
           fv3jedi_name = 'bcphilic', units = 'kg kg-1', staggerloc = center, arraypointer = self%bcphilic, &
           tracer = .true.)
    case("ocphobic","OCPHOBIC","oc1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_hydrophobic_organic_carbon', &
           fv3jedi_name = 'ocphobic', units = 'kg kg-1', staggerloc = center, arraypointer = self%ocphobic, &
           tracer = .true.)
    case("ocphilic","OCPHILIC","oc2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_hydrophilic_organic_carbon', &
           fv3jedi_name = 'ocphilic', units = 'kg kg-1', staggerloc = center, arraypointer = self%ocphilic, &
           tracer = .true.)
    case("no3an1","NO3AN1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_nitrate_size_bin_1', &
           fv3jedi_name = 'no3an1', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an1, &
           tracer = .true.)
    case("no3an2","NO3AN2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_nitrate_size_bin_2', &
           fv3jedi_name = 'no3an2', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an2, &
           tracer = .true.)
    case("no3an3","NO3AN3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_nitrate_size_bin_3', &
           fv3jedi_name = 'no3an3', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an3, &
           tracer = .true.)
    case("so4","SO4","sulf")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%fldnames(var), long_name = 'increment_of_sulphate_aerosol', &
           fv3jedi_name = 'so4', units = 'kg kg-1', staggerloc = center, arraypointer = self%so4, &
           tracer = .true.)
    !Not found
    case default
      call abor1_ftn("fv3jedi_increment_mod.create: unknown variable "//trim(vars%fldnames(var)))

  end select

enddo

if (vcount .ne. self%nf) &
call abor1_ftn("fv3jedi_increment_mod.create: vcount does not equal self%nf")

self%hydrostatic = .true.
if (associated(self%w) .and. associated(self%delz)) self%hydrostatic = .false.

! Initialize all domain arrays to zero
call zeros(self)

! For convenience
self%isc    = geom%isc
self%iec    = geom%iec
self%jsc    = geom%jsc
self%jec    = geom%jec
self%isd    = geom%isd
self%ied    = geom%ied
self%jsd    = geom%jsd
self%jed    = geom%jed
self%npx    = geom%npx
self%npy    = geom%npy
self%npz    = geom%npz
self%ntile  = geom%ntile
self%ntiles = geom%ntiles

! Pointer to fv3jedi communicator
self%f_comm = fckit_mpi_comm()

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

!Deallocate fields
do var = 1, self%nf
  call self%fields(var)%deallocate_field()
enddo
deallocate(self%fields)

!Nullify pointers
if (associated(self%ud  )) nullify(self%ud  )
if (associated(self%vd  )) nullify(self%vd  )
if (associated(self%ua  )) nullify(self%ua  )
if (associated(self%va  )) nullify(self%va  )
if (associated(self%t   )) nullify(self%t   )
if (associated(self%ps  )) nullify(self%ps  )
if (associated(self%delp)) nullify(self%delp)
if (associated(self%w   )) nullify(self%w   )
if (associated(self%delz)) nullify(self%delz)
if (associated(self%q   )) nullify(self%q   )
if (associated(self%qi  )) nullify(self%qi  )
if (associated(self%ql  )) nullify(self%ql  )
if (associated(self%qr  )) nullify(self%qr  )
if (associated(self%qs  )) nullify(self%qs  )
if (associated(self%o3  )) nullify(self%o3  )
if (associated(self%psi )) nullify(self%psi )
if (associated(self%chi )) nullify(self%chi )
if (associated(self%tv  )) nullify(self%tv  )
if (associated(self%rh  )) nullify(self%rh  )
!Aerosols
if (associated(self%du001    )) nullify(self%du001    )
if (associated(self%du002    )) nullify(self%du002    )
if (associated(self%du003    )) nullify(self%du003    )
if (associated(self%du004    )) nullify(self%du004    )
if (associated(self%du005    )) nullify(self%du005    )
if (associated(self%ss001    )) nullify(self%ss001    )
if (associated(self%ss002    )) nullify(self%ss002    )
if (associated(self%ss003    )) nullify(self%ss003    )
if (associated(self%ss004    )) nullify(self%ss004    )
if (associated(self%ss005    )) nullify(self%ss005    )
if (associated(self%no3an1   )) nullify(self%no3an1   )
if (associated(self%no3an2   )) nullify(self%no3an2   )
if (associated(self%no3an3   )) nullify(self%no3an3   )
if (associated(self%so4      )) nullify(self%so4      )
if (associated(self%bcphobic )) nullify(self%bcphobic )
if (associated(self%bcphilic )) nullify(self%bcphilic )
if (associated(self%ocphobic )) nullify(self%ocphobic )
if (associated(self%ocphilic )) nullify(self%ocphilic )

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 0.0_kind_real
enddo

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine ones(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var

do var = 1,self%nf
  self%fields(var)%array = 1.0_kind_real
enddo

end subroutine ones

! ------------------------------------------------------------------------------

subroutine random(self)

implicit none
type(fv3jedi_increment), intent(inout) :: self

integer :: var
integer, parameter :: rseed = 7

do var = 1,self%nf
  call normal_distribution(self%fields(var)%array, 0.0_kind_real, 1.0_kind_real, rseed)
enddo

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

logical :: found
integer :: self_var, rhs_var

self%isc            = rhs%isc
self%iec            = rhs%iec
self%jsc            = rhs%jsc
self%jec            = rhs%jec
self%isd            = rhs%isd
self%ied            = rhs%ied
self%jsd            = rhs%jsd
self%jed            = rhs%jed
self%npx            = rhs%npx
self%npy            = rhs%npy
self%npz            = rhs%npz
self%ntiles         = rhs%ntiles
self%ntile          = rhs%ntile
self%hydrostatic    = rhs%hydrostatic
self%tladphystrj    = rhs%tladphystrj
self%calendar_type  = rhs%calendar_type
self%date_init      = rhs%date_init

!Copy the individual fields
if (.not.allocated(self%fields)) then

  !Direct copy of one increment to another
  self%nf = rhs%nf
  allocate(self%fields(self%nf))

  do self_var = 1, self%nf

    !Field copy
    self%fields(self_var) = rhs%fields(self_var)

  enddo

else

  !Increment copy, potentialy with differnt fields
  do self_var = 1, self%nf
    found = .false.
    do rhs_var = 1, rhs%nf
      if (trim(self%fields(self_var)%fv3jedi_name) == trim(rhs%fields(rhs_var)%fv3jedi_name)) then
        self%fields(self_var) = rhs%fields(rhs_var)
        found = .true.
        exit
      endif
    enddo
    if (.not.found) call abor1_ftn("fv3jedi_increment_mod.copy: Field "//&
                    trim(self%fields(self_var)%fv3jedi_name)//" not found in increment being copied from." )

  enddo

endif

!Set pointers
do self_var = 1, self%nf
  select case (trim(self%fields(self_var)%fv3jedi_name))
    case("u","ud")
      call self%fields(self_var)%array_pointer(self%ud)
    case("v","vd")
      call self%fields(self_var)%array_pointer(self%vd)
    case("ua")
      call self%fields(self_var)%array_pointer(self%ua)
    case("va")
      call self%fields(self_var)%array_pointer(self%va)
    case("t","T")
      call self%fields(self_var)%array_pointer(self%t)
    case("delp","DELP")
      call self%fields(self_var)%array_pointer(self%delp)
    case("ps")
      call self%fields(self_var)%array_pointer(self%ps)
    case("w","W")
      call self%fields(self_var)%array_pointer(self%w)
    case("delz","DZ")
      call self%fields(self_var)%array_pointer(self%delz)
    case("q","sphum")
      call self%fields(self_var)%array_pointer(self%q)
    case("qi","ice_wat")
      call self%fields(self_var)%array_pointer(self%qi)
    case("ql","liq_wat")
      call self%fields(self_var)%array_pointer(self%ql)
    case("qs","snowwat")
      call self%fields(self_var)%array_pointer(self%qs)
    case("qr","rainwat")
      call self%fields(self_var)%array_pointer(self%qr)
    case("o3","o3mr")
      call self%fields(self_var)%array_pointer(self%o3)
    case("psi")
      call self%fields(self_var)%array_pointer(self%psi)
    case("chi")
      call self%fields(self_var)%array_pointer(self%chi)
    case("tv")
      call self%fields(self_var)%array_pointer(self%tv)
    case("rh")
      call self%fields(self_var)%array_pointer(self%rh)
    !Aerosols
    case("du001","DU001","dust1")
      call self%fields(self_var)%array_pointer(self%du001)
    case("du002","DU002","dust2")
      call self%fields(self_var)%array_pointer(self%du002)
    case("du003","DU003","dust3")
      call self%fields(self_var)%array_pointer(self%du003)
    case("du004","DU004","dust4")
      call self%fields(self_var)%array_pointer(self%du004)
    case("du005","DU005","dust5")
      call self%fields(self_var)%array_pointer(self%du005)
    case("ss001","SS001","seas1")
      call self%fields(self_var)%array_pointer(self%ss001)
    case("ss002","SS002","seas2")
      call self%fields(self_var)%array_pointer(self%ss002)
    case("ss003","SS003","seas3")
      call self%fields(self_var)%array_pointer(self%ss003)
    case("ss004","SS004","seas4")
      call self%fields(self_var)%array_pointer(self%ss004)
    case("ss005","SS005","seas5")
      call self%fields(self_var)%array_pointer(self%ss005)
    case("no3an1","NO3AN1")
      call self%fields(self_var)%array_pointer(self%no3an1)
    case("no3an2","NO3AN2")
      call self%fields(self_var)%array_pointer(self%no3an2)
    case("no3an3","NO3AN3")
      call self%fields(self_var)%array_pointer(self%no3an3)
    case("so4","SO4","sulf")
      call self%fields(self_var)%array_pointer(self%so4)
    case("bcphobic","BCPHOBIC","bc1")
      call self%fields(self_var)%array_pointer(self%bcphobic)
    case("bcphilic","BCPHILIC","bc2")
      call self%fields(self_var)%array_pointer(self%bcphilic)
    case("ocphobic","OCPHOBIC","oc1")
      call self%fields(self_var)%array_pointer(self%ocphobic)
    case("ocphilic","OCPHILIC","oc2")
      call self%fields(self_var)%array_pointer(self%ocphilic)
    case default
      !Not found
      call abor1_ftn("fv3jedi_increment_mod.copy: unknown variable "//trim(self%fields(self_var)%fv3jedi_name))
  end select
enddo

end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self,rhs)

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + rhs%fields(var)%array
enddo

end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self,rhs)

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array * rhs%fields(var)%array
enddo

end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self,rhs)

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array - rhs%fields(var)%array
enddo

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz

integer :: var

do var = 1,self%nf
  self%fields(var)%array = zz * self%fields(var)%array
enddo

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy_inc(self,zz,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var

call checksame(self,rhs)

do var = 1,self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine axpy_inc

! ------------------------------------------------------------------------------

subroutine axpy_state(self,zz,rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
real(kind=kind_real),    intent(in)    :: zz
type(fv3jedi_state),     intent(in)    :: rhs

integer :: var
type(fv3jedi_field), pointer :: field_pointer

do var = 1,self%nf

  !Get pointer to equivalent state
  call get_field(rhs%nf,rhs%fields,self%fields(var)%fv3jedi_name,field_pointer)

  !inc = inc + zz state
  self%fields(var)%array = self%fields(var)%array + zz * field_pointer%array

  !Clear pointer
  nullify(field_pointer)

enddo

end subroutine axpy_state

! ------------------------------------------------------------------------------

subroutine dot_prod(self,other,zprod)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(fv3jedi_increment), intent(in)    :: other
real(kind=kind_real),    intent(inout) :: zprod

real(kind=kind_real) :: zp
integer :: i,j,k
integer :: var

call checksame(self,other)

zp=0.0_kind_real
do var = 1,self%nf
  do k = 1,self%fields(var)%npz
    do j = self%jsc,self%jec
      do i = self%isc,self%iec
        zp = zp + self%fields(var)%array(i,j,k) * other%fields(var)%array(i,j,k)
      enddo
    enddo
  enddo
enddo

!Get global dot product
call self%f_comm%allreduce(zp,zprod,fckit_mpi_sum())

!For debugging print result:
if (self%f_comm%rank() == 0) print*, "Dot product test result: ", zprod

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine diff_incr(self,x1,x2,geom)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_state),     intent(in)    :: x1
type(fv3jedi_state),     intent(in)    :: x2
type(fv3jedi_geom),      intent(inout) :: geom

integer :: var
type(fv3jedi_field), pointer :: x1p, x2p
real(kind=kind_real), allocatable :: x1_ua(:,:,:), x1_va(:,:,:)
real(kind=kind_real), allocatable :: x2_ua(:,:,:), x2_va(:,:,:)
real(kind=kind_real), allocatable :: x1_ps(:,:,:), x2_ps(:,:,:)

call zeros(self)

!D-Grid to A-Grid (if needed)
if (associated(self%ua) .and. associated(self%va)) then
  allocate(x1_ua(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  allocate(x1_va(x1%isc:x1%iec,x1%jsc:x1%jec,1:x1%npz))
  if (associated(x1%ud) .and. associated(x1%vd)) then
    if (.not.associated(x1%ua) .and. .not. associated(x1%va)) then
      call d2a(geom, x1%ud, x1%vd, x1_ua, x1_va)
    else
      x1_ua = x1%ua
      x1_va = x1%va
    endif
  else
    x1_ua = 0.0_kind_real
    x1_va = 0.0_kind_real
  endif
  allocate(x2_ua(x2%isc:x2%iec,x2%jsc:x2%jec,1:x2%npz))
  allocate(x2_va(x2%isc:x2%iec,x2%jsc:x2%jec,1:x2%npz))
  if (associated(x2%ud) .and. associated(x2%vd)) then
    if (.not.associated(x2%ua) .and. .not. associated(x2%va)) then
      call d2a(geom, x2%ud, x2%vd, x2_ua, x2_va)
    else
      x2_ua = x2%ua
      x2_va = x2%va
    endif
  else
    x2_ua = 0.0_kind_real
    x2_va = 0.0_kind_real
  endif
endif

!delp to ps
if (associated(self%ps)) then

  allocate(x1_ps(x1%isc:x1%iec,x1%jsc:x1%jec,1))
  allocate(x2_ps(x2%isc:x2%iec,x2%jsc:x2%jec,1))

  if (associated(x1%delp)) then
    x1_ps(:,:,1) = sum(x1%delp,3)
  elseif (associated(x1%ps)) then
    x1_ps = x1%ps
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x1")
  endif

  if (associated(x2%delp)) then
    x2_ps(:,:,1) = sum(x2%delp,3)
  elseif (associated(x2%ps)) then
    x2_ps = x2%ps
  else
    call abor1_ftn("fv3jedi_increment_mod.diff_incr: problem getting ps from state x2")
  endif

endif

do var = 1,self%nf

  !A-Grid winds can be a special case
  if (self%fields(var)%fv3jedi_name == 'ua') then

    self%ua = x1_ua - x2_ua

  elseif (self%fields(var)%fv3jedi_name == 'va') then

    self%va = x1_va - x2_va

  !Ps can be a special case
  elseif (self%fields(var)%fv3jedi_name == 'ps') then

    self%ps = x1_ps - x2_ps

  else

    !Get pointer to states
    call get_field(x1%nf,x1%fields,self%fields(var)%fv3jedi_name,x1p)
    call get_field(x2%nf,x2%fields,self%fields(var)%fv3jedi_name,x2p)

    !inc = state - state
    self%fields(var)%array = x1p%array - x2p%array

    !Nullify pointers
    nullify(x1p,x2p)

  endif

enddo

if (allocated(x1_ua)) deallocate(x1_ua)
if (allocated(x1_va)) deallocate(x1_va)
if (allocated(x2_ua)) deallocate(x2_ua)
if (allocated(x2_va)) deallocate(x2_va)
if (allocated(x1_ps)) deallocate(x1_ps)
if (allocated(x2_ps)) deallocate(x2_ps)

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(self,geom,rhs,geom_rhs)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_increment), intent(in)    :: rhs
type(fv3jedi_geom),      intent(in)    :: geom_rhs

integer :: check

check = (rhs%iec-rhs%isc+1) - (self%iec-self%isc+1)

call checksame(self,rhs)

if (check==0) then

  !Resolution is the same so copy
  call copy(self, rhs)

else

  !Interpolate from rhs to self resolution
  call bilinear_bump_interp(self%nf, geom_rhs, rhs%fields, geom, self%fields)

endif

end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(geom, self, c_conf, vdate)

  implicit none
  type(fv3jedi_geom),      intent(inout) :: geom     !< Geometry
  type(fv3jedi_increment), intent(inout) :: self     !< Increment
  type(c_ptr),             intent(in)    :: c_conf   !< Configuration
  type(datetime),          intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype
  character(len=255) :: filename

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%read_meta(geom, vdate, self%calendar_type, self%date_init)
    call gfs%read_fields(geom, self%fields)
  elseif (trim(filetype) == 'geos') then
    filename = config_get_string(c_conf,len(filename),"filename")
    call geos%create(geom, 'read', filename)
    call geos%read_time(vdate)
    call geos%read_fields(geom, self%fields)
    call geos%delete()
  else
     call abor1_ftn("fv3jedi_increment_mod.read: restart type not supported")
  endif

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, self, c_conf, vdate)

  implicit none

  type(fv3jedi_geom),      intent(inout) :: geom     !< Geometry
  type(fv3jedi_increment), intent(in)    :: self     !< Increment
  type(c_ptr),             intent(in)    :: c_conf   !< Configuration
  type(datetime),          intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%write_all(geom, self%fields, vdate, self%calendar_type, self%date_init)
  elseif (trim(filetype) == 'geos') then
    call geos%create(geom, 'write')
    call geos%write_all(geom, self%fields, c_conf, vdate)
    call geos%delete()
  else
     call abor1_ftn("fv3jedi_increment_mod.write: restart type not supported")
  endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine increment_print(self)

implicit none
type(fv3jedi_increment), intent(in) :: self

call fields_print(self%nf, self%fields, "Increment")

end subroutine increment_print

! ------------------------------------------------------------------------------

subroutine gpnorm(self, nf, pstat)

implicit none
type(fv3jedi_increment), intent(in)    :: self
integer,                 intent(in)    :: nf
real(kind=kind_real),    intent(inout) :: pstat(3, nf)

if (nf .ne. self%nf) then
  call abor1_ftn("fv3jedi_increment.gpnorm: nf passed in does not match expeted nf")
endif

call fields_gpnorm(nf, self%fields, pstat)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_increment), intent(in)  :: self
real(kind=kind_real),    intent(out) :: prms

call fields_rms(self%nf,self%fields,prms)

end subroutine rms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf, geom)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(c_ptr),             intent(in)    :: c_conf
type(fv3jedi_geom),      intent(in)    :: geom

integer :: ndir,idir,ildir,itiledir,var
integer,allocatable :: ixdir(:),iydir(:)
character(len=3) :: idirchar
character(len=32) :: ifdir
logical :: found

! Get Diracs positions
ndir = config_get_int(c_conf,"ndir")
allocate(ixdir(ndir))
allocate(iydir(ndir))

do idir=1,ndir
   write(idirchar,'(i3)') idir
   ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
   iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
end do
ildir = config_get_int(c_conf,"ildir")
ifdir = config_get_string(c_conf,len(ifdir),"ifdir")
itiledir = config_get_int(c_conf,"itiledir")

! Check
if (ndir<1) call abor1_ftn("fv3jedi_increment_mod.dirac: non-positive ndir")
if (any(ixdir<1).or.any(ixdir>self%npx)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid ixdir")
endif
if (any(iydir<1).or.any(iydir>geom%size_cubic_grid)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid iydir")
endif
if ((ildir<1).or.(ildir>self%npz)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac invalid ildir")
endif
if ((itiledir<1).or.(itiledir>6)) then
   call abor1_ftn("fv3jedi_increment_mod.dirac: invalid itiledir")
endif

! Setup Diracs
call zeros(self)

! only u, v, T, ps and tracers allowed
do idir=1,ndir

   ! is specified grid point, tile number on this processor
   if (geom%ntile == itiledir .and. &
       ixdir(idir) >= self%isc .and. ixdir(idir) <= self%iec .and. &
       iydir(idir) >= self%jsc .and. iydir(idir) <= self%jec) then
       ! If so, perturb desired increment and level
       found = .false.
       do var = 1,self%nf
         if (trim(self%fields(var)%fv3jedi_name) == trim(ifdir)) then
           found = .true.
           self%fields(var)%array(ixdir(idir),iydir(idir),ildir) = 1.0
         endif
       enddo
       if (.not.found) call abor1_ftn("fv3jedi_increment_mod.dirac:: field "&
                                       //trim(ifdir)//" not found")
   endif
end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine ug_size(self, ug)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug

! Set number of grids
if (ug%colocated==1) then
   ! Colocatd
   ug%ngrid = 1
else
   ! Not colocated
   ug%ngrid = 1
   call abor1_ftn("fv3jedi_increment_mod.ug_size: Uncolocated grids not coded yet")
end if

! Allocate grid instances
if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

! Set local number of points
ug%grid(1)%nmga = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1)

! Set number of levels
ug%grid(1)%nl0 = self%npz

! Set number of variables
ug%grid(1)%nv = self%nf

! Set number of timeslots
ug%grid(1)%nts = 1

end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine ug_coord(self, ug, geom)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug
type(fv3jedi_geom),      intent(in)    :: geom

integer :: imga,jx,jy,jl
real(kind=kind_real) :: sigmaup,sigmadn

! Define size
call ug_size(self, ug)

! Allocate unstructured grid coordinates
call allocate_unstructured_grid_coord(ug)

if (ug%colocated==1) then
  imga = 0
  do jy=self%jsc,self%jec
    do jx=self%isc,self%iec
      imga = imga+1
      ug%grid(1)%lon(imga) = rad2deg*geom%grid_lon(jx,jy)
      ug%grid(1)%lat(imga) = rad2deg*geom%grid_lat(jx,jy)
      ug%grid(1)%area(imga) = geom%area(jx,jy)
      do jl=1,self%npz
        sigmaup = geom%ak(jl+1)/101300.0+geom%bk(jl+1) ! si are now sigmas
        sigmadn = geom%ak(jl  )/101300.0+geom%bk(jl  )
        ug%grid(1)%vunit(imga,jl) = 0.5*(sigmaup+sigmadn) ! 'fake' sigma coordinates
        ug%grid(1)%lmask(imga,jl) = .true.
      enddo
    enddo
  enddo
endif

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine increment_to_ug(self, ug, its)

implicit none
type(fv3jedi_increment), intent(in)    :: self
type(unstructured_grid), intent(inout) :: ug
integer,                 intent(in)    :: its

integer :: var,imga,jx,jy,jl

! Define size
call ug_size(self, ug)

! Allocate unstructured grid increment
call allocate_unstructured_grid_field(ug)

! Copy increment

ug%grid(1)%fld = 0.0_kind_real

if (ug%colocated==1) then
  do var = 1,self%nf
    imga = 0
    do jy=self%jsc,self%jec
      do jx=self%isc,self%iec
        imga = imga+1
        do jl=1,self%fields(var)%npz
          ug%grid(1)%fld(imga,jl,var,its) = self%fields(var)%array(jx,jy,jl)
        enddo
      enddo
    enddo
  enddo
endif

end subroutine increment_to_ug

! -----------------------------------------------------------------------------

subroutine increment_from_ug(self, ug, its)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(unstructured_grid), intent(in)    :: ug
integer,                 intent(in)    :: its

integer :: imga,jx,jy,jl,var

! Copy increment

if (ug%colocated==1) then
  do var = 1,self%nf
    imga = 0
    do jy=self%jsc,self%jec
      do jx=self%isc,self%iec
        imga = imga+1
        do jl=1,self%fields(var)%npz
          self%fields(var)%array(jx,jy,jl) = ug%grid(1)%fld(imga,jl,var,its)
        enddo
      enddo
    enddo
  enddo
endif

end subroutine increment_from_ug

! ------------------------------------------------------------------------------

subroutine jnormgrad(self,geom,ref,c_conf)

implicit none
type(fv3jedi_increment), intent(inout) :: self
type(fv3jedi_geom),      intent(in)    :: geom
type(fv3jedi_state),     intent(in)    :: ref !To linearize around if nl
type(c_ptr),             intent(in)    :: c_conf

integer :: i,j,k
integer :: isc,iec,jsc,jec,npz
real(kind=kind_real), allocatable :: cellweight(:,:,:), ref_ps(:,:)

real(kind=kind_real) :: global_area

real(kind=kind_real) :: Ufac
real(kind=kind_real) :: Tfac, Tref
real(kind=kind_real) :: qfac, qeps
real(kind=kind_real) :: pfac, pref

!Code to compute a vector norm for an increment, e.g. the energy norm for FSOI

isc = self%isc
iec = self%iec
jsc = self%jsc
jec = self%jec
npz = self%npz

! Constants
! ---------
tref = config_get_real(c_conf,"Tref")
qeps = config_get_real(c_conf,"qepsilon")
pref = config_get_real(c_conf,"pref")

Ufac = 0.5_kind_real
Tfac = 0.5_kind_real*cp/tref
qfac = 0.5_kind_real*qeps*alhl*alhl/(cp*tref)
pfac = 0.5_kind_real*Rgas*tref/pref**2

! Compute grid weighting based on volume
! --------------------------------------

global_area = mpp_global_sum(geom%domain, geom%area, flags=bitwise_efp_sum)

allocate(ref_ps(isc:iec,jsc:jec))
ref_ps = sum(ref%delp,3)

allocate(cellweight(isc:iec,jsc:jec,1:npz))
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      cellweight(i,j,k) = (ref%delp(i,j,k)/ref_ps(i,j)) * geom%area(i,j)/global_area
    enddo
  enddo
enddo

!ua
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%ua(i,j,k) = Ufac * 2.0_kind_real * ref%ua(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!va
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%va(i,j,k) = Ufac * 2.0_kind_real * ref%va(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!T
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%t(i,j,k) = Tfac * 2.0_kind_real * ref%t(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!q
do k = 1, npz
  do j = jsc,jec
    do i = isc,iec
      self%q(i,j,k) = qfac * 2.0_kind_real * ref%q(i,j,k) * cellweight(i,j,k)
    enddo
  enddo
enddo

!ps
if (associated(self%ps)) then
  do j = jsc,jec
    do i = isc,iec
      self%ps(i,j,1) = pfac * 2.0_kind_real * ref_ps (i,j) * cellweight(i,j,npz) &
                                          / (ref%delp(i,j,npz)/ref_ps(i,j))
    enddo
  enddo
else
  call abor1_ftn("fv3jedi_increment_mod.jnormgrad: Ps must be in the increment")
endif

deallocate(cellweight)
deallocate(ref_ps)

end subroutine jnormgrad

! ------------------------------------------------------------------------------

subroutine checksame(self,other)

implicit none
type(fv3jedi_increment), intent(in) :: self
type(fv3jedi_increment), intent(in) :: other

integer :: var

if (self%nf .ne. other%nf) then
  call abor1_ftn("fv3jedi_increment_mod.checksame: Different number of fields")
endif

do var = 1,self%nf
  if (self%fields(var)%fv3jedi_name .ne. other%fields(var)%fv3jedi_name) then
      call abor1_ftn("fv3jedi_increment_mod.checksame: field "//trim(self%fields(var)%fv3jedi_name)//&
                     " not in the equivalent position in the right hand side")
  endif
enddo

end subroutine checksame

! ------------------------------------------------------------------------------

end module fv3jedi_increment_mod
