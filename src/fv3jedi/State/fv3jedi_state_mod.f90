! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_state_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use fckit_mpi_module
use oops_variables_mod

use fv3jedi_field_mod
use fv3jedi_constants_mod,       only: rad2deg, constoz
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_interpolation_mod,   only: field2field_interp
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use fv3jedi_state_utils_mod,     only: fv3jedi_state
use fv3jedi_getvalues_mod,       only: getvalues

use wind_vt_mod, only: a2d

use mpp_domains_mod,             only: east, north, center

implicit none
private
public :: fv3jedi_state, create, delete, zeros, copy, axpy, add_incr, &
          read_file, write_file, gpnorm, rms, &
          change_resol, getvalues, analytic_IC, state_print

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_state),  intent(inout) :: self
type(fv3jedi_geom),   intent(in)    :: geom
type(oops_variables), intent(in)    :: vars

integer :: var, vcount, nlev

! Total fields
! ------------
self%nf = vars%nvars()

! Allocate fields structure
! -------------------------
allocate(self%fields(self%nf))

! Loop through and allocate main state fields
! -------------------------------------------
vcount = 0
do var = 1, vars%nvars()
   select case (trim(vars%variable(var)))

     ! In the below the variable names are designed to cover the potential names encountered in
     ! GEOS and GFS restart and history files. E.g. U for GEOS, u for GFS and ud for history.
     ! Within fv3-jedi variables can be accessed using the fv3jedi_name.

     case("ud","u","U")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'eastward_wind_on_native_D-Grid', &
            fv3jedi_name = 'ud', units = 'm s-1', staggerloc = north, arraypointer = self%ud )
     case("vd","v","V")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'northward_wind_on_native_D-Grid', &
            fv3jedi_name = 'vd', units = 'm s-1', staggerloc = east, arraypointer = self%vd )
     case("ua")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'eastward_wind', &
            fv3jedi_name = 'ua', units = 'm s-1', staggerloc = center, arraypointer = self%ua)
     case("va")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'northward_wind', &
            fv3jedi_name = 'va', units = 'm s-1', staggerloc = center, arraypointer = self%va )
     case("t","T")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'air_temperature', &
            fv3jedi_name = 't', units = 'K', staggerloc = center, arraypointer = self%t)
     case("tv","TV")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'virtual_temperature', &
            fv3jedi_name = 'tv', units = 'K', staggerloc = center, arraypointer = self%tv)
     case("pt","PT")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'potential_temperature', &
            fv3jedi_name = 'pt', units = 'K', staggerloc = center, arraypointer = self%pt)
     case("delp","DELP")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'pressure_thickness', &
            fv3jedi_name = 'delp', units = 'Pa', staggerloc = center, arraypointer = self%delp)
     case("pkz","PKZ")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'pressure_to_kappa', &
            fv3jedi_name = 'pkz', units = 'Pa', staggerloc = center, arraypointer = self%pkz)
     case("pe","PE")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz+1, &
            short_name = vars%variable(var), long_name = 'pressure', &
            fv3jedi_name = 'pe', units = 'Pa', staggerloc = center, arraypointer = self%pe)
     case("ps")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_pressure', &
            fv3jedi_name = 'ps', units = 'Pa', staggerloc = center, arraypointer = self%ps)
     case("q","sphum","Q")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'specific_humidity', &
            fv3jedi_name = 'q', units = 'kg kg-1', staggerloc = center, arraypointer = self%q, &
            tracer = .true.)
     case("rh","RH")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'relative_humidity', &
            fv3jedi_name = 'rh', units = '1', staggerloc = center, arraypointer = self%rh, &
            tracer = .true.)
     case("qi","ice_wat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'cloud_liquid_ice', &
            fv3jedi_name = 'qi', units = 'kg kg-1', staggerloc = center, arraypointer = self%qi, &
            tracer = .true.)
     case("ql","liq_wat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'cloud_liquid_ice_water', &
            fv3jedi_name = 'ql', units = 'kg kg-1', staggerloc = center, arraypointer = self%ql, &
            tracer = .true.)
     case("qils","QILS")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'mass_fraction_of_large_scale_cloud_ice_water', &
            fv3jedi_name = 'qils', units = 'kg kg-1', staggerloc = center, arraypointer = self%qils, &
            tracer = .true.)


     case("qlls","QLLS")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'mass_fraction_of_large_scale_cloud_liquid_water', &
            fv3jedi_name = 'qlls', units = 'kg kg-1', staggerloc = center, arraypointer = self%qlls, &
            tracer = .true.)
     case("qicn","QICN")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'mass_fraction_of_convective_cloud_ice_water', &
            fv3jedi_name = 'qicn', units = 'kg kg-1', staggerloc = center, arraypointer = self%qicn, &
            tracer = .true.)
     case("qlcn","QLCN")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'mass_fraction_of_convective_cloud_liquid_water', &
            fv3jedi_name = 'qlcn', units = 'kg kg-1', staggerloc = center, arraypointer = self%qlcn, &
            tracer = .true.)
     case("qicnf","QICNF")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
             short_name = vars%variable(var), long_name = 'fraction_of_large_scale_cloud_ice_water_in_total', &
             fv3jedi_name = 'qicnf', units = '1', staggerloc = center, arraypointer = self%qicnf, &
             tracer = .true.)
     case("qilsf","QILSF")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
             short_name = vars%variable(var), long_name = 'fraction_of_large_scale_cloud_ice_water_in_total', &
             fv3jedi_name = 'qilsf', units = '1', staggerloc = center, arraypointer = self%qilsf, &
             tracer = .true.)
     case("qs","snowwat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'snow_water', &
            fv3jedi_name = 'qs', units = 'kg kg-1', staggerloc = center, arraypointer = self%qs, &
            tracer = .true.)
     case("qr","rainwat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'rain_water', &
            fv3jedi_name = 'qr', units = 'kg kg-1', staggerloc = center, arraypointer = self%qr, &
            tracer = .true.)
     case("graupel")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'graupel', &
            fv3jedi_name = 'graupel', units = 'kg kg-1', staggerloc = center, arraypointer = self%gr, &
            tracer = .true.)
     case("cld_amt")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'cld_amt', &
            fv3jedi_name = 'cld_amt', units = 'kg kg-1', staggerloc = center, arraypointer = self%ca, &
            tracer = .true.)
     case("o3","o3mr")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'ozone_mass_mixing_ratio', &
            fv3jedi_name = 'o3', units = 'kg kg-1', staggerloc = center, arraypointer = self%o3, &
            tracer = .true.)
     case("ox","OX")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'odd_oxygen_mixing_ratio', &
            fv3jedi_name = 'ox', units = 'kg kg-1', staggerloc = center, arraypointer = self%ox, &
            tracer = .true.)
     case("w","W")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'vertical_wind', &
            fv3jedi_name = 'w', units = 'm s-1', staggerloc = center, arraypointer = self%w)
     case("delz","DZ")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'layer_thickness', &
            fv3jedi_name = 'delz', units = 'm', staggerloc = center, arraypointer = self%delz)
     case("phis","PHIS")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_geopotential_height', &
            fv3jedi_name = 'phis', units = 'm', staggerloc = center, arraypointer = self%phis)
     case("psi")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'stream_function', &
            fv3jedi_name = 'psi', units = 'm+2 s', staggerloc = center, arraypointer = self%psi)
     case("chi")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'velocity_potential', &
            fv3jedi_name = 'chi', units = 'm+2 s', staggerloc = center, arraypointer = self%chi)
     case("vort")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'vorticity', &
            fv3jedi_name = 'vort', units = 'm+2 s', staggerloc = center, arraypointer = self%vort)
     case("divg")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'divergence', &
            fv3jedi_name = 'divg', units = 'm+2 s', staggerloc = center, arraypointer = self%divg)
     !CRTM
     case("slmsk")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'slmsk', &
            fv3jedi_name = 'slmsk', units = 'none', staggerloc = center, arraypointer = self%slmsk)
     case("sheleg")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'sheleg', &
            fv3jedi_name = 'sheleg', units = 'none', staggerloc = center, arraypointer = self%sheleg)
     case("tsea","ts")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'tsea', &
            fv3jedi_name = 'tsea', units = 'none', staggerloc = center, arraypointer = self%tsea)
     case("vtype")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'vtype', &
            fv3jedi_name = 'vtype', units = 'none', staggerloc = center, arraypointer = self%vtype, &
            integerfield = .true.)
     case("stype")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'stype', &
            fv3jedi_name = 'stype', units = 'none', staggerloc = center, arraypointer = self%stype, &
            integerfield = .true.)
     case("vfrac")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'vfrac', &
            fv3jedi_name = 'vfrac', units = 'none', staggerloc = center, arraypointer = self%vfrac)
     case("stc","soilt")
       vcount=vcount+1;
       if (trim(vars%variable(var)) == "soilt") then
         nlev = 1 !geos
       else
         nlev = 4 !gfs
       endif
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,nlev, &
            short_name = vars%variable(var), long_name = 'stc', &
            fv3jedi_name = 'stc', units = 'none', staggerloc = center, arraypointer = self%stc)
     case("smc","soilm")
       vcount=vcount+1;
       if (trim(vars%variable(var)) == "soilm") then
         nlev = 1 !geos
       else
         nlev = 4 !gfs
       endif
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,nlev, &
            short_name = vars%variable(var), long_name = 'smc', &
            fv3jedi_name = 'smc', units = 'none', staggerloc = center, arraypointer = self%smc)
     case("snwdph")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'snwdph', &
            fv3jedi_name = 'snwdph', units = 'none', staggerloc = center, arraypointer = self%snwdph)
     case("u_srf","u10m")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'u_srf', &
            fv3jedi_name = 'u_srf', units = 'none', staggerloc = center, arraypointer = self%u_srf)
     case("v_srf","v10m")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'v_srf', &
            fv3jedi_name = 'v_srf', units = 'none', staggerloc = center, arraypointer = self%v_srf)
     case("f10m")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'f10m', &
            fv3jedi_name = 'f10m', units = 'none', staggerloc = center, arraypointer = self%f10m)
     !TL/AD trajectory
     case("qls")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'initial_mass_fraction_of_large_scale_cloud_condensate', &
            fv3jedi_name = 'qls', units = 'kg kg-1', staggerloc = center, arraypointer = self%qls)
     case("qcn")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'initial_mass_fraction_of_convective_cloud_condensate', &
            fv3jedi_name = 'qcn', units = 'kg kg-1', staggerloc = center, arraypointer = self%qcn)
     case("cfcn")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%variable(var), long_name = 'convective_cloud_area_fraction', &
            fv3jedi_name = 'cfcn', units = '1', staggerloc = center, arraypointer = self%cfcn)
     case("frocean")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'fraction_of_ocean', &
            fv3jedi_name = 'frocean', units = '1', staggerloc = center, arraypointer = self%frocean)
     case("frland")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'fraction_of_land', &
            fv3jedi_name = 'frland', units = '1', staggerloc = center, arraypointer = self%frland)
     case("frlandice")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'fraction_of_landice', &
            fv3jedi_name = 'frlandice', units = '1', staggerloc = center, arraypointer = self%frlandice)
     case("frlake")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'fraction_of_lake', &
            fv3jedi_name = 'frlake', units = '1', staggerloc = center, arraypointer = self%frlake)
     case("frseaice")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'fraction_of_ice', &
            fv3jedi_name = 'frseaice', units = '1', staggerloc = center, arraypointer = self%frseaice)
     case("varflt")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'isotropic_variance_of_filtered_topography', &
            fv3jedi_name = 'varflt', units = 'm+2', staggerloc = center, arraypointer = self%varflt)
     case("ustar")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_velocity_scale', &
            fv3jedi_name = 'ustar', units = 'm s-1', staggerloc = center, arraypointer = self%ustar)
     case("bstar")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_bouyancy_scale', &
            fv3jedi_name = 'bstar', units = 'm s-2', staggerloc = center, arraypointer = self%bstar)
     case("zpbl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'planetary_boundary_layer_height', &
            fv3jedi_name = 'zpbl', units = 'm', staggerloc = center, arraypointer = self%zpbl)
     case("cm")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_exchange_coefficient_for_momentum', &
            fv3jedi_name = 'cm', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%cm)
     case("ct")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_exchange_coefficient_for_heat', &
            fv3jedi_name = 'ct', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%ct)
     case("cq")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_exchange_coefficient_for_moisture', &
            fv3jedi_name = 'cq', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%cq)
     case("kcbl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'KCBL_before_moist', &
            fv3jedi_name = 'kcbl', units = '1', staggerloc = center, arraypointer = self%kcbl)
     case("tsm")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'surface_temp_before_moist', &
            fv3jedi_name = 'tsm', units = 'K', staggerloc = center, arraypointer = self%tsm)
     case("khl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'lower_index_where_Kh_greater_than_2', &
            fv3jedi_name = 'khl', units = '1', staggerloc = center, arraypointer = self%khl)
     case("khu")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%variable(var), long_name = 'upper_index_where_Kh_greater_than_2', &
            fv3jedi_name = 'khu', units = '1', staggerloc = center, arraypointer = self%khu)
    !Aerosols
    case("du001","DU001","dust1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_dust001_in_air', &
           fv3jedi_name = 'du001', units = 'kg kg-1', staggerloc = center, arraypointer = self%du001, &
           tracer = .true.)
    case("du002","DU002","dust2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_dust002_in_air', &
           fv3jedi_name = 'du002', units = 'kg kg-1', staggerloc = center, arraypointer = self%du002, &
           tracer = .true.)
    case("du003","DU003","dust3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_dust003_in_air', &
           fv3jedi_name = 'du003', units = 'kg kg-1', staggerloc = center, arraypointer = self%du003, &
           tracer = .true.)
    case("du004","DU004","dust4")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_dust004_in_air', &
           fv3jedi_name = 'du004', units = 'kg kg-1', staggerloc = center, arraypointer = self%du004, &
           tracer = .true.)
    case("du005","DU005","dust5")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_dust005_in_air', &
           fv3jedi_name = 'du005', units = 'kg kg-1', staggerloc = center, arraypointer = self%du005, &
           tracer = .true.)
    case("ss001","SS001","seas1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sea_salt001_in_air', &
           fv3jedi_name = 'ss001', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss001)
    case("ss002","SS002","seas2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sea_salt002_in_air', &
           fv3jedi_name = 'ss002', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss002, &
           tracer = .true.)
    case("ss003","SS003","seas3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sea_salt003_in_air', &
           fv3jedi_name = 'ss003', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss003, &
           tracer = .true.)
    case("ss004","SS004","seas4")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sea_salt004_in_air', &
           fv3jedi_name = 'ss004', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss004, &
           tracer = .true.)
    case("ss005","SS005","seas5")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sea_salt005_in_air', &
           fv3jedi_name = 'ss005', units = 'kg kg-1', staggerloc = center, arraypointer = self%ss005, &
           tracer = .true.)
    case("bcphobic","BCPHOBIC","bc1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_hydrophobic_black_carbon_in_air', &
           fv3jedi_name = 'bcphobic', units = 'kg kg-1', staggerloc = center, arraypointer = self%bcphobic, &
           tracer = .true.)
    case("bcphilic","BCPHILIC","bc2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_hydrophilic_black_carbon_in_air', &
           fv3jedi_name = 'bcphilic', units = 'kg kg-1', staggerloc = center, arraypointer = self%bcphilic, &
           tracer = .true.)
    case("ocphobic","OCPHOBIC","oc1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_hydrophobic_organic_carbon_in_air', &
           fv3jedi_name = 'ocphobic', units = 'kg kg-1', staggerloc = center, arraypointer = self%ocphobic, &
           tracer = .true.)
    case("ocphilic","OCPHILIC","oc2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_hydrophilic_organic_carbon_in_air', &
           fv3jedi_name = 'ocphilic', units = 'kg kg-1', staggerloc = center, arraypointer = self%ocphilic, &
           tracer = .true.)
    case("no3an1","NO3AN1")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_nitrate001_in_air', &
           fv3jedi_name = 'no3an1', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an1, &
           tracer = .true.)
    case("no3an2","NO3AN2")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_nitrate002_in_air', &
           fv3jedi_name = 'no3an2', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an2, &
           tracer = .true.)
    case("no3an3","NO3AN3")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_nitrate003_in_air', &
           fv3jedi_name = 'no3an3', units = 'kg kg-1', staggerloc = center, arraypointer = self%no3an3, &
           tracer = .true.)
    case("so4","SO4","sulf")
      vcount=vcount+1;
      call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
           short_name = vars%variable(var), long_name = 'mass_fraction_of_sulfate_in_air', &
           fv3jedi_name = 'so4', units = 'kg kg-1', staggerloc = center, arraypointer = self%so4, &
           tracer = .true.)
     !Not found
     case default
       call abor1_ftn("Create: unknown variable "//trim(vars%variable(var)))
   end select
enddo

if (vcount .ne. self%nf) &
call abor1_ftn("fv3jedi_state_mod create: vcount does not equal self%nf")

self%hydrostatic = .true.
if (associated(self%w) .and. associated(self%delz)) self%hydrostatic = .false.

! Initialize all arrays to zero
call zeros(self)

! Copy some geometry for convenience
self%isc    = geom%isc
self%iec    = geom%iec
self%jsc    = geom%jsc
self%jec    = geom%jec
self%npx    = geom%npx
self%npy    = geom%npy
self%npz    = geom%npz
self%ntile  = geom%ntile
self%ntiles = geom%ntiles

! Pointer to fv3jedi communicator
self%f_comm = geom%f_comm

! Check winds
if (associated(self%ua) .and. .not.associated(self%va)) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid u but not v")
if (.not.associated(self%ua) .and. associated(self%va)) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid v but not u")
if (associated(self%ud) .and. .not.associated(self%vd)) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid u but not v")
if (.not.associated(self%ud) .and. associated(self%vd)) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid v but not u")

self%have_agrid = .false.
self%have_dgrid = .false.
if (associated(self%ua)) self%have_agrid = .true.
if (associated(self%ud)) self%have_dgrid = .true.

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_state), intent(inout) :: self
integer :: var

if (associated(self%ud       )) nullify(self%ud       )
if (associated(self%vd       )) nullify(self%vd       )
if (associated(self%ua       )) nullify(self%ua       )
if (associated(self%va       )) nullify(self%va       )
if (associated(self%t        )) nullify(self%t        )
if (associated(self%tv       )) nullify(self%tv       )
if (associated(self%pt       )) nullify(self%pt       )
if (associated(self%delp     )) nullify(self%delp     )
if (associated(self%pe       )) nullify(self%pe       )
if (associated(self%pkz      )) nullify(self%pkz      )
if (associated(self%ps       )) nullify(self%ps       )
if (associated(self%q        )) nullify(self%q        )
if (associated(self%rh       )) nullify(self%rh       )
if (associated(self%qi       )) nullify(self%qi       )
if (associated(self%ql       )) nullify(self%ql       )
if (associated(self%qils     )) nullify(self%qils     )
if (associated(self%qlls     )) nullify(self%qlls     )
if (associated(self%qicn     )) nullify(self%qicn     )
if (associated(self%qlcn     )) nullify(self%qlcn     )
if (associated(self%qilsf    )) nullify(self%qilsf    )
if (associated(self%qicnf    )) nullify(self%qicnf    )
if (associated(self%qs       )) nullify(self%qs       )
if (associated(self%qr       )) nullify(self%qr       )
if (associated(self%gr       )) nullify(self%gr       )
if (associated(self%ca       )) nullify(self%ca       )
if (associated(self%o3       )) nullify(self%o3       )
if (associated(self%ox       )) nullify(self%ox       )
if (associated(self%w        )) nullify(self%w        )
if (associated(self%delz     )) nullify(self%delz     )
if (associated(self%phis     )) nullify(self%phis     )
if (associated(self%psi      )) nullify(self%psi      )
if (associated(self%chi      )) nullify(self%chi      )
if (associated(self%vort     )) nullify(self%vort     )
if (associated(self%divg     )) nullify(self%divg     )
if (associated(self%slmsk    )) nullify(self%slmsk    )
if (associated(self%sheleg   )) nullify(self%sheleg   )
if (associated(self%tsea     )) nullify(self%tsea     )
if (associated(self%vtype    )) nullify(self%vtype    )
if (associated(self%stype    )) nullify(self%stype    )
if (associated(self%vfrac    )) nullify(self%vfrac    )
if (associated(self%stc      )) nullify(self%stc      )
if (associated(self%smc      )) nullify(self%smc      )
if (associated(self%snwdph   )) nullify(self%snwdph   )
if (associated(self%u_srf    )) nullify(self%u_srf    )
if (associated(self%v_srf    )) nullify(self%v_srf    )
if (associated(self%f10m     )) nullify(self%f10m     )
if (associated(self%qls      )) nullify(self%qls      )
if (associated(self%qcn      )) nullify(self%qcn      )
if (associated(self%cfcn     )) nullify(self%cfcn     )
if (associated(self%frocean  )) nullify(self%frocean  )
if (associated(self%frland   )) nullify(self%frland   )
if (associated(self%frlandice)) nullify(self%frlandice)
if (associated(self%frlake   )) nullify(self%frlake   )
if (associated(self%frseaice )) nullify(self%frseaice )
if (associated(self%varflt   )) nullify(self%varflt   )
if (associated(self%ustar    )) nullify(self%ustar    )
if (associated(self%bstar    )) nullify(self%bstar    )
if (associated(self%zpbl     )) nullify(self%zpbl     )
if (associated(self%cm       )) nullify(self%cm       )
if (associated(self%ct       )) nullify(self%ct       )
if (associated(self%cq       )) nullify(self%cq       )
if (associated(self%kcbl     )) nullify(self%kcbl     )
if (associated(self%tsm      )) nullify(self%tsm      )
if (associated(self%khl      )) nullify(self%khl      )
if (associated(self%khu      )) nullify(self%khu      )
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

do var = 1, self%nf
  call self%fields(var)%deallocate_field()
enddo
deallocate(self%fields)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)

implicit none
type(fv3jedi_state), intent(inout) :: self
integer :: var

do var = 1, self%nf
  self%fields(var)%array = 0.0_kind_real
enddo

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)

implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.copy")

do var = 1, self%nf
  self%fields(var) = rhs%fields(var)
enddo

self%calendar_type = rhs%calendar_type
self%date_init = rhs%date_init

end subroutine copy

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

implicit none
type(fv3jedi_state),  intent(inout) :: self
real(kind=kind_real), intent(in)    :: zz
type(fv3jedi_state),  intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.axpy")

do var = 1, self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine add_incr(geom,self,rhs)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_state),     intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var,i,j,k
logical :: found_neg
type(fv3jedi_field), pointer :: field_pointer

real(kind=kind_real), allocatable :: rhs_ud(:,:,:), rhs_vd(:,:,:)
real(kind=kind_real), allocatable :: rhs_delp(:,:,:)
real(kind=kind_real), allocatable :: rhs_pt(:,:,:)


! Handle special cases, e.g. D grid winds to A grid winds
! -------------------------------------------------------

! Get D-Grid winds if necessary
if (rhs%have_agrid) then !A-Grid in increment
  if (.not.rhs%have_dgrid) then !D-Grid not in increment
    if (self%have_dgrid) then !D-grid in state
      allocate(rhs_ud(rhs%isc:rhs%iec  ,rhs%jsc:rhs%jec+1,1:rhs%npz))
      allocate(rhs_vd(rhs%isc:rhs%iec+1,rhs%jsc:rhs%jec  ,1:rhs%npz))
      call a2d(geom, rhs%ua, rhs%va, rhs_ud, rhs_vd) !Linear
    endif
  endif
endif


! Convert ps to delp if necessary
! -------------------------------
if (associated(rhs%ps)) then !ps in increment
  if (.not.associated(rhs%delp)) then !delp not in increment
    if (associated(self%delp)) then !delp in state
      allocate(rhs_delp(rhs%isc:rhs%iec,rhs%jsc:rhs%jec,1:rhs%npz))
      do k = 1,rhs%npz
        rhs_delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*rhs%ps(:,:,1) !TLM
      enddo
    endif
  endif
endif

!Fields to add determined from increment
do var = 1,rhs%nf

  !Winds are a special case
  if (rhs%fields(var)%fv3jedi_name == 'ua') then

    if (associated(self%ua)) self%ua = self%ua + rhs%ua
    if (associated(self%ud) .and. .not.associated(rhs%ud)) self%ud = self%ud + rhs_ud

  elseif (rhs%fields(var)%fv3jedi_name == 'va') then

    if (associated(self%va)) self%va = self%va + rhs%va
    if (associated(self%vd) .and. .not.associated(rhs%vd)) self%vd = self%vd + rhs_vd

  elseif (rhs%fields(var)%fv3jedi_name == 't') then

    if (associated(self%t)) self%t = self%t + rhs%t
    if (associated(self%pt)) self%pt = self%pt + rhs_pt

  elseif (rhs%fields(var)%fv3jedi_name == 'ps') then

    if (associated(self%ps)) self%ps = self%ps + rhs%ps
    if (associated(self%delp) .and. .not.associated(rhs%delp)) self%delp = self%delp + rhs_delp

  else

    !Get pointer to state
    call get_field(self%nf,self%fields,rhs%fields(var)%fv3jedi_name,field_pointer)

    !Add increment to state
    field_pointer%array = field_pointer%array + rhs%fields(var)%array

    !Nullify pointer
    nullify(field_pointer)

  endif

enddo

if (allocated(rhs_ud)) deallocate(rhs_ud)
if (allocated(rhs_vd)) deallocate(rhs_vd)
if (allocated(rhs_pt)) deallocate(rhs_pt)
if (allocated(rhs_delp)) deallocate(rhs_delp)

!Check for negative tracers and increase to 0.0
do var = 1,self%nf

  if (self%fields(var)%tracer) then

    found_neg = .false.

    do k = 1,self%fields(var)%npz
      do j = geom%jsc,geom%jec
        do i = geom%isc,geom%iec
          if (self%fields(var)%array(i,j,k) < 0.0_kind_real) then
            found_neg = .true.
            self%fields(var)%array(i,j,k) = 0.0_kind_real
          endif
        enddo
      enddo
    enddo

    !Print message warning about negative tracer removal
    if (found_neg .and. self%f_comm%rank() == 0) print*, &
             'fv3jedi_state_mod.add_incr: Removed negative values for '&
             //trim(self%fields(var)%fv3jedi_name)

  endif

enddo

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine change_resol(self,geom,rhs,geom_rhs)

implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_geom),  intent(inout) :: geom
type(fv3jedi_state), intent(in)    :: rhs
type(fv3jedi_geom),  intent(inout) :: geom_rhs

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.change_resol")

if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1) == 0) then
  call copy(self, rhs)
else
  call field2field_interp(self%nf, geom_rhs, rhs%fields, geom, self%fields)
  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init
endif

end subroutine change_resol

! ------------------------------------------------------------------------------
!> Analytic Initialization for the FV3 Model
!!
!! \details **analytic_IC()** initializes the FV3JEDI state and State objects using one of
!! several alternative idealized analytic models.  This is intended to facilitate testing by
!! eliminating the need to read in the initial state from a file and by providing exact expressions
!! to test interpolations.  This function is activated by setting the "analytic_init" state in the
!! "initial" or "StateFile" section of the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!! * dcmip-test-3-1: Non-hydrostatic gravity wave
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \author M. Miesch (adapted from a pre-existing call to invent_state)
!! \date March, 2018: Created
!! \date May, 2018: Added dcmip-test-3-1
!! \date June, 2018: Added dcmip-test-4-0
!!
!! \warning This routine initializes the fv3jedi_state object.  However, since the fv_atmos_type
!! component of fv3jedi_state is a subset of the corresponding object in the fv3 model,
!! this initialization routine is not sufficient to comprehensively define the full fv3 state.
!! So, this intitialization can be used for interpolation and other tests within JEDI but it is
!! cannot currently be used to initiate a forecast with gfs.
!!
!! \warning This routine does not initialize the fv3jedi_interp member of the fv3jedi_state object
!!
!! \warning Though an input state file is not required for these analytic initialization routines,
!! some grid information (in particular the hybrid vertical grid coefficients ak and bk)
!! is still read in from an input file when creating the geometry object that is a required
!! member of fv3jedi_state; see c_fv3jedi_geo_setup() in fv3jedi_geom_mod.F90.
!!
!! \warning It's unclear whether the pt member of the fv_atmos_type structure is potential temperature
!! or temperature.  This routine assumes the latter.  If this is not correct, then we will need to
!! implement a conversion
!!
subroutine analytic_IC(self, geom, c_conf, vdate)

  use fv3jedi_kinds_mod
  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log
  use constants_mod, only: pi=>pi_8
  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

  !FV3 Test Cases
  use fv_arrays_nlm_mod,  only: fv_atmos_type, deallocate_fv_atmos_type
  use test_cases_nlm_mod, only: init_case, test_case
  use fv_control_nlm_mod, only: fv_init, pelist_all

  implicit none

  type(fv3jedi_state), intent(inout) :: self    !< State
  type(fv3jedi_geom),  intent(inout) :: geom    !< Geometry
  type(c_ptr), intent(in)            :: c_conf  !< Configuration
  type(datetime), intent(inout)      :: vdate   !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: i,j,k
  real(kind=kind_real) :: rlat, rlon
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  type(fv_atmos_type), allocatable :: FV_AtmIC(:)
  real(kind=kind_real)             :: DTdummy = 900.0
  logical, allocatable             :: grids_on_this_pe(:)
  integer                          :: p_split = 1

  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  ! Fortran configuration
  ! ---------------------
  f_conf = fckit_configuration(c_conf)

  If (f_conf%has("analytic_init")) Then
     call f_conf%get_or_die("analytic_init",str)
     IC = str
     deallocate(str)
  EndIf

  call log%warning("fv3jedi_state:analytic_init: "//IC)
  call f_conf%get_or_die("date",str)
  sdate = str
  deallocate(str)
  WRITE(buf,*) 'validity date is: '//sdate
  call log%info(buf)
  call datetime_set(sdate, vdate)

  !===================================================================
  int_option: Select Case (IC)

     Case("fv3_init_case")

        !Initialize temporary FV_Atm fv3 construct
        call fv_init(FV_AtmIC, DTdummy, grids_on_this_pe, p_split)
        deallocate(pelist_all)

        !Test case to run, see fv3: /tools/test_cases.F90 for possibilities
        call f_conf%get_or_die("fv3_test_case",test_case)

        call init_case( FV_AtmIC(1)%u,FV_AtmIC(1)%v,FV_AtmIC(1)%w,FV_AtmIC(1)%pt,FV_AtmIC(1)%delp,FV_AtmIC(1)%q, &
                        FV_AtmIC(1)%phis, FV_AtmIC(1)%ps,FV_AtmIC(1)%pe, FV_AtmIC(1)%peln,FV_AtmIC(1)%pk,FV_AtmIC(1)%pkz, &
                        FV_AtmIC(1)%uc,FV_AtmIC(1)%vc, FV_AtmIC(1)%ua,FV_AtmIC(1)%va,        &
                        FV_AtmIC(1)%ak, FV_AtmIC(1)%bk, FV_AtmIC(1)%gridstruct, FV_AtmIC(1)%flagstruct,&
                        FV_AtmIC(1)%npx, FV_AtmIC(1)%npy, FV_AtmIC(1)%npz, FV_AtmIC(1)%ng, &
                        FV_AtmIC(1)%flagstruct%ncnst, FV_AtmIC(1)%flagstruct%nwat,  &
                        FV_AtmIC(1)%flagstruct%ndims, FV_AtmIC(1)%flagstruct%ntiles, &
                        FV_AtmIC(1)%flagstruct%dry_mass, &
                        FV_AtmIC(1)%flagstruct%mountain,       &
                        FV_AtmIC(1)%flagstruct%moist_phys, FV_AtmIC(1)%flagstruct%hydrostatic, &
                        FV_AtmIC(1)%flagstruct%hybrid_z, FV_AtmIC(1)%delz, FV_AtmIC(1)%ze0, &
                        FV_AtmIC(1)%flagstruct%adiabatic, FV_AtmIC(1)%ks, FV_AtmIC(1)%neststruct%npx_global, &
                        FV_AtmIC(1)%ptop, FV_AtmIC(1)%domain, FV_AtmIC(1)%tile, FV_AtmIC(1)%bd )

        !Copy from temporary structure into state
        self%ud = FV_AtmIC(1)%u
        self%vd = FV_AtmIC(1)%v
        self%t = FV_AtmIC(1)%pt
        self%delp = FV_AtmIC(1)%delp
        self%q = FV_AtmIC(1)%q(:,:,:,1)
        self%phis(:,:,1) = FV_AtmIC(1)%phis
        geom%ak = FV_AtmIC(1)%ak
        geom%ak = FV_AtmIC(1)%ak
        geom%ptop = FV_AtmIC(1)%ptop
        if (.not. self%hydrostatic) then
           self%w = FV_AtmIC(1)%w
           self%delz = FV_AtmIC(1)%delz
        endif

        !Deallocate temporary FV_Atm fv3 structure
        call deallocate_fv_atmos_type(FV_AtmIC(1))
        deallocate(FV_AtmIC)
        deallocate(grids_on_this_pe)

     Case ("dcmip-test-1-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
                                               phis0,ps,rho0,hum0,q1,q2,q3,q4)

              self%phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                                                  phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 self%ud(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 self%vd(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
                 If (.not.self%hydrostatic) self%w(i,j,k) = w0
                 self%t(i,j,k) = t0
                 self%delp(i,j,k) = pe2-pe1
                 self%q (i,j,k) = hum0
                 self%qi(i,j,k) = q1
                 self%ql(i,j,k) = q2
                 self%o3(i,j,k) = q3

              enddo
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                          t0,phis0,ps,rho0,hum0,q1)

              self%phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                             t0,phis0,ps,rho0,hum0,q1)

                 self%ud(i,j,k) = u0 !ATTN comment above
                 self%vd(i,j,k) = v0
                 If (.not.self%hydrostatic) self%w(i,j,k) = w0
                 self%t(i,j,k) = t0
                 self%delp(i,j,k) = pe2-pe1
                 self%q(i,j,k) = hum0
                 self%qi(i,j,k) = q1

              enddo
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                      t0,phis0,ps,rho0,hum0)

              self%phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0)

                 self%ud(i,j,k) = u0 !ATTN comment above
                 self%vd(i,j,k) = v0
                 If (.not.self%hydrostatic) self%w(i,j,k) = w0
                 self%t(i,j,k) = t0
                 self%delp(i,j,k) = pe2-pe1
                 self%q(i,j,k) = hum0

              enddo
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

              self%phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

                 self%ud(i,j,k) = u0 !ATTN comment above
                 self%vd(i,j,k) = v0
                 If (.not.self%hydrostatic) self%w(i,j,k) = w0
                 self%t(i,j,k) = t0
                 self%delp(i,j,k) = pe2-pe1
                 self%q(i,j,k) = hum0

              enddo
           enddo
        enddo

     Case Default

        call abor1_ftn("fv3jedi_state analytic_IC: provide analytic_init")

     End Select int_option

end subroutine analytic_IC

! ------------------------------------------------------------------------------

subroutine read_file(geom, self, c_conf, vdate)
use string_utils

implicit none

type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
type(fv3jedi_state), intent(inout) :: self     !< State
type(c_ptr),         intent(in)    :: c_conf   !< Configuration
type(datetime),      intent(inout) :: vdate    !< DateTime

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos

character(len=10) :: filetype
character(len=255) :: filename
integer :: flipvert
type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)

call f_conf%get_or_die("filetype",str)
filetype = str
deallocate(str)


if (trim(filetype) == 'gfs') then

  call gfs%setup(f_conf)
  call gfs%read_meta(geom, vdate, self%calendar_type, self%date_init)
  call gfs%read_fields(geom, self%fields)

  flipvert = 0
  if (f_conf%has("flip_vertically")) then
     call f_conf%get_or_die("flip_vertically",flipvert)
  endif
  if (flipvert==1) call flip_array_vertical(self%nf, self%fields)

elseif (trim(filetype) == 'geos') then

  call geos%setup(geom, self%fields, vdate, 'read', f_conf)
  call geos%read_meta(geom, vdate, self%calendar_type, self%date_init)
  call geos%read_fields(geom, self%fields)
  call geos%delete()

else

  call abor1_ftn("fv3jedi_state_mod.read: restart type not supported")

endif

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(geom, self, c_conf, vdate)

  use fv3jedi_io_latlon_mod

  implicit none

  type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
  type(fv3jedi_state), intent(inout) :: self     !< State
  type(c_ptr),         intent(in)    :: c_conf   !< Configuration
  type(datetime),      intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype
  integer :: flipvert
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str


  ! Fortran configuration
  ! ---------------------
  f_conf = fckit_configuration(c_conf)


  call f_conf%get_or_die("filetype",str)
  filetype = str
  deallocate(str)

  if (trim(filetype) == 'gfs') then

    flipvert = 0
    if (f_conf%has("flip_vertically")) then
      call f_conf%get_or_die("flip_vertically",flipvert)
    endif
    if (flipvert==1) call flip_array_vertical(self%nf, self%fields)

    call gfs%setup(f_conf)
    call gfs%write_all(geom, self%fields, vdate, self%calendar_type, self%date_init)

    if (flipvert==1) call flip_array_vertical(self%nf, self%fields)

  elseif (trim(filetype) == 'geos') then

    call geos%setup(geom, self%fields, vdate, 'write', f_conf)
    call geos%write_all(geom, self%fields, vdate)
    call geos%delete()

  else

     call abor1_ftn("fv3jedi_state_mod.write: restart type not supported")

  endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine state_print(self)

implicit none
type(fv3jedi_state), intent(in) :: self

call fields_print(self%nf, self%fields, "State", self%f_comm)

end subroutine state_print

! ------------------------------------------------------------------------------

subroutine gpnorm(self, nf, pstat)

implicit none
type(fv3jedi_state),  intent(in)    :: self
integer,              intent(in)    :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)

if (nf .ne. self%nf) then
  call abor1_ftn("fv3jedi_state: gpnorm | nf passed in does not match expeted nf")
endif

call fields_gpnorm(nf, self%fields, pstat, self%f_comm)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_state),  intent(in)  :: self
real(kind=kind_real), intent(out) :: prms

call fields_rms(self%nf, self%fields, prms, self%f_comm)

end subroutine rms

! ------------------------------------------------------------------------------

end module fv3jedi_state_mod
