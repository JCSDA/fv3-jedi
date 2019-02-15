! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_state_mod

use iso_c_binding
use config_mod
use datetime_mod
use fckit_mpi_module

use fv3jedi_field_mod,           only: fv3jedi_field, fields_rms, fields_rms, fields_gpnorm, fields_print, get_field
use fv3jedi_constants_mod,       only: rad2deg, constoz
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs 
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos 
use fv3jedi_state_utils_mod,     only: fv3jedi_state
use fv3jedi_vars_mod,            only: fv3jedi_vars
use fv3jedi_getvalues_mod,       only: getvalues

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
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_geom),  intent(in)    :: geom
type(fv3jedi_vars),  intent(in)    :: vars

integer :: var, vcount

! Total fields
! ------------
self%nf = vars%nv

! Allocate fields structure
! -------------------------
allocate(self%fields(self%nf))

! Loop through and allocate main state fields
! -------------------------------------------
vcount = 0
do var = 1, vars%nv
   select case (trim(vars%fldnames(var)))
     case("ud","u")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'eastward_wind_on_native_D-Grid', &
            fv3jedi_name = 'ud', units = 'm s-1', staggerloc = north, arraypointer = self%ud )
     case("vd","v")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'northward_wind_on_native_D-Grid', &
            fv3jedi_name = 'vd', units = 'm s-1', staggerloc = east, arraypointer = self%vd )
     case("ua")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'eastward_wind', &
            fv3jedi_name = 'ua', units = 'm s-1', staggerloc = center, arraypointer = self%ua)
     case("va")
       vcount=vcount+1; 
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'northward_wind', &
            fv3jedi_name = 'va', units = 'm s-1', staggerloc = center, arraypointer = self%va )
     case("t","T")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'air_temperature', &
            fv3jedi_name = 't', units = 'K', staggerloc = center, arraypointer = self%t)
     case("delp","DELP")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'pressure_thickness', &
            fv3jedi_name = 'delp', units = 'Pa', staggerloc = center, arraypointer = self%delp)
     case("ps")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_pressure', &
            fv3jedi_name = 'ps', units = 'Pa', staggerloc = center, arraypointer = self%ps)
     case("q","sphum")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'specific_humidity', &
            fv3jedi_name = 'q', units = 'kg kg-1', staggerloc = center, arraypointer = self%q, &
            tracer = .true.)
     case("qi","ice_wat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'cloud_liquid_ice', &
            fv3jedi_name = 'qi', units = 'kg kg-1', staggerloc = center, arraypointer = self%qi, &
            tracer = .true.)
     case("ql","liq_wat")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'cloud_liquid_ice_water', &
            fv3jedi_name = 'ql', units = 'kg kg-1', staggerloc = center, arraypointer = self%ql, &
            tracer = .true.)
     case("o3","o3mr")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'ozone_mass_mixing_ratio', &
            fv3jedi_name = 'o3', units = 'kg kg-1', staggerloc = center, arraypointer = self%o3, &
            tracer = .true.)
     case("w","W")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'vertical_wind', &
            fv3jedi_name = 'w', units = 'm s-1', staggerloc = center, arraypointer = self%w)
     case("delz","DZ")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'layer_thickness', &
            fv3jedi_name = 'delz', units = 'm', staggerloc = center, arraypointer = self%delz)
     case("phis")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_geopotential_height', &
            fv3jedi_name = 'phis', units = 'm', staggerloc = center, arraypointer = self%phis)
     !CRTM
     case("slmsk")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'slmsk', &
            fv3jedi_name = 'slmsk', units = 'none', staggerloc = center, arraypointer = self%slmsk)
     case("sheleg")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'sheleg', &
            fv3jedi_name = 'sheleg', units = 'none', staggerloc = center, arraypointer = self%sheleg)
     case("tsea")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'tsea', &
            fv3jedi_name = 'tsea', units = 'none', staggerloc = center, arraypointer = self%tsea)
     case("vtype")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'vtype', &
            fv3jedi_name = 'vtype', units = 'none', staggerloc = center, arraypointer = self%vtype)
     case("stype")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'stype', &
            fv3jedi_name = 'stype', units = 'none', staggerloc = center, arraypointer = self%stype)
     case("vfrac")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'vfrac', &
            fv3jedi_name = 'vfrac', units = 'none', staggerloc = center, arraypointer = self%vfrac)
     case("stc")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,4, &
            short_name = vars%fldnames(var), long_name = 'stc', &
            fv3jedi_name = 'stc', units = 'none', staggerloc = center, arraypointer = self%stc)
     case("smc")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,4, &
            short_name = vars%fldnames(var), long_name = 'smc', &
            fv3jedi_name = 'smc', units = 'none', staggerloc = center, arraypointer = self%smc)
     case("snwdph")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'snwdph', &
            fv3jedi_name = 'snwdph', units = 'none', staggerloc = center, arraypointer = self%snwdph)
     case("u_srf")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'u_srf', &
            fv3jedi_name = 'u_srf', units = 'none', staggerloc = center, arraypointer = self%u_srf)
     case("v_srf")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'v_srf', &
            fv3jedi_name = 'v_srf', units = 'none', staggerloc = center, arraypointer = self%v_srf)
     case("f10m")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'f10m', &
            fv3jedi_name = 'f10m', units = 'none', staggerloc = center, arraypointer = self%f10m)
     !TL/AD trajectory
     case("qls")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'initial_mass_fraction_of_large_scale_cloud_condensate', &
            fv3jedi_name = 'qls', units = 'kg kg-1', staggerloc = center, arraypointer = self%qls)
     case("qcn")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'initial_mass_fraction_of_convective_cloud_condensate', &
            fv3jedi_name = 'qcn', units = 'kg kg-1', staggerloc = center, arraypointer = self%qcn)
     case("cfcn")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,geom%npz, &
            short_name = vars%fldnames(var), long_name = 'convective_cloud_area_fraction', &
            fv3jedi_name = 'cfcn', units = '1', staggerloc = center, arraypointer = self%cfcn)
     case("frocean")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'fraction_of_ocean', &
            fv3jedi_name = 'frocean', units = '1', staggerloc = center, arraypointer = self%frocean)
     case("frland")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'fraction_of_land', &
            fv3jedi_name = 'frland', units = '1', staggerloc = center, arraypointer = self%frland)
     case("varflt")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'isotropic_variance_of_filtered_topography', &
            fv3jedi_name = 'varflt', units = 'm+2', staggerloc = center, arraypointer = self%varflt)
     case("ustar")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_velocity_scale', &
            fv3jedi_name = 'ustar', units = 'm s-1', staggerloc = center, arraypointer = self%ustar)
     case("bstar")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_bouyancy_scale', &
            fv3jedi_name = 'bstar', units = 'm s-2', staggerloc = center, arraypointer = self%bstar)
     case("zpbl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'planetary_boundary_layer_height', &
            fv3jedi_name = 'zpbl', units = 'm', staggerloc = center, arraypointer = self%zpbl)
     case("cm")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_exchange_coefficient_for_momentum', &
            fv3jedi_name = 'cm', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%cm)
     case("ct")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_exchange_coefficient_for_heat', &
            fv3jedi_name = 'ct', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%ct)
     case("cq")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_exchange_coefficient_for_moisture', &
            fv3jedi_name = 'cq', units = 'kg m-2 s-1', staggerloc = center, arraypointer = self%cq)
     case("kcbl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'KCBL_before_moist', &
            fv3jedi_name = 'kcbl', units = '1', staggerloc = center, arraypointer = self%kcbl)
     case("ts")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'surface_temp_before_moist', &
            fv3jedi_name = 'ts', units = 'K', staggerloc = center, arraypointer = self%ts)
     case("khl")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'lower_index_where_Kh_greater_than_2', &
            fv3jedi_name = 'khl', units = '1', staggerloc = center, arraypointer = self%khl)
     case("khu")
       vcount=vcount+1;
       call self%fields(vcount)%allocate_field(geom%isc,geom%iec,geom%jsc,geom%jec,1, &
            short_name = vars%fldnames(var), long_name = 'upper_index_where_Kh_greater_than_2', &
            fv3jedi_name = 'khu', units = '1', staggerloc = center, arraypointer = self%khu)
     case default 
       call abor1_ftn("Create: unknown variable "//trim(vars%fldnames(var)))
   end select
enddo

if (vcount .ne. self%nf) &
call abor1_ftn("fv3jedi_state_mod create: vcount does not equal self%nf")

self%hydrostatic = .true.
if (associated(self%w) .and. associated(self%delz)) self%hydrostatic = .false.

self%tladphystrj = .false.
if (associated(self%khu)) self%tladphystrj = .true.

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
self%f_comm = fckit_mpi_comm()

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_state), intent(inout) :: self
integer :: var

do var = 1, self%nf
  call self%fields(var)%deallocate_field()
enddo
deallocate(self%fields)

if (associated(self%ud     )) nullify(self%ud     )
if (associated(self%vd     )) nullify(self%vd     )
if (associated(self%ua     )) nullify(self%ua     )
if (associated(self%va     )) nullify(self%va     )
if (associated(self%t      )) nullify(self%t      )
if (associated(self%delp   )) nullify(self%delp   )
if (associated(self%ps     )) nullify(self%ps     )
if (associated(self%q      )) nullify(self%q      )
if (associated(self%qi     )) nullify(self%qi     )
if (associated(self%ql     )) nullify(self%ql     )
if (associated(self%o3     )) nullify(self%o3     )
if (associated(self%w      )) nullify(self%w      )
if (associated(self%delz   )) nullify(self%delz   )
if (associated(self%phis   )) nullify(self%phis   )
if (associated(self%slmsk  )) nullify(self%slmsk  )
if (associated(self%sheleg )) nullify(self%sheleg )
if (associated(self%tsea   )) nullify(self%tsea   )
if (associated(self%vtype  )) nullify(self%vtype  )
if (associated(self%stype  )) nullify(self%stype  )
if (associated(self%vfrac  )) nullify(self%vfrac  )
if (associated(self%stc    )) nullify(self%stc    )
if (associated(self%smc    )) nullify(self%smc    )
if (associated(self%snwdph )) nullify(self%snwdph )
if (associated(self%u_srf  )) nullify(self%u_srf  )
if (associated(self%v_srf  )) nullify(self%v_srf  )
if (associated(self%f10m   )) nullify(self%f10m   )
if (associated(self%qls    )) nullify(self%qls    )
if (associated(self%qcn    )) nullify(self%qcn    )
if (associated(self%cfcn   )) nullify(self%cfcn   )
if (associated(self%frocean)) nullify(self%frocean)
if (associated(self%frland )) nullify(self%frland )
if (associated(self%varflt )) nullify(self%varflt )
if (associated(self%ustar  )) nullify(self%ustar  )
if (associated(self%bstar  )) nullify(self%bstar  )
if (associated(self%zpbl   )) nullify(self%zpbl   )
if (associated(self%cm     )) nullify(self%cm     )
if (associated(self%ct     )) nullify(self%ct     )
if (associated(self%cq     )) nullify(self%cq     )
if (associated(self%kcbl   )) nullify(self%kcbl   )
if (associated(self%ts     )) nullify(self%ts     )
if (associated(self%khl    )) nullify(self%khl    )
if (associated(self%khu    )) nullify(self%khu    )

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

logical :: found
integer :: self_var, rhs_var

self%isc            = rhs%isc
self%iec            = rhs%iec
self%jsc            = rhs%jsc
self%jec            = rhs%jec
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

  !Direct copy of one state to another
  self%nf = rhs%nf
  allocate(self%fields(self%nf))
  do self_var = 1, self%nf
    self%fields(self_var) = rhs%fields(self_var)
  enddo

else

  !State copy, potentialy with differnt fields
  do self_var = 1, self%nf
    found = .false.
    do rhs_var = 1, rhs%nf
      if (trim(self%fields(self_var)%fv3jedi_name) == trim(rhs%fields(rhs_var)%fv3jedi_name)) then
        self%fields(self_var) = rhs%fields(rhs_var)
        found = .true.
        exit
      endif
    enddo
    if (.not.found) call abor1_ftn("fv3jedi_state: Error in state copy, field "//&
                    trim(self%fields(self_var)%fv3jedi_name)//" not found in state being copied from." )
  enddo

endif  

! Set pointers
do self_var = 1, self%nf
  select case (trim(self%fields(self_var)%fv3jedi_name))
    case("ud","u")
      call self%fields(self_var)%array_pointer(self%ud)
    case("vd","v")
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
    case("q","sphum")
      call self%fields(self_var)%array_pointer(self%q)
    case("qi","ice_wat")
      call self%fields(self_var)%array_pointer(self%qi)
    case("ql","liq_wat")
      call self%fields(self_var)%array_pointer(self%ql)
    case("o3","o3mr")
      call self%fields(self_var)%array_pointer(self%o3)
    case("w","W")
      call self%fields(self_var)%array_pointer(self%w)
    case("delz","DZ")
      call self%fields(self_var)%array_pointer(self%delz)
    case("phis")
      call self%fields(self_var)%array_pointer(self%phis)
    case("slmsk")
      call self%fields(self_var)%array_pointer(self%slmsk)
    case("sheleg")
      call self%fields(self_var)%array_pointer(self%sheleg)
    case("tsea")
      call self%fields(self_var)%array_pointer(self%tsea)
    case("vtype")
      call self%fields(self_var)%array_pointer(self%vtype)
    case("stype")
      call self%fields(self_var)%array_pointer(self%stype)
    case("vfrac")
      call self%fields(self_var)%array_pointer(self%vfrac)
    case("stc")
      call self%fields(self_var)%array_pointer(self%stc)
    case("smc")
      call self%fields(self_var)%array_pointer(self%smc)
    case("snwdph")
      call self%fields(self_var)%array_pointer(self%snwdph)
    case("u_srf")
      call self%fields(self_var)%array_pointer(self%u_srf)
    case("v_srf")
      call self%fields(self_var)%array_pointer(self%v_srf)
    case("f10m")
      call self%fields(self_var)%array_pointer(self%f10m)
    case("qls")
      call self%fields(self_var)%array_pointer(self%qls)
    case("qcn")
      call self%fields(self_var)%array_pointer(self%qcn)
    case("cfcn")
      call self%fields(self_var)%array_pointer(self%cfcn)
    case("frocean")
      call self%fields(self_var)%array_pointer(self%frocean)
    case("frland")
      call self%fields(self_var)%array_pointer(self%frland)
    case("varflt")
      call self%fields(self_var)%array_pointer(self%varflt)
    case("ustar")
      call self%fields(self_var)%array_pointer(self%ustar)
    case("bstar")
      call self%fields(self_var)%array_pointer(self%bstar)
    case("zpbl")
      call self%fields(self_var)%array_pointer(self%zpbl)
    case("cm")
      call self%fields(self_var)%array_pointer(self%cm)
    case("ct")
      call self%fields(self_var)%array_pointer(self%ct)
    case("cq")
      call self%fields(self_var)%array_pointer(self%cq)
    case("kcbl")
      call self%fields(self_var)%array_pointer(self%kcbl)
    case("ts")
      call self%fields(self_var)%array_pointer(self%ts)
    case("khl")
      call self%fields(self_var)%array_pointer(self%khl)
    case("khu")
      call self%fields(self_var)%array_pointer(self%khu)
    case default 
      call abor1_ftn("fv3jedi_state_mod.copy: unknown variable "//trim(self%fields(self_var)%fv3jedi_name))
  end select
enddo

end subroutine copy

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

implicit none
type(fv3jedi_state),  intent(inout) :: self
real(kind=kind_real), intent(in)    :: zz
type(fv3jedi_state),  intent(in)    :: rhs

integer :: var

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

!Check for matching resolution between state and increment
if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1)==0) then

  do var = 1,self%nf
  
    !Get pointer to increment
    call get_field(rhs%nf,rhs%fields,self%fields(var)%fv3jedi_name,field_pointer)

    !Add increment to state
    self%fields(var)%array = self%fields(var)%array + field_pointer%array

    !Nullify pointer
    nullify(field_pointer)

  enddo

else

   call abor1_ftn("fv3jedi state:  add_incr not implemented for low res increment yet")

endif

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

subroutine change_resol(self,rhs)
implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: rhs

integer :: check
check = (rhs%iec-rhs%isc+1) - (self%iec-self%isc+1)

if (check==0) then
   call copy(self, rhs)
else
   call abor1_ftn("fv3jedi_state: change_resol not implmeneted yet")
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

  If (config_element_exists(c_conf,"analytic_init")) Then
     IC = Trim(config_get_string(c_conf,len(IC),"analytic_init"))
  EndIf

  call log%warning("fv3jedi_state:analytic_init: "//IC)
  sdate = config_get_string(c_conf,len(sdate),"date")
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
        test_case = config_get_int(c_conf,"fv3_test_case")

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

  implicit none

  type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
  type(fv3jedi_state), intent(inout) :: self     !< State
  type(c_ptr),         intent(in)    :: c_conf   !< Configuration
  type(datetime),      intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype
  character(len=255) :: filename

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%write(geom, self%fields, vdate, self%calendar_type, self%date_init)
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

  use fv3jedi_io_latlon_mod

  implicit none

  type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
  type(fv3jedi_state), intent(in)    :: self     !< State
  type(c_ptr),         intent(in)    :: c_conf   !< Configuration
  type(datetime),      intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos

  character(len=10) :: filetype

  filetype = config_get_string(c_conf,len(filetype),"filetype")

  if (trim(filetype) == 'gfs') then
    call gfs%setup(c_conf)
    call gfs%write(geom, self%fields, vdate, self%calendar_type, self%date_init)
  elseif (trim(filetype) == 'geos') then
    call geos%create(geom, 'write')
    call geos%write_all(geom, self%fields, c_conf, vdate)
    call geos%delete()
  else
     call abor1_ftn("fv3jedi_increment_mod.write: restart type not supported")
  endif

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine state_print(self)

implicit none
type(fv3jedi_state), intent(in) :: self

call fields_print(self%nf, self%fields, "State")

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

call fields_gpnorm(nf, self%fields, pstat)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_state),  intent(in)  :: self
real(kind=kind_real), intent(out) :: prms

call fields_rms(self%nf,self%fields,prms)

end subroutine rms

! ------------------------------------------------------------------------------

end module fv3jedi_state_mod
