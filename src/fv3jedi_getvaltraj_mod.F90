
!> Fortran module handling interpolation trajectory for the FV3 model

module fv3jedi_getvaltraj_mod

!General JEDI uses
use kinds
use iso_c_binding
use type_bump, only: bump_type

implicit none
private

public fv3jedi_getvaltraj
public fv3jedi_getvaltraj_registry
public c_fv3jedi_getvaltraj_setup, c_fv3jedi_getvaltraj_delete

type :: fv3jedi_getvaltraj
 integer :: nobs, ngrid
 real(kind=kind_real), allocatable :: pt(:,:,:)
 real(kind=kind_real), allocatable :: q(:,:,:,:)
 type(bump_type), allocatable :: bump(:)
 logical :: lalloc = .false.
end type fv3jedi_getvaltraj

#define LISTED_TYPE fv3jedi_getvaltraj

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_getvaltraj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_getvaltraj_setup(c_key_self) bind(c,name='fv3jedi_getvaltraj_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_getvaltraj), pointer :: self

! Init, add and get key
! ---------------------
call fv3jedi_getvaltraj_registry%init()
call fv3jedi_getvaltraj_registry%add(c_key_self)
call fv3jedi_getvaltraj_registry%get(c_key_self,self)

print*, 'dh: getvaltraj_setup', c_key_self

self%lalloc = .false.
self%nobs = 0
self%ngrid = 0
allocate(self%bump(1))

end subroutine c_fv3jedi_getvaltraj_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_getvaltraj_delete(c_key_self) bind(c,name='fv3jedi_getvaltraj_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_getvaltraj), pointer :: self

! Get key
call fv3jedi_getvaltraj_registry%get(c_key_self, self)

if (self%lalloc) then
  self%nobs = 0
  self%ngrid = 0
  if (allocated(self%pt)) deallocate(self%pt)
  if (allocated(self%q)) deallocate(self%q)
  deallocate(self%bump)
  self%lalloc = .false.
endif

! Remove key
call fv3jedi_getvaltraj_registry%remove(c_key_self)

end subroutine c_fv3jedi_getvaltraj_delete

! ------------------------------------------------------------------------------

end module fv3jedi_getvaltraj_mod
