
!> Fortran module handling interpolation trajectory for the FV3 model

module fv3jedi_getvalues_traj_mod

!General JEDI uses
use fv3jedi_kinds_mod, only: kind_real
use iso_c_binding
use type_bump, only: bump_type

implicit none
private

public fv3jedi_getvalues_traj
public fv3jedi_getvalues_traj_registry
public c_fv3jedi_getvalues_traj_setup, c_fv3jedi_getvalues_traj_delete

type :: fv3jedi_getvalues_traj
 integer :: nobs, ngrid
 real(kind=kind_real), allocatable :: t(:,:,:)
 real(kind=kind_real), allocatable :: q(:,:,:)
 type(bump_type) :: bump
 logical :: lalloc = .false.
end type fv3jedi_getvalues_traj

#define LISTED_TYPE fv3jedi_getvalues_traj

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_getvalues_traj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_getvalues_traj_setup(c_key_self) bind(c,name='fv3jedi_getvalues_traj_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_getvalues_traj), pointer :: self

! Init, add and get key
! ---------------------
call fv3jedi_getvalues_traj_registry%init()
call fv3jedi_getvalues_traj_registry%add(c_key_self)
call fv3jedi_getvalues_traj_registry%get(c_key_self,self)

self%lalloc = .false.
self%nobs = 0
self%ngrid = 0

end subroutine c_fv3jedi_getvalues_traj_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_getvalues_traj_delete(c_key_self) bind(c,name='fv3jedi_getvalues_traj_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_getvalues_traj), pointer :: self

! Get key
call fv3jedi_getvalues_traj_registry%get(c_key_self, self)

if (self%lalloc) then
  self%nobs = 0
  self%ngrid = 0
  if (allocated(self%t)) deallocate(self%t)
  if (allocated(self%q)) deallocate(self%q)
  call self%bump%dealloc
  self%lalloc = .false.
endif

! Remove key
call fv3jedi_getvalues_traj_registry%remove(c_key_self)

end subroutine c_fv3jedi_getvalues_traj_delete

! ------------------------------------------------------------------------------

end module fv3jedi_getvalues_traj_mod
