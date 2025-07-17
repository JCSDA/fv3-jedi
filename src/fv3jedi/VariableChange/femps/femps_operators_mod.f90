! (C) Copyright 2019 UCAR and 2011-2018 John Thuburn, University of Exeter, UK
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module femps_operators_mod

use netcdf

use femps_kinds_mod
use femps_grid_mod
use femps_utils_mod

implicit none
private
public fempsoprs

! Type to hold the pre-built operators
! ------------------------------------
type fempsoprs

  logical :: poisson   ! Some opertors are not needed if just used in Poisson problem

  ! Dimensions
  ! ----------
  integer :: nmsmx
  integer :: nxmisx

  integer, private :: nlsmx
  integer, private :: njsmx
  integer, private :: nhsmx
  integer, private :: nrsmx
  integer, private :: nrxsmx
  integer, private :: nwsmx
  integer, private :: ntsmx
  integer, private :: ninjmx

  real(kind=kind_real), allocatable, dimension(:,:) :: varea !area of dual faces on each grid

  integer, allocatable, dimension(:,:,:) :: eoffin !Indicates whether the normal at the j'th edge is inward or outward relative to face f
  integer, allocatable, dimension(:,:,:) :: eofvin !*Indicates whether the tangent at the j'th edge is inward or outward relative to vertex v.

  ! Hodge star, mass matrix, and related operators
  ! ----------------------------------------------
  integer, allocatable, dimension(:,:) :: nlsten     ! number of faces in stencil for L mass matrix
  integer, allocatable, dimension(:,:) :: nmsten     ! number of faces in stencil for M mass matrix
  integer, allocatable, dimension(:,:) :: njsten     ! number of vertices in stencil for J operator
  integer, allocatable, dimension(:,:) :: nhsten     ! number of edges in stencil for H operator
  integer, allocatable, dimension(:,:) :: nrsten     ! number of vertices in stencil for R operator (= self%neoff)
  integer, allocatable, dimension(:,:) :: nrxsten    ! number of faces in stencil for R transpose operator (= self%neofv)
  integer, allocatable, dimension(:,:) :: nwsten     ! number of edges in stencil for W operator
  integer, allocatable, dimension(:,:) :: ntsten     ! number of edges in stencel for T operator
  integer, allocatable, dimension(:,:) :: nxminvsten ! number of edges in stencil for approximate inverse of M

  integer, allocatable, dimension(:,:,:) :: lsten     ! stencil for L mass matrix
  integer, allocatable, dimension(:,:,:) :: msten     ! stencil for M mass matrix
  integer, allocatable, dimension(:,:,:) :: jsten     ! stencil for J operator
  integer, allocatable, dimension(:,:,:) :: hsten     ! stencil for H operator
  integer, allocatable, dimension(:,:,:) :: rsten     ! stencil for R operator (= voff)
  integer, allocatable, dimension(:,:,:) :: rxsten    ! stencil for R transpose operator (= fofv)
  integer, allocatable, dimension(:,:,:) :: wsten     ! stencil for W operator
  integer, allocatable, dimension(:,:,:) :: tsten     ! stencil for T operator
  integer, allocatable, dimension(:,:,:) :: xminvsten ! stencil for approximate inverse of M

  real(kind=kind_real), allocatable, dimension(:,:) :: jlump ! coefficients of lumped J matrix
  real(kind=kind_real), allocatable, dimension(:,:) :: mlump ! coefficients of lumped M matrix
  real(kind=kind_real), allocatable, dimension(:,:) :: hlump ! coefficients of lumped H matrix

  real(kind=kind_real), allocatable, dimension(:,:,:)   :: lmass    ! coefficients for L mass matrix
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: mmass    ! coefficients for M mass matrix
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: jstar    ! coefficients for J operator
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: hstar    ! coefficients for H operator
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: rcoeff   ! coefficients for R operator
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: rxcoeff  ! coefficients for R transpose operator
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: wcoeff   ! coefficients for W operator
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: xminv    ! coefficients for approximate inverse of M
  real(kind=kind_real), allocatable, dimension(:,:,:,:) :: tcoeff   ! coefficients for T operator

  real(kind=kind_real), allocatable, dimension(:,:) :: elong
  real(kind=kind_real), allocatable, dimension(:,:) :: elat

  ! Restriction and prolongation operators for multigrid
  ! ----------------------------------------------------
  integer,              allocatable, dimension(:,:)   :: ninj    ! number of faces in stencil for restriction operator
  integer,              allocatable, dimension(:,:,:) :: injsten ! stencil for restriction operator
  real(kind=kind_real), allocatable, dimension(:,:,:) :: injwgt  ! weights for restriction operator

  real(kind=kind_real), allocatable, dimension(:,:) :: lapdiag
  real(kind=kind_real), allocatable, dimension(:)   :: underrel

  ! Information defining compound elements
  ! --------------------------------------
  integer, private :: ncvpmx, ncspmx, ncepmx, ncvdmx, ncsdmx

  integer, private, allocatable, dimension(:) :: ncvp ! Number of internal dofs to define a compound P0 element in space Vp.
  integer, private, allocatable, dimension(:) :: ncsp ! Number of internal dofs to define a compound RT0 element in space Sp.
  integer, private, allocatable, dimension(:) :: ncep ! Number of internal dofs to define a compound P1 element in space Ep.
  integer, private, allocatable, dimension(:) :: ncvd ! Number of internal dofs to define a compound P0 element in space Vd.
  integer, private, allocatable, dimension(:) :: ncsd ! Number of internal dofs to define a compound N0 element in space Sd.

  real(kind=kind_real), private, allocatable, dimension(:,:) :: cvp ! Dofs to define a compound element in space Vp.
  real(kind=kind_real), private, allocatable, dimension(:,:) :: csp ! Dofs to define a compound element in space Sp.
  real(kind=kind_real), private, allocatable, dimension(:,:) :: cep ! Dofs to define a compound element in space Ep.
  real(kind=kind_real), private, allocatable, dimension(:,:) :: cvd ! Dofs to define a compound element in space Vp.
  real(kind=kind_real), private, allocatable, dimension(:,:) :: csd ! Dofs to define a compound element in space Sd.

  contains

   ! Public methods
   procedure, public :: setup
   procedure, public :: pdelete
   procedure, public :: delete
   procedure, public :: writeoperators
   procedure, public :: build

   ! Internal methods
   procedure, private :: ordering
   procedure, private :: buildgeom
   procedure, private :: buildcomp
   procedure, private :: buildmat
   procedure, private :: buildvp
   procedure, private :: buildsp
   procedure, private :: buildep
   procedure, private :: buildvd
   procedure, private :: buildsd
   procedure, private :: buildL
   procedure, private :: buildM
   procedure, private :: buildR
   procedure, private :: buildW
   procedure, private :: buildJ
   procedure, private :: buildH
   procedure, private :: buildinj

end type fempsoprs

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine setup(self,grid,poisson)

implicit none
class(fempsoprs),  intent(inout) :: self
type(fempsgrid),   intent(in)    :: grid
logical, optional, intent(in)    :: poisson

! These operators can be used in the more general Helmholtz solver
! some opertors are not needed for the Poisson problem, set true to avoid allocation
self%poisson = .true.
if (present(poisson)) self%poisson = poisson

self%nlsmx  = 1
self%nmsmx  = 2*grid%nefmx - 1
self%njsmx  = grid%nevmx*(grid%nefmx - 2) + 1
self%nhsmx  = 2*(grid%nevmx - 1)*(grid%nefmx - 1) - 1
self%nrsmx  = grid%nefmx
self%nrxsmx = grid%nevmx
self%nwsmx  = 2*(grid%nefmx - 1)
self%ntsmx  = grid%nefmx

allocate(self%varea(grid%nvertx,grid%ngrids))
allocate(self%eoffin(grid%nfacex,grid%nefmx,grid%ngrids))

! Hodge star, mass matrix, and related operators
! ----------------------------------------------

allocate(self%nlsten    (grid%nfacex,grid%ngrids))
allocate(self%nmsten    (grid%nedgex,grid%ngrids))
allocate(self%njsten    (grid%nvertx,grid%ngrids))
allocate(self%nhsten    (grid%nedgex,grid%ngrids))
allocate(self%nxminvsten(grid%nedgex,grid%ngrids))

allocate(self%lsten  (grid%nfacex,self%nlsmx,grid%ngrids))
allocate(self%msten  (grid%nedgex,self%nmsmx,grid%ngrids))
allocate(self%jsten  (grid%nvertx,self%njsmx,grid%ngrids))
allocate(self%hsten  (grid%nedgex,self%nhsmx,grid%ngrids))
!self%xminvsten deferred until buildxminv

allocate(self%jlump(grid%nvertx,grid%ngrids))
allocate(self%mlump(grid%nedgex,grid%ngrids))
allocate(self%hlump(grid%nedgex,grid%ngrids))

allocate(self%lmass  (grid%nfacex,self%nlsmx,grid%ngrids))
allocate(self%mmass  (grid%nedgex,self%nmsmx,grid%ngrids))
allocate(self%jstar  (grid%nvertx,self%njsmx,grid%ngrids))
allocate(self%hstar  (grid%nedgex,self%nhsmx,grid%ngrids))
!self%xminv deferred until buildxminv

allocate(self%elong(grid%nedgex,grid%ngrids))
allocate(self%elat (grid%nedgex,grid%ngrids))

allocate(self%ninj(grid%nfacex,grid%ngrids-1))
!self%injsten deferred until buildinj
!self%injwgt deferred until buildinj

allocate(self%lapdiag(grid%nfacex,grid%ngrids))
allocate(self%underrel(grid%ngrids))

! Some operators are not used for the Poisson problem
if (.not.self%poisson) then

  allocate(self%eofvin (grid%nvertx,grid%nevmx,grid%ngrids))
  allocate(self%nrsten (grid%nfacex,grid%ngrids))
  allocate(self%nrxsten(grid%nvertx,grid%ngrids))
  allocate(self%nwsten (grid%nedgex,grid%ngrids))
  allocate(self%ntsten (grid%nfacex,grid%ngrids))
  allocate(self%rsten  (grid%nfacex,self%nrsmx,grid%ngrids))
  allocate(self%rxsten (grid%nvertx,self%nrxsmx,grid%ngrids))
  allocate(self%wsten  (grid%nedgex,self%nwsmx,grid%ngrids))
  allocate(self%tsten  (grid%nfacex,self%ntsmx,grid%ngrids))
  allocate(self%rcoeff (grid%nfacex,self%nrsmx,grid%ngrids))
  allocate(self%rxcoeff(grid%nvertx,self%nrxsmx,grid%ngrids))
  allocate(self%wcoeff (grid%nedgex,self%nwsmx,grid%ngrids))
  allocate(self%tcoeff (grid%nfacex,self%ntsmx,self%ntsmx,grid%ngrids))

endif

! Information defining compound elements
! --------------------------------------
self%ncvpmx = 2*grid%nefmx
self%ncspmx = 2 + 4*grid%nefmx
self%ncepmx = 2*grid%nefmx
self%ncvdmx = 2*grid%nefmx
self%ncsdmx = 2 + 4*grid%nevmx

allocate(self%ncvp(grid%nfacex))
allocate(self%ncsp(grid%nedgex))
allocate(self%ncep(grid%nvertx))
allocate(self%ncvd(grid%nvertx))
allocate(self%ncsd(grid%nedgex))

allocate(self%cvp (grid%nfacex,self%ncvpmx))
allocate(self%csp (grid%nedgex,self%ncspmx))
allocate(self%cep (grid%nvertx,self%ncepmx))
allocate(self%cvd (grid%nvertx,self%ncvdmx))
allocate(self%csd (grid%nedgex,self%ncsdmx))

end subroutine setup

! --------------------------------------------------------------------------------------------------

subroutine pdelete(self)

implicit none
class(fempsoprs), intent(inout) :: self

! Partial deallocate after build has been run

if(allocated(self%eofvin    )) deallocate(self%eofvin    )
if(allocated(self%nrsten    )) deallocate(self%nrsten    )
if(allocated(self%nrxsten   )) deallocate(self%nrxsten   )
if(allocated(self%nwsten    )) deallocate(self%nwsten    )
if(allocated(self%ntsten    )) deallocate(self%ntsten    )
if(allocated(self%rsten     )) deallocate(self%rsten     )
if(allocated(self%rxsten    )) deallocate(self%rxsten    )
if(allocated(self%wsten     )) deallocate(self%wsten     )
if(allocated(self%tsten     )) deallocate(self%tsten     )
if(allocated(self%rcoeff    )) deallocate(self%rcoeff    )
if(allocated(self%rxcoeff   )) deallocate(self%rxcoeff   )
if(allocated(self%wcoeff    )) deallocate(self%wcoeff    )
if(allocated(self%tcoeff    )) deallocate(self%tcoeff    )
if(allocated(self%ncvp      )) deallocate(self%ncvp      )
if(allocated(self%ncsp      )) deallocate(self%ncsp      )
if(allocated(self%ncep      )) deallocate(self%ncep      )
if(allocated(self%ncvd      )) deallocate(self%ncvd      )
if(allocated(self%ncsd      )) deallocate(self%ncsd      )
if(allocated(self%cvp       )) deallocate(self%cvp       )
if(allocated(self%csp       )) deallocate(self%csp       )
if(allocated(self%cep       )) deallocate(self%cep       )
if(allocated(self%cvd       )) deallocate(self%cvd       )
if(allocated(self%csd       )) deallocate(self%csd       )

end subroutine pdelete

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fempsoprs), intent(inout) :: self

if(allocated(self%varea     )) deallocate(self%varea     )
if(allocated(self%eoffin    )) deallocate(self%eoffin    )
if(allocated(self%eofvin    )) deallocate(self%eofvin    )
if(allocated(self%nlsten    )) deallocate(self%nlsten    )
if(allocated(self%nmsten    )) deallocate(self%nmsten    )
if(allocated(self%njsten    )) deallocate(self%njsten    )
if(allocated(self%nhsten    )) deallocate(self%nhsten    )
if(allocated(self%nrsten    )) deallocate(self%nrsten    )
if(allocated(self%nrxsten   )) deallocate(self%nrxsten   )
if(allocated(self%nwsten    )) deallocate(self%nwsten    )
if(allocated(self%ntsten    )) deallocate(self%ntsten    )
if(allocated(self%nxminvsten)) deallocate(self%nxminvsten)
if(allocated(self%lsten     )) deallocate(self%lsten     )
if(allocated(self%msten     )) deallocate(self%msten     )
if(allocated(self%jsten     )) deallocate(self%jsten     )
if(allocated(self%hsten     )) deallocate(self%hsten     )
if(allocated(self%rsten     )) deallocate(self%rsten     )
if(allocated(self%rxsten    )) deallocate(self%rxsten    )
if(allocated(self%wsten     )) deallocate(self%wsten     )
if(allocated(self%tsten     )) deallocate(self%tsten     )
if(allocated(self%xminvsten )) deallocate(self%xminvsten )
if(allocated(self%jlump     )) deallocate(self%jlump     )
if(allocated(self%mlump     )) deallocate(self%mlump     )
if(allocated(self%hlump     )) deallocate(self%hlump     )
if(allocated(self%lmass     )) deallocate(self%lmass     )
if(allocated(self%mmass     )) deallocate(self%mmass     )
if(allocated(self%jstar     )) deallocate(self%jstar     )
if(allocated(self%hstar     )) deallocate(self%hstar     )
if(allocated(self%rcoeff    )) deallocate(self%rcoeff    )
if(allocated(self%rxcoeff   )) deallocate(self%rxcoeff   )
if(allocated(self%wcoeff    )) deallocate(self%wcoeff    )
if(allocated(self%xminv     )) deallocate(self%xminv     )
if(allocated(self%tcoeff    )) deallocate(self%tcoeff    )
if(allocated(self%elong     )) deallocate(self%elong     )
if(allocated(self%elat      )) deallocate(self%elat      )
if(allocated(self%ninj      )) deallocate(self%ninj      )
if(allocated(self%injsten   )) deallocate(self%injsten   )
if(allocated(self%injwgt    )) deallocate(self%injwgt    )
if(allocated(self%lapdiag   )) deallocate(self%lapdiag   )
if(allocated(self%underrel  )) deallocate(self%underrel  )
if(allocated(self%ncvp      )) deallocate(self%ncvp      )
if(allocated(self%ncsp      )) deallocate(self%ncsp      )
if(allocated(self%ncep      )) deallocate(self%ncep      )
if(allocated(self%ncvd      )) deallocate(self%ncvd      )
if(allocated(self%ncsd      )) deallocate(self%ncsd      )
if(allocated(self%cvp       )) deallocate(self%cvp       )
if(allocated(self%csp       )) deallocate(self%csp       )
if(allocated(self%cep       )) deallocate(self%cep       )
if(allocated(self%cvd       )) deallocate(self%cvd       )
if(allocated(self%csd       )) deallocate(self%csd       )

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine writeoperators(self,grid,filename)

implicit none
class(fempsoprs), intent(in) :: self
type(fempsgrid),  intent(in) :: grid
character(len=*), intent(in) :: filename

integer :: ncid, vc, varid(1000)
integer :: nvertx_dimid, ngrids_dimid, nfacex_dimid, nedgex_dimid,    nefmx_dimid, &
           nevmx_dimid,   nlsmx_dimid,  nmsmx_dimid,  njsmx_dimid,    nhsmx_dimid,  nrsmx_dimid, &
           nrxsmx_dimid,  nwsmx_dimid,  ntsmx_dimid, nxmisx_dimid, ngridsm1_dimid, ninjmx_dimid, &
           ncvpmx_dimid, ncspmx_dimid, ncepmx_dimid, ncvdmx_dimid,   ncsdmx_dimid

! Create file
! -----------
call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_CLOBBER), ncid), "nf90_create" )

! Define dimensions
! -----------------
call nccheck( nf90_def_dim(ncid, "nvertx"  , grid%nvertx  ,   nvertx_dimid), "nf90_def_dim nvertx"   )
call nccheck( nf90_def_dim(ncid, "ngrids"  , grid%ngrids  ,   ngrids_dimid), "nf90_def_dim ngrids"   )
call nccheck( nf90_def_dim(ncid, "nfacex"  , grid%nfacex  ,   nfacex_dimid), "nf90_def_dim nfacex"   )
call nccheck( nf90_def_dim(ncid, "nedgex"  , grid%nedgex  ,   nedgex_dimid), "nf90_def_dim nedgex"   )
call nccheck( nf90_def_dim(ncid, "nefmx"   , grid%nefmx   ,    nefmx_dimid), "nf90_def_dim nefmx"    )
call nccheck( nf90_def_dim(ncid, "nevmx"   , grid%nevmx   ,    nevmx_dimid), "nf90_def_dim nevmx"    )
call nccheck( nf90_def_dim(ncid, "nlsmx"   , self%nlsmx   ,    nlsmx_dimid), "nf90_def_dim nlsmx"    )
call nccheck( nf90_def_dim(ncid, "nmsmx"   , self%nmsmx   ,    nmsmx_dimid), "nf90_def_dim nmsmx"    )
call nccheck( nf90_def_dim(ncid, "njsmx"   , self%njsmx   ,    njsmx_dimid), "nf90_def_dim njsmx"    )
call nccheck( nf90_def_dim(ncid, "nhsmx"   , self%nhsmx   ,    nhsmx_dimid), "nf90_def_dim nhsmx"    )
call nccheck( nf90_def_dim(ncid, "nrsmx"   , self%nrsmx   ,    nrsmx_dimid), "nf90_def_dim nrsmx"    )
call nccheck( nf90_def_dim(ncid, "nrxsmx"  , self%nrxsmx  ,   nrxsmx_dimid), "nf90_def_dim nrxsmx"   )
call nccheck( nf90_def_dim(ncid, "nwsmx"   , self%nwsmx   ,    nwsmx_dimid), "nf90_def_dim nwsmx"    )
call nccheck( nf90_def_dim(ncid, "ntsmx"   , self%ntsmx   ,    ntsmx_dimid), "nf90_def_dim ntsmx"    )
call nccheck( nf90_def_dim(ncid, "nxmisx"  , self%nxmisx  ,   nxmisx_dimid), "nf90_def_dim nxmisx"   )
call nccheck( nf90_def_dim(ncid, "ngridsm1", grid%ngrids-1, ngridsm1_dimid), "nf90_def_dim ngridsm1" )
call nccheck( nf90_def_dim(ncid, "ninjmx"  , self%ninjmx  ,   ninjmx_dimid), "nf90_def_dim ninjmx"   )
call nccheck( nf90_def_dim(ncid, "ncvpmx"  , self%ncvpmx  ,   ncvpmx_dimid), "nf90_def_dim ncvpmx"   )
call nccheck( nf90_def_dim(ncid, "ncspmx"  , self%ncspmx  ,   ncspmx_dimid), "nf90_def_dim ncspmx"   )
call nccheck( nf90_def_dim(ncid, "ncepmx"  , self%ncepmx  ,   ncepmx_dimid), "nf90_def_dim ncepmx"   )
call nccheck( nf90_def_dim(ncid, "ncvdmx"  , self%ncvdmx  ,   ncvdmx_dimid), "nf90_def_dim ncvdmx"   )
call nccheck( nf90_def_dim(ncid, "ncsdmx"  , self%ncsdmx  ,   ncsdmx_dimid), "nf90_def_dim ncsdmx"   )

! Define variables
! ----------------
vc = 1

if (allocated(self%varea     )) call nccheck( nf90_def_var(ncid, "varea",      NF90_FLOAT, &
    (/nvertx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var varea" ); vc=vc+1
if (allocated(self%eoffin    )) call nccheck( nf90_def_var(ncid, "eoffin",     NF90_INT  , &
    (/nfacex_dimid,nefmx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var eoffin" ); vc=vc+1
if (allocated(self%eofvin    )) call nccheck( nf90_def_var(ncid, "eofvin",     NF90_INT  , &
    (/nvertx_dimid,nevmx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var eofvin" ); vc=vc+1
if (allocated(self%nlsten    )) call nccheck( nf90_def_var(ncid, "nlsten",     NF90_INT  , &
    (/nfacex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nlsten" ); vc=vc+1
if (allocated(self%nmsten    )) call nccheck( nf90_def_var(ncid, "nmsten",     NF90_INT  , &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nmsten" ); vc=vc+1
if (allocated(self%njsten    )) call nccheck( nf90_def_var(ncid, "njsten",     NF90_INT  , &
    (/nvertx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var njsten" ); vc=vc+1
if (allocated(self%nhsten    )) call nccheck( nf90_def_var(ncid, "nhsten",     NF90_INT  , &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nhsten" ); vc=vc+1
if (allocated(self%nrsten    )) call nccheck( nf90_def_var(ncid, "nrsten",     NF90_INT  , &
    (/nfacex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nrsten" ); vc=vc+1
if (allocated(self%nrxsten   )) call nccheck( nf90_def_var(ncid, "nrxsten",    NF90_INT  , &
    (/nvertx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nrxsten" ); vc=vc+1
if (allocated(self%nwsten    )) call nccheck( nf90_def_var(ncid, "nwsten",     NF90_INT  , &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nwsten" ); vc=vc+1
if (allocated(self%ntsten    )) call nccheck( nf90_def_var(ncid, "ntsten",     NF90_INT  , &
    (/nfacex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var ntsten" ); vc=vc+1
if (allocated(self%nxminvsten)) call nccheck( nf90_def_var(ncid, "nxminvsten", NF90_INT  , &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var nxminvsten" ); vc=vc+1
if (allocated(self%lsten     )) call nccheck( nf90_def_var(ncid, "lsten",      NF90_INT  , &
    (/nfacex_dimid,nlsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var lsten" ); vc=vc+1
if (allocated(self%msten     )) call nccheck( nf90_def_var(ncid, "msten",      NF90_INT  , &
    (/nedgex_dimid,nmsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var msten" ); vc=vc+1
if (allocated(self%jsten     )) call nccheck( nf90_def_var(ncid, "jsten",      NF90_INT  , &
    (/nvertx_dimid,njsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var jsten" ); vc=vc+1
if (allocated(self%hsten     )) call nccheck( nf90_def_var(ncid, "hsten",      NF90_INT  , &
    (/nedgex_dimid,nhsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var hsten" ); vc=vc+1
if (allocated(self%rsten     )) call nccheck( nf90_def_var(ncid, "rsten",      NF90_INT  , &
    (/nfacex_dimid,nrsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var rsten" ); vc=vc+1
if (allocated(self%rxsten    )) call nccheck( nf90_def_var(ncid, "rxsten",     NF90_INT  , &
    (/nvertx_dimid,nrxsmx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var rxsten" ); vc=vc+1
if (allocated(self%wsten     )) call nccheck( nf90_def_var(ncid, "wsten",      NF90_INT  , &
    (/nedgex_dimid,nwsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var wsten" ); vc=vc+1
if (allocated(self%tsten     )) call nccheck( nf90_def_var(ncid, "tsten",      NF90_INT  , &
    (/nfacex_dimid,ntsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var tsten" ); vc=vc+1
if (allocated(self%xminvsten )) call nccheck( nf90_def_var(ncid, "xminvsten",  NF90_INT  , &
    (/nedgex_dimid,nxmisx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var xminvsten" ); vc=vc+1
if (allocated(self%jlump     )) call nccheck( nf90_def_var(ncid, "jlump",      NF90_FLOAT, &
    (/nvertx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var jlump" ); vc=vc+1
if (allocated(self%mlump     )) call nccheck( nf90_def_var(ncid, "mlump",      NF90_FLOAT, &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var mlump" ); vc=vc+1
if (allocated(self%hlump     )) call nccheck( nf90_def_var(ncid, "hlump",      NF90_FLOAT, &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var hlump" ); vc=vc+1
if (allocated(self%lmass     )) call nccheck( nf90_def_var(ncid, "lmass",      NF90_FLOAT, &
    (/nfacex_dimid,nlsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var lmass" ); vc=vc+1
if (allocated(self%mmass     )) call nccheck( nf90_def_var(ncid, "mmass",      NF90_FLOAT, &
    (/nedgex_dimid,nmsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var mmass" ); vc=vc+1
if (allocated(self%jstar     )) call nccheck( nf90_def_var(ncid, "jstar",      NF90_FLOAT, &
    (/nvertx_dimid,njsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var jstar" ); vc=vc+1
if (allocated(self%hstar     )) call nccheck( nf90_def_var(ncid, "hstar",      NF90_FLOAT, &
    (/nedgex_dimid,nhsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var hstar" ); vc=vc+1
if (allocated(self%rcoeff    )) call nccheck( nf90_def_var(ncid, "rcoeff",     NF90_FLOAT, &
    (/nfacex_dimid,nrsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var rcoeff" ); vc=vc+1
if (allocated(self%rxcoeff   )) call nccheck( nf90_def_var(ncid, "rxcoeff",    NF90_FLOAT, &
    (/nvertx_dimid,nrxsmx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var rxcoeff" ); vc=vc+1
if (allocated(self%wcoeff    )) call nccheck( nf90_def_var(ncid, "wcoeff",     NF90_FLOAT, &
    (/nedgex_dimid,nwsmx_dimid, ngrids_dimid/), varid(vc)), "nf90_def_var wcoeff" ); vc=vc+1
if (allocated(self%xminv     )) call nccheck( nf90_def_var(ncid, "xminv",      NF90_FLOAT, &
    (/nedgex_dimid,nxmisx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var xminv" ); vc=vc+1
if (allocated(self%tcoeff    )) call nccheck( nf90_def_var(ncid, "tcoeff",     NF90_FLOAT, &
    (/nfacex_dimid,ntsmx_dimid, ntsmx_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var tcoeff" ); vc=vc+1
if (allocated(self%elong     )) call nccheck( nf90_def_var(ncid, "elong",      NF90_FLOAT, &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var elong" ); vc=vc+1
if (allocated(self%elat      )) call nccheck( nf90_def_var(ncid, "elat",       NF90_FLOAT, &
    (/nedgex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var elat" ); vc=vc+1
if (allocated(self%ninj      )) call nccheck( nf90_def_var(ncid, "ninj",       NF90_INT  , &
    (/nfacex_dimid,ngridsm1_dimid/), varid(vc)), "nf90_def_var ninj" ); vc=vc+1
if (allocated(self%injsten   )) call nccheck( nf90_def_var(ncid, "injsten",    NF90_INT  , &
    (/nfacex_dimid,ninjmx_dimid,ngridsm1_dimid/), varid(vc)), "nf90_def_var injsten" ); vc=vc+1
if (allocated(self%injwgt    )) call nccheck( nf90_def_var(ncid, "injwgt",     NF90_FLOAT, &
    (/nfacex_dimid,ninjmx_dimid,ngridsm1_dimid/), varid(vc)), "nf90_def_var injwgt" ); vc=vc+1
if (allocated(self%lapdiag   )) call nccheck( nf90_def_var(ncid, "lapdiag",    NF90_FLOAT, &
    (/nfacex_dimid,ngrids_dimid/), varid(vc)), "nf90_def_var lapdiag" ); vc=vc+1
if (allocated(self%underrel  )) call nccheck( nf90_def_var(ncid, "underrel",   NF90_FLOAT, &
    (/ngrids_dimid/), varid(vc)), "nf90_def_var underrel" ); vc=vc+1
if (allocated(self%ncvp      )) call nccheck( nf90_def_var(ncid, "ncvp",       NF90_INT  , &
    (/nfacex_dimid/), varid(vc)), "nf90_def_var ncvp" ); vc=vc+1
if (allocated(self%ncsp      )) call nccheck( nf90_def_var(ncid, "ncsp",       NF90_INT  , &
    (/nedgex_dimid/), varid(vc)), "nf90_def_var ncsp" ); vc=vc+1
if (allocated(self%ncep      )) call nccheck( nf90_def_var(ncid, "ncep",       NF90_INT  , &
    (/nvertx_dimid/), varid(vc)), "nf90_def_var ncep" ); vc=vc+1
if (allocated(self%ncvd      )) call nccheck( nf90_def_var(ncid, "ncvd",       NF90_INT  , &
    (/nvertx_dimid/), varid(vc)), "nf90_def_var ncvd" ); vc=vc+1
if (allocated(self%ncsd      )) call nccheck( nf90_def_var(ncid, "ncsd",       NF90_INT  , &
    (/nedgex_dimid/), varid(vc)), "nf90_def_var ncsd" ); vc=vc+1
if (allocated(self%cvp       )) call nccheck( nf90_def_var(ncid, "cvp",        NF90_FLOAT, &
    (/nfacex_dimid,ncvpmx_dimid/), varid(vc)), "nf90_def_var cvp" ); vc=vc+1
if (allocated(self%csp       )) call nccheck( nf90_def_var(ncid, "csp",        NF90_FLOAT, &
    (/nedgex_dimid,ncspmx_dimid/), varid(vc)), "nf90_def_var csp" ); vc=vc+1
if (allocated(self%cep       )) call nccheck( nf90_def_var(ncid, "cep",        NF90_FLOAT, &
    (/nvertx_dimid,ncepmx_dimid/), varid(vc)), "nf90_def_var cep" ); vc=vc+1
if (allocated(self%cvd       )) call nccheck( nf90_def_var(ncid, "cvd",        NF90_FLOAT, &
    (/nvertx_dimid,ncvdmx_dimid/), varid(vc)), "nf90_def_var cvd" ); vc=vc+1
if (allocated(self%csd       )) call nccheck( nf90_def_var(ncid, "csd",        NF90_FLOAT, &
    (/nedgex_dimid,ncsdmx_dimid/), varid(vc)), "nf90_def_var csd" ); vc=vc+1

! Write variables
! ---------------
vc = 1

if (allocated(self%varea     )) call nccheck( nf90_put_var( ncid, varid(vc), self%varea     ), "nf90_put_var varea"      ); vc=vc+1
if (allocated(self%eoffin    )) call nccheck( nf90_put_var( ncid, varid(vc), self%eoffin    ), "nf90_put_var eoffin"     ); vc=vc+1
if (allocated(self%eofvin    )) call nccheck( nf90_put_var( ncid, varid(vc), self%eofvin    ), "nf90_put_var eofvin"     ); vc=vc+1
if (allocated(self%nlsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%nlsten    ), "nf90_put_var nlsten"     ); vc=vc+1
if (allocated(self%nmsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%nmsten    ), "nf90_put_var nmsten"     ); vc=vc+1
if (allocated(self%njsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%njsten    ), "nf90_put_var njsten"     ); vc=vc+1
if (allocated(self%nhsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%nhsten    ), "nf90_put_var nhsten"     ); vc=vc+1
if (allocated(self%nrsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%nrsten    ), "nf90_put_var nrsten"     ); vc=vc+1
if (allocated(self%nrxsten   )) call nccheck( nf90_put_var( ncid, varid(vc), self%nrxsten   ), "nf90_put_var nrxsten"    ); vc=vc+1
if (allocated(self%nwsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%nwsten    ), "nf90_put_var nwsten"     ); vc=vc+1
if (allocated(self%ntsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%ntsten    ), "nf90_put_var ntsten"     ); vc=vc+1
if (allocated(self%nxminvsten)) call nccheck( nf90_put_var( ncid, varid(vc), self%nxminvsten), "nf90_put_var nxminvsten" ); vc=vc+1
if (allocated(self%lsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%lsten     ), "nf90_put_var lsten"      ); vc=vc+1
if (allocated(self%msten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%msten     ), "nf90_put_var msten"      ); vc=vc+1
if (allocated(self%jsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%jsten     ), "nf90_put_var jsten"      ); vc=vc+1
if (allocated(self%hsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%hsten     ), "nf90_put_var hsten"      ); vc=vc+1
if (allocated(self%rsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%rsten     ), "nf90_put_var rsten"      ); vc=vc+1
if (allocated(self%rxsten    )) call nccheck( nf90_put_var( ncid, varid(vc), self%rxsten    ), "nf90_put_var rxsten"     ); vc=vc+1
if (allocated(self%wsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%wsten     ), "nf90_put_var wsten"      ); vc=vc+1
if (allocated(self%tsten     )) call nccheck( nf90_put_var( ncid, varid(vc), self%tsten     ), "nf90_put_var tsten"      ); vc=vc+1
if (allocated(self%xminvsten )) call nccheck( nf90_put_var( ncid, varid(vc), self%xminvsten ), "nf90_put_var xminvsten"  ); vc=vc+1
if (allocated(self%jlump     )) call nccheck( nf90_put_var( ncid, varid(vc), self%jlump     ), "nf90_put_var jlump"      ); vc=vc+1
if (allocated(self%mlump     )) call nccheck( nf90_put_var( ncid, varid(vc), self%mlump     ), "nf90_put_var mlump"      ); vc=vc+1
if (allocated(self%hlump     )) call nccheck( nf90_put_var( ncid, varid(vc), self%hlump     ), "nf90_put_var hlump"      ); vc=vc+1
if (allocated(self%lmass     )) call nccheck( nf90_put_var( ncid, varid(vc), self%lmass     ), "nf90_put_var lmass"      ); vc=vc+1
if (allocated(self%mmass     )) call nccheck( nf90_put_var( ncid, varid(vc), self%mmass     ), "nf90_put_var mmass"      ); vc=vc+1
if (allocated(self%jstar     )) call nccheck( nf90_put_var( ncid, varid(vc), self%jstar     ), "nf90_put_var jstar"      ); vc=vc+1
if (allocated(self%hstar     )) call nccheck( nf90_put_var( ncid, varid(vc), self%hstar     ), "nf90_put_var hstar"      ); vc=vc+1
if (allocated(self%rcoeff    )) call nccheck( nf90_put_var( ncid, varid(vc), self%rcoeff    ), "nf90_put_var rcoeff"     ); vc=vc+1
if (allocated(self%rxcoeff   )) call nccheck( nf90_put_var( ncid, varid(vc), self%rxcoeff   ), "nf90_put_var rxcoeff"    ); vc=vc+1
if (allocated(self%wcoeff    )) call nccheck( nf90_put_var( ncid, varid(vc), self%wcoeff    ), "nf90_put_var wcoeff"     ); vc=vc+1
if (allocated(self%xminv     )) call nccheck( nf90_put_var( ncid, varid(vc), self%xminv     ), "nf90_put_var xminv"      ); vc=vc+1
if (allocated(self%tcoeff    )) call nccheck( nf90_put_var( ncid, varid(vc), self%tcoeff    ), "nf90_put_var tcoeff"     ); vc=vc+1
if (allocated(self%elong     )) call nccheck( nf90_put_var( ncid, varid(vc), self%elong     ), "nf90_put_var elong"      ); vc=vc+1
if (allocated(self%elat      )) call nccheck( nf90_put_var( ncid, varid(vc), self%elat      ), "nf90_put_var elat"       ); vc=vc+1
if (allocated(self%ninj      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ninj      ), "nf90_put_var ninj"       ); vc=vc+1
if (allocated(self%injsten   )) call nccheck( nf90_put_var( ncid, varid(vc), self%injsten   ), "nf90_put_var injsten"    ); vc=vc+1
if (allocated(self%injwgt    )) call nccheck( nf90_put_var( ncid, varid(vc), self%injwgt    ), "nf90_put_var injwgt"     ); vc=vc+1
if (allocated(self%lapdiag   )) call nccheck( nf90_put_var( ncid, varid(vc), self%lapdiag   ), "nf90_put_var lapdiag"    ); vc=vc+1
if (allocated(self%underrel  )) call nccheck( nf90_put_var( ncid, varid(vc), self%underrel  ), "nf90_put_var underrel"   ); vc=vc+1
if (allocated(self%ncvp      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ncvp      ), "nf90_put_var ncvp"       ); vc=vc+1
if (allocated(self%ncsp      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ncsp      ), "nf90_put_var ncsp"       ); vc=vc+1
if (allocated(self%ncep      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ncep      ), "nf90_put_var ncep"       ); vc=vc+1
if (allocated(self%ncvd      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ncvd      ), "nf90_put_var ncvd"       ); vc=vc+1
if (allocated(self%ncsd      )) call nccheck( nf90_put_var( ncid, varid(vc), self%ncsd      ), "nf90_put_var ncsd"       ); vc=vc+1
if (allocated(self%cvp       )) call nccheck( nf90_put_var( ncid, varid(vc), self%cvp       ), "nf90_put_var cvp"        ); vc=vc+1
if (allocated(self%csp       )) call nccheck( nf90_put_var( ncid, varid(vc), self%csp       ), "nf90_put_var csp"        ); vc=vc+1
if (allocated(self%cep       )) call nccheck( nf90_put_var( ncid, varid(vc), self%cep       ), "nf90_put_var cep"        ); vc=vc+1
if (allocated(self%cvd       )) call nccheck( nf90_put_var( ncid, varid(vc), self%cvd       ), "nf90_put_var cvd"        ); vc=vc+1
if (allocated(self%csd       )) call nccheck( nf90_put_var( ncid, varid(vc), self%csd       ), "nf90_put_var csd"        ); vc=vc+1

! Close file
! ----------
call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine writeoperators

! --------------------------------------------------------------------------------------------------

subroutine build(self,grid)

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(inout) :: grid

integer :: igrid
character(len=255) :: infomessage

! Construct some basic geometrical information
call buildgeom(self,grid)
call message('Done buildgeom',trace)

! Information is assumed to be in a certain order in arrays
! defining connectivity. this may already hold for the output
! of the grid generator, but check and correct if necessary.
call ordering(self,grid)
call message('Done ordering',trace)

! Loop over hierarchy of grids
do igrid = 1, grid%ngrids

  call message(' ',trace)
  write(infomessage,*) 'Grid ',igrid
  call message(infomessage,trace)

  ! Build compound elements on grid igrid
  call buildcomp(self,grid,igrid)
  call message('Done buildcomp',trace)

  ! Build matrix operators on grid igrid
  call buildmat(self,grid,igrid)
  call message('Done buildmat',trace)

enddo
call message(' ',trace)

! Build restriction operator
call buildinj(self,grid)
call message('Multigrid restriction operator built',trace)

end subroutine build

! --------------------------------------------------------------------------------------------------

subroutine ordering(self,grid)

! Some of the later routines assume the grid connectivity data
! obeys certain ordering conventions. Some of these may be
! imposed by the grid generation program, but check and fix
! here if necessary.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(inout) :: grid

integer :: igrid, if0, ne, ix1, if1, ixmin, ifmin, ix2, if2, &
           if3, iv0, ie1, if21, if22, iv11, iv12, iv21, iv22, &
           ie2, ie3, iemin, if11, if12, ie0, iv1, iv2
real(kind=kind_real) :: pi, long, lat, x0(3), x1(3), d1x, d1y, d1z, thetamin, &
                        x2(3), d2x, d2y, d2z, cs, sn, theta
logical :: lfound

character(len=255) :: errormessage

! Note: this method updates:
! grid%fnxtf
! grid%eoff
! grid%voff
! grid%eofv
! grid%fofv
! grid%vofe

pi = 4.0_kind_real*atan(1.0_kind_real)

! Loop over all grids in the hierarchy
do igrid = 1, grid%ngrids

  ! Loop over all faces on this grid
  do if0 = 1, grid%nface(igrid)

    ! Coordinates of face if0
    long = grid%flong(if0,igrid)
    lat = grid%flat(if0,igrid)
    call ll2xyz(long,lat,x0)

    ! Number of edges/neighbours
    ne = grid%neoff(if0,igrid)

    ! Sort FNXTF into anticlockwise order
    do ix1 = 1, ne - 2

      ! Coordinates of IX1'th neighbour
      if1 = grid%fnxtf(if0,ix1,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      call ll2xyz(long,lat,x1)
      d1x = x1(1) - x0(1)
      d1y = x1(2) - x0(2)
      d1z = x1(3) - x0(3)

      ! Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0

      do ix2 = ix1 + 1, ne

        ! Coordinates of IX2'th neighbour
        if2 = grid%fnxtf(if0,ix2,igrid)
        long = grid%flong(if2,igrid)
        lat = grid%flat(if2,igrid)
        call ll2xyz(long,lat,x2)
        d2x=x2(1) - x0(1)
        d2y=x2(2) - x0(2)
        d2z=x2(3) - x0(3)
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0(1)*(d1y*d2z - d1z*d2y) &
           + x0(2)*(d1z*d2x - d1x*d2z) &
           + x0(3)*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        if ((theta < thetamin) .and. (theta > 0.0_kind_real)) then
          ixmin = ix2
          ifmin = if2
          thetamin = theta
        endif
      enddo

!     The face in position IXMIN belongs in position IX1+1 so swap them
      if3 = grid%fnxtf(if0,ix1+1,igrid)
      grid%fnxtf(if0,ix1+1,igrid) = ifmin
      grid%fnxtf(if0,ixmin,igrid) = if3

    enddo

    ! Now sort EOFF to correspond to FNXTF
    do ix1 = 1, ne

      if1 = grid%fnxtf(if0,ix1,igrid)
      ix2 = ix1 - 1
      lfound = .false.
      do while (.not. lfound)
         ix2 = ix2 + 1
         ie1 = grid%eoff(if0,ix2,igrid)
         if21 = grid%fnxte(ie1,1,igrid)
         if22 = grid%fnxte(ie1,2,igrid)
         if ((if21 + if22) == (if0 + if1)) lfound = .true.
      enddo

      ! Edge IE2 corresponds to face IF1
      grid%eoff(if0,ix2,igrid) = grid%eoff(if0,ix1,igrid)
      grid%eoff(if0,ix1,igrid) = ie1

    enddo

    ! Order VOFF so that the k'th vertex is between the
    ! k'th and (k+1)'th edges in EOFF
    do ix1 = 1, ne
      ix2 = ix1 + 1
      if (ix2 > ne) ix2 = 1
      ie1 = grid%eoff(if0,ix1,igrid)
      ie2 = grid%eoff(if0,ix2,igrid)
      ! Find the common vertex of IE1 and IE2
      iv11 = grid%vofe(ie1,1,igrid)
      iv12 = grid%vofe(ie1,2,igrid)
      iv21 = grid%vofe(ie2,1,igrid)
      iv22 = grid%vofe(ie2,2,igrid)
      if ((iv11 == iv21) .or. (iv11 == iv22)) then
        iv0 = iv11
      elseif ((iv12 == iv21) .or. (iv12 == iv22)) then
        iv0 = iv12
      else
        call message('common vertex not found',fatal)
      endif
      grid%voff(if0,ix1,igrid) = iv0
    enddo

  enddo

  ! Loop over all vertices on this grid
  do iv0 = 1, grid%nvert(igrid)

    ! Coordinates of vertex iv0
    long = grid%vlong(iv0,igrid)
    lat = grid%vlat(iv0,igrid)
    call ll2xyz(long,lat,x0)

    ! Number of edges / adjacent faces
    ne = grid%neofv(iv0,igrid)

    ! Sort EOFV into anticlockwise order
    do ix1 = 1, ne - 2

      ! Coordinates of IX1'th edge
      ie1 = grid%eofv(iv0,ix1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x1)
      d1x = x1(1) - x0(1)
      d1y = x1(2) - x0(2)
      d1z = x1(3) - x0(3)

!     Find next neighbour (anticlockwise)
      thetamin = pi
      ixmin = 0
      ifmin = 0

      do ix2 = ix1 + 1, ne

        ! Coordinates of IX2'th neighbour
        ie2 = grid%eofv(iv0,ix2,igrid)
        long = self%elong(ie2,igrid)
        lat = self%elat(ie2,igrid)
        call ll2xyz(long,lat,x2)
        d2x=x2(1) - x0(1)
        d2y=x2(2) - x0(2)
        d2z=x2(3) - x0(3)
        cs = d1x*d2x + d1y*d2y + d1z*d2z
        sn = x0(1)*(d1y*d2z - d1z*d2y) &
           + x0(2)*(d1z*d2x - d1x*d2z) &
           + x0(3)*(d1x*d2y - d1y*d2x)
        theta = ATAN2(sn,cs)
        if ((theta < thetamin) .AND. (theta > 0.0_kind_real)) THEN
          ixmin = ix2
          iemin = ie2
          thetamin = theta
        endif
      enddo

!     The edge in position IXMIN belongs in position IX1+1 so swap them
      ie3 = grid%eofv(iv0,ix1+1,igrid)
      grid%eofv(iv0,ix1+1,igrid) = iemin
      grid%eofv(iv0,ixmin,igrid) = ie3

    enddo

    ! Order FOFV so that the k'th face is between the
    ! k'th and (k+1)'th edges in EOFV
    do ix1 = 1, ne
      ix2 = ix1 + 1
      if (ix2 > ne) ix2 = 1
      ie1 = grid%eofv(iv0,ix1,igrid)
      ie2 = grid%eofv(iv0,ix2,igrid)
      ! Find the common face of IE1 and IE2
      if11 = grid%fnxte(ie1,1,igrid)
      if12 = grid%fnxte(ie1,2,igrid)
      if21 = grid%fnxte(ie2,1,igrid)
      if22 = grid%fnxte(ie2,2,igrid)
      if ((if11 == if21) .or. (if11 == if22)) then
        if0 = if11
      elseif ((if12 == if21) .or. (if12 == if22)) then
        if0 = if12
      else
        write(errormessage,*) 'grid ',igrid,' vertex ',iv0
        call message('common face not found')
        call message(errormessage,fatal)
      endif
      grid%fofv(iv0,ix1,igrid) = if0
    enddo

  enddo

  ! Loop over all edges on this grid
  do ie0 = 1, grid%nedge(igrid)

    ! Sort VOFE so that VOFE(1) -> VOFE(2) (tangent vector)
    ! is 90 degrees anticlockwise of FNXTE(1) -> FNXTE(2) (normal vector)
    if1 = grid%fnxte(ie0,1,igrid)
    if2 = grid%fnxte(ie0,2,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x0)
    long = grid%flong(if2,igrid)
    lat = grid%flat(if2,igrid)
    call ll2xyz(long,lat,x1)
    d1x = x1(1) - x0(1)
    d1y = x1(2) - x0(2)
    d1z = x1(3) - x0(3)
    iv1 = grid%vofe(ie0,1,igrid)
    iv2 = grid%vofe(ie0,2,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x0)
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    call ll2xyz(long,lat,x1)
    d2x = x1(1) - x0(1)
    d2y = x1(2) - x0(2)
    d2z = x1(3) - x0(3)
    sn = x0(1)*(d1y*d2z - d1z*d2y) &
       + x0(2)*(d1z*d2x - d1x*d2z) &
       + x0(3)*(d1x*d2y - d1y*d2x)
    if (sn < 0.0_kind_real) THEN
      ! Swap the two vertices
      grid%vofe(ie0,1,igrid) = iv2
      grid%vofe(ie0,2,igrid) = iv1
    endif

  enddo

enddo

end subroutine ordering

! --------------------------------------------------------------------------------------------------

subroutine buildgeom(self,grid)

! Build some basic geometrical information:
! locations of edge crossing points;
! lengths of primal and dual edges;
! cell areas.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(inout) :: grid

integer :: igrid, ie0, if1, if2, iv1, iv2, if0, ie1, ix
real(kind=kind_real) :: long, lat, x1(3), x2(3), y1(3), y2(3), &
                        n1(3), n2(3), r1(3), mag, l1sq, l2sq, l3sq, &
                        area, a

! Note: this method updates:
! grid%ddist
! grid%ldist
! grid%farea

! Loop over grids
do igrid = 1, grid%ngrids

  ! Loop over edges
  do ie0 = 1, grid%nedge(igrid)

    ! Locate cell centres either side (i.e dual vertices)
    if1 = grid%fnxte(ie0,1,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)
    if2 = grid%fnxte(ie0,2,igrid)
    long = grid%flong(if2,igrid)
    lat = grid%flat(if2,igrid)
    call ll2xyz(long,lat,x2)

    ! And ends of edge (i.e. primal vertices)
    iv1 = grid%vofe(ie0,1,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,y1)
    iv2 = grid%vofe(ie0,2,igrid)
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    call ll2xyz(long,lat,y2)

    ! Normal to plane of dual edge
    n1(1) = x1(2)*x2(3) - x1(3)*x2(2)
    n1(2) = x1(3)*x2(1) - x1(1)*x2(3)
    n1(3) = x1(1)*x2(2) - x1(2)*x2(1)
    ! Normal to plane of primal edge
    n2(1) = y1(2)*y2(3) - y1(3)*y2(2)
    n2(2) = y1(3)*y2(1) - y1(1)*y2(3)
    n2(3) = y1(1)*y2(2) - y1(2)*y2(1)
    ! Hence radial vector of crossing point
    r1(1) = n1(2)*n2(3) - n1(3)*n2(2)
    r1(2) = n1(3)*n2(1) - n1(1)*n2(3)
    r1(3) = n1(1)*n2(2) - n1(2)*n2(1)
    ! Normalize to unit sphere
    mag = SQRT(r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3))
    r1 = r1/mag
    ! Convert back to lat long and save
    call xyz2ll(r1,long,lat)
    self%elong(ie0,igrid) = long
    self%elat(ie0,igrid) = lat

    ! Dual edge is now composed of two straight edge pieces
    l1sq = (r1(1) - x1(1))**2 + (r1(2) - x1(2))**2 + (r1(3) - x1(3))**2
    l2sq = (r1(1) - x2(1))**2 + (r1(2) - x2(2))**2 + (r1(3) - x2(3))**2
    grid%ddist(ie0,igrid) = SQRT(l1sq) + SQRT(l2sq)

    ! Primal edge is now composed of two straight edge pieces
    l1sq = (r1(1) - y1(1))**2 + (r1(2) - y1(2))**2 + (r1(3) - y1(3))**2
    l2sq = (r1(1) - y2(1))**2 + (r1(2) - y2(2))**2 + (r1(3) - y2(3))**2
    grid%ldist(ie0,igrid) = SQRT(l1sq) + SQRT(l2sq)

  enddo

  ! Intialize dual cell areas to zero and
  ! accumulate as we touch each vertex
  self%varea(:,igrid) = 0.0_kind_real
  ! Loop over primal cells
  do if0 = 1, grid%nface(igrid)

    ! Initialize area to zero and locate call centre
    area = 0.0_kind_real
    long = grid%flong(if0,igrid)
    lat = grid%flat(if0,igrid)
    call ll2xyz(long,lat,r1)
    ! Accumulate contributions to cell area from supermesh cells
    ! Loop over primal cell edges
    do ix = 1, grid%neoff(if0,igrid)
      ! There are two supermesh cells incident on each
      ! primal cell edge
      ie1 = grid%eoff(if0,ix,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x1)
      ! First supermesh cell
      iv1 = grid%vofe(ie1,1,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(r1,x1,x2,l1sq,l2sq,l3sq,a)
      area = area + a
      self%varea(iv1,igrid) = self%varea(iv1,igrid) + a
      ! Second supermesh cell
      iv1 = grid%vofe(ie1,2,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(r1,x1,x2,l1sq,l2sq,l3sq,a)
      area = area + a
      self%varea(iv1,igrid) = self%varea(iv1,igrid) + a
    enddo
    grid%farea(if0,igrid) = area

  enddo

enddo

! Construct the tables eoffin and eofvin
do igrid = 1, grid%ngrids

  do if1 = 1, grid%nface(igrid)
    do ix = 1, grid%neoff(if1,igrid)
      ie1 = grid%eoff(if1,ix,igrid)
      if2 = grid%fnxte(ie1,1,igrid)
      if (if1 == if2) THEN
        ! This edge points out of face if1
        self%eoffin(if1,ix,igrid) = -1.0_kind_real
      else
        ! This edge points into face if1
        self%eoffin(if1,ix,igrid) = 1.0_kind_real
      endif
    enddo
  enddo

  if (.not.self%poisson) then
    do iv1 = 1, grid%nvert(igrid)
      do ix = 1, grid%neofv(iv1,igrid)
        ie1 = grid%eofv(iv1,ix,igrid)
        iv2 = grid%vofe(ie1,1,igrid)
        if (iv1 == iv2) THEN
          ! This edge points away from vertex iv1
          self%eofvin(iv1,ix,igrid) = -1.0_kind_real
        else
          ! This edge points towards vertex iv1
          self%eofvin(iv1,ix,igrid) = 1.0_kind_real
        endif
      enddo
    enddo
  endif

enddo

end subroutine buildgeom

! --------------------------------------------------------------------------------------------------

subroutine buildcomp(self,grid,igrid)

! Construct the information about `internal' degrees of freedom
! needed to define the compound elements on grid igrid.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

! Compound P0 elements on primal cells (space Vp)
call buildvp(self,grid,igrid)
call message('  Done buildvp',trace)

! Compound RT0 elements on primal cells (space Sp)
call buildsp(self,grid,igrid)
call message('  Done buildsp',trace)

! Compound P1 elements on primal cells (space Ep)
call buildep(self,grid,igrid)
call message('  Done buildep',trace)

! Compound P0 elements on dual cells (space Vd)
call buildvd(self,grid,igrid)
call message('  Done buildvd',trace)

! Compound N0 elements on dual cells (space Sd)
call buildsd(self,grid,igrid)
call message('  Done buildsd',trace)

end subroutine buildcomp

! --------------------------------------------------------------------------------------------------

subroutine buildmat(self,grid,igrid)

! Build the matrix representations, in terms of stencils
! and ceofficients, of the various operators and mass matrices
! needed.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid), intent(in) :: grid
integer,         intent(in) :: igrid

! Mass matrix L
call buildL(self,grid,igrid)
call message('  Done buildL',trace)

! Mass matrix M
call buildM(self,grid,igrid)
call message('  Done buildM',trace)

! Vp to Ep transfer operator R
if (.not.self%poisson) then
  call buildR(self,grid,igrid)
  call message('  Done buildR',trace)
endif

! W operator for constructing perp of vectors
! (maps E_p tp E_p)
if (.not.self%poisson) then
  call buildW(self,grid,igrid)
  call message('  Done buildW',trace)
endif

! Vd to Ep transfer operator J
call buildJ(self,grid,igrid)
call message('  Done buildJ',trace)

! Sd to Sp transfer operator H
call buildH(self,grid,igrid)
call message('  Done buildH',trace)

end subroutine buildmat

! --------------------------------------------------------------------------------------------------

subroutine buildvp(self,grid,igrid)

! Find the internal degrees of freedom that define the
! compound P0 elements on primal cells (space Vp).
!
! The compound element is constant over each cell with value
! equal to the inverse of the cell area. The micro elements
! are themselves constant over each supermesh cell with value
! equal to the inverse of the supermesh cell area. The coefficient
! of the micro element is thus the supermesh cell area divided by
! the compound cell area.
!
! There are 2*neoff micro elements in one compound element,
! and they are assumed to be ordered anticlockwise, starting with
! those incident on the first edge of the compound element.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nf, if0, ie1, iv1, ix, ixm
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), l1sq, l2sq, l3sq, a

! Number of faces on this grid
nf = grid%nface(igrid)

! Number of coefficients to define element - depends on
! shape of primal cell.
self%ncvp(1:nf) = 2*grid%neoff(1:nf,igrid)

! Loop over primal cells
do if0 = 1, nf

  ! Cell centre
  long = grid%flong(if0,igrid)
  lat = grid%flat(if0,igrid)
  call ll2xyz(long,lat,x0)

  ! Loop over edges of cell
  do ix = 1, grid%neoff(if0,igrid)
    ixm = ix - 1
    if (ixm < 1) ixm = grid%neoff(if0,igrid)
    ie1 = grid%eoff(if0,ix,igrid)

    ! Edge crossing point
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x1)

    ! First micro element
    iv1 = grid%voff(if0,ixm,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    self%cvp(if0,2*ix - 1) = a/grid%farea(if0,igrid)

    ! Second micro element
    iv1 = grid%voff(if0,ix,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    self%cvp(if0,2*ix) = a/grid%farea(if0,igrid)

  enddo

enddo

end subroutine buildvp

! --------------------------------------------------------------------------------------------------

subroutine buildsp(self,grid,igrid)

! Find the internal degrees of freedom that define the
! compound RT0 elements on primal cells (space Sp).
!
! Each compound element is associated with a unit normal flux
! at some primal edge. The first two degrees of freedom are
! the fluxes at the two micro edges that make up the primal edge.
! The sign convention is the same as the sign convention for
! the primal edge.
! The remaining degrees of freedom are the micro edge fluxes
! internal to the primal cells either side of the primal edge.
! These dofs are assumed to be ordered anticlockwise starting
! with the microedge that touches the middle of the primal edge,
! for the first primal cell incident on the primal edge, then for
! the second primal cell incident on the primal edge. The sign convention
! is that clockwise circulation about a primal cell centre is positive.
!
! The compound element has constant divergence in the primal cell
! each side of the primal edge. This constraint fixes all but two
! degrees of freedom. These last two are fixed by demanding that
! a certain measure of the vorticity in each compound cell should
! vanish (this is the usual mixed finite element vorticity defined
! by integration by parts against a test function on the supermesh).

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: ne, ie0, if1, if2, ixf, ixv, iv1, nef1, ixc, &
           ie1, ixe, ix, ixv1, ixv2, ixe1, offset(2), &
           i1, i2, ixcp
real(kind=kind_real) :: long, lat, x0(3), x1(3), l1, sgnflx, &
                        l1sq, l2sq, l3sq, a, atot, x2(3), x3(3), &
                        sum1, sum2, c00, sg, contrib

! Number of edges on this grid
ne = grid%nedge(igrid)

! Number of coefficients to define element - depends on
! shape of cells either side of edge. Computed below.


! Loop over edges
do ie0 = 1, ne

  ! Find number of coefficients for this edge
  if1 = grid%fnxte(ie0,1,igrid)
  if2 = grid%fnxte(ie0,2,igrid)
  self%ncsp(ie0) = 2*(grid%neoff(if1,igrid) + grid%neoff(if2,igrid) + 1)

  ! Save related offsets
  offset(1) = 2
  offset(2) = 2 + 2*grid%neoff(if1,igrid)

  ! Find coordinates of edge ie0
  long = self%elong(ie0,igrid)
  lat = self%elat(ie0,igrid)
  call ll2xyz(long,lat,x0)

  ! Find the lengths of the two micro edges that make up
  ! edge ie0, and hence determine the first two coefficients
  do ixv = 1, 2
    iv1 = grid%vofe(ie0,ixv,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x1)
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    self%csp(ie0,ixv) = l1/grid%ldist(ie0,igrid)
  enddo

  ! Loop over faces next to edge ie0
  do ixf = 1, 2
    if1 = grid%fnxte(ie0,ixf,igrid)

    ! How many edges does this cell have
    nef1 = grid%neoff(if1,igrid)

    ! Is edge flux in or out of cell if1 ?
    if (ixf == 1) THEN
      sgnflx =  1.0_kind_real ! outward
    else
      sgnflx = -1.0_kind_real ! inward
    endif

    ! Area of cell if1
    atot = grid%farea(if1,igrid)

    ! Find coordinates of centre of face if1
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)

    ! Which edge of cell if1 is ie0 ?
    ixe = 1
    do while (ie0 .NE. grid%eoff(if1,ixe,igrid))
      ixe = ixe + 1
    enddo

    ! Reset pointer to dofs
    ixc = offset(ixf)

    ! Internal fluxes
    ! First guess to get the divergence right
    ! First two internal fluxes
    ixc = ixc + 1
    self%csp(ie0,ixc) = 0.0_kind_real
    ixv = 3 - ixf
    iv1 = grid%vofe(ie0,ixv,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
    ixc = ixc + 1
    self%csp(ie0,ixc) = sgnflx*(self%csp(ie0,ixv) - a/atot)

    ! Loop to find remaining internal fluxes
    ixe1 = ixe
    do ix = 1, nef1 - 1
      ixv1 = ixe1
      ixe1 = ixe1 + 1
      if (ixe1 > nef1) ixe1 = 1
      ixv2 = ixe1
      ie1 = grid%eoff(if1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      iv1 = grid%voff(if1,ixv1,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      self%csp(ie0,ixc) = self%csp(ie0,ixc - 1) - sgnflx*a/atot
      iv1 = grid%voff(if1,ixv2,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      self%csp(ie0,ixc) = self%csp(ie0,ixc - 1) - sgnflx*a/atot
    enddo

    ! Now correct to make vorticity vanish
    ! Accumulate integrals over micro elements
    sum1 = 0.0_kind_real
    sum2 = 0.0_kind_real
    ! First the contributions from the flux across edge ie0
    sg = 1.0_kind_real
    do ixv = 1, 2
      iv1 = grid%vofe(ie0,ixv,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
      sum2 = sum2 + sg*self%csp(ie0,ixv)*(l2sq - l3sq)/(12.0_kind_real*a)
      sg = -sg
    enddo
    ! Now the contributions from all the internal fluxes
    ! Loop over vertices with two micro elements on each vertex
    ixe1 = ixe
    ixv1 = ixe
    ixc = offset(ixf)
    do ix = 1, nef1

      ! Find the vertex
      iv1 = grid%voff(if1,ixv1,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = grid%eoff(if1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      sum1 = sum1 + l1sq/(4.0_kind_real*a)
      contrib = ( self%csp(ie0,ixc )*(3.0_kind_real*l1sq + l3sq - l2sq)   &
                + self%csp(ie0,ixcp)*(3.0_kind_real*l1sq + l2sq - l3sq) ) &
                                               /(24.0_kind_real*a)
      sum2 = sum2 + contrib
      ! Find the following edge
      ixe1 = ixe1 + 1
      if (ixe1 > nef1) ixe1 = 1
      ie1 = grid%eoff(if1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      if (ix == nef1) ixcp = offset(ixf) + 1
      sum1 = sum1 + l1sq/(4.0_kind_real*a)
      contrib = ( self%csp(ie0,ixcp)*(3.0_kind_real*l1sq + l3sq - l2sq)   &
                + self%csp(ie0,ixc )*(3.0_kind_real*l1sq + l2sq - l3sq) ) &
                                               /(24.0_kind_real*a)
      sum2 = sum2 + contrib

      ixv1 = ixv1 + 1
      if (ixv1 > nef1) ixv1 = 1
    enddo

    ! Correction to remove vorticity
    c00 = - sum2/sum1
    i1 = offset(ixf) + 1
    i2 = offset(ixf) + 2*nef1
    self%csp(ie0,i1:i2) = self%csp(ie0,i1:i2) + c00

  enddo

enddo

end subroutine buildsp

! --------------------------------------------------------------------------------------------------

subroutine buildep(self,grid,igrid)

! Find the internal degrees of freedom that define the
! compound P1 elements on primal cells (space Ep).
!
! The compound element is piecewise linear. It takes the value
! 1 at the central primal grid vertex, zero at all other
! primal grid vertices, and varies linearly along primal edges.
! To define the compound element we need to compute
! (i) the values at the supermesh vertices adjacent to the
! vertex with value 1, and
! (ii) the values at the supermesh vertices at the centres of
! of the primal cells incident on the central primal vertex.
! Thus there are 2*neofv degrees of freedom to be determined.
!
! The degrees of freedom (i) are easily obtained by linear
! interpolation. The degrees of freedom (ii) are fixed by
! demanding that k X grad of the compound element has zero
! vorticity within each face, where the vorticity is defined
! integration by parts, as for the Sp elements.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nv, iv0, ix, ixp, ie1, if1, ixv, ixvp, iv1, iv2, &
           nef1, nev0
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), x3(3), l1, &
                        l1sq, l2sq, l3sq, a, sum1, sum2

! Number of vertices on this grid
nv = grid%nvert(igrid)

! Number of coefficients to define element - depends on
! shape of primal cell.
self%ncep(1:nv) = 2*grid%neofv(1:nv,igrid)


! Loop over primal vertices
do iv0 = 1, nv

  ! Locate the vertex
  long = grid%vlong(iv0,igrid)
  lat = grid%vlat(iv0,igrid)
  call ll2xyz(long,lat,x0)

  ! Loop over edges around this vertex and set the
  ! first neofv degres of freedom
  nev0 = grid%neofv(iv0,igrid)
  do ix = 1, nev0

    ! Locate the supermesh vertex on this edge
    ie1 = grid%eofv(iv0,ix,igrid)
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x1)
    ! calculate distance of supermesh vertex from central vertex
    ! and hence the coefficient
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    self%cep(iv0,ix) = 1.0_kind_real - l1/grid%ldist(ie1,igrid)

  enddo

  ! Now loop over the faces around this vertex and build the compound
  ! element in that face
  do ix = 1, nev0

    ! Locate face centre
    if1 = grid%fofv(iv0,ix,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)

    ! Two edges of this face (four supermesh edges/cells)
    ! contribute to an integral
    ! First edge
    ie1 = grid%eofv(iv0,ix,igrid)
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x2)
    ! Find the vertex of this edge that isn't iv0
    iv2 = grid%vofe(ie1,1,igrid) + grid%vofe(ie1,2,igrid) - iv0
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    call ll2xyz(long,lat,x3)
    call triangle(x2,x3,x1,l1sq,l2sq,l3sq,a)
    sum2 = self%cep(iv0,ix)*(l1sq - l2sq + l3sq)/(8.0_kind_real*a)
    call triangle(x0,x2,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + (l1sq - l2sq + l3sq)/(8.0_kind_real*a) &
                + self%cep(iv0,ix)*(-l1sq + l2sq + l3sq)/(8.0_kind_real*a)
    ! Second edge
    ixp = ix + 1
    if (ixp > nev0) ixp = 1
    ie1 = grid%eofv(iv0,ixp,igrid)
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x2)
    ! Find the vertex of this edge that isn't iv0
    iv2 = grid%vofe(ie1,1,igrid) + grid%vofe(ie1,2,igrid) - iv0
    long = grid%vlong(iv2,igrid)
    lat = grid%vlat(iv2,igrid)
    call ll2xyz(long,lat,x3)
    call triangle(x3,x2,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + self%cep(iv0,ixp)*(-l1sq + l2sq + l3sq)/(8.0_kind_real*a)
    call triangle(x2,x0,x1,l1sq,l2sq,l3sq,a)
    sum2 = sum2 + (-l1sq + l2sq + l3sq)/(8.0_kind_real*a) &
                + self%cep(iv0,ixp)*(l1sq - l2sq + l3sq)/(8.0_kind_real*a)

    ! Now construct a second integral summing over all microelements
    sum1 = 0.0_kind_real
    nef1 = grid%neoff(if1,igrid)
    do ixv = 1, nef1

      ! Find the vertex
      iv1 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = grid%eoff(if1,ixv,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + l1sq/(4.0_kind_real*a)
      ! Find the following edge
      ixvp = ixv + 1
      if (ixvp > nef1) ixvp = 1
      ie1 = grid%eoff(if1,ixvp,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + l1sq/(4.0_kind_real*a)

    enddo

    ! Finish calculation of degree of freedom
    self%cep(iv0,ix + nev0) = sum2/sum1

  enddo

enddo

end subroutine buildep

! --------------------------------------------------------------------------------------------------

subroutine buildvd(self,grid,igrid)

! Find the internal degrees of freedom that define the
! compound P0 elements on dual cells (space Vd).
!
! The compound element is constant over each dual cell with value
! equal to the inverse of the dual cell area. The micro elements
! are themselves constant over each supermesh cell with value
! equal to the inverse of the supermesh cell area. The coefficient
! of the micro element is thus the supermesh cell area divided by
! the compound dual cell area.
!
! There are 2*neofv micro elements in one compound element,
! and they are assumed to be ordered anticlockwise, starting with
! those incident on the first edge of the compound element.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nv, iv0, ie1, if1, ix, ixm
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), l1sq, l2sq, l3sq, a

! Number of dual cells (vertices) on this grid
nv = grid%nvert(igrid)

! Number of coefficients to define element - depends on
! shape of dual cell.
self%ncvd(1:nv) = 2*grid%neofv(1:nv,igrid)


! Loop over dual cells
do iv0 = 1, nv

  ! Dual cell centre
  long = grid%vlong(iv0,igrid)
  lat = grid%vlat(iv0,igrid)
  call ll2xyz(long,lat,x0)

  ! Loop over edges of dual cell
  do ix = 1, grid%neofv(iv0,igrid)
    ixm = ix - 1
    if (ixm < 1) ixm = grid%neofv(iv0,igrid)
    ie1 = grid%eofv(iv0,ix,igrid)

    ! Edge crossing point
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x1)

    ! First micro element
    if1 = grid%fofv(iv0,ixm,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    self%cvd(iv0,2*ix - 1) = a/self%varea(iv0,igrid)

    ! Second micro element
    if1 = grid%fofv(iv0,ix,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    self%cvd(iv0,2*ix) = a/self%varea(iv0,igrid)

  enddo

enddo

end subroutine buildvd

! --------------------------------------------------------------------------------------------------

subroutine buildsd(self,grid,igrid)

! Find the internal degrees of freedom that define the
! compound N0 elements on dual cells (space Sd).
!
! Each compound element is associated with a unit tangential circulation
! at some dual edge. The first two degrees of freedom are
! the circulations at the two micro edges that make up the dual edge.
! The sign convention is the same as the sign convention for
! the dual edge.
! The remaining degrees of freedom are the micro edge circulations
! internal to the dual cells either side of the dual edge.
! These dofs are assumed to be ordered anticlockwise starting
! with the microedge that touches the middle of the dual edge,
! for the first dual cell incident on the dual edge, then for
! the second dual cell incident on the dual edge. The sign convention
! is that outwards from the dual cell centre is positive.
!
! The compound element has constant vorticity in the dual cell
! each side of the dual edge. This constraint fixes all but two
! degrees of freedom. These last two are fixed by demanding that
! a certain measure of the divergence in each compound cell should
! vanish (this is the usual mixed finite element divergence defined
! by integration by parts against a test function on the supermesh).

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: ne, ie0, iv1, iv2, offset(2), ixf, if1, ixv, nev1, &
           ixe, ixc, ixe1, ixf1, ixf2, ixcp, ix, i1, i2, ie1
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), x3(3), l1, sgncrc, atot, &
                        l1sq, l2sq, l3sq, a, sum1, sum2, contrib, c00, sg

! Number of edges on this grid
ne = grid%nedge(igrid)

! Number of coefficients to define element - depends on
! shape of cells either side of edge. Computed below.


! Loop over edges
do ie0 = 1, ne

  ! Find number of coefficients for this edge
  iv1 = grid%vofe(ie0,1,igrid)
  iv2 = grid%vofe(ie0,2,igrid)
  self%ncsd(ie0) = 2*(grid%neofv(iv1,igrid) + grid%neofv(iv2,igrid) + 1)

  ! Save related offsets
  offset(1) = 2
  offset(2) = 2 + 2*grid%neofv(iv1,igrid)

  ! Find coordinates of edge ie0
  long = self%elong(ie0,igrid)
  lat = self%elat(ie0,igrid)
  call ll2xyz(long,lat,x0)

  ! Find the lengths of the two micro edges that make up
  ! edge ie0, and hence determine the first two coefficients
  do ixf = 1, 2
    if1 = grid%fnxte(ie0,ixf,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)
    l1 = SQRT((x1(1) - x0(1))**2 + (x1(2) - x0(2))**2 + (x1(3) - x0(3))**2)
    self%csd(ie0,ixf) = l1/grid%ddist(ie0,igrid)
  enddo

  ! Loop over dual cells next to edge ie0
  do ixv = 1, 2
    iv1 = grid%vofe(ie0,ixv,igrid)

    ! How many edges does this dual cell have
    nev1 = grid%neofv(iv1,igrid)

    ! Is circulation positive or negative for dual cell iv1 ?
    if (ixv == 1) THEN
      sgncrc = -1.0_kind_real
    else
      sgncrc = 1.0_kind_real
    endif

    ! Area of dual cell iv1
    atot = self%varea(iv1,igrid)

    ! Find coordinates of centre of dual cell iv1
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x1)

    ! Which edge of dual cell iv1 is ie0 ?
    ixe = 1
    do while (ie0 .NE. grid%eofv(iv1,ixe,igrid))
      ixe = ixe + 1
    enddo

    ! Reset pointer to dofs
    ixc = offset(ixv)

    ! Internal circulations
    ! First guess to get the vorticity right
    ! First two internal circulations
    ixc = ixc + 1
    self%csd(ie0,ixc) = 0.0_kind_real
    ixf = ixv
    if1 = grid%fnxte(ie0,ixf,igrid)
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
    ixc = ixc + 1
    self%csd(ie0,ixc) = sgncrc*(self%csd(ie0,ixf) - a/atot)

    ! Loop to find remaining internal circulations
    ixe1 = ixe
    do ix = 1, nev1 - 1
      ixf1 = ixe1
      ixe1 = ixe1 + 1
      if (ixe1 > nev1) ixe1 = 1
      ixf2 = ixe1
      ie1 = grid%eofv(iv1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      if1 = grid%fofv(iv1,ixf1,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      self%csd(ie0,ixc) = self%csd(ie0,ixc - 1) - sgncrc*a/atot
      if1 = grid%fofv(iv1,ixf2,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      self%csd(ie0,ixc) = self%csd(ie0,ixc - 1) - sgncrc*a/atot
    enddo

    ! Now correct to make vorticity vanish
    ! Accumulate integrals over micro elements
    sum1 = 0.0_kind_real
    sum2 = 0.0_kind_real
    ! First the contributions from the circulations at edge ie0
    sg = -1.0_kind_real
    do ixf = 1, 2
      if1 = grid%fnxte(ie0,ixf,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
!      Restore this erroneous line for backward compatibility
!      sum2 = sum2 + sg*self%csd(ie0,ixf)*(l2sq - l3sq)/(12.0_kind_real*a)
!     This is the correct calculation
      sum2 = sum2 + sg*self%csd(ie0,ixf)*(l3sq - l2sq)/(12.0_kind_real*a)
      sg = -sg
    enddo

    ! Now the contributions from all the internal circulations
    ! Loop over dual vertices with two micro elements on each vertex
    ixe1 = ixe
    ixf1 = ixe
    ixc = offset(ixv)
    do ix = 1, nev1

      ! Find the dual vertex
      if1 = grid%fofv(iv1,ixf1,igrid)
      long = grid%flong(if1,igrid)
      lat = grid%flat(if1,igrid)
      call ll2xyz(long,lat,x2)
      ! Find the preceding edge
      ie1 = grid%eofv(iv1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      sum1 = sum1 + l1sq/(4.0_kind_real*a)
      contrib = ( self%csd(ie0,ixc )*(3.0_kind_real*l1sq + l3sq - l2sq)   &
                + self%csd(ie0,ixcp)*(3.0_kind_real*l1sq + l2sq - l3sq) ) &
                                               /(24.0_kind_real*a)
      sum2 = sum2 + contrib
      ! Find the following edge
      ixe1 = ixe1 + 1
      if (ixe1 > nev1) ixe1 = 1
      ie1 = grid%eofv(iv1,ixe1,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      ixc = ixc + 1
      ixcp = ixc + 1
      if (ix == nev1) ixcp = offset(ixv) + 1
      sum1 = sum1 + l1sq/(4.0_kind_real*a)
      contrib = ( self%csd(ie0,ixcp)*(3.0_kind_real*l1sq + l3sq - l2sq)   &
                + self%csd(ie0,ixc )*(3.0_kind_real*l1sq + l2sq - l3sq) ) &
                                               /(24.0_kind_real*a)
      sum2 = sum2 + contrib

      ixf1 = ixf1 + 1
      if (ixf1 > nev1) ixf1 = 1
    enddo

    ! Correction to remove vorticity
    c00 = - sum2/sum1
    i1 = offset(ixv) + 1
    i2 = offset(ixv) + 2*nev1
    self%csd(ie0,i1:i2) = self%csd(ie0,i1:i2) + c00

  enddo

enddo

end subroutine buildsd

! --------------------------------------------------------------------------------------------------

subroutine buildL(self,grid,igrid)

! To build the stencil and coefficients for the Vp mass matrix L
!
! The stencil is a single cell. The coefficient is just the reciprocal
! of the cell area, but let's go through the motions of summing the
! integrals over microelements to check.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nf, if0, ix, ixm, ie1, iv1
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), c, &
                        l1sq, l2sq, l3sq, a, sum1

nf = grid%nface(igrid)

! Loop over primal cells
do if0 = 1, nf

  ! Stencil is just the cell itself
  self%nlsten(if0,igrid) = 1
  self%lsten(if0,1,igrid) = if0

  sum1 = 0.0_kind_real

  ! Cell centre
  long = grid%flong(if0,igrid)
  lat = grid%flat(if0,igrid)
  call ll2xyz(long,lat,x0)

  ! Loop over edges of cell
  do ix = 1, grid%neoff(if0,igrid)
    ixm = ix - 1
    if (ixm < 1) ixm = grid%neoff(if0,igrid)
    ie1 = grid%eoff(if0,ix,igrid)

    ! Edge crossing point
    long = self%elong(ie1,igrid)
    lat = self%elat(ie1,igrid)
    call ll2xyz(long,lat,x1)

    ! First micro element
    iv1 = grid%voff(if0,ixm,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    c = self%cvp(if0,2*ix - 1)
    sum1 = sum1 + c*c/a

    ! Second micro element
    iv1 = grid%voff(if0,ix,igrid)
    long = grid%vlong(iv1,igrid)
    lat = grid%vlat(iv1,igrid)
    call ll2xyz(long,lat,x2)
    call triangle(x0,x1,x2,l1sq,l2sq,l3sq,a)
    c = self%cvp(if0,2*ix)
    sum1 = sum1 + c*c/a
  enddo
  self%lmass(if0,1,igrid) = sum1

enddo

end subroutine buildL

! --------------------------------------------------------------------------------------------------

subroutine buildM(self,grid,igrid)

! To build the stencil and coefficients for the Sp mass matrix M
!
! The stencil for edge ie0 comprises all the edges of the two
! primal cells either side of edge ie0
!
! Also build the stencil and coefficients for the T operator
!
! The stencil for cell if0 is all the edges of cell if0

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: ne, ie0, ixsten, ixf, if1, ix0, ixe, nef1, off1, &
           ixf2, ie2, ix3, ie3, ixm, ixc, ixp, ixm2, ixc2, ixp2, &
           disp, iv2, ixv
real(kind=kind_real) :: long, lat, x0(3), x1(3), x2(3), x3(3), &
                        l1sq, l2sq, l3sq, a, sum1, &
                        a1(2), a2(2), b1(2*grid%nefmx), b2(2*grid%nefmx)

if (.not.self%poisson) then
  self%ntsten(:,igrid) = grid%neoff(:,igrid)
  self%tsten(:,:,igrid) = grid%eoff(:,:,igrid)
endif

ne = grid%nedge(igrid)

! Loop over edges
do ie0 = 1, ne

  ! Locate the edge crossing point
  long = self%elong(ie0,igrid)
  lat = self%elat(ie0,igrid)
  call ll2xyz(long,lat,x0)

  ixsten = 1

  ! Make sure first stencil edge is ie0 itself
  self%msten(ie0,1,igrid) = ie0
  self%mmass(ie0,1,igrid) = 0.0_kind_real

  ! Loop over the two primal cells either side
  do ixf = 1, 2
    if1 = grid%fnxte(ie0,ixf,igrid)
    nef1 = grid%neoff(if1,igrid)

    ! Extract coefficients of ie0 basis function in this cell
    if (ixf == 1) THEN
      a1(1) = self%csp(ie0,1)
      a1(2) = self%csp(ie0,2)
      b1(1:2*nef1) = self%csp(ie0,3:2*(nef1 + 1))
    else
      off1 = 2*grid%neoff(grid%fnxte(ie0,1,igrid),igrid)
      a1(1) = -self%csp(ie0,2)
      a1(2) = -self%csp(ie0,1)
      b1(1:2*nef1) = self%csp(ie0,3 + off1:2*(nef1 + 1) + off1)
    endif

    ! Locate the primal cell centre
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)

    ! Which edge of if1 is ie0 ?
    ix0 = 1
    do while (grid%eoff(if1,ix0,igrid) .NE. ie0)
      ix0 = ix0 + 1
    enddo

    ! Loop over edges of if1
    do ixe = 1, nef1
      ie2 = grid%eoff(if1,ixe,igrid)

      ! Which face of ie2 is if1 ?
      ixf2 = 1
      if (grid%fnxte(ie2,2,igrid) == if1) ixf2 = 2

      ! Locate the edge crossing point
      long = self%elong(ie2,igrid)
      lat = self%elat(ie2,igrid)
      call ll2xyz(long,lat,x3)

      ! Extract coefficients of ie2 basis function in this cell
      if (ixf2 == 1) THEN
        a2(1) = self%csp(ie2,1)
        a2(2) = self%csp(ie2,2)
        b2(1:2*nef1) = self%csp(ie2,3:2*(nef1 + 1))
      else
        off1 = 2*grid%neoff(grid%fnxte(ie2,1,igrid),igrid)
        a2(1) = -self%csp(ie2,2)
        a2(2) = -self%csp(ie2,1)
        b2(1:2*nef1) = self%csp(ie2,3 + off1:2*(nef1 + 1) + off1)
      endif

      ! Displacement of edge ie2 from edge ie0 around cell if1
      disp = ixe - ix0
      if (disp < 0) disp = disp + nef1

      ! Integrate the product of the two basis functions

      if (disp > 0) THEN
        ! Initialize sum to zero
        sum1 = 0.0_kind_real
      else
        ! We need the product of the exterior flux contributions
        ixv = ix0 - 1
        if (ixv < 1) ixv = nef1
        iv2 = grid%voff(if1,ixv,igrid)
        long = grid%vlong(iv2,igrid)
        lat = grid%vlat(iv2,igrid)
        call ll2xyz(long,lat,x2)
        call triangle(x1,x2,x0,l1sq,l2sq,l3sq,a)
        sum1 = a1(1)*a2(1)*(3.0_kind_real*(l2sq + l3sq) - l1sq)/(48.0_kind_real*a)
        ixv = ix0
        iv2 = grid%voff(if1,ixv,igrid)
        long = grid%vlong(iv2,igrid)
        lat = grid%vlat(iv2,igrid)
        call ll2xyz(long,lat,x2)
        call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + a1(2)*a2(2)*(3.0_kind_real*(l2sq + l3sq) - l1sq)/(48.0_kind_real*a)
      endif

      ! Now cross products of exterior and interior fluxes
      ! Indices of edges of micro elements on ie0 relative to ie2
      ixm = 2*(nef1 - disp)
      ixc = ixm + 1
      if (ixc > 2*nef1) ixc = 1
      ixp = ixc + 1
      ixv = ix0 - 1
      if (ixv < 1) ixv = nef1
      iv2 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv2,igrid)
      lat = grid%vlat(iv2,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x0,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a1(1)*( - b2(ixc)*(l1sq + l2sq - 3.0_kind_real*l3sq)    &
                            + b2(ixm)*(l1sq + l3sq - 3.0_kind_real*l2sq) )  &
                                                     /(48.0_kind_real*a)
      ixv = ix0
      iv2 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv2,igrid)
      lat = grid%vlat(iv2,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x0,x2,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a1(2)*( - b2(ixp)*(l1sq + l2sq - 3.0_kind_real*l3sq)    &
                            + b2(ixc)*(l1sq + l3sq - 3.0_kind_real*l2sq) )  &
                                                     /(48.0_kind_real*a)
      ! Indices of edges of micro elements on ie2 relative to ie0
      ixm = 2*disp
      ixc = ixm + 1
      ixp = ixc + 1
      if (ixm == 0) ixm = 2*nef1
      ixv = ixe - 1
      if (ixv < 1) ixv = nef1
      iv2 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv2,igrid)
      lat = grid%vlat(iv2,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a2(1)*( - b1(ixc)*(l1sq + l2sq - 3.0_kind_real*l3sq)    &
                            + b1(ixm)*(l1sq + l3sq - 3.0_kind_real*l2sq) )  &
                                                     /(48.0_kind_real*a)
      ixv = ixe
      iv2 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv2,igrid)
      lat = grid%vlat(iv2,igrid)
      call ll2xyz(long,lat,x2)
      call triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      sum1 = sum1 + a2(2)*( - b1(ixp)*(l1sq + l2sq - 3.0_kind_real*l3sq)    &
                            + b1(ixc)*(l1sq + l3sq - 3.0_kind_real*l2sq) )  &
                                                     /(48.0_kind_real*a)

      ! And finally products of interior flux contributions
      ! Loop over edges of if1
      do ix3 = 1, nef1

        ! Indices of elements on edge ie3 relative to ie0
        disp = ix3 - ix0
        if (disp < 0) disp = disp + nef1
        ixm = 2*disp
        ixc = ixm + 1
        ixp = ixc + 1
        if (ixm == 0) ixm = 2*nef1
        ! Indices of elements on edge ie3 relative to ie2
        disp = ix3 - ixe
        if (disp < 0) disp = disp + nef1
        ixm2 = 2*disp
        ixc2 = ixm2 + 1
        ixp2 = ixc2 + 1
        if (ixm2 == 0) ixm2 = 2*nef1

        ! Locate edge ie3
        ie3 = grid%eoff(if1,ix3,igrid)
        long = self%elong(ie3,igrid)
        lat = self%elat(ie3,igrid)
        call ll2xyz(long,lat,x3)

        ! Two micro elements incident on edge ie3
        ixv = ix3 - 1
        if (ixv < 1) ixv = nef1
        iv2 = grid%voff(if1,ixv,igrid)
        long = grid%vlong(iv2,igrid)
        lat = grid%vlat(iv2,igrid)
        call ll2xyz(long,lat,x2)
        call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + ( b1(ixm)*b2(ixm2)*(3.0_kind_real*(l1sq + l2sq) - l3sq)        &
                      + b1(ixc)*b2(ixc2)*(3.0_kind_real*(l1sq + l3sq) - l2sq)        &
                     - (b1(ixm)*b2(ixc2) + b1(ixc)*b2(ixm2))                 &
                                        *(l2sq + l3sq - 3.0_kind_real*l1sq)    )     &
                                                       /(48.0_kind_real*a)
        ixv = ix3
        iv2 = grid%voff(if1,ixv,igrid)
        long = grid%vlong(iv2,igrid)
        lat = grid%vlat(iv2,igrid)
        call ll2xyz(long,lat,x2)
        call triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
        sum1 = sum1 + ( b1(ixc)*b2(ixc2)*(3.0_kind_real*(l1sq + l2sq) - l3sq)        &
                      + b1(ixp)*b2(ixp2)*(3.0_kind_real*(l1sq + l3sq) - l2sq)        &
                     - (b1(ixc)*b2(ixp2) + b1(ixp)*b2(ixc2))                 &
                                        *(l2sq + l3sq - 3.0_kind_real*l1sq)    )     &
                                                       /(48.0_kind_real*a)

      enddo

      ! Finally save the result
      if (ie2 == ie0) THEN
        self%mmass(ie0,1,igrid) = self%mmass(ie0,1,igrid) + sum1
      else
        ixsten = ixsten + 1
        self%msten(ie0,ixsten,igrid) = ie2
        self%mmass(ie0,ixsten,igrid) = sum1
      endif
      if (.not.self%poisson) then
        self%tcoeff(if1,ix0,ixe,igrid) = sum1
      endif

    enddo

  enddo
  self%nmsten(ie0,igrid) = ixsten

enddo

end subroutine buildM

! --------------------------------------------------------------------------------------------------

subroutine buildR(self,grid,igrid)

! To build the R operator to transfer from Vp to Ep
!
! The operator is built in two forms, one suitable for
! looping over faces (which can be used for the TRiSK
! implementation of W), one suitable for looping over
! vertices (avoids reduce_all).

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nv, nf, iv0, nev1, if1, ixv0, ixf, ixfp, nef1, &
           ixv, ixvm, iv1, ixe, ie1, ne1
real(kind=kind_real) :: v1, v2, vc, long, lat, x1(3), x2(3), x3(3), &
                        sum1, l1sq, l2sq, l3sq, a, aby3

nv = grid%nvert(igrid)

self%nrsten(:,igrid) = grid%neoff(:,igrid)

! Loop over vertices
do iv0 = 1, nv

  nev1 = grid%neofv(iv0,igrid)

  ! Loop over the faces around iv0
  do ixf = 1, nev1

    if1 = grid%fofv(iv0,ixf,igrid)

    ! Locate the face centre
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)

    ! Pick out degrees of freedom defining Ep basis element
    ! in this cell
    ixfp = ixf + 1
    if (ixfp > nev1) ixfp = 1
    v1 = self%cep(iv0,ixfp)
    v2 = self%cep(iv0,ixf)
    vc = self%cep(iv0,ixf + nev1)

    ! Which vertex of if1 is iv0 ?
    ixv0 = 1
    do while (grid%voff(if1,ixv0,igrid) .NE. iv0)
      ixv0 = ixv0 + 1
    enddo

    ! Now loop over the micro elements in if1, integrating
    ! the product of the basis functions
    nef1 = grid%neoff(if1,igrid)
    sum1 = 0.0_kind_real
    do ixe = 1, nef1
      ie1 = grid%eoff(if1,ixe,igrid)
      ! Locate edge crossing point
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x2)
      ! Two micro elements on this edge
      ixv = ixe
      ixvm = ixv - 1
      if (ixvm < 1) ixvm = nef1
      ! First one
      iv1 = grid%voff(if1,ixvm,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0_kind_real
      sum1 = sum1 + vc*aby3
      if (ixv  == ixv0) sum1 = sum1 + v1*aby3
      if (ixvm == ixv0) sum1 = sum1 + (1.0_kind_real + v2)*aby3
      ! Second one
      iv1 = grid%voff(if1,ixv,igrid)
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0_kind_real
      sum1 = sum1 + vc*aby3
      if (ixv  == ixv0) sum1 = sum1 + (1.0_kind_real + v1)*aby3
      if (ixvm == ixv0) sum1 = sum1 + v2*aby3
    enddo
    self%rsten(if1,ixv0,igrid) = iv0
    self%rcoeff(if1,ixv0,igrid) = sum1/grid%farea(if1,igrid)

  enddo

enddo


! Re-tabulate the coefficients of the R operator in a form
! suitable for looping over vertices rather than faces

nf = grid%nface(igrid)

self%nrxsten(:,igrid) = grid%neofv(:,igrid)

do if1 = 1, nf

  ne1 = grid%neoff(if1,igrid)
  do ixv = 1, ne1

    iv1 = grid%voff(if1,ixv,igrid)

    ! What is the index of face if1 relative to vertex iv1 ?
    ixf = 1
    do while (grid%fofv(iv1,ixf,igrid) .NE. if1)
      ixf = ixf + 1
    enddo

    ! Save the stencil face and coefficient
    self%rxsten(iv1,ixf,igrid) = if1
    self%rxcoeff(iv1,ixf,igrid) = self%rcoeff(if1,ixv,igrid)

  enddo

enddo

end subroutine buildR

! --------------------------------------------------------------------------------------------------

subroutine buildW(self,grid,igrid)

! To construct the stencil and coefficients for the W operator,
! which maps from Ep to Ep.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: ie0, ixf, if1, ne1, ix1, ix, ixv, ix2, ie2, &
           isten, ne
real(kind=kind_real) :: ss, w

! The W operator is built from the R operator a la TRiSK
! (avoids having to do more integrals, which would be another
! way of doing it)

ne = grid%nedge(igrid)

! Loop over edges
do ie0 = 1, ne

  isten = 0

  ! Loop over the faces either side of this edge
  do ixf = 1, 2
    if1 = grid%fnxte(ie0,ixf,igrid)
    ne1 = grid%neoff(if1,igrid)
    ss = -0.5_kind_real

    ! Which edge of face if1 is edge ie0?
    ix1 = 1
    do while (grid%eoff(if1,ix1,igrid) .NE. ie0)
      ix1 = ix1 + 1
    enddo

    ! Find the contribution from every other
    ! edge of this face
    do ix = 0, ne1 - 2
      ixv = MOD(ix1 + ix - 1,ne1) + 1
      ix2 = MOD(ix1 + ix,ne1) + 1
      ie2 = grid%eoff(if1,ix2,igrid)
      ss = ss + self%rcoeff(if1,ixv,igrid)
      w = -ss*self%eoffin(if1,ix1,igrid)*self%eoffin(if1,ix2,igrid)
      isten = isten + 1
      self%wsten(ie0,isten,igrid) = ie2
      self%wcoeff(ie0,isten,igrid) = w
    enddo
  enddo
  self%nwsten(ie0,igrid) = isten
enddo

end subroutine buildW

! --------------------------------------------------------------------------------------------------

subroutine buildJ(self,grid,igrid)

! To construct the stencil and coefficients for the J operator,
! which maps from Vd to Ep.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: nv, iv0, nev1, ixf, nef1, ns, ixfp, ix, iv1, &
           ixv0, ixvm, ixvp, ixe, ixep, ie1, ixv, if1, temp(self%njsmx)
real(kind=kind_real) :: v1, v2, vc, long, lat, x1(3), x2(3), x3(3), &
                        l1sq, l2sq, l3sq, a, aby3, sum1

nv = grid%nvert(igrid)

! Loop over vertices
do iv0 = 1, nv

  nev1 = grid%neofv(iv0,igrid)

  ! Initialize to make sure iv0 itself is first in the stencil
  self%jsten(iv0,:,igrid) = 0
  self%jsten(iv0,1,igrid) = iv0
  ns = 1
  self%jstar(iv0,:,igrid) = 0.0_kind_real

  ! Loop over the faces comprising the element in Ep
  do ixf = 1, nev1

    if1 = grid%fofv(iv0,ixf,igrid)
    nef1 = grid%neoff(if1,igrid)

    ! Locate the face centre
    long = grid%flong(if1,igrid)
    lat = grid%flat(if1,igrid)
    call ll2xyz(long,lat,x1)

    ! Which vertex of if1 is iv0 ?
    ixv0 = 1
    do while (grid%voff(if1,ixv0,igrid) .NE. iv0)
      ixv0 = ixv0 + 1
    enddo
    ixvm = ixv0 - 1
    if (ixvm < 1) ixvm = nef1
    ixvp = ixv0 + 1
    if (ixvp > nef1) ixvp = 1

    ! Pick out degrees of freedom defining Ep basis element
    ! in this cell
    ixfp = ixf + 1
    if (ixfp > nev1) ixfp = 1
    v1 = self%cep(iv0,ixfp)
    v2 = self%cep(iv0,ixf)
    vc = self%cep(iv0,ixf + nev1)

    ! Loop over vertices of face if1
    do ixv = 1, nef1
      ixe = ixv
      ixep = ixe + 1
      if (ixep > nef1) ixep = 1
      iv1 = grid%voff(if1,ixv,igrid)
      ! Locate the vertex iv1
      long = grid%vlong(iv1,igrid)
      lat = grid%vlat(iv1,igrid)
      call ll2xyz(long,lat,x2)
      ! Two micro elements are incident on vertex iv1 in cell if1
      ! First one
      ie1 = grid%eoff(if1,ixe,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x3,x2,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0_kind_real
      sum1 = vc*aby3
      if (ixv == ixv0) sum1 = sum1 + (1.0_kind_real + v1)*aby3
      if (ixv == ixvp) sum1 = sum1 + v2*aby3
      ! Second one
      ie1 = grid%eoff(if1,ixep,igrid)
      long = self%elong(ie1,igrid)
      lat = self%elat(ie1,igrid)
      call ll2xyz(long,lat,x3)
      call triangle(x1,x2,x3,l1sq,l2sq,l3sq,a)
      aby3 = a/3.0_kind_real
      sum1 = sum1 + vc*aby3
      if (ixv == ixvm) sum1 = sum1 + v1*aby3
      if (ixv == ixv0) sum1 = sum1 + (1.0_kind_real + v2)*aby3
      ! Now include this contribution in the coefficient
      temp = self%jsten(iv0,:,igrid)
      call findinlist(iv1,temp,self%njsmx,ix)
      self%jsten(iv0,:,igrid) = temp
      ns = MAX(ns,ix)
      self%jstar(iv0,ix,igrid) = self%jstar(iv0,ix,igrid) + sum1/self%varea(iv1,igrid)
    enddo

  enddo
  self%njsten(iv0,igrid) = ns

enddo

end subroutine buildJ

! --------------------------------------------------------------------------------------------------

subroutine buildH(self,grid,igrid)

! To construct the stencil and coefficients for the H operator,
! which maps from Sd to Sp.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid
integer,          intent(in)    :: igrid

integer :: ne, ie0, ns, ixf, if1, nef1, iv1, nev1, &
           ie2, ixf2, ixe2, offset, offset2, ivfrst, d, &
           ip3, ip4, ip5, iq3, iq4, iq5, ixe0, ixv, ix, temp(self%nhsmx)
real(kind=kind_real) :: flx1, flx2, flx3, flx4, flx5, crc1, crc2, crc3, crc4, crc5, &
                        extf1, extf2, extc1, extc2, sum1

ne = grid%nedge(igrid)

! Loop over all edges
do ie0 = 1, ne

  ! Initialize to make sure ie0 itself is first in the stencil
  self%hsten(ie0,:,igrid) = 0
  self%hsten(ie0,1,igrid) = ie0
  ns = 1
  self%hstar(ie0,:,igrid) = 0.0_kind_real

  ! Two faces either side of edge ie0
  do ixf = 1, 2

    if1 = grid%fnxte(ie0,ixf,igrid)
    nef1 = grid%neoff(if1,igrid)

    ! Which edge of if1 is ie0 ?
    ixe0 = 1
    do while (grid%eoff(if1,ixe0,igrid) .NE. ie0)
      ixe0 = ixe0 + 1
    enddo

    ! Offset to coefficients for internal degrees of freedom,
    ! and orient external fluxes
    if (ixf == 1) THEN
      offset = 2
      extf1 = self%csp(ie0,1)
      extf2 = self%csp(ie0,2)
    else
      offset = 2 + 2*grid%neoff(grid%fnxte(ie0,1,igrid),igrid)
      extf1 = -self%csp(ie0,2)
      extf2 = -self%csp(ie0,1)
    endif

    ! Loop over vertices of face if1
    do ixv = 1, nef1

      iv1 = grid%voff(if1,ixv,igrid)

      ! Pick out coefficients defining the part of the Sp
      ! element in this corner
      d = ixv - ixe0
      if (d < 0) d = d + nef1
      if (d == 0) THEN
        flx1 = extf2
      else
        flx1 = 0.0_kind_real
      endif
      if (d == nef1 - 1) THEN
        flx2 = extf1
      else
        flx2 = 0.0_kind_real
      endif
      ip3 = 2*d + 1
      flx3 = self%csp(ie0,offset + ip3)
      ip4 = ip3 + 1
      flx4 = self%csp(ie0,offset + ip4)
      ip5 = ip4 + 1
      if (ip5 > 2*nef1) ip5 = 1
      flx5 = self%csp(ie0,offset + ip5)

      ! Which face of iv1 is if1 ?
      ixf2 = 1
      do while (grid%fofv(iv1,ixf2,igrid) .NE. if1)
        ixf2 = ixf2 + 1
      enddo

      ! Loop over all edges around vertex iv1 - these
      ! are the stencil edges
      nev1 = grid%neofv(iv1,igrid)
      do ixe2 = 1, nev1

        ie2 = grid%eofv(iv1,ixe2,igrid)

        ! Is iv1 the first or second vertex of edge ie2 ?
        ivfrst = grid%vofe(ie2,1,igrid)
        if (iv1 == ivfrst) THEN
          offset2 = 2
          extc1 = self%csd(ie2,1)
          extc2 = self%csd(ie2,2)
        else
          offset2 = 2 + 2*grid%neofv(ivfrst,igrid)
          extc1 = -self%csd(ie2,2)
          extc2 = -self%csd(ie2,1)
        endif

        ! Pick out coefficients defining the part of the Sd
        ! element in this corner
        d = ixf2 - ixe2
        if (d < 0) d = d + nev1
        if (d == 0) THEN
          crc2 = -extc1
        else
          crc2 = 0.0_kind_real
        endif
        if (d == nev1 - 1) THEN
          crc1 = -extc2
        else
          crc1 = 0.0_kind_real
        endif
        iq3 = 2*d + 1
        crc3 = self%csd(ie2,offset2 + iq3)
        iq4 = iq3 + 1
        crc4 = self%csd(ie2,offset2 + iq4)
        iq5 = iq4 + 1
        if (iq5 > 2*nev1) iq5 = 1
        crc5 = self%csd(ie2,offset2 + iq5)

        ! Integral over the two overlapping microelements
        sum1 = (-flx1*crc4 + flx1*crc1 + flx4*crc5 + flx4*crc1   &
               + flx3*crc5 + flx3*crc4 - flx2*crc2 - flx2*crc4   &
               - flx5*crc3 - flx5*crc4 - flx4*crc3 + flx4*crc2)  &
                                             / 6.0_kind_real

        ! Update stencil and coefficients
        temp = self%hsten(ie0,:,igrid)
        call findinlist(ie2,temp,self%nhsmx,ix)
        self%hsten(ie0,:,igrid) = temp

        ns = MAX(ns,ix)
        self%hstar(ie0,ix,igrid) = self%hstar(ie0,ix,igrid) + sum1

      enddo

    enddo

  enddo
  self%nhsten(ie0,igrid) = ns

enddo

end subroutine buildH

! --------------------------------------------------------------------------------------------------

subroutine buildinj(self,grid)

! To build restriction operator for multigrid
! Note this code assumes the same grid cell numbering as
! the gengrid_hex.f and gengrid_cube.f90 grid generation programs.

implicit none
class(fempsoprs), intent(inout) :: self
type(fempsgrid),  intent(in)    :: grid

integer :: if0, if1, igrid, igridp, nf, ix, ixp, iv1, iv2, &
           n2, n, n2p, np, i, j, &
           t1, s1, s2, s3, s4, p
character(len=255) :: errormessage

! if ((grid%nedgex == 3*(grid%nfacex - 2)) .AND. (grid%nefmx == 6) .AND. (grid%nevmx == 3)) THEN
if (grid%gtype == 'ih') THEN

  ! It looks like we're on a hexagonal icoshedral grid, so build the
  ! appropriate multigrid restriction operator

  ! Stencil is 6 or 7 cells on standard buckyball grid
  self%ninjmx = 7

  ! Allocate arrays for operator stencils and coefficients
  allocate(self%injsten(grid%nfacex,self%ninjmx,grid%ngrids-1))
  allocate(self%injwgt(grid%nfacex,self%ninjmx,grid%ngrids-1))

  ! Initialize to zero
  self%injsten = 0
  self%injwgt = 0.0_kind_real

  ! And define stencil and coefficients
  do igrid = 1, grid%ngrids-1
    igridp = igrid + 1
    do if0 = 1, grid%nface(igrid)
      ! Face if0 exists on grid igrid and grid igrid+1 and is
      ! the centre of the stencil
      self%injsten(if0,1,igrid) = if0
      self%injwgt(if0,1,igrid) = 1.0_kind_real
      ! The neighbours of if0 on grid igrid+1 are the other members
      ! of the stencil
      nf = grid%neoff(if0,igridp)
      self%ninj(if0,igrid) = 1 + nf
      do ix = 1, nf
        ixp = ix + 1
        if1 = grid%fnxtf(if0,ix,igridp)
        self%injsten(if0,ixp,igrid) = if1
        self%injwgt(if0,ixp,igrid) = 0.5_kind_real
      enddo
    enddo
  enddo

!elseif ((grid%nedgex == 2*grid%nfacex) .AND. (grid%nefmx == 4) .AND. (grid%nevmx == 4)) THEN
elseif (grid%gtype == 'cs') THEN

  ! It looks like we're on a cubed sphere grid, so build the
  ! appropriate multigrid restriction operator

  ! Stencil is always 4 cells
  self%ninjmx = 4

  ! Allocate arrays for operator stencils and coefficients
  allocate(self%injsten(grid%nfacex,self%ninjmx,grid%ngrids-1))
  allocate(self%injwgt(grid%nfacex,self%ninjmx,grid%ngrids-1))

  ! Initialize to zero
  self%injsten = 0
  self%injwgt = 0.0_kind_real

  ! And define stencil and coefficients
  do igrid = 1, grid%ngrids-1
    igridp = igrid + 1
    ! Number of cells per face and cells per edge on grid igrid
    n2 = grid%nface(igrid)/6
    n = NINT(SQRT(1.0_kind_real*n2))
    ! And on grid igridp
    n2p = grid%nface(igridp)/6
    np = NINT(SQRT(1.0_kind_real*n2p))
    ! Loop over cells on a panel
    do j = 1, n
      do i = 1, n
        t1 = (j-1)*n + i
        s1 = (2*j - 2)*np + (2*i - 1)
        s2 = s1 + 1
        s3 = s2 + np
        s4 = s3 - 1
        ! Loop over panels
        do p = 1, 6
          if0 = t1 + (p - 1)*n2
          self%ninj(if0,igrid) = 4
          self%injsten(if0,1,igrid) = s1 + (p - 1)*n2p
          self%injsten(if0,2,igrid) = s2 + (p - 1)*n2p
          self%injsten(if0,3,igrid) = s3 + (p - 1)*n2p
          self%injsten(if0,4,igrid) = s4 + (p - 1)*n2p
          self%injwgt(if0,1:4,igrid) = 1.0_kind_real
        enddo
      enddo
    enddo
  enddo

!elseif ((grid%nedgex == 2*grid%nfacex) .AND. (grid%nefmx == 4) .AND. (grid%nevmx == 4)) THEN
elseif (grid%gtype == 'di') THEN

  ! It looks like we're on a diamond grid, so build the
  ! appropriate multigrid restriction operator

  ! Stencil is always 4 cells
  self%ninjmx = 4

  ! Allocate arrays for operator stencils and coefficients
  allocate(self%injsten(grid%nfacex,self%ninjmx,grid%ngrids-1))
  allocate(self%injwgt(grid%nfacex,self%ninjmx,grid%ngrids-1))

  ! Initialize to zero
  self%injsten = 0
  self%injwgt = 0.0_kind_real

  ! And define stencil and coefficients
  do igrid = 1, grid%ngrids-1
    igridp = igrid + 1
    ! Length (cells) of cubed sphere panel
    n = NINT(SQRT(grid%nface(igrid)/12.0_kind_real))
    ! Diamond cells on the coarse grid correspond to cubed sphere edges.
    ! Match these to fine cubed sphere vertices at the centres of the edges.
    ! These vertices are also vertices of the diamond grid and their neighboring
    ! faces comprise the restriction stencil.
    if0 = 0
    do j = 1, n
      do i = 1, n

        ! Odd cube edges / diamond faces
        if0 = if0 + 1
        iv1 = (2*j - 1)*2*n + (2*i - 1)
        do  p = 1, 6
          if1 = if0 + (p - 1)*2*n*n
          iv2 = iv1 + (p - 1)*4*n*n
          self%injsten(if1,:,igrid) = grid%fofv(iv2,:,igrid)
          self%injwgt(if1,1:4,igrid) = 1.0_kind_real
        enddo

        ! Even cube edges / diamond faces
        if0 = if0 + 1
        iv1 = 2*(j - 1)*2*n + 2*i
        do  p = 1, 6
          if1 = if0 + (p - 1)*2*n*n
          iv2 = iv1 + (p - 1)*4*n*n
          self%injsten(if1,:,igrid) = grid%fofv(iv2,:,igrid)
          self%injwgt(if1,1:4,igrid) = 1.0_kind_real
        enddo

      enddo
    enddo
  enddo

else

  write(errormessage,*) 'grid%gtype = ', grid%gtype
  call message(errormessage)
  call message('This grid is not recognized. Please programme the multigrid restriction operator for this grid.', fatal)

endif

end subroutine buildinj

! --------------------------------------------------------------------------------------------------

end module femps_operators_mod
