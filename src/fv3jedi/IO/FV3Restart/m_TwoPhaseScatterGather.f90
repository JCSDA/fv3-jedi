module TwoPhaseScatterGather
  ! Inspired by a post by Jonathan Dursi at
  ! https://stackoverflow.com/questions/29325513/scatter-matrix-blocks-of-different-sizes-using-mpi
  use mpi
  use fv3jedi_geom_mod,             only: fv3jedi_geom
  implicit none

  ! GENERIC INTERFACES
  ! ------------------
  interface TwoPhaseScatter_Phase1
    procedure TwoPhaseScatter_Phase1_r4, TwoPhaseScatter_Phase1_r8
  end interface TwoPhaseScatter_Phase1

  interface TwoPhaseScatter_Phase2
    procedure TwoPhaseScatter_Phase2_r4, TwoPhaseScatter_Phase2_r8
  end interface TwoPhaseScatter_Phase2

  interface TwoPhaseGather_Phase1
    procedure TwoPhaseGather_Phase1_r4, TwoPhaseGather_Phase1_r8
  end interface TwoPhaseGather_Phase1

  interface TwoPhaseGather_Phase2
    procedure TwoPhaseGather_Phase2_r4, TwoPhaseGather_Phase2_r8
  end interface TwoPhaseGather_Phase2

  ! CUSTOM TYPES
  ! ------------
  type :: scatter_t
    logical :: lalloc = .false.
    integer, allocatable :: sendcounts_phase1(:), senddispls_phase1(:)
    integer, allocatable :: sendcounts_phase2(:), senddispls_phase2(:)
    integer :: vec_r4, localvec_r4
    integer :: vec_r8, localvec_r8
  endtype scatter_t

  type :: gather_t
    logical :: lalloc = .false.
    integer, allocatable :: recvcounts_phase1(:), recvdispls_phase1(:)
    integer, allocatable :: recvcounts_phase2(:), recvdispls_phase2(:)
    integer :: vec_r4, localvec_r4
    integer :: vec_r8, localvec_r8
  endtype gather_t

  ! COMMUNICATION STATE OF ONE READ OR ONE WRITE
  ! --------------------------------------------
  ! Everything the two-phase scatter/gather needs for ONE geometry.  Each read or write call
  ! owns one instance: it initializes it (TwoPhaseScatter_init or TwoPhaseGather_init), passes
  ! it to the Phase1/Phase2 routines, and frees it with the matching delete routine before
  ! returning.  Nothing is shared between calls, so a read and a write (or two grids) can no
  ! longer see each other's sizes, communicators or workspaces.
  !
  ! Declare instances TARGET, ASYNCHRONOUS.  The workspaces and the per-owner count and
  ! displacement arrays are handed to nonblocking collectives that are completed later by an
  ! MPI_Waitall in the caller.  ASYNCHRONOUS cannot be given to components, so it goes on the
  ! variable (and on every dummy argument that receives it); its subobjects inherit it.  Do not
  ! copy or modify an instance while requests are outstanding.
  type :: twophase_decomp
    ! Per-owner counts, displacements and MPI datatypes (created lazily in Phase1)
    type(scatter_t), allocatable :: TwoPhaseScatter(:)
    type(gather_t), allocatable  :: TwoPhaseGather(:)

    ! Dedicated Scatter Workspaces
    real(kind=4), allocatable :: scatter_workspace_r4(:,:,:)
    real(kind=8), allocatable :: scatter_workspace_r8(:,:,:)

    ! Dedicated Gather Workspaces
    real(kind=4), allocatable :: gather_workspace_r4(:,:,:)
    real(kind=8), allocatable :: gather_workspace_r8(:,:,:)

    ! Workspaces for local type conversion
    real(kind=4), allocatable :: scatter_recv_cast_workspace_r4(:,:,:)
    real(kind=4), allocatable :: gather_send_cast_workspace_r4(:,:,:)

    ! Geometrical arrays
    integer :: EWindex = -1, NSindex = -1, rowrank = -1, colrank = -1
    integer :: colComm = MPI_COMM_NULL, rowComm = MPI_COMM_NULL
    integer(kind=4), allocatable :: NumColsPerRank(:), NumRowsPerRank(:)
    integer(kind=4), allocatable :: MyRowGlobal(:), MyColGlobal(:)
    integer(kind=4), allocatable :: MyRankInRowComm(:), MyRankInColComm(:)
    integer, allocatable :: ibegin(:), iend(:), jbegin(:), jend(:)
    integer :: layout(2) = 0, globalsizes(2) = 0, localsizes(2) = 0
  end type twophase_decomp

contains

  ! SUBROUTINES
  ! -----------
  subroutine TwoPhaseScatter_init(decomp, geom, npes, batch_size)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    type(fv3jedi_geom), intent(in):: geom
    integer, intent(in) :: npes, batch_size
    integer :: r, ierr
    integer :: rowSize, colSize
    integer :: myRowRank, myColRank
    integer :: mpicomm, wrank, wsize

    mpicomm = geom%f_comm%communicator()
    call MPI_Comm_rank(mpicomm, wrank, ierr)
    call MPI_Comm_size(mpicomm, wsize, ierr)

    ! A decomp is initialized once and freed by the matching delete routine
    if (decomp%rowComm /= MPI_COMM_NULL .or. decomp%colComm /= MPI_COMM_NULL) then
      write(6,'("TwoPhaseScatter_init: decomp is already initialized; call TwoPhaseScatter_delete first")')
      call flush(6)
      call MPI_Abort(mpicomm, 50, ierr)
    endif

    decomp%layout = geom%layout

    decomp%EWindex = modulo(wrank,decomp%layout(1))
    decomp%NSindex = (wrank/decomp%layout(1))

    call MPI_Comm_split(mpicomm, decomp%NSindex, wrank, decomp%rowComm, ierr)
    call MPI_Comm_split(mpicomm, decomp%EWindex, wrank, decomp%colComm, ierr)

    call MPI_Comm_size(decomp%rowComm, rowSize, ierr)
    call MPI_Comm_size(decomp%colComm, colSize, ierr)

    ! Allocate based on actual comm sizes
    allocate(decomp%ibegin(0:rowSize-1), decomp%iend(0:rowSize-1))
    allocate(decomp%jbegin(0:colSize-1), decomp%jend(0:colSize-1))

    call MPI_AllGather(geom%isc,1,MPI_Integer,decomp%ibegin(0:),1,MPI_Integer, decomp%rowComm, ierr)
    call MPI_AllGather(geom%iec,1,MPI_Integer,decomp%iend(0:)  ,1,MPI_Integer, decomp%rowComm, ierr)
    call MPI_AllGather(geom%jsc,1,MPI_Integer,decomp%jbegin(0:),1,MPI_Integer, decomp%colComm, ierr)
    call MPI_AllGather(geom%jec,1,MPI_Integer,decomp%jend(0:)  ,1,MPI_Integer, decomp%colComm, ierr)

    ! Let other ranks know my row and column index
    allocate(decomp%MyRowGlobal(0:wsize-1), decomp%MyColGlobal(0:wsize-1))
    decomp%MyRowGlobal=-999; decomp%MyColGlobal=-999
    call MPI_AllGather(decomp%NSindex,1,MPI_Integer,decomp%MyRowGlobal,1,MPI_Integer, mpicomm, ierr)
    call MPI_AllGather(decomp%EWindex,1,MPI_Integer,decomp%MyColGlobal,1,MPI_Integer, mpicomm, ierr)

    ! Let other ranks know my rank in the row and column communicators
    call MPI_Comm_rank(decomp%rowComm, decomp%rowrank, ierr)
    call MPI_Comm_rank(decomp%colComm, decomp%colrank, ierr)
    allocate(decomp%MyRankInRowComm(0:wsize-1), decomp%MyRankInColComm(0:wsize-1))
    decomp%MyRankInRowComm=-999; decomp%MyRankInColComm=-999
    call MPI_AllGather(decomp%rowrank,1,MPI_Integer,decomp%MyRankInRowComm,1,MPI_Integer, mpicomm, ierr)
    call MPI_AllGather(decomp%colrank,1,MPI_Integer,decomp%MyRankInColComm,1,MPI_Integer, mpicomm, ierr)

    ! dimensions of my subdomain
    call MPI_Comm_rank(decomp%rowComm, myRowRank, ierr)  ! 0..rowSize-1
    call MPI_Comm_rank(decomp%colComm, myColRank, ierr)  ! 0..colSize-1
    decomp%localsizes(1) = decomp%iend(myRowRank) - decomp%ibegin(myRowRank) + 1
    decomp%localsizes(2) = decomp%jend(myColRank) - decomp%jbegin(myColRank) + 1

    ! Let other ranks in my row and column know how many rows and columns I have in my subdomain
    allocate(decomp%NumColsPerRank(0:rowSize-1))
    allocate(decomp%NumRowsPerRank(0:colSize-1))
    decomp%NumColsPerRank=-999; decomp%NumRowsPerRank=-999
    call MPI_Allgather(decomp%localsizes(1), 1, MPI_Integer, decomp%NumColsPerRank, 1, MPI_Integer, decomp%rowComm, ierr)
    call MPI_Allgather(decomp%localsizes(2), 1, MPI_Integer, decomp%NumRowsPerRank, 1, MPI_Integer, decomp%colComm, ierr)

    ! Horizontal dimensions of the files being read
    decomp%globalsizes(1) = geom%npx-1
    decomp%globalsizes(2) = geom%npy-1

    ! Allocate the structure
    ! ----------------------
    allocate(decomp%TwoPhaseScatter(0:npes-1))

    ! Initialize all MPI Datatype handles to NULL
    ! (Guarantees safe datatype creation and teardown checks)
    ! -------------------------------------------------------
    do r = 0, npes-1
      ! Scatter handles
      decomp%TwoPhaseScatter(r)%localvec_r4 = MPI_DATATYPE_NULL
      decomp%TwoPhaseScatter(r)%localvec_r8 = MPI_DATATYPE_NULL
      decomp%TwoPhaseScatter(r)%vec_r4      = MPI_DATATYPE_NULL
      decomp%TwoPhaseScatter(r)%vec_r8      = MPI_DATATYPE_NULL
      decomp%TwoPhaseScatter(r)%lalloc      = .false.
    enddo

    ! Allocate the workspaces
    allocate(decomp%scatter_workspace_r4(decomp%globalsizes(1), decomp%localsizes(2), batch_size))
    allocate(decomp%scatter_workspace_r8(decomp%globalsizes(1), decomp%localsizes(2), batch_size))
    allocate(decomp%scatter_recv_cast_workspace_r4(decomp%localsizes(1), decomp%localsizes(2), batch_size))

  end subroutine TwoPhaseScatter_init

  subroutine TwoPhaseGather_init(decomp, geom, npes, batch_size)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    type(fv3jedi_geom), intent(in):: geom
    integer, intent(in) :: npes, batch_size
    integer :: r, ierr
    integer :: rowSize, colSize
    integer :: myRowRank, myColRank
    integer :: mpicomm, wrank, wsize

    mpicomm = geom%f_comm%communicator()
    call MPI_Comm_rank(mpicomm, wrank, ierr)
    call MPI_Comm_size(mpicomm, wsize, ierr)

    ! A decomp is initialized once and freed by the matching delete routine
    if (decomp%rowComm /= MPI_COMM_NULL .or. decomp%colComm /= MPI_COMM_NULL) then
      write(6,'("TwoPhaseGather_init: decomp is already initialized; call TwoPhaseGather_delete first")')
      call flush(6)
      call MPI_Abort(mpicomm, 51, ierr)
    endif

    decomp%layout = geom%layout

    decomp%EWindex = modulo(wrank,decomp%layout(1))
    decomp%NSindex = (wrank/decomp%layout(1))

    call MPI_Comm_split(mpicomm, decomp%NSindex, wrank, decomp%rowComm, ierr)
    call MPI_Comm_split(mpicomm, decomp%EWindex, wrank, decomp%colComm, ierr)

    call MPI_Comm_size(decomp%rowComm, rowSize, ierr)
    call MPI_Comm_size(decomp%colComm, colSize, ierr)

    ! Allocate based on actual comm sizes
    allocate(decomp%ibegin(0:rowSize-1), decomp%iend(0:rowSize-1))
    allocate(decomp%jbegin(0:colSize-1), decomp%jend(0:colSize-1))

    call MPI_AllGather(geom%isc,1,MPI_Integer,decomp%ibegin(0:),1,MPI_Integer, decomp%rowComm, ierr)
    call MPI_AllGather(geom%iec,1,MPI_Integer,decomp%iend(0:)  ,1,MPI_Integer, decomp%rowComm, ierr)
    call MPI_AllGather(geom%jsc,1,MPI_Integer,decomp%jbegin(0:),1,MPI_Integer, decomp%colComm, ierr)
    call MPI_AllGather(geom%jec,1,MPI_Integer,decomp%jend(0:)  ,1,MPI_Integer, decomp%colComm, ierr)

    ! Let other ranks know my row and column index
    allocate(decomp%MyRowGlobal(0:wsize-1), decomp%MyColGlobal(0:wsize-1))
    decomp%MyRowGlobal=-999; decomp%MyColGlobal=-999
    call MPI_AllGather(decomp%NSindex,1,MPI_Integer,decomp%MyRowGlobal,1,MPI_Integer, mpicomm, ierr)
    call MPI_AllGather(decomp%EWindex,1,MPI_Integer,decomp%MyColGlobal,1,MPI_Integer, mpicomm, ierr)

    ! Let other ranks know my rank in the row and column communicators
    call MPI_Comm_rank(decomp%rowComm, decomp%rowrank, ierr)
    call MPI_Comm_rank(decomp%colComm, decomp%colrank, ierr)
    allocate(decomp%MyRankInRowComm(0:wsize-1), decomp%MyRankInColComm(0:wsize-1))
    decomp%MyRankInRowComm=-999; decomp%MyRankInColComm=-999
    call MPI_AllGather(decomp%rowrank,1,MPI_Integer,decomp%MyRankInRowComm,1,MPI_Integer, mpicomm, ierr)
    call MPI_AllGather(decomp%colrank,1,MPI_Integer,decomp%MyRankInColComm,1,MPI_Integer, mpicomm, ierr)

    ! dimensions of my subdomain
    call MPI_Comm_rank(decomp%rowComm, myRowRank, ierr)  ! 0..rowSize-1
    call MPI_Comm_rank(decomp%colComm, myColRank, ierr)  ! 0..colSize-1
    decomp%localsizes(1) = decomp%iend(myRowRank) - decomp%ibegin(myRowRank) + 1
    decomp%localsizes(2) = decomp%jend(myColRank) - decomp%jbegin(myColRank) + 1

    ! Let other ranks in my row and column know how many rows and columns I have in my subdomain
    allocate(decomp%NumColsPerRank(0:rowSize-1))
    allocate(decomp%NumRowsPerRank(0:colSize-1))
    decomp%NumColsPerRank=-999; decomp%NumRowsPerRank=-999
    call MPI_Allgather(decomp%localsizes(1), 1, MPI_Integer, decomp%NumColsPerRank, 1, MPI_Integer, decomp%rowComm, ierr)
    call MPI_Allgather(decomp%localsizes(2), 1, MPI_Integer, decomp%NumRowsPerRank, 1, MPI_Integer, decomp%colComm, ierr)

    ! Horizontal dimensions of the files being written
    decomp%globalsizes(1) = geom%npx-1
    decomp%globalsizes(2) = geom%npy-1

    ! Allocate the structure
    ! ----------------------
    allocate(decomp%TwoPhaseGather(0:npes-1))

    ! Initialize all MPI Datatype handles to NULL
    ! (Guarantees safe datatype creation and teardown checks)
    ! -------------------------------------------------------
    do r = 0, npes-1
      ! Gather handles
      decomp%TwoPhaseGather(r)%localvec_r4  = MPI_DATATYPE_NULL
      decomp%TwoPhaseGather(r)%localvec_r8  = MPI_DATATYPE_NULL
      decomp%TwoPhaseGather(r)%vec_r4       = MPI_DATATYPE_NULL
      decomp%TwoPhaseGather(r)%vec_r8       = MPI_DATATYPE_NULL
      decomp%TwoPhaseGather(r)%lalloc       = .false.
    enddo

    ! Allocate the workspaces
    allocate(decomp%gather_workspace_r4(decomp%globalsizes(1), decomp%localsizes(2), batch_size))
    allocate(decomp%gather_workspace_r8(decomp%globalsizes(1), decomp%localsizes(2), batch_size))
    allocate(decomp%gather_send_cast_workspace_r4(decomp%localsizes(1), decomp%localsizes(2), batch_size))
  end subroutine TwoPhaseGather_init

  ! Free everything TwoPhaseScatter_init created.  Call after the last MPI_Waitall.
  ! The geometry communicator is borrowed and is not freed.
  subroutine TwoPhaseScatter_delete(decomp)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer :: r, ierr

    if (allocated(decomp%TwoPhaseScatter)) then
      do r = lbound(decomp%TwoPhaseScatter, 1), ubound(decomp%TwoPhaseScatter, 1)
        associate(ts => decomp%TwoPhaseScatter(r))
        if (allocated(ts%sendcounts_phase1)) deallocate(ts%sendcounts_phase1)
        if (allocated(ts%senddispls_phase1)) deallocate(ts%senddispls_phase1)
        if (allocated(ts%sendcounts_phase2)) deallocate(ts%sendcounts_phase2)
        if (allocated(ts%senddispls_phase2)) deallocate(ts%senddispls_phase2)

        ! Safely free Datatypes ONLY if they were actually created
        if (ts%localvec_r4 /= MPI_DATATYPE_NULL) call MPI_Type_free(ts%localvec_r4, ierr)
        if (ts%localvec_r8 /= MPI_DATATYPE_NULL) call MPI_Type_free(ts%localvec_r8, ierr)

        if (ts%vec_r4 /= MPI_DATATYPE_NULL) call MPI_Type_free(ts%vec_r4, ierr)
        if (ts%vec_r8 /= MPI_DATATYPE_NULL) call MPI_Type_free(ts%vec_r8, ierr)
        end associate
      enddo
      deallocate(decomp%TwoPhaseScatter)
    endif
    call free_geometry(decomp)
    if (allocated(decomp%scatter_workspace_r4)) deallocate(decomp%scatter_workspace_r4)
    if (allocated(decomp%scatter_workspace_r8)) deallocate(decomp%scatter_workspace_r8)
    if (allocated(decomp%scatter_recv_cast_workspace_r4)) deallocate(decomp%scatter_recv_cast_workspace_r4)
  end subroutine TwoPhaseScatter_delete

  ! Free everything TwoPhaseGather_init created.  Call after the last MPI_Waitall.
  ! The geometry communicator is borrowed and is not freed.
  subroutine TwoPhaseGather_delete(decomp)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer :: r, ierr

    if (allocated(decomp%TwoPhaseGather)) then
      do r = lbound(decomp%TwoPhaseGather, 1), ubound(decomp%TwoPhaseGather, 1)
        associate(tg => decomp%TwoPhaseGather(r))
        if (allocated(tg%recvcounts_phase1)) deallocate(tg%recvcounts_phase1)
        if (allocated(tg%recvdispls_phase1)) deallocate(tg%recvdispls_phase1)
        if (allocated(tg%recvcounts_phase2)) deallocate(tg%recvcounts_phase2)
        if (allocated(tg%recvdispls_phase2)) deallocate(tg%recvdispls_phase2)

        ! Safely free Datatypes ONLY if they were actually created
        if (tg%localvec_r4 /= MPI_DATATYPE_NULL) call MPI_Type_free(tg%localvec_r4, ierr)
        if (tg%localvec_r8 /= MPI_DATATYPE_NULL) call MPI_Type_free(tg%localvec_r8, ierr)

        if (tg%vec_r4 /= MPI_DATATYPE_NULL) call MPI_Type_free(tg%vec_r4, ierr)
        if (tg%vec_r8 /= MPI_DATATYPE_NULL) call MPI_Type_free(tg%vec_r8, ierr)
        end associate
      enddo
      deallocate(decomp%TwoPhaseGather)
    endif
    call free_geometry(decomp)
    if (allocated(decomp%gather_workspace_r4))  deallocate(decomp%gather_workspace_r4)
    if (allocated(decomp%gather_workspace_r8))  deallocate(decomp%gather_workspace_r8)
    if (allocated(decomp%gather_send_cast_workspace_r4)) deallocate(decomp%gather_send_cast_workspace_r4)
  end subroutine TwoPhaseGather_delete

  ! Free the geometrical arrays and the row/column communicators (shared by both deletes)
  subroutine free_geometry(decomp)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer :: ierr

    if(allocated(decomp%ibegin)) deallocate(decomp%ibegin)
    if(allocated(decomp%iend)) deallocate(decomp%iend)
    if(allocated(decomp%jbegin)) deallocate(decomp%jbegin)
    if(allocated(decomp%jend)) deallocate(decomp%jend)
    if(allocated(decomp%MyRowGlobal)) deallocate(decomp%MyRowGlobal)
    if(allocated(decomp%MyColGlobal)) deallocate(decomp%MyColGlobal)
    if(allocated(decomp%MyRankInRowComm)) deallocate(decomp%MyRankInRowComm)
    if(allocated(decomp%MyRankInColComm)) deallocate(decomp%MyRankInColComm)
    if(allocated(decomp%NumColsPerRank)) deallocate(decomp%NumColsPerRank)
    if(allocated(decomp%NumRowsPerRank)) deallocate(decomp%NumRowsPerRank)

    ! MPI_Comm_free resets the handle to MPI_COMM_NULL
    if (decomp%rowComm /= MPI_COMM_NULL) call MPI_Comm_free(decomp%rowComm, ierr)
    if (decomp%colComm /= MPI_COMM_NULL) call MPI_Comm_free(decomp%colComm, ierr)
  end subroutine free_geometry

  subroutine TwoPhaseScatter_Phase1_r4(decomp, owner, rank, sendbuf, b_ind, req_p1)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=4), contiguous, target, asynchronous, intent(inout) :: sendbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p1

    integer :: row, col, ierr, temptype
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent
    integer :: myrow, mycol

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(ts => decomp%TwoPhaseScatter(owner))
    if (.not. ts%lalloc) then
      if (mycol == decomp%MyColGlobal(owner)) then
        allocate(ts%sendcounts_phase1(0:decomp%layout(2)-1))
        allocate(ts%senddispls_phase1(0:decomp%layout(2)-1))
        ts%senddispls_phase1(0) = 0
        ts%sendcounts_phase1(0) = decomp%globalsizes(1) * decomp%NumRowsPerRank(0)
        do row=1, decomp%layout(2)-1
          ts%sendcounts_phase1(row) = decomp%globalsizes(1) * decomp%NumRowsPerRank(row)
          ts%senddispls_phase1(row) = ts%senddispls_phase1(row-1) + ts%sendcounts_phase1(row-1)
        enddo
      endif
      allocate(ts%sendcounts_phase2(0:decomp%layout(1)-1))
      allocate(ts%senddispls_phase2(0:decomp%layout(1)-1))
      ts%senddispls_phase2(0) = 0
      ts%sendcounts_phase2(0) = decomp%NumColsPerRank(0)
      do col = 1, decomp%layout(1)-1
        ts%sendcounts_phase2(col) = decomp%NumColsPerRank(col)
        ts%senddispls_phase2(col) = ts%senddispls_phase2(col-1) + ts%sendcounts_phase2(col-1)
      enddo
      ts%lalloc = .true.
    endif

    if (ts%localvec_r4 == MPI_DATATYPE_NULL) then
      lb=0; extent=4
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%NumColsPerRank(mycol), MPI_REAL, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, ts%localvec_r4, ierr)
      call MPI_Type_commit(ts%localvec_r4, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (mycol == decomp%MyColGlobal(owner) .and. ts%vec_r4 == MPI_DATATYPE_NULL) then
      lb=0; extent=4
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%globalsizes(1), MPI_REAL, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, ts%vec_r4, ierr)
      call MPI_Type_commit(ts%vec_r4, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (mycol == decomp%MyColGlobal(owner)) then
      if (rank == owner) then
        call MPI_Iscatterv(sendbuf(1,1), ts%sendcounts_phase1, ts%senddispls_phase1, &
                           MPI_REAL, decomp%scatter_workspace_r4(1,1,b_ind), ts%sendcounts_phase1(myrow), &
                           MPI_REAL, decomp%MyRankInColComm(owner), decomp%colComm, req_p1, ierr)
      else
        call MPI_Iscatterv(MPI_BOTTOM, ts%sendcounts_phase1, ts%senddispls_phase1, &
                           MPI_REAL, decomp%scatter_workspace_r4(1,1,b_ind), ts%sendcounts_phase1(myrow), &
                           MPI_REAL, decomp%MyRankInColComm(owner), decomp%colComm, req_p1, ierr)
      endif
    else
      req_p1 = MPI_REQUEST_NULL
    endif
    end associate
  end subroutine TwoPhaseScatter_Phase1_r4


  subroutine TwoPhaseScatter_Phase1_r8(decomp, owner, rank, sendbuf, b_ind, req_p1)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=8), contiguous, target, asynchronous, intent(inout) :: sendbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p1

    integer :: row, col, ierr, temptype
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent
    integer :: myrow, mycol

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(ts => decomp%TwoPhaseScatter(owner))
    if (.not. ts%lalloc) then
      if (mycol == decomp%MyColGlobal(owner)) then
        allocate(ts%sendcounts_phase1(0:decomp%layout(2)-1))
        allocate(ts%senddispls_phase1(0:decomp%layout(2)-1))
        ts%senddispls_phase1(0) = 0
        ts%sendcounts_phase1(0) = decomp%globalsizes(1) * decomp%NumRowsPerRank(0)
        do row=1, decomp%layout(2)-1
          ts%sendcounts_phase1(row) = decomp%globalsizes(1) * decomp%NumRowsPerRank(row)
          ts%senddispls_phase1(row) = ts%senddispls_phase1(row-1) + ts%sendcounts_phase1(row-1)
        enddo
      endif

      allocate(ts%sendcounts_phase2(0:decomp%layout(1)-1))
      allocate(ts%senddispls_phase2(0:decomp%layout(1)-1))
      ts%senddispls_phase2(0) = 0
      ts%sendcounts_phase2(0) = decomp%NumColsPerRank(0)
      do col = 1, decomp%layout(1)-1
        ts%sendcounts_phase2(col) = decomp%NumColsPerRank(col)
        ts%senddispls_phase2(col) = ts%senddispls_phase2(col-1) + ts%sendcounts_phase2(col-1)
      enddo
      ts%lalloc = .true.
    endif

    if (ts%localvec_r8 == MPI_DATATYPE_NULL) then
      lb=0; extent=8
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%NumColsPerRank(mycol), MPI_DOUBLE_PRECISION, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, ts%localvec_r8, ierr)
      call MPI_Type_commit(ts%localvec_r8, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (mycol == decomp%MyColGlobal(owner) .and. ts%vec_r8 == MPI_DATATYPE_NULL) then
      lb=0; extent=8
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%globalsizes(1), MPI_DOUBLE_PRECISION, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, ts%vec_r8, ierr)
      call MPI_Type_commit(ts%vec_r8, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (mycol == decomp%MyColGlobal(owner)) then
      if (rank == owner) then
        call MPI_Iscatterv(sendbuf(1,1), ts%sendcounts_phase1, ts%senddispls_phase1, &
                           MPI_DOUBLE_PRECISION, decomp%scatter_workspace_r8(1,1,b_ind), ts%sendcounts_phase1(myrow), &
                           MPI_DOUBLE_PRECISION, decomp%MyRankInColComm(owner), decomp%colComm, req_p1, ierr)
      else
        call MPI_Iscatterv(MPI_BOTTOM, ts%sendcounts_phase1, ts%senddispls_phase1, &
                           MPI_DOUBLE_PRECISION, decomp%scatter_workspace_r8(1,1,b_ind), ts%sendcounts_phase1(myrow), &
                           MPI_DOUBLE_PRECISION, decomp%MyRankInColComm(owner), decomp%colComm, req_p1, ierr)
      endif
    else
      req_p1 = MPI_REQUEST_NULL
    endif
    end associate
  end subroutine TwoPhaseScatter_Phase1_r8


  subroutine TwoPhaseScatter_Phase2_r4(decomp, owner, rank, recvbuf, b_ind, req_p2)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=4), contiguous, target, asynchronous, intent(inout) :: recvbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p2

    integer :: ierr, vec
    integer :: mycol

    mycol = decomp%EWindex

    associate(ts => decomp%TwoPhaseScatter(owner))
    if (mycol == decomp%MyColGlobal(owner)) then
      vec = ts%vec_r4
      call MPI_Iscatterv(decomp%scatter_workspace_r4(1,1,b_ind), ts%sendcounts_phase2, &
                         ts%senddispls_phase2, vec, recvbuf(1,1), &
                         ts%sendcounts_phase2(mycol), ts%localvec_r4, &
                         decomp%MyRankInRowComm(owner), decomp%rowComm, req_p2, ierr)
    else
      call MPI_Iscatterv(MPI_BOTTOM, ts%sendcounts_phase2, ts%senddispls_phase2, &
                         MPI_REAL, recvbuf(1,1), ts%sendcounts_phase2(mycol), &
                         ts%localvec_r4, decomp%MyRankInRowComm(owner), decomp%rowComm, req_p2, ierr)
    endif
    end associate
  end subroutine TwoPhaseScatter_Phase2_r4


  subroutine TwoPhaseScatter_Phase2_r8(decomp, owner, rank, recvbuf, b_ind, req_p2)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=8), contiguous, target, asynchronous, intent(inout) :: recvbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p2

    integer :: ierr, vec
    integer :: mycol

    mycol = decomp%EWindex

    associate(ts => decomp%TwoPhaseScatter(owner))
    if (mycol == decomp%MyColGlobal(owner)) then
      vec = ts%vec_r8
      call MPI_Iscatterv(decomp%scatter_workspace_r8(1,1,b_ind), ts%sendcounts_phase2, &
                         ts%senddispls_phase2, vec, recvbuf(1,1), &
                         ts%sendcounts_phase2(mycol), ts%localvec_r8, &
                         decomp%MyRankInRowComm(owner), decomp%rowComm, req_p2, ierr)
    else
      call MPI_Iscatterv(MPI_BOTTOM, ts%sendcounts_phase2, ts%senddispls_phase2, &
                         MPI_DOUBLE_PRECISION, recvbuf(1,1), ts%sendcounts_phase2(mycol), &
                         ts%localvec_r8, decomp%MyRankInRowComm(owner), decomp%rowComm, req_p2, ierr)
    endif
    end associate
  end subroutine TwoPhaseScatter_Phase2_r8

  subroutine TwoPhaseGather_Phase1_r4(decomp, owner, rank, sendbuf, b_ind, req_p1)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=4), contiguous, target, asynchronous, intent(inout) :: sendbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p1

    integer :: col, row, ierr, temptype, vec, localvec
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent
    integer :: myrow, mycol

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(tg => decomp%TwoPhaseGather(owner))
    if (.not. tg%lalloc) then
      allocate(tg%recvcounts_phase1(0:decomp%layout(1)-1))
      allocate(tg%recvdispls_phase1(0:decomp%layout(1)-1))
      tg%recvdispls_phase1(0) = 0
      tg%recvcounts_phase1(0) = decomp%NumColsPerRank(0)
      do col = 1, decomp%layout(1)-1
        tg%recvcounts_phase1(col) = decomp%NumColsPerRank(col)
        tg%recvdispls_phase1(col) = tg%recvdispls_phase1(col-1) + tg%recvcounts_phase1(col-1)
      enddo

      if (mycol == decomp%MyColGlobal(owner)) then
        allocate(tg%recvcounts_phase2(0:decomp%layout(2)-1))
        allocate(tg%recvdispls_phase2(0:decomp%layout(2)-1))
        tg%recvdispls_phase2(0) = 0
        tg%recvcounts_phase2(0) = decomp%globalsizes(1) * decomp%NumRowsPerRank(0)
        do row=1, decomp%layout(2)-1
          tg%recvcounts_phase2(row) = decomp%globalsizes(1) * decomp%NumRowsPerRank(row)
          tg%recvdispls_phase2(row) = tg%recvdispls_phase2(row-1) + tg%recvcounts_phase2(row-1)
        enddo
      endif
      tg%lalloc = .true.
    endif

    if (tg%vec_r4 == MPI_DATATYPE_NULL) then
      lb=0; extent=4
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%globalsizes(1), MPI_REAL, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, tg%vec_r4, ierr)
      call MPI_Type_commit(tg%vec_r4, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (tg%localvec_r4 == MPI_DATATYPE_NULL) then
      lb=0; extent=4
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%NumColsPerRank(mycol), MPI_REAL, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, tg%localvec_r4, ierr)
      call MPI_Type_commit(tg%localvec_r4, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    vec = tg%vec_r4
    localvec = tg%localvec_r4

    if (mycol == decomp%MyColGlobal(owner)) then
      call MPI_Igatherv(sendbuf(1,1), tg%recvcounts_phase1(mycol), localvec, &
                        decomp%gather_workspace_r4(1,1,b_ind), tg%recvcounts_phase1, &
                        tg%recvdispls_phase1, vec, decomp%MyRankInRowComm(owner), &
                        decomp%rowComm, req_p1, ierr)
    else
      call MPI_Igatherv(sendbuf(1,1), tg%recvcounts_phase1(mycol), localvec, &
                        MPI_BOTTOM, tg%recvcounts_phase1, &
                        tg%recvdispls_phase1, vec, decomp%MyRankInRowComm(owner), &
                        decomp%rowComm, req_p1, ierr)
    endif
    end associate
  end subroutine TwoPhaseGather_Phase1_r4

  subroutine TwoPhaseGather_Phase1_r8(decomp, owner, rank, sendbuf, b_ind, req_p1)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=8), contiguous, target, asynchronous, intent(inout) :: sendbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p1

    integer :: col, row, ierr, temptype, vec, localvec
    integer(kind=MPI_ADDRESS_KIND) :: lb, extent
    integer :: myrow, mycol

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(tg => decomp%TwoPhaseGather(owner))
    if (.not. tg%lalloc) then
      allocate(tg%recvcounts_phase1(0:decomp%layout(1)-1))
      allocate(tg%recvdispls_phase1(0:decomp%layout(1)-1))
      tg%recvdispls_phase1(0) = 0
      tg%recvcounts_phase1(0) = decomp%NumColsPerRank(0)
      do col = 1, decomp%layout(1)-1
        tg%recvcounts_phase1(col) = decomp%NumColsPerRank(col)
        tg%recvdispls_phase1(col) = tg%recvdispls_phase1(col-1) + tg%recvcounts_phase1(col-1)
      enddo

      if (mycol == decomp%MyColGlobal(owner)) then
        allocate(tg%recvcounts_phase2(0:decomp%layout(2)-1))
        allocate(tg%recvdispls_phase2(0:decomp%layout(2)-1))
        tg%recvdispls_phase2(0) = 0
        tg%recvcounts_phase2(0) = decomp%globalsizes(1) * decomp%NumRowsPerRank(0)
        do row=1, decomp%layout(2)-1
          tg%recvcounts_phase2(row) = decomp%globalsizes(1) * decomp%NumRowsPerRank(row)
          tg%recvdispls_phase2(row) = tg%recvdispls_phase2(row-1) + tg%recvcounts_phase2(row-1)
        enddo
      endif
      tg%lalloc = .true.
    endif

    if (tg%vec_r8 == MPI_DATATYPE_NULL) then
      lb=0; extent=8
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%globalsizes(1), MPI_DOUBLE_PRECISION, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, tg%vec_r8, ierr)
      call MPI_Type_commit(tg%vec_r8, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    if (tg%localvec_r8 == MPI_DATATYPE_NULL) then
      lb=0; extent=8
      call MPI_Type_vector(decomp%localsizes(2), 1, decomp%NumColsPerRank(mycol), MPI_DOUBLE_PRECISION, temptype, ierr)
      call MPI_Type_create_resized(temptype, lb, extent, tg%localvec_r8, ierr)
      call MPI_Type_commit(tg%localvec_r8, ierr)
      call MPI_Type_free(temptype, ierr)
    endif

    vec = tg%vec_r8
    localvec = tg%localvec_r8

    if (mycol == decomp%MyColGlobal(owner)) then
      call MPI_Igatherv(sendbuf(1,1), tg%recvcounts_phase1(mycol), localvec, &
                        decomp%gather_workspace_r8(1,1,b_ind), tg%recvcounts_phase1, &
                        tg%recvdispls_phase1, vec, decomp%MyRankInRowComm(owner), &
                        decomp%rowComm, req_p1, ierr)
    else
      call MPI_Igatherv(sendbuf(1,1), tg%recvcounts_phase1(mycol), localvec, &
                        MPI_BOTTOM, tg%recvcounts_phase1, &
                        tg%recvdispls_phase1, vec, decomp%MyRankInRowComm(owner), &
                        decomp%rowComm, req_p1, ierr)
    endif
    end associate
  end subroutine TwoPhaseGather_Phase1_r8


  subroutine TwoPhaseGather_Phase2_r4(decomp, owner, rank, recvbuf, b_ind, req_p2)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=4), contiguous, target, asynchronous, intent(inout) :: recvbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p2

    integer :: myrow, mycol, ierr

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(tg => decomp%TwoPhaseGather(owner))
    if (mycol == decomp%MyColGlobal(owner)) then
      if (rank == owner) then
        call MPI_Igatherv(decomp%gather_workspace_r4(1,1,b_ind), tg%recvcounts_phase2(myrow), MPI_REAL, &
                          recvbuf(1,1), tg%recvcounts_phase2, &
                          tg%recvdispls_phase2, MPI_REAL, decomp%MyRankInColComm(owner), &
                          decomp%colComm, req_p2, ierr)
      else
        call MPI_Igatherv(decomp%gather_workspace_r4(1,1,b_ind), tg%recvcounts_phase2(myrow), MPI_REAL, &
                          MPI_BOTTOM, tg%recvcounts_phase2, &
                          tg%recvdispls_phase2, MPI_REAL, decomp%MyRankInColComm(owner), &
                          decomp%colComm, req_p2, ierr)
      endif
    else
      req_p2 = MPI_REQUEST_NULL
    endif
    end associate
  end subroutine TwoPhaseGather_Phase2_r4


  subroutine TwoPhaseGather_Phase2_r8(decomp, owner, rank, recvbuf, b_ind, req_p2)
    implicit none
    type(twophase_decomp), target, asynchronous, intent(inout) :: decomp
    integer, intent(in) :: owner, rank
    real(kind=8), contiguous, target, asynchronous, intent(inout) :: recvbuf(:,:)
    integer, intent(in)  :: b_ind
    integer, intent(out) :: req_p2

    integer :: myrow, mycol, ierr

    myrow = decomp%NSindex
    mycol = decomp%EWindex

    associate(tg => decomp%TwoPhaseGather(owner))
    if (mycol == decomp%MyColGlobal(owner)) then
      if (rank == owner) then
        call MPI_Igatherv(decomp%gather_workspace_r8(1,1,b_ind), tg%recvcounts_phase2(myrow), MPI_DOUBLE_PRECISION, &
                          recvbuf(1,1), tg%recvcounts_phase2, &
                          tg%recvdispls_phase2, MPI_DOUBLE_PRECISION, decomp%MyRankInColComm(owner), &
                          decomp%colComm, req_p2, ierr)
      else
        call MPI_Igatherv(decomp%gather_workspace_r8(1,1,b_ind), tg%recvcounts_phase2(myrow), MPI_DOUBLE_PRECISION, &
                          MPI_BOTTOM, tg%recvcounts_phase2, &
                          tg%recvdispls_phase2, MPI_DOUBLE_PRECISION, decomp%MyRankInColComm(owner), &
                          decomp%colComm, req_p2, ierr)
      endif
    else
      req_p2 = MPI_REQUEST_NULL
    endif
    end associate
  end subroutine TwoPhaseGather_Phase2_r8

end module TwoPhaseScatterGather
