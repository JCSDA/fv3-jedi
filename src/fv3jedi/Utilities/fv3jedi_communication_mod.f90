module fv3jedi_communication_mod

use mpi
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real

implicit none

private
public gather_field, scatter_field

!---------------------------------------------------------------------------------------------------

contains

!---------------------------------------------------------------------------------------------------

subroutine gather_field(geom,comm,gproc,field_in,field_out)

 implicit none

 type(fv3jedi_geom),   intent(in)  :: geom  !fv3-jedi geom
 integer,              intent(in)  :: comm  !MPI communicator
 integer,              intent(in)  :: gproc !Processor on which to gather
 real(kind=kind_real), intent(in)  :: field_in (geom%isc:geom%iec,geom%jsc:geom%jec)
 real(kind=kind_real), intent(out) :: field_out(1:geom%npx-1,1:geom%npy-1,6)

 integer :: comm_size, comm_rank, ierr, n, ji, jj, jc, npx_l, npy_l, npx_g, npy_g
 integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:), til_l(:)
 integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)
 real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)

 ! Get comm size
 ! -------------
 call mpi_comm_size(comm, comm_size, ierr)
 call mpi_comm_rank(comm, comm_rank, ierr)

 npx_l = geom%iec-geom%isc+1
 npy_l = geom%jec-geom%jsc+1
 npx_g = geom%npx-1
 npy_g = geom%npy-1

 ! Array of counts and displacement
 ! --------------------------------
 allocate(counts(comm_size), displs(comm_size))
 do n = 1,comm_size
    displs(n) = n-1
    counts(n) = 1
 enddo

 ! Gather local dimensions
 ! -----------------------
 allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size), til_l(comm_size))
 call mpi_allgatherv(geom%isc,   1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%iec,   1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%jsc,   1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%jec,   1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%ntile, 1, mpi_int, til_l, counts, displs, mpi_int, comm, ierr)
 deallocate(counts,displs)

 ! Gather counts and displacement and pack local array into vector
 ! ---------------------------------------------------------------
allocate(vectorcounts(comm_size), vectordispls(comm_size))

 n = 0
 do jc = 1,comm_size
   vectordispls(jc) = n
   do jj = jsc_l(jc),jec_l(jc)
     do ji = isc_l(jc),iec_l(jc)
       n = n+1
     enddo
   enddo
   vectorcounts(jc) = n - vectordispls(jc)
 enddo

 allocate(vector_l(npx_l*npy_l))
 n = 0
 do jj = geom%jsc,geom%jec
   do ji = geom%isc,geom%iec
     n = n+1
     vector_l(n) = field_in(ji,jj)
   enddo
 enddo

 ! Gather the full field
 allocate(vector_g(npx_g*npy_g*6))
 call mpi_gatherv( vector_l, npx_l*npy_l, mpi_double_precision, &
                   vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                   0, comm, ierr)
 deallocate(vector_l,vectorcounts,vectordispls)


 ! Unpack global vector into array
 ! -------------------------------

 if (comm_rank == 0) then
   n = 0
   do jc = 1,comm_size
     do jj = jsc_l(jc),jec_l(jc)
       do ji = isc_l(jc),iec_l(jc)
         n = n+1
         field_out(ji,jj,til_l(jc)) = vector_g(n)
       enddo
     enddo
   enddo
 endif
 deallocate(isc_l, iec_l, jsc_l, jec_l)

 deallocate(vector_g)

end subroutine gather_field

!---------------------------------------------------------------------------------------------------

subroutine scatter_field(geom,comm,gproc,field_in,field_out)

 implicit none

 type(fv3jedi_geom),   intent(in)  :: geom  !fv3-jedi geom
 integer,              intent(in)  :: comm  !MPI communicator
 integer,              intent(in)  :: gproc !Processor on which to gather
 real(kind=kind_real), intent(in)  :: field_in (1:geom%npx-1,1:geom%npy-1,6)
 real(kind=kind_real), intent(out) :: field_out(geom%isc:geom%iec,geom%jsc:geom%jec)

 integer :: comm_size, comm_rank, ierr, n, ji, jj, jc, npx_l, npy_l, npx_g, npy_g
 integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:), til_l(:)
 integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)
 real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)

 ! Get comm size
 ! -------------
 call mpi_comm_size(comm, comm_size, ierr)
 call mpi_comm_rank(comm, comm_rank, ierr)

 npx_l = geom%iec-geom%isc+1
 npy_l = geom%jec-geom%jsc+1
 npx_g = geom%npx-1
 npy_g = geom%npy-1

 ! Array of counts and displacement
 ! --------------------------------
 allocate(counts(comm_size), displs(comm_size))
 do n = 1,comm_size
    displs(n) = n-1
    counts(n) = 1
 enddo

 ! Gather local dimensions
 ! -----------------------
 allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size), til_l(comm_size))
 call mpi_allgatherv(geom%isc,   1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%iec,   1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%jsc,   1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%jec,   1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
 call mpi_allgatherv(geom%ntile, 1, mpi_int, til_l, counts, displs, mpi_int, comm, ierr)
 deallocate(counts,displs)


 ! Pack whole tile array into vector
 ! ---------------------------------
 allocate(vector_g(npx_g*npy_g*6))
 allocate(vectorcounts(comm_size), vectordispls(comm_size))

 n = 0
 do jc = 1,comm_size
   vectordispls(jc) = n
   do jj = jsc_l(jc),jec_l(jc)
     do ji = isc_l(jc),iec_l(jc)
       n = n+1
     enddo
   enddo
   vectorcounts(jc) = n - vectordispls(jc)
 enddo

 if (comm_rank == 0) then
   n = 0
   do jc = 1,comm_size
     do jj = jsc_l(jc),jec_l(jc)
       do ji = isc_l(jc),iec_l(jc)
         n = n+1
         vector_g(n) = field_in(ji,jj,til_l(jc))
       enddo
     enddo
   enddo
 endif

 deallocate(isc_l, iec_l, jsc_l, jec_l)

 ! Scatter tile array to processors
 allocate(vector_l(npx_l*npy_l))

 call mpi_scatterv( vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                    vector_l, npx_l*npy_l, mpi_double_precision, &
                    0, comm, ierr )

  deallocate(vector_g,vectorcounts,vectordispls)

  ! Unpack local vector into array
  ! ------------------------------
  n = 0
  do jj = geom%jsc,geom%jec
    do ji = geom%isc,geom%iec
      n = n+1
      field_out(ji,jj) = vector_l(n)
    enddo
  enddo

  deallocate(vector_l)

end subroutine scatter_field

!---------------------------------------------------------------------------------------------------

end module fv3jedi_communication_mod
