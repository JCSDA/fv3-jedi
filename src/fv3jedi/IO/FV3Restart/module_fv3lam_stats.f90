module module_ncfile_stat
!
!   PRGMMR: Ming Hu          ORG: GSL        DATE: 2022-03-08
!
! ABSTRACT:
!     This module read fv3lam fields and figure out dimension of each fields
!
! PROGRAM HISTORY LOG:
!
!   variable list
!
! USAGE:
!   INPUT FILES:
!   OUTPUT FILES:
!
! REMARKS:
!
! ATTRIBUTES:
!
!$$$
!
!_____________________________________________________________________

  implicit none

  integer,parameter :: max_varname_length=72
!
! Rset default to private
!

  private

  public :: ncfile_stat
  public :: inquire_var_shape

  type :: ncfile_stat

      integer :: numfiles
      character(len=512),allocatable :: filename(:)   ! full paths; 120 truncated long paths
      integer,allocatable :: numvarfile(:)
      integer :: numvar
      character(len=max_varname_length),allocatable :: list_varname(:)
      integer,allocatable :: num_dim(:)
      integer,allocatable :: dim_1(:)
      integer,allocatable :: dim_2(:)
      integer,allocatable :: dim_3(:)
      integer,allocatable :: vartype(:)
      integer :: num_totalvl

    contains
      procedure :: init
      procedure :: fill_dims
      procedure :: close
  end type ncfile_stat
!
! constants
!
contains

  subroutine init(this,numfiles,filein,numvar,varlist)
!                .      .    .                                       .
! subprogram:
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    implicit none

    integer, intent(in)          :: numfiles
    character(len=*), intent(in) :: filein(numfiles)
    integer, intent(in)          :: numvar(numfiles)
    character(len=*), intent(in) :: varlist(numfiles)
    class(ncfile_stat) :: this
!
    integer :: i,n1,n2
    character(len=500) :: varlistlocal
!
    this%numfiles=numfiles
    this%numvar=0
    if(this%numfiles>0) then
      allocate(this%filename(this%numfiles))
      allocate(this%numvarfile(this%numfiles))
      do i=1,numfiles
       this%filename(i)=trim(filein(i))
       this%numvarfile(i)=numvar(i)
       this%numvar=this%numvar+numvar(i)
      enddo
    endif
!
    if( this%numvar>0 ) then
       allocate(this%list_varname(this%numvar))
       allocate(this%dim_1(this%numvar))
       allocate(this%dim_2(this%numvar))
       allocate(this%dim_3(this%numvar))
       allocate(this%num_dim(this%numvar))
       allocate(this%vartype(this%numvar))
       this%dim_1=1
       this%dim_2=1
       this%dim_3=1
       this%num_dim=1
       this%vartype=0
       n1=1
       do i = 1, numfiles
          n2=n1+this%numvarfile(i)-1
          varlistlocal=trim(varlist(i))
          read(varlistlocal,*) this%list_varname(n1:n2)
          n1=n2+1
       end do
    endif
    this%num_totalvl=0

    do i=1,numfiles
      write(6,*) 'process file: ',trim(this%filename(i))
      write(6,*) 'variable numbers in this file: ',this%numvarfile(i)
    enddo
    write(6,*) 'variable will be process:',this%numvar
    write(6,*) this%list_varname
!
  end subroutine init

  subroutine close(this)
    implicit none
    class(ncfile_stat) :: this

    if(this%numfiles>0) then
      deallocate(this%filename)
      deallocate(this%numvarfile)
    endif

    if(this%numvar>0) then
       deallocate(this%list_varname)
       deallocate(this%dim_1)
       deallocate(this%dim_2)
       deallocate(this%dim_3)
       deallocate(this%num_dim)
       deallocate(this%vartype)
    endif
    this%numfiles=0
    this%numvar=0
      this%num_totalvl=0

  end subroutine close

  subroutine fill_dims(this)

    use netcdf, only: nf90_open,nf90_close,nf90_noerr,nf90_nowrite,nf90_strerror

    implicit none
    class(ncfile_stat) :: this

    integer :: ncid
    integer :: k
    character(len=max_varname_length) :: name
    character(len=256) :: errmsg
    integer :: iret
    logical :: ifexist
    integer :: n,lb,kk
!
!
    if(this%numfiles<=0) return

    lb=1
    do n=1,this%numfiles
       if(this%numvarfile(n) <=0) cycle

       ! A missing file is fatal: skipping it would pair every later variable with the
       ! wrong entry of list_varname
       inquire(file=trim(this%filename(n)),exist=ifexist )
       if(.not.ifexist) then
         write(6,*) 'fill_dims: file does not exist ',trim(this%filename(n))
         call fatal_error()
       else
         write(6,*) ' find dimensions for file ',trim(this%filename(n))
       endif
!
       iret=nf90_open(trim(this%filename(n)),nf90_nowrite,ncid)
       if(iret /= nf90_noerr) then
         write(6,*) 'fill_dims: cannot open ',trim(this%filename(n)),': ',trim(nf90_strerror(iret))
         call fatal_error()
       endif

       do kk=1,this%numvarfile(n)
          k=lb+kk-1
          name=this%list_varname(k)
          call inquire_var_shape(ncid, name, this%dim_1(k), this%dim_2(k), this%dim_3(k), &
                                 this%num_dim(k), this%vartype(k), iret, errmsg)
          if(iret /= 0) then
            write(6,'("fill_dims: ",3A)') trim(errmsg),' in file ',trim(this%filename(n))
            call fatal_error()
          endif
       enddo   ! end of loop for var list in a file

       iret=nf90_close(ncid)
       lb=lb+this%numvarfile(n)

    enddo  ! end of file loop

    this%num_totalvl=0
    do k=1,this%numvar
       this%num_totalvl=this%num_totalvl+this%dim_3(k)
    enddo

  end subroutine fill_dims

  ! Shape and storage type of one variable of an open file, interpreted the way the regional
  ! restart I/O stores fields (x, y, then levels):
  !   (x, y, z, t) with t = 1   -> dim_3 = z, num_dim = 3
  !   (x, y, z)    with z > 1   -> dim_3 = z, num_dim = 3
  !   (x, y, t)    with t = 1   -> dim_3 = 1, num_dim = 2
  !   (x, y)                    -> dim_3 = 1, num_dim = 2
  ! Returns status /= 0 and a message, without aborting, if the variable is missing, has
  ! another shape, or any netCDF inquiry fails, so each caller can abort on the
  ! communicator it runs on.
  subroutine inquire_var_shape(ncid, name, dim_1, dim_2, dim_3, num_dim, xtype, status, errmsg)
    use netcdf, only: nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension, &
                      nf90_noerr, nf90_strerror
    implicit none
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: dim_1, dim_2, dim_3, num_dim, xtype, status
    character(len=*), intent(out) :: errmsg

    integer :: varid, ndim, k, iret
    integer :: dim_id(4), dim_len(4)

    dim_1 = 1; dim_2 = 1; dim_3 = 1; num_dim = 1; xtype = 0
    status = 0
    errmsg = ''

    iret = nf90_inq_varid(ncid, trim(name), varid)
    if (iret /= nf90_noerr) then
      status = 1
      errmsg = 'variable '//trim(name)//' not found: '//trim(nf90_strerror(iret))
      return
    endif

    iret = nf90_inquire_variable(ncid, varid, ndims=ndim, xtype=xtype)
    if (iret /= nf90_noerr) then
      status = 4
      errmsg = 'cannot inquire variable '//trim(name)//': '//trim(nf90_strerror(iret))
      return
    endif
    if (ndim < 2 .or. ndim > 4) then
      status = 2
      write(errmsg,'("variable ",a," has ",i0," dimensions; expected 2, 3 or 4")') trim(name), ndim
      return
    endif

    iret = nf90_inquire_variable(ncid, varid, dimids=dim_id(1:ndim))
    if (iret /= nf90_noerr) then
      status = 4
      errmsg = 'cannot inquire the dimensions of variable '//trim(name)//': '//trim(nf90_strerror(iret))
      return
    endif
    do k = 1, ndim
      iret = nf90_inquire_dimension(ncid, dim_id(k), len=dim_len(k))
      if (iret /= nf90_noerr) then
        status = 4
        write(errmsg,'("cannot inquire dimension ",i0," of variable ",a,": ",a)') &
              k, trim(name), trim(nf90_strerror(iret))
        return
      endif
    enddo

    dim_1 = dim_len(1)
    dim_2 = dim_len(2)
    select case (ndim)
    case (4)
      if (dim_len(4) /= 1) then
        status = 3
        write(errmsg,'("variable ",a,": 4th dimension must be 1, is ",i0)') trim(name), dim_len(4)
        return
      endif
      dim_3 = dim_len(3)
      num_dim = 3
    case (3)
      if (dim_len(3) == 1) then
        dim_3 = 1
        num_dim = 2
      else
        dim_3 = dim_len(3)
        num_dim = 3
      endif
    case (2)
      dim_3 = 1
      num_dim = 2
    end select
  end subroutine inquire_var_shape

  ! fill_dims runs on one rank only while the others wait in the next collective, so a
  ! plain STOP would leave them hanging: abort the whole job instead.
  subroutine fatal_error()
    use mpi, only: MPI_COMM_WORLD
    implicit none
    integer :: ierr
    call flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 123, ierr)
  end subroutine fatal_error

end module module_ncfile_stat
