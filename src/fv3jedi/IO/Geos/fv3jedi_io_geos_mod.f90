! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_geos_mod

! libs
use mpi
use netcdf

! fckit
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module

! oops
use datetime_mod

! fv3-jedi
use fv3jedi_constants_mod,    only: rad2deg
use fv3jedi_geom_mod,         only: fv3jedi_geom
use fv3jedi_field_mod,        only: fv3jedi_field
use fv3jedi_io_utils_mod
use fv3jedi_kinds_mod,        only: kind_real
use fv3jedi_netcdf_utils_mod, only: nccheck

implicit none
private
public fv3jedi_io_geos

integer, parameter :: numfiles = 5 ! May change as more restarts are used

type fv3jedi_io_geos
 ! File names and paths
 integer :: ncid(numfiles)
 integer, allocatable :: ncid_forfield(:)
 logical :: ncid_isneeded(numfiles)
 ! Whether to expect/use the tile dimension
 logical :: tiledim(numfiles) = .false.
 logical:: restart(numfiles) = .false.
 ! Filename
 character(len=maxstring) :: datapath = ''
 character(len=maxstring) :: filenames(numfiles) = ''
 character(len=maxstring) :: filenames_conf(numfiles) = ''
 character(len=maxstring) :: filenames_default(numfiles) = ''
 ! Write extra meta data needed for GEOS ingest
 logical :: geosingestmeta = .false.
 ! Comms
 logical :: iam_io_proc
 type(fckit_mpi_comm) :: ccomm
 integer :: tcomm, ocomm       !Communicator for each tile and for output
 integer :: trank, tsize       !Tile come info
 integer :: crank, csize       !Component comm info
 integer :: orank, osize       !Output comm info
 ! IO ranges
 integer :: is_r3_tile(5), ic_r3_tile(5)
 integer :: is_r2_tile(4), ic_r2_tile(4)
 integer :: is_r3_noti(4), ic_r3_noti(4)
 integer :: is_r2_noti(3), ic_r2_noti(3)
 integer :: vindex_tile = 4
 integer :: vindex_noti = 3
 ! Clobber option
 logical :: clobber = .true.
 integer :: x_dimid, y_dimid, n_dimid, z_dimid, e_dimid, t_dimid, f_dimid, c_dimid, o_dimid
 ! Ps in file option
 logical :: ps_in_file = .false.
 ! Geometry copies
 integer :: isc, iec, jsc, jec
 integer :: npx, npy, npz, ntiles
 real(kind=kind_real), allocatable :: grid_lat(:,:), grid_lon(:,:)
 logical :: input_is_date_templated
 contains
  procedure :: create
  procedure :: delete
  procedure :: read
  procedure :: write
  final     :: dummy_final
end type fv3jedi_io_geos

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_io_geos),    intent(inout) :: self
type(fv3jedi_geom),        intent(in)    :: geom
type(fckit_configuration), intent(in)    :: conf

integer :: ierr, n, var

integer :: tileoffset, dt_in_name
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss

! Default file names
! ------------------
n = 0
n = n+1; self%filenames_default(n) = 'bkg'
n = n+1; self%filenames_default(n) = 'crtmsrf'
n = n+1; self%filenames_default(n) = 'fvcore_internal_rst'
n = n+1; self%filenames_default(n) = 'moist_internal_rst'
n = n+1; self%filenames_default(n) = 'surf_import_rst'

if (n .ne. numfiles) &
  call abor1_ftn("fv3jedi_io_geos_mod.setup: number of potential restart files &
                  does not match numfiles")

! Set default tile dim and restart flags
! --------------------------------------
self%tiledim(1) = .true.  ! History defult tile dim
self%tiledim(2) = .true.  ! History defult tile dim
self%tiledim(3) = .false. ! Restarts do not use tiledim
self%tiledim(4) = .false. ! Restarts do not use tiledim
self%tiledim(5) = .false. ! Restarts do not use tiledim

self%restart(1) = .false. !Is not a restart file
self%restart(2) = .false. !Is not a restart file
self%restart(3) = .true.  !Is a restart file
self%restart(4) = .true.  !Is a restart file
self%restart(5) = .true.  !Is a restart file

! Get configuration
! -----------------
call get_conf(self, conf)

! Config filenames to filenames
! -----------------------------
do n = 1,numfiles
  self%filenames(n) = trim(self%filenames_conf(n))
enddo

! Component communicator / get the main communicator from fv3 geometry
! --------------------------------------------------------------------
self%ccomm = geom%f_comm
self%csize = self%ccomm%size()
self%crank = self%ccomm%rank()

! Split comm and set flag for IO procs
! ------------------------------------
self%iam_io_proc = .true.

if (self%csize > 6) then

  ! Communicator for all procs on each tile. To communicate to tiles write proc
  call mpi_comm_split(self%ccomm%communicator(), geom%ntile, self%ccomm%rank(), self%tcomm, ierr)
  call mpi_comm_rank(self%tcomm, self%trank, ierr)
  call mpi_comm_size(self%tcomm, self%tsize, ierr)

  ! Communicator for procs that will write, one per tile
  call mpi_comm_split(self%ccomm%communicator(), self%trank, geom%ntile, self%ocomm, ierr)
  call mpi_comm_rank(self%ocomm, self%orank, ierr)
  call mpi_comm_size(self%ocomm, self%osize, ierr)

  if (self%trank .ne. 0) self%iam_io_proc = .false.

else

  ! For 6 cores output comm is componenent comm
  call mpi_comm_dup(self%ccomm%communicator(), self%ocomm, ierr)

endif

! Create local to this proc start/count
! -------------------------------------
if (self%iam_io_proc) then

  ! Starts/counts with tile dimension
  self%is_r3_tile(1) = 1;           self%ic_r3_tile(1) = geom%npx-1  !X
  self%is_r3_tile(2) = 1;           self%ic_r3_tile(2) = geom%npy-1  !Y
  self%is_r3_tile(3) = geom%ntile;  self%ic_r3_tile(3) = 1           !Tile
  self%is_r3_tile(4) = 1;           self%ic_r3_tile(4) = 1           !Lev
  self%is_r3_tile(5) = 1;           self%ic_r3_tile(5) = 1           !Time
  self%is_r2_tile(1) = 1;           self%ic_r2_tile(1) = geom%npx-1
  self%is_r2_tile(2) = 1;           self%ic_r2_tile(2) = geom%npy-1
  self%is_r2_tile(3) = geom%ntile;  self%ic_r2_tile(3) = 1
  self%is_r2_tile(4) = 1;           self%ic_r2_tile(4) = 1

  ! Starts/counts with no tile dimension
  tileoffset = (geom%ntile-1)*(6*(geom%npy-1)/geom%ntiles)
  self%is_r3_noti(1) = 1;              self%ic_r3_noti(1) = geom%npx-1
  self%is_r3_noti(2) = tileoffset+1;   self%ic_r3_noti(2) = geom%npy-1
  self%is_r3_noti(3) = 1;              self%ic_r3_noti(3) = 1
  self%is_r3_noti(4) = 1;              self%ic_r3_noti(4) = 1
  self%is_r2_noti(1) = 1;              self%ic_r2_noti(1) = geom%npx-1
  self%is_r2_noti(2) = tileoffset+1;   self%ic_r2_noti(2) = geom%npy-1
  self%is_r2_noti(3) = 1;              self%ic_r2_noti(3) = 1

endif

! Copy some geometry for later use
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npx = geom%npx
self%npy = geom%npy
self%npz = geom%npz
self%ntiles = geom%ntiles
allocate(self%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec))
allocate(self%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec))
self%grid_lat = rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec)
self%grid_lon = rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec)

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_io_geos), intent(inout) :: self

integer :: ierr, n

! Deallocate
! ----------
if (allocated(self%ncid_forfield)) deallocate(self%ncid_forfield)

if (allocated(self%grid_lat)) deallocate(self%grid_lat)
if (allocated(self%grid_lon)) deallocate(self%grid_lon)

! Release split comms
! -------------------
if (self%csize > 6) call MPI_Comm_free(self%tcomm, ierr)
call MPI_Comm_free(self%ocomm, ierr)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine read(self, vdate, calendar_type, date_init, fields)

class(fv3jedi_io_geos), intent(inout) :: self
type(datetime),         intent(inout) :: vdate
integer,                intent(inout) :: calendar_type
integer,                intent(inout) :: date_init(6)
type(fv3jedi_field),    intent(inout) :: fields(:)

! Overwrite any datetime templates in the file names
! --------------------------------------------------
if (self%input_is_date_templated) call setup_date(self, vdate)

! Read meta data
! --------------
call read_meta(self, vdate, calendar_type, date_init, fields)

! Read fields
! -----------
call read_fields(self, fields)

end subroutine read

! ------------------------------------------------------------------------------

subroutine write(self, fields, vdate)

class(fv3jedi_io_geos), intent(inout) :: self
type(fv3jedi_field),    intent(in)    :: fields(:)     !< Fields to be written
type(datetime),         intent(in)    :: vdate         !< DateTime

! Overwrite any datetime templates in the file names
! --------------------------------------------------
call setup_date(self, vdate)

! Set filenames
! -------------
call set_file_names(self, fields)

! Open/create files
! -----------------
call create_files(self)

! Write meta data
! ---------------
if (self%clobber) call write_meta(self, fields, vdate)

! Write fields
! ------------
call write_fields(self, fields, vdate)

! Close files
! -----------
call close_files(self)

end subroutine write

! ------------------------------------------------------------------------------

subroutine setup_date(self, vdate)

type(fv3jedi_io_geos), intent(inout) :: self
type(datetime),        intent(in)    :: vdate

integer :: n
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss

! Datetime to strings
! -------------------
call vdate_to_datestring(vdate, yyyy=yyyy, mm=mm, dd=dd, hh=hh, min=min, ss=ss)

do n = 1, numfiles

  ! Config filenames to filenames
  ! -----------------------------
  self%filenames(n) = trim(self%filenames_conf(n))

  ! Swap out datetime templates if needed
  if (index(self%filenames(n),"%yyyy") > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%yyyy',yyyy)
  if (index(self%filenames(n),"%mm"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%mm'  ,mm  )
  if (index(self%filenames(n),"%dd"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%dd'  ,dd  )
  if (index(self%filenames(n),"%hh"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%hh'  ,hh  )
  if (index(self%filenames(n),"%MM"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%MM'  ,min )
  if (index(self%filenames(n),"%ss"  ) > 0) &
    self%filenames(n) = replace_text(self%filenames(n),'%ss'  ,ss  )

enddo

end subroutine setup_date

! ------------------------------------------------------------------------------

subroutine get_conf(self, conf)

type(fv3jedi_io_geos),     intent(inout) :: self
type(fckit_configuration), intent(in)    :: conf

integer :: n

! Path where files are read from or saved to
! ------------------------------------------
call string_from_conf(conf,"datapath",self%datapath,'Data',memberswap=.true.)

! User can ask for extra meta data needed for GEOS ingest
! -------------------------------------------------------
if (.not. conf%get('geosingestmeta',self%geosingestmeta)) self%geosingestmeta = .false.

! Whether to expect/use the tile dimenstion in the file
! -----------------------------------------------------
if (.not. conf%get('clobber',self%clobber)) self%clobber = .true.

! Whether to expect/use the tile dimenstion in the file
! -----------------------------------------------------
if (.not. conf%get('tiledim',self%tiledim(1))) self%tiledim(1) = .true.
if (.not. conf%get('tiledim',self%tiledim(2))) self%tiledim(2) = .true.

! User can optionally specify the file names
! ------------------------------------------
n = 0
n = n+1; call string_from_conf(conf,"filename_bkgd",self%filenames_conf(1), &
                               self%filenames_default(1),memberswap=.true.)
n = n+1; call string_from_conf(conf,"filename_crtm",self%filenames_conf(2), &
                               self%filenames_default(2),memberswap=.true.)
n = n+1; call string_from_conf(conf,"filename_core",self%filenames_conf(3), &
                               self%filenames_default(3),memberswap=.true.)
n = n+1; call string_from_conf(conf,"filename_mois",self%filenames_conf(4), &
                               self%filenames_default(4),memberswap=.true.)
n = n+1; call string_from_conf(conf,"filename_surf",self%filenames_conf(5), &
                               self%filenames_default(5),memberswap=.true.)

! Optionally the file name to be read is datetime templated
! ---------------------------------------------------------
if (conf%has("filename is datetime templated")) then
  call conf%get_or_die("filename is datetime templated", self%input_is_date_templated)
else
  self%input_is_date_templated = .false.
endif

! Sanity check
! ------------
if (n .ne. numfiles) &
  call abor1_ftn("fv3jedi_io_geos_mod.get_conf: number of potential restart files &
                  does not match numfiles")

! Option to allow for ps infile
if (conf%has("psinfile")) call conf%get_or_die("psinfile",self%ps_in_file)

end subroutine get_conf

! ------------------------------------------------------------------------------

subroutine read_meta(self, vdate, calendar_type, date_init, fields)

type(fv3jedi_io_geos), intent(inout) :: self
type(datetime),        intent(inout) :: vdate         !< DateTime
integer,               intent(inout) :: calendar_type !< Calendar type
integer,               intent(inout) :: date_init(6)  !< Date intialized
type(fv3jedi_field),   intent(in)    :: fields(:)

integer :: varid, date(6), intdate, inttime, idate, isecs
character(len=8) :: cdate
character(len=6) :: ctime
integer :: cidate, cisecs, n, df

! Set filenames
! -------------
call set_file_names(self, fields)

! Open files
! ----------
call open_files(self)

calendar_type = -1
date_init = 0

idate = 0
isecs = 0

if (self%iam_io_proc) then

  ! Get time attributes
  do n = 1,numfiles
    if (self%ncid_isneeded(n)) then
      df = n
      exit
    endif
  enddo

  call nccheck ( nf90_inq_varid(self%ncid(df), "time", varid), "nf90_inq_varid time" )
  call nccheck ( nf90_get_att(self%ncid(df), varid, "begin_date", intdate), &
                 "nf90_get_att begin_date" )
  call nccheck ( nf90_get_att(self%ncid(df), varid, "begin_time", inttime), &
                 "nf90_get_att begin_time" )

  ! Pad with leading zeros if need be
  write(cdate,"(I0.8)") intdate
  write(ctime,"(I0.6)") inttime

  ! Back to integer
  read(cdate(1:4),*) date(1)
  read(cdate(5:6),*) date(2)
  read(cdate(7:8),*) date(3)
  read(ctime(1:2),*) date(4)
  read(ctime(3:4),*) date(5)
  read(ctime(5:6),*) date(6)

  ! To idate/isecs for Jedi
  idate = date(1)*10000 + date(2)*100 + date(3)
  isecs = date(4)*3600  + date(5)*60  + date(6)

endif

! Set the object date from the date of the file
call self%ccomm%broadcast(idate,0)
call self%ccomm%broadcast(isecs,0)
cidate = idate
cisecs = isecs
call datetime_from_ifs(vdate, cidate, cisecs)

! Close files
! -----------
call close_files(self)

end subroutine read_meta

! ------------------------------------------------------------------------------

subroutine read_fields(self, fields)

type(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_field),           intent(inout) :: fields(:)

integer :: varid, var, lev, ncid
logical :: tiledim
integer, pointer :: istart(:), icount(:)
integer, allocatable, target :: is_r3_tile(:), is_r3_noti(:)
real(kind=kind_real), allocatable :: arrayg(:,:), delp(:,:,:)

! Set filenames
! -------------
call set_file_names(self, fields)

! Open files
! ----------
call open_files(self)

! Local copy of starts for rank 3 in order to do one level at a time
! ------------------------------------------------------------------
allocate(is_r3_tile(size(self%is_r3_tile)))
is_r3_tile = self%is_r3_tile
allocate(is_r3_noti(size(self%is_r3_noti)))
is_r3_noti = self%is_r3_noti

! Array for level of whole tile
! -----------------------------
allocate(arrayg(1:self%npx-1,1:self%npy-1))


! Loop over fields
! ----------------
do var = 1,size(fields)

  ! ncid for this variable
  ncid = self%ncid(self%ncid_forfield(var))

  !If ps then switch to delp
  if (.not. self%ps_in_file) then
    if (trim(fields(var)%fv3jedi_name) == 'ps') then
      fields(var)%short_name = 'delp'
      fields(var)%npz = self%npz
      allocate(delp(self%isc:self%iec,self%jsc:self%jec,1:self%npz))
    endif
  endif

  ! Set pointers to the appropriate array ranges
  ! --------------------------------------------
  if (self%iam_io_proc) then

    tiledim = self%tiledim(self%ncid_forfield(var))

    if (associated(istart)) nullify(istart)
    if (associated(icount)) nullify(icount)

    if (fields(var)%npz == 1) then
      if (tiledim) then
        istart => self%is_r2_tile
        icount => self%ic_r2_tile
      else
        istart => self%is_r2_noti
        icount => self%ic_r2_noti
      endif
    elseif (fields(var)%npz > 1) then
      if (tiledim) then
        istart => is_r3_tile;
        icount => self%ic_r3_tile
      else
        istart => is_r3_noti
        icount => self%ic_r3_noti
      endif
    endif
  endif

  if (self%iam_io_proc) then
    call nccheck ( nf90_inq_varid (ncid, trim(fields(var)%short_name), varid), &
                  "nf90_inq_varid "//trim(fields(var)%short_name) )
  endif

  ! Loop over level and read the data
  ! ---------------------------------
  do lev = 1,fields(var)%npz

    arrayg = 0.0_kind_real

    if (self%iam_io_proc) then

      !Set start to current level
      is_r3_tile(self%vindex_tile) = lev
      is_r3_noti(self%vindex_noti) = lev

      ! Read the level
      call nccheck ( nf90_get_var( ncid, varid, arrayg, istart, icount), &
                    "nf90_get_var "//trim(fields(var)%short_name) )
    endif


    ! Scatter the field to all processors on the tile
    ! -----------------------------------------------
    if (.not. self%ps_in_file .and. trim(fields(var)%fv3jedi_name) == 'ps') then
      if (self%csize > 6) then
        call scatter_tile(self%isc, self%iec, self%jsc, self%jec, self%npx, self%npy, self%tcomm, &
                          1, arrayg, delp(self%isc:self%iec,self%jsc:self%jec,lev))
      else
        delp(self%isc:self%iec,self%jsc:self%jec,lev) = arrayg(self%isc:self%iec,self%jsc:self%jec)
      endif
    else
      if (self%csize > 6) then
        call scatter_tile(self%isc, self%iec, self%jsc, self%jec, self%npx, self%npy, self%tcomm, &
                          1, arrayg, fields(var)%array(self%isc:self%iec,self%jsc:self%jec,lev))
      else
        fields(var)%array(self%isc:self%iec,self%jsc:self%jec,lev) = &
                                                         arrayg(self%isc:self%iec,self%jsc:self%jec)
      endif
    endif

  enddo

  if (.not. self%ps_in_file .and. trim(fields(var)%fv3jedi_name) == 'ps') then
    fields(var)%short_name = 'ps'
    fields(var)%npz = 1
    fields(var)%array(:,:,1) = sum(delp,3)
    deallocate(delp)
  endif

enddo


! Deallocate locals
! -----------------
deallocate(is_r3_tile,is_r3_noti)
deallocate(arrayg)

! Close files
! -----------
call close_files(self)

end subroutine read_fields

! ------------------------------------------------------------------------------

subroutine write_meta(self, fields, vdate)

type(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_field),           intent(in)    :: fields(:)     !< Fields to be written
type(datetime),                intent(in)    :: vdate         !< DateTime

integer :: var, ymult, k, n, vc
character(len=15)  :: datefile
integer :: date(6), date8, time6
character(len=8)   :: date8s, cubesize
character(len=6)   :: time6s
character(len=4) :: XdimVar, YdimVar
integer :: varid(10000)

integer, pointer :: istart(:), icount(:)
integer, allocatable :: dimidsg(:), tiles(:), levels(:)
real(kind=kind_real), allocatable :: latg(:,:), long(:,:), xdimydim(:)


! Gathered lats/lons
! ------------------
allocate(latg(1:self%npx-1,1:self%npy-1))
allocate(long(1:self%npx-1,1:self%npy-1))

if (self%csize > 6) then
  call gather_tile(self%isc, self%iec, self%jsc, self%jec, self%npx, self%npy, self%tcomm, 1, &
                   self%grid_lat, latg)
  call gather_tile(self%isc, self%iec, self%jsc, self%jec, self%npx, self%npy, self%tcomm, 1, &
                   self%grid_lon, long)
else
  latg = self%grid_lat
  long = self%grid_lon
endif

write(cubesize,'(I8)') self%npx-1

! IO processors write the metadata
! --------------------------------
if (self%iam_io_proc) then

  ! Get datetime information ready to write
  ! ---------------------------------------
  call vdate_to_datestring(vdate,datest=datefile,date=date)
  write(date8s,'(I4,I0.2,I0.2)')   date(1),date(2),date(3)
  write(time6s,'(I0.2,I0.2,I0.2)') date(4),date(5),date(6)
  read(date8s,*) date8
  read(time6s,*) time6


  ! Xdim, Ydim arrays
  ! -----------------
  allocate(xdimydim(self%npx-1))
  do k = 1,self%npx-1
    xdimydim(k) = real(k,kind_real)
  enddo


  ! Tile array
  ! ----------
  allocate(tiles(6))
  do k = 1,6
    tiles(k) = k
  enddo


  ! Level and edge arrays
  ! ---------------------
  allocate(levels(self%npz+1))
  do k = 1,self%npz+1
    levels(k) = k
  enddo


  ! Loop over all files to be created/written to
  ! --------------------------------------------
  do n = 1, numfiles

    if (self%ncid_isneeded(n)) then

      ! Multiplication factor when no tile dimension
      ! --------------------------------------------
      ymult = 1
      if (.not. self%tiledim(n)) ymult = 6

      if ( self%tiledim(n) ) then
        XdimVar = 'Xdim'
        YdimVar = 'Ydim'
      else
        XdimVar = 'lon'
        YdimVar = 'lat'
      endif

      ! Main dimensions
      call nccheck ( nf90_def_dim(self%ncid(n), trim(XdimVar), self%npx-1,         self%x_dimid), &
                     "nf90_def_dim "//trim(XdimVar) )
      call nccheck ( nf90_def_dim(self%ncid(n), trim(YdimVar), ymult*(self%npy-1), self%y_dimid), &
                     "nf90_def_dim "//trim(YdimVar) )
      if (self%tiledim(n)) &
        call nccheck ( nf90_def_dim(self%ncid(n), "n",  self%ntiles, self%n_dimid), "nf90_def_dim n"    )
      call nccheck ( nf90_def_dim(self%ncid(n), "lev",  self%npz,    self%z_dimid), "nf90_def_dim lev"  )
      call nccheck ( nf90_def_dim(self%ncid(n), "edge", self%npz+1,  self%e_dimid), "nf90_def_dim edge" )
      call nccheck ( nf90_def_dim(self%ncid(n), "time", 1,           self%t_dimid), "nf90_def_dim time" )

      ! In case the four level surface fields need to be written
      do var = 1,size(fields)
        if (fields(var)%npz == 4) then
          call nccheck ( nf90_def_dim(self%ncid(n), "lev4", 4, self%f_dimid), "nf90_def_dim lev"  )
          exit
        endif
      enddo

      !Needed by GEOS for ingesting cube sphere field
      call nccheck ( nf90_def_dim(self%ncid(n), "ncontact",          4, self%c_dimid), &
                     "nf90_def_dim ncontact" )
      call nccheck ( nf90_def_dim(self%ncid(n), "orientationStrLen", 5, self%o_dimid), &
                     "nf90_def_dim orientationStrLend" )

      ! Dimension ID array for lat/lon arrays
      if (allocated(dimidsg)) deallocate(dimidsg)
      if ( self%tiledim(n) ) then
        allocate(dimidsg(3))
        dimidsg(:) = (/ self%x_dimid, self%y_dimid, self%n_dimid /)
      else
        allocate(dimidsg(2))
        dimidsg(:) = (/ self%x_dimid, self%y_dimid /)
      endif

      ! Define fields to be written
      vc=0;

      if (self%tiledim(n)) then
        vc=vc+1;
        call nccheck( nf90_def_var(self%ncid(n), "n", NF90_INT, self%n_dimid, varid(vc)), "nf90_def_var n" )
        call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "cubed-sphere face") )
        call nccheck( nf90_put_att(self%ncid(n), varid(vc), "axis", "e") )
        call nccheck( nf90_put_att(self%ncid(n), varid(vc), "grads_dim", "e") )
      endif

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), trim(XdimVar), NF90_DOUBLE, self%x_dimid, varid(vc)), &
                    "nf90_def_var "//trim(XdimVar) )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "Fake Longitude for GrADS Compatibility") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "degrees_east") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), trim(YdimVar), NF90_DOUBLE, self%y_dimid, varid(vc)), &
                    "nf90_def_var "//trim(YdimVar) )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "Fake Latitude for GrADS Compatibility") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "degrees_north") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), "lons", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lons" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "longitude") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "degrees_east") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), "lats", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lats" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "latitude") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "degrees_north") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), "lev", NF90_DOUBLE, self%z_dimid, varid(vc)), "nf90_def_var lev" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "vertical level") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "layer") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "positive", "down") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "coordinate", "eta") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "standard_name", "model_layers") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), "edge", NF90_DOUBLE, self%e_dimid, varid(vc)), "nf90_def_var edge" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "vertical level edges") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "units", "layer") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "positive", "down") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "coordinate", "eta") )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "standard_name", "model_layers") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(n), "time", NF90_INT, self%t_dimid, varid(vc)), "nf90_def_var time" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "long_name", "time"), "nf90_def_var time long_name" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "begin_date", date8), "nf90_def_var time begin_date" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "begin_time", time6), "nf90_def_var time begin_time" )

      vc=vc+1; !(Needed by GEOS to ingest cube sphere analysis)
      call nccheck( nf90_def_var(self%ncid(n), "cubed_sphere", NF90_CHAR, varid(vc)), &
                    "nf90_def_var cubed_sphere" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "grid_mapping_name", "gnomonic cubed-sphere"), &
                    "nf90_def_var time grid_mapping_name" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "file_format_version", "2.90"), &
                    "nf90_def_var time file_format_version" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "additional_vars", "contacts,orientation,anchor"), &
                    "nf90_def_var time additional_vars" )
      call nccheck( nf90_put_att(self%ncid(n), varid(vc), "gridspec_file", "C"//trim(cubesize)//"_gridspec.nc4"), &
                    "nf90_def_var gridspec_file" )

      !vc=vc+1; !(Needed by GEOS to ingest cube sphere analysis)
      !call nccheck( nf90_def_var(self%ncid(n), "ncontact", NF90_INT, varid(vc)), "nf90_def_var ncontact" )


      ! End define mode
      ! ---------------
      call nccheck( nf90_enddef(self%ncid(n)), "nf90_enddef" )


      ! Reset counter
      ! -------------
      vc=0


      ! Write metadata arrays
      ! ---------------------

      ! Tiles
      if (self%tiledim(n)) then
        vc=vc+1
        call nccheck( nf90_put_var( self%ncid(n), varid(vc), tiles ), "nf90_put_var n" )
      endif

      ! Xdim & Ydim arrays
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), xdimydim ), "nf90_put_var "//trim(XdimVar) )
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), xdimydim ), "nf90_put_var "//trim(YdimVar) )

      ! Start/counts
      if (associated(istart)) nullify(istart)
      if (associated(icount)) nullify(icount)
      if (self%tiledim(n)) then
        istart => self%is_r2_tile(1:3); icount => self%ic_r2_tile(1:3)
      else
        istart => self%is_r2_noti(1:2); icount => self%ic_r2_noti(1:2)
      endif

      ! Lat/lon arrays
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), long, &
                                                  start = istart, &
                                                  count = icount ), &
                                                  "nf90_put_var lons" )

      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), latg, &
                                                  start = istart, &
                                                  count = icount ), &
                                                  "nf90_put_var lats" )

      ! Write model levels & time
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), levels(1:self%npz) ), "nf90_put_var lev" )
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), levels ), "nf90_put_var edge" )

      ! Time
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(n), varid(vc), 0 ), "nf90_put_var time" )

    endif

  enddo

endif

end subroutine write_meta

! ------------------------------------------------------------------------------

subroutine write_fields(self, fields, vdate)

type(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_field),           intent(in)    :: fields(:)     !< Fields to be written
type(datetime),                intent(in)    :: vdate         !< DateTime

integer :: var, lev, n, ncid, filei, varid
integer, target :: dimids2_tile(4), dimids3_tile(5), dimidse_tile(5), dimids4_tile(5)
integer, target :: dimids2_noti(3), dimids3_noti(4), dimidse_noti(4), dimids4_noti(4)
integer, pointer :: dimids2(:), dimids3(:), dimidse(:), dimids4(:), dimids(:)
integer, allocatable, target :: is_r3_tile(:), is_r3_noti(:)
real(kind=kind_real), allocatable :: arrayg(:,:)
integer, pointer :: istart(:), icount(:)


! Whole level of tile array
! -------------------------
allocate(arrayg(1:self%npx-1,1:self%npy-1))


! Local counts for 3 to allow changing start point
! ------------------------------------------------
allocate(is_r3_tile(size(self%is_r3_tile)))
is_r3_tile = self%is_r3_tile
allocate(is_r3_noti(size(self%is_r3_noti)))
is_r3_noti = self%is_r3_noti


! Dimension ID arrays for the various fields with and without tile dimension
! --------------------------------------------------------------------------
dimids2_tile = (/ self%x_dimid, self%y_dimid, self%n_dimid,               self%t_dimid /)
dimids3_tile = (/ self%x_dimid, self%y_dimid, self%n_dimid, self%z_dimid, self%t_dimid /)
dimidse_tile = (/ self%x_dimid, self%y_dimid, self%n_dimid, self%e_dimid, self%t_dimid /)
dimids4_tile = (/ self%x_dimid, self%y_dimid, self%n_dimid, self%f_dimid, self%t_dimid /)

dimids2_noti = (/ self%x_dimid, self%y_dimid,                             self%t_dimid /)
dimids3_noti = (/ self%x_dimid, self%y_dimid,               self%z_dimid, self%t_dimid /)
dimidse_noti = (/ self%x_dimid, self%y_dimid,               self%e_dimid, self%t_dimid /)
dimids4_noti = (/ self%x_dimid, self%y_dimid,               self%f_dimid, self%t_dimid /)


! Loop over the fields
! --------------------
do var = 1,size(fields)

  if (self%iam_io_proc) then

    ! ncid for this field
    ! -------------------
    filei = self%ncid_forfield(var)
    ncid = self%ncid(filei)

    ! Redefine
    ! --------
    if (self%clobber) then

      call nccheck( nf90_redef(ncid), "nf90_enddef" )

      ! Dimension IDs for this field
      ! ----------------------------
      if (associated(dimids2)) nullify(dimids2)
      if (associated(dimids3)) nullify(dimids3)
      if (associated(dimidse)) nullify(dimidse)
      if (associated(dimids4)) nullify(dimids4)

      if (self%tiledim(filei)) then
        dimids2 => dimids2_tile
        dimids3 => dimids3_tile
        dimidse => dimidse_tile
        dimids4 => dimids4_tile
      else
        dimids2 => dimids2_noti
        dimids3 => dimids3_noti
        dimidse => dimidse_noti
        dimids4 => dimids4_noti
      endif

      if (associated(dimids)) nullify (dimids)

      if (fields(var)%npz == 1) then
        dimids => dimids2
      elseif (fields(var)%npz == self%npz) then
        dimids => dimids3
      elseif (fields(var)%npz == self%npz+1) then
        dimids => dimidse
      elseif (fields(var)%npz == 4) then
        dimids => dimids4
      else
        call abor1_ftn("write_geos: vertical dimension not supported")
      endif

      ! Define field
      call nccheck( nf90_def_var(ncid, trim(fields(var)%short_name), NF90_DOUBLE, dimids, varid), &
                    "nf90_def_var"//trim(fields(var)%short_name))

      ! Write attributes if clobbering
      if (self%clobber) then

        ! Long name and units
        call nccheck( nf90_put_att(ncid, varid, "long_name"    , trim(fields(var)%long_name) ), "nf90_put_att" )
        call nccheck( nf90_put_att(ncid, varid, "units"        , trim(fields(var)%units)     ), "nf90_put_att" )

        ! Additional attributes for history and or plotting compatibility
        if (.not.self%restart(filei)) then
          call nccheck( nf90_put_att(ncid, varid, "standard_name", trim(fields(var)%long_name) ), "nf90_put_att" )
          call nccheck( nf90_put_att(ncid, varid, "coordinates"  , "lons lats"                 ), "nf90_put_att" )
          call nccheck( nf90_put_att(ncid, varid, "grid_mapping" , "cubed_sphere"              ), "nf90_put_att" )
        endif

      endif

      ! End define mode
      call nccheck( nf90_enddef(ncid), "nf90_enddef" )

    else

      ! Get existing variable id to write to
      call nccheck ( nf90_inq_varid (ncid, trim(fields(var)%short_name), varid), &
                    "nf90_inq_varid "//trim(fields(var)%short_name) )

    endif

    ! Set starts and counts based on levels and tiledim flag
    ! ------------------------------------------------------

    if (associated(istart)) nullify(istart)
    if (associated(icount)) nullify(icount)

    if (fields(var)%npz == 1) then
      if (self%tiledim(filei)) then
        istart => self%is_r2_tile; icount => self%ic_r2_tile
      else
        istart => self%is_r2_noti; icount => self%ic_r2_noti
      endif
    elseif (fields(var)%npz > 1) then
      if (self%tiledim(filei)) then
        istart => is_r3_tile; icount => self%ic_r3_tile
      else
        istart => is_r3_noti; icount => self%ic_r3_noti
      endif
    endif

  endif

  do lev = 1,fields(var)%npz

    if (self%csize > 6) then
      call gather_tile(self%isc, self%iec, self%jsc, self%jec, self%npx, self%npy, self%tcomm, 1, &
                       fields(var)%array(self%isc:self%iec,self%jsc:self%jec,lev), arrayg)
    else
      arrayg = fields(var)%array(self%isc:self%iec,self%jsc:self%jec,lev)
    endif

    if (self%iam_io_proc) then

      is_r3_tile(self%vindex_tile) = lev
      is_r3_noti(self%vindex_noti) = lev

      call nccheck( nf90_put_var( ncid, varid, arrayg, start = istart, count = icount ), &
                                  "nf90_put_var "//trim(fields(var)%short_name) )

    endif

  enddo

enddo

! Deallocate locals
! -----------------
deallocate(is_r3_tile,is_r3_noti)
deallocate(arrayg)

end subroutine write_fields

! ------------------------------------------------------------------------------

subroutine gather_tile(isc, iec, jsc, jec, npx, npy, comm, nlev, array_l, array_g)

integer,              intent(in)    :: isc, iec, jsc, jec, npx, npy
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_l(isc:iec,jsc:jec,1:nlev)  ! Local array
real(kind=kind_real), intent(inout) :: array_g(1:npx-1,1:npy-1,1:nlev)            ! Gathered array (only valid on root)

real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)
integer :: comm_size, ierr
integer :: ji, jj, jk, jc, n
integer :: npx_g, npy_g, npx_l, npy_l
integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:)
integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)

!Get comm size
call mpi_comm_size(comm, comm_size, ierr)

!Array of counts and displacement
allocate(counts(comm_size), displs(comm_size))
do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

!Horizontal size for global and local
npx_g = npx-1
npy_g = npy-1
npx_l = iec-isc+1
npy_l = jec-jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
deallocate(counts,displs)

! Pack whole tile array into vector
allocate(vectorcounts(comm_size), vectordispls(comm_size))

!Gather counts and displacement
n = 0
do jc = 1,comm_size
  vectordispls(jc) = n
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
      enddo
    enddo
  enddo
  vectorcounts(jc) = n - vectordispls(jc)
enddo

! Pack local array into vector
allocate(vector_l(npx_l*npy_l*nlev))
n = 0
do jk = 1,nlev
  do jj = jsc,jec
    do ji = isc,iec
      n = n+1
      vector_l(n) = array_l(ji,jj,jk)
    enddo
  enddo
enddo

! Gather the full field
allocate(vector_g(npx_g*npy_g*nlev))
call mpi_gatherv( vector_l, npx_l*npy_l, mpi_double_precision, &
                  vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                  0, comm, ierr)
deallocate(vector_l,vectorcounts,vectordispls)

!Unpack global vector into array
n = 0
do jc = 1,comm_size
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
        array_g(ji,jj,jk) = vector_g(n)
      enddo
    enddo
  enddo
enddo
deallocate(isc_l, iec_l, jsc_l, jec_l)

deallocate(vector_g)

end subroutine gather_tile

! ------------------------------------------------------------------------------

subroutine scatter_tile(isc, iec, jsc, jec, npx, npy, comm, nlev, array_g, array_l)

integer,              intent(in)    :: isc, iec, jsc, jec, npx, npy
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_g(1:npx-1,1:npy-1,nlev)            ! Gathered array (only valid on root)
real(kind=kind_real), intent(inout) :: array_l(isc:iec,jsc:jec,nlev)  ! Local array

real(kind=kind_real), allocatable :: vector_g(:), vector_l(:)
integer :: comm_size, ierr
integer :: ji, jj, jk, jc, n
integer :: npx_g, npy_g, npx_l, npy_l
integer, allocatable :: isc_l(:), iec_l(:), jsc_l(:), jec_l(:)
integer, allocatable :: counts(:), displs(:), vectorcounts(:), vectordispls(:)

!Get comm size
call mpi_comm_size(comm, comm_size, ierr)

!Array of counts and displacement
allocate(counts(comm_size), displs(comm_size))
do n = 1,comm_size
   displs(n) = n-1
   counts(n) = 1
enddo

!Horizontal size for global and local
npx_g = npx-1
npy_g = npy-1
npx_l = iec-isc+1
npy_l = jec-jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
deallocate(counts,displs)

! Pack whole tile array into vector
allocate(vector_g(npx_g*npy_g*nlev))
allocate(vectorcounts(comm_size), vectordispls(comm_size))

n = 0
do jc = 1,comm_size
  vectordispls(jc) = n
  do jk = 1,nlev
    do jj = jsc_l(jc),jec_l(jc)
      do ji = isc_l(jc),iec_l(jc)
        n = n+1
        vector_g(n) = array_g(ji,jj,jk)
      enddo
    enddo
  enddo
  vectorcounts(jc) = n - vectordispls(jc)
enddo
deallocate(isc_l, iec_l, jsc_l, jec_l)

! Scatter tile array to processors
allocate(vector_l(npx_l*npy_l*nlev))

call mpi_scatterv( vector_g, vectorcounts, vectordispls, mpi_double_precision, &
                   vector_l, npx_l*npy_l, mpi_double_precision, &
                   0, comm, ierr )

deallocate(vector_g,vectorcounts,vectordispls)

! Unpack local vector into array
n = 0
do jk = 1,nlev
  do jj = jsc,jec
    do ji = isc,iec
      n = n+1
      array_l(ji,jj,jk) = vector_l(n)
    enddo
  enddo
enddo
deallocate(vector_l)

end subroutine scatter_tile

! --------------------------------------------------------------------------------------------------

! Not really needed but prevents gnu compiler bug
subroutine dummy_final(self)
type(fv3jedi_io_geos), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

subroutine open_files(self)

type(fv3jedi_io_geos), intent(inout) :: self

integer :: n

self%ncid = -1

! Open files for reading
do n = 1,numfiles
  if (self%ncid_isneeded(n)) then
    call nccheck ( nf90_open( trim(self%datapath)//'/'//trim(self%filenames(n)), NF90_NOWRITE, &
                   self%ncid(n) ), "nf90_open"//trim(self%filenames(n)) )
  endif
enddo

end subroutine open_files

! --------------------------------------------------------------------------------------------------

subroutine create_files(self)

type(fv3jedi_io_geos), intent(inout) :: self

integer :: n, fileopts

self%ncid = -1

if (self%iam_io_proc) then

  if (self%clobber) then

    fileopts = ior(NF90_NETCDF4, NF90_MPIIO)

    ! Create/open files for writing
    do n = 1,numfiles
      if (self%ncid_isneeded(n)) then
        call nccheck( nf90_create( trim(self%datapath)//'/'//trim(self%filenames(n)), fileopts, &
                                   self%ncid(n), comm = self%ocomm, info = MPI_INFO_NULL), &
                                   "nf90_create"//trim(self%filenames(n)) )
      endif
    enddo

  else

    ! Open files for writing
    do n = 1,numfiles

      if (self%ncid_isneeded(n)) then
        call nccheck ( nf90_open( trim(self%datapath)//'/'//trim(self%filenames(n)), NF90_WRITE, &
                       self%ncid(n) ), "nf90_open"//trim(self%filenames(n)) )
      endif
    enddo

  endif

endif

end subroutine create_files

! --------------------------------------------------------------------------------------------------

subroutine close_files(self)

type(fv3jedi_io_geos), intent(inout) :: self

integer :: n

! Close the files
! ---------------
if (self%iam_io_proc) then
  do n = 1,numfiles
    if (self%ncid_isneeded(n)) then
      call nccheck ( nf90_close(self%ncid(n)), "nf90_close" )
    endif
  enddo
endif

end subroutine close_files

! --------------------------------------------------------------------------------------------------

subroutine set_file_names(self, fields)

type(fv3jedi_io_geos), intent(inout) :: self
type(fv3jedi_field),   intent(in)    :: fields(:)

integer :: var, grp, indfile

if (allocated(self%ncid_forfield)) deallocate(self%ncid_forfield)
allocate(self%ncid_forfield(size(fields)))
self%ncid_forfield = -1
self%ncid_isneeded = .false.

! Set files for the fields
! ------------------------
do var = 1,size(fields)

  select case (trim(fields(var)%short_name))

  case("vtype","stype","vfrac")                       ! CRTM surface
    grp = 2
  case("U","V","W","PT","PKZ","PE","DZ")              ! GEOS fv_core restart
    grp = 3
  case("Q","QILS","QICN","QLLS","QLCN","CLLS","CLCN") ! GEOS moist restart
    grp = 4
  case("PHIS")                                        ! GEOS surface restart
    grp = 5
  case default                                        ! Background file
    grp = 1

  endselect

  self%ncid_forfield(var) = grp
  self%ncid_isneeded(grp) = .true.

enddo

end subroutine set_file_names

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_geos_mod
