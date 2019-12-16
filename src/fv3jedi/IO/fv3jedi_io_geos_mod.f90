! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_geos_mod

use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use iso_c_binding

use string_utils, only: swap_name_member

use fckit_mpi_module

use fv3jedi_constants_mod,    only: rad2deg
use fv3jedi_geom_mod,         only: fv3jedi_geom
use fv3jedi_field_mod,        only: fv3jedi_field
use fv3jedi_kinds_mod,        only: kind_real
use fv3jedi_netcdf_utils_mod, only: nccheck

use mpi
use netcdf

implicit none
private
public fv3jedi_io_geos

integer, parameter :: maxstring = 2048

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
 logical :: clobber = .false.
 integer :: x_dimid, y_dimid, n_dimid, z_dimid, e_dimid, t_dimid, f_dimid, c_dimid, o_dimid
 contains
  procedure :: setup
  procedure :: delete
  procedure :: read_meta
  procedure :: read_fields
  procedure :: write_all
  final     :: dummy_final
end type fv3jedi_io_geos

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine setup(self, geom, fields, vdate, readorwrite, f_conf)

class(fv3jedi_io_geos),              intent(inout) :: self
type(fv3jedi_geom),                  intent(in)    :: geom
type(fv3jedi_field),                 intent(in)    :: fields(:)
type(datetime),                      intent(in)    :: vdate
character(len=*),                    intent(in)    :: readorwrite
type(fckit_configuration), optional, intent(in)    :: f_conf

integer :: ierr, nf, var, fileopts

integer :: tileoffset, dt_in_name
character(len=4) :: yyyy
character(len=2) :: mm, dd, hh, min, ss


! Allocatable arrays
! ------------------
allocate(self%ncid_forfield(size(fields)))

self%ncid = -1
self%ncid_isneeded = .false.
self%ncid_forfield = -1


! Default file names
! ------------------
nf = 0
nf = nf+1; self%filenames_default(nf) = 'bkg'
nf = nf+1; self%filenames_default(nf) = 'crtmsrf'
nf = nf+1; self%filenames_default(nf) = 'fvcore_internal_rst'
nf = nf+1; self%filenames_default(nf) = 'moist_internal_rst'
nf = nf+1; self%filenames_default(nf) = 'surf_import_rst'

if (nf .ne. numfiles) &
  call abor1_ftn("fv3jedi_io_geos_mod.setup: number of potential restart files &
                  does not match numfiles")


! Set default tile dim and restart flags
! --------------------------------------
self%tiledim(1) = .true.  ! History defult tile dim
self%tiledim(2) = .true.  ! History defult tile dim
self%restart(3) = .false. ! Restarts do not use tiledim
self%restart(4) = .false. ! Restarts do not use tiledim
self%restart(5) = .false. ! Restarts do not use tiledim

self%restart(1) = .false. !Is not a restart file
self%restart(2) = .false. !Is not a restart file
self%restart(3) = .true.  !Is a restart file
self%restart(4) = .true.  !Is a restart file
self%restart(5) = .true.  !Is a restart file


! Get configuration
! -----------------
if (present(f_conf)) call get_conf(self,f_conf)


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


! Formatting of the files
! -----------------------

! Check if filename needs datetime information
do nf = 1, numfiles
  dt_in_name = index(self%filenames(nf),"%yyyy") + index(self%filenames(nf),"%mm") + index(self%filenames(nf),"%dd") + &
               index(self%filenames(nf),"%hh"  ) + index(self%filenames(nf),"%MM") + index(self%filenames(nf),"%ss")

  ! If needed in filename get datetime
  if (dt_in_name > 0) then
    call vdate_to_datestring(vdate, yyyy=yyyy, mm=mm, dd=dd, hh=hh, min=min, ss=ss)
    exit
  endif

end do

do nf = 1, numfiles

  ! Replace %yyyy, %mm, %dd, %hh, %MM, %ss with current date if user requested
  if (index(self%filenames(nf),"%yyyy") > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%yyyy',yyyy)
  if (index(self%filenames(nf),"%mm"  ) > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%mm'  ,mm  )
  if (index(self%filenames(nf),"%dd"  ) > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%dd'  ,dd  )
  if (index(self%filenames(nf),"%hh"  ) > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%hh'  ,hh  )
  if (index(self%filenames(nf),"%MM"  ) > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%MM'  ,min )
  if (index(self%filenames(nf),"%ss"  ) > 0) self%filenames(nf) = replace_text(self%filenames(nf),'%ss'  ,ss  )

  ! Prepend filenames with path
  self%filenames(nf) = trim(self%datapath)//'/'//trim(self%filenames(nf))

enddo


! Set files for the fields
! ------------------------
do var = 1,size(fields)

  select case (trim(fields(var)%short_name))

    ! Standard background history file
    case("ud","vd","u","v","ua","va","t","q","delp","qi","ql","qs","qr","o3mr","qls","qcn","cfcn","ps","phis",&
         "hs_stdv","frland","frlandice","frlake","frocean","frseaice","kcbl","tsm","khl","khu",&
         "varflt","ustar","bstar","zpbl","cm","ct","cq","u10m","v10m","ts","sheleg","soilt","soilm",&
         "DU001","DU002","DU003","DU004","DU005","SS001","SS002","SS003","SS004","SS005",&
         "BCPHOBIC","BCPHILIC","OCPHOBIC","OCPHILIC","NO3AN1","NO3AN2","NO3AN3","SO4",&
         "T","DELP","sphum","ice_wat","liq_wat")
      call set_file_names(self,var,1)

    ! CRTM surface quantities, usually from GFS output
    case("vtype","stype","vfrac")
      call set_file_names(self,var,2)

    ! GEOS fv_core restart
    case("U","V","W","PT","PKZ","PE","DZ")
      call set_file_names(self,var,3)

    ! GEOS moist restart
    case("Q","QILS","QICN","QLLS","QLCN","CLLS","CLCN")
      call set_file_names(self,var,4)

    ! GEOS surface restart
    case("PHIS")
      call set_file_names(self,var,5)

    ! Default to abort
    case default
      call abor1_ftn("fv3jedi_io_geos_mod.create: geos restart file for "//trim(fields(var)%short_name)//" not defined")

  endselect

enddo


! Open/create the files
! ---------------------
self%ncid = -1

if (self%iam_io_proc) then

  if (trim(readorwrite) == 'read') then

    fileopts = NF90_NOWRITE

    ! Open files for reading
    do nf = 1,numfiles
      if (self%ncid_isneeded(nf)) then
        call nccheck ( nf90_open( trim(self%filenames(nf)), fileopts, self%ncid(nf) ), &
                       "nf90_open"//trim(self%filenames(nf)) )
      endif
    enddo

  elseif (trim(readorwrite) == 'write') then

    fileopts = ior(NF90_NETCDF4, NF90_MPIIO)
    if (self%clobber) fileopts = ior(fileopts, NF90_CLOBBER)

    ! Create/open files for writing
    do nf = 1,numfiles
      if (self%ncid_isneeded(nf)) then
        call nccheck( nf90_create( trim(self%filenames(nf)), fileopts, self%ncid(nf), &
                                   comm = self%ocomm, info = MPI_INFO_NULL), &
                                   "nf90_create"//trim(self%filenames(nf)) )
      endif
    enddo

  else

    call abor1_ftn("fv3jedi_io_geos_mod.setup: readorwrite must be read or write not "//readorwrite)

  endif

endif

end subroutine setup

! ------------------------------------------------------------------------------

subroutine get_conf(self,f_conf)

implicit none
class(fv3jedi_io_geos),    intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_conf

integer :: n, nf

! Path where files are read from or saved to
! ------------------------------------------
call string_from_conf(f_conf,"datapath",self%datapath,'Data',memberswap=.true.)

! User can ask for extra meta data needed for GEOS ingest
! -------------------------------------------------------
if (.not. f_conf%get('geosingestmeta',self%geosingestmeta)) self%geosingestmeta = .false.

! Whether to expect/use the tile dimenstion in the file
! -----------------------------------------------------
if (.not. f_conf%get('clobber',self%clobber)) self%clobber = .false.

! Whether to expect/use the tile dimenstion in the file
! -----------------------------------------------------
if (.not. f_conf%get('tiledim',self%tiledim(1))) self%tiledim(1) = .true.
if (.not. f_conf%get('tiledim',self%tiledim(2))) self%tiledim(2) = .true.

! User can optionally specify the file names
! -------------------------------------------
nf = 0
nf = nf+1; call string_from_conf(f_conf,"filename_bkgd",self%filenames(1),self%filenames_default(1),memberswap=.true.)
nf = nf+1; call string_from_conf(f_conf,"filename_crtm",self%filenames(2),self%filenames_default(2),memberswap=.true.)
nf = nf+1; call string_from_conf(f_conf,"filename_core",self%filenames(3),self%filenames_default(3),memberswap=.true.)
nf = nf+1; call string_from_conf(f_conf,"filename_mois",self%filenames(4),self%filenames_default(4),memberswap=.true.)
nf = nf+1; call string_from_conf(f_conf,"filename_surf",self%filenames(5),self%filenames_default(5),memberswap=.true.)

! Sanity check
! ------------
if (nf .ne. numfiles) &
  call abor1_ftn("fv3jedi_io_geos_mod.get_conf: number of potential restart files &
                  does not match numfiles")

end subroutine get_conf

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3jedi_io_geos), intent(inout) :: self

integer :: ierr, nf

! Close the files
! ---------------
if (self%iam_io_proc) then
  do nf = 1,numfiles
    if (self%ncid_isneeded(nf)) then
      call nccheck ( nf90_close(self%ncid(nf)), "nf90_close" )
    endif
  enddo
endif

! Deallocate
! ----------
if (allocated(self%ncid_forfield)) deallocate(self%ncid_forfield)

! Release split comms
! -------------------
if (self%csize > 6) call MPI_Comm_free(self%tcomm, ierr)
call MPI_Comm_free(self%ocomm, ierr)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine read_meta(self, geom, vdate, calendar_type, date_init)

implicit none
class(fv3jedi_io_geos), intent(inout) :: self
type(fv3jedi_geom),     intent(inout) :: geom          !< Geometry
type(datetime),         intent(inout) :: vdate         !< DateTime
integer,                intent(inout) :: calendar_type !< Calendar type
integer,                intent(inout) :: date_init(6)  !< Date intialized

integer :: varid, date(6), intdate, inttime, idate, isecs
character(len=8) :: cdate
character(len=6) :: ctime
integer(kind=c_int) :: cidate, cisecs, nf, df

calendar_type = -1
date_init = 0

idate = 0
isecs = 0

if (self%iam_io_proc) then

  ! Get time attributes
  do nf = 1,numfiles
    if (self%ncid_isneeded(nf)) then
      df = nf
      exit
    endif
  enddo

  call nccheck ( nf90_inq_varid(self%ncid(df), "time", varid), "nf90_inq_varid time" )
  call nccheck ( nf90_get_att(self%ncid(df), varid, "begin_date", intdate), "nf90_get_att begin_date" )
  call nccheck ( nf90_get_att(self%ncid(df), varid, "begin_time", inttime), "nf90_get_att begin_time" )

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

end subroutine read_meta

! ------------------------------------------------------------------------------

subroutine read_fields(self, geom, fields)

implicit none
class(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_geom),             intent(inout) :: geom
type(fv3jedi_field),            intent(inout) :: fields(:)

integer :: varid, var, lev, ncid
logical :: tiledim
integer, pointer :: istart(:), icount(:)
integer, allocatable, target :: is_r3_tile(:), is_r3_noti(:)
real(kind=kind_real), allocatable :: arrayg(:,:), delp(:,:,:)


! Local copy of starts for rank 3 in order to do one level at a time
! ------------------------------------------------------------------
allocate(is_r3_tile(size(self%is_r3_tile)))
is_r3_tile = self%is_r3_tile
allocate(is_r3_noti(size(self%is_r3_noti)))
is_r3_noti = self%is_r3_noti

! Array for level of whole tile
! -----------------------------
allocate(arrayg(1:geom%npx-1,1:geom%npy-1))


! Loop over fields
! ----------------
do var = 1,size(fields)

  ! ncid for this variable
  ncid = self%ncid(self%ncid_forfield(var))

  !If ps then switch to delp
  if (trim(fields(var)%fv3jedi_name) == 'ps') then
    fields(var)%short_name = 'delp'
    fields(var)%npz = geom%npz
    allocate(delp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
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
    if (trim(fields(var)%fv3jedi_name) .ne. 'ps') then
      if (self%csize > 6) then
        call scatter_tile(geom, self%tcomm, 1, arrayg, fields(var)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev))
      else
        fields(var)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev) = arrayg(geom%isc:geom%iec,geom%jsc:geom%jec)
      endif
    else
      if (self%csize > 6) then
        call scatter_tile(geom, self%tcomm, 1, arrayg, delp(geom%isc:geom%iec,geom%jsc:geom%jec,lev))
      else
        delp(geom%isc:geom%iec,geom%jsc:geom%jec,lev) = arrayg(geom%isc:geom%iec,geom%jsc:geom%jec)
      endif
    endif

  enddo

  if (trim(fields(var)%fv3jedi_name) == 'ps') then
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

end subroutine read_fields

! ------------------------------------------------------------------------------

subroutine write_all(self, geom, fields, vdate)

implicit none
class(fv3jedi_io_geos), intent(inout) :: self
type(fv3jedi_geom),     intent(inout) :: geom          !< Geom
type(fv3jedi_field),    intent(in)    :: fields(:)     !< Fields to be written
type(datetime),         intent(in)    :: vdate         !< DateTime

! Write meta data
! ---------------
if (.not. self%clobber) call write_meta(self, geom, fields, vdate)

! Write fields
! ------------
call write_fields(self, geom, fields, vdate)

end subroutine write_all

! ------------------------------------------------------------------------------

subroutine write_meta(self, geom, fields, vdate)

implicit none
class(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_geom),             intent(inout) :: geom          !< Geom
type(fv3jedi_field),            intent(in)    :: fields(:)     !< Fields to be written
type(datetime),                 intent(in)    :: vdate         !< DateTime

integer :: var, ymult, k, nf, vc
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
allocate(latg(1:geom%npx-1,1:geom%npy-1))
allocate(long(1:geom%npx-1,1:geom%npy-1))

if (self%csize > 6) then
  call gather_tile(geom, self%tcomm, 1, rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec), latg)
  call gather_tile(geom, self%tcomm, 1, rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec), long)
else
  latg = rad2deg*geom%grid_lat(geom%isc:geom%iec,geom%jsc:geom%jec)
  long = rad2deg*geom%grid_lon(geom%isc:geom%iec,geom%jsc:geom%jec)
endif

write(cubesize,'(I8)') geom%npx-1

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
  allocate(xdimydim(geom%npx-1))
  do k = 1,geom%npx-1
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
  allocate(levels(geom%npz+1))
  do k = 1,geom%npz+1
    levels(k) = k
  enddo


  ! Loop over all files to be created/written to
  ! --------------------------------------------
  do nf = 1, numfiles

    if (self%ncid_isneeded(nf)) then

      ! Multiplication factor when no tile dimension
      ! --------------------------------------------
      ymult = 1
      if (.not. self%tiledim(nf)) ymult = 6

      if ( self%tiledim(nf) ) then
        XdimVar = 'Xdim'
        YdimVar = 'Ydim'
      else
        XdimVar = 'lon'
        YdimVar = 'lat'
      endif

      ! Main dimensions
      call nccheck ( nf90_def_dim(self%ncid(nf), trim(XdimVar), geom%npx-1,         self%x_dimid), "nf90_def_dim "//trim(XdimVar) )
      call nccheck ( nf90_def_dim(self%ncid(nf), trim(YdimVar), ymult*(geom%npy-1), self%y_dimid), "nf90_def_dim "//trim(YdimVar) )
      if (self%tiledim(nf)) &
      call nccheck ( nf90_def_dim(self%ncid(nf), "nf",   geom%ntiles, self%n_dimid),  "nf90_def_dim nf"   )
      call nccheck ( nf90_def_dim(self%ncid(nf), "lev",  geom%npz,    self%z_dimid),  "nf90_def_dim lev"  )
      call nccheck ( nf90_def_dim(self%ncid(nf), "edge", geom%npz+1,  self%e_dimid), "nf90_def_dim edge" )
      call nccheck ( nf90_def_dim(self%ncid(nf), "time", 1,           self%t_dimid),  "nf90_def_dim time" )

      ! In case the four level surface fields need to be written
      do var = 1,size(fields)
        if (fields(var)%npz == 4) then
          call nccheck ( nf90_def_dim(self%ncid(nf), "lev4", 4, self%f_dimid), "nf90_def_dim lev"  )
          exit
        endif
      enddo

      !Needed by GEOS for ingesting cube sphere field
      call nccheck ( nf90_def_dim(self%ncid(nf), "ncontact",          4, self%c_dimid), "nf90_def_dim ncontact" )
      call nccheck ( nf90_def_dim(self%ncid(nf), "orientationStrLen", 5, self%o_dimid), "nf90_def_dim orientationStrLend" )

      ! Dimension ID array for lat/lon arrays
      if (allocated(dimidsg)) deallocate(dimidsg)
      if ( self%tiledim(nf) ) then
        allocate(dimidsg(3))
        dimidsg(:) = (/ self%x_dimid, self%y_dimid, self%n_dimid /)
      else
        allocate(dimidsg(2))
        dimidsg(:) = (/ self%x_dimid, self%y_dimid /)
      endif

      ! Define fields to be written (geom)
      vc=0;

      if (self%tiledim(nf)) then
        vc=vc+1;
        call nccheck( nf90_def_var(self%ncid(nf), "nf", NF90_INT, self%n_dimid, varid(vc)), "nf90_def_var nf" )
        call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "cubed-sphere face") )
        call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "axis", "e") )
        call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "grads_dim", "e") )
      endif

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), trim(XdimVar), NF90_DOUBLE, self%x_dimid, varid(vc)), "nf90_def_var "//trim(XdimVar) )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "Fake Longitude for GrADS Compatibility") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "degrees_east") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), trim(YdimVar), NF90_DOUBLE, self%y_dimid, varid(vc)), "nf90_def_var "//trim(YdimVar) )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "Fake Latitude for GrADS Compatibility") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "degrees_north") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), "lons", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lons" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "longitude") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "degrees_east") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), "lats", NF90_DOUBLE, dimidsg, varid(vc)), "nf90_def_var lats" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "latitude") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "degrees_north") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), "lev", NF90_DOUBLE, self%z_dimid, varid(vc)), "nf90_def_var lev" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "vertical level") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "layer") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "positive", "down") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "coordinate", "eta") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "standard_name", "model_layers") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), "edge", NF90_DOUBLE, self%e_dimid, varid(vc)), "nf90_def_var edge" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "vertical level edges") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "units", "layer") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "positive", "down") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "coordinate", "eta") )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "standard_name", "model_layers") )

      vc=vc+1;
      call nccheck( nf90_def_var(self%ncid(nf), "time", NF90_INT, self%t_dimid, varid(vc)), "nf90_def_var time" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "long_name", "time"), "nf90_def_var time long_name" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "begin_date", date8), "nf90_def_var time begin_date" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "begin_time", time6), "nf90_def_var time begin_time" )

      vc=vc+1; !(Needed by GEOS to ingest cube sphere analysis)
      call nccheck( nf90_def_var(self%ncid(nf), "cubed_sphere", NF90_CHAR, varid(vc)), "nf90_def_var cubed_sphere" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "grid_mapping_name", "gnomonic cubed-sphere"), "nf90_def_var time grid_mapping_name" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "file_format_version", "2.90"), "nf90_def_var time file_format_version" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "additional_vars", "contacts,orientation,anchor"), "nf90_def_var time additional_vars" )
      call nccheck( nf90_put_att(self%ncid(nf), varid(vc), "gridspec_file", "C"//trim(cubesize)//"_gridspec.nc4"), "nf90_def_var gridspec_file" )

      !vc=vc+1; !(Needed by GEOS to ingest cube sphere analysis)
      !call nccheck( nf90_def_var(self%ncid(nf), "ncontact", NF90_INT, varid(vc)), "nf90_def_var ncontact" )


      ! End define mode
      ! ---------------
      call nccheck( nf90_enddef(self%ncid(nf)), "nf90_enddef" )


      ! Reset counter
      ! -------------
      vc=0


      ! Write metadata arrays
      ! ---------------------

      ! Tiles
      if (self%tiledim(nf)) &
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), tiles ), "nf90_put_var nf" )

      ! Xdim & Ydim arrays
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), xdimydim ), "nf90_put_var "//trim(XdimVar) )
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), xdimydim ), "nf90_put_var "//trim(YdimVar) )

      ! Start/counts
      if (associated(istart)) nullify(istart)
      if (associated(icount)) nullify(icount)
      if (self%tiledim(nf)) then
        istart => self%is_r2_tile(1:3); icount => self%ic_r2_tile(1:3)
      else
        istart => self%is_r2_noti(1:2); icount => self%ic_r2_noti(1:2)
      endif

      ! Lat/lon arrays
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), long, &
                                                  start = istart, &
                                                  count = icount ), &
                                                  "nf90_put_var lons" )

      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), latg, &
                                                  start = istart, &
                                                  count = icount ), &
                                                  "nf90_put_var lats" )

      ! Write model levels & time
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), levels(1:geom%npz) ), "nf90_put_var lev" )
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), levels ), "nf90_put_var edge" )

      ! Time
      vc=vc+1;call nccheck( nf90_put_var( self%ncid(nf), varid(vc), 0 ), "nf90_put_var time" )

    endif

  enddo

endif

end subroutine write_meta

! ------------------------------------------------------------------------------

subroutine write_fields(self, geom, fields, vdate)

implicit none
class(fv3jedi_io_geos), target, intent(inout) :: self
type(fv3jedi_geom),             intent(inout) :: geom          !< Geom
type(fv3jedi_field),            intent(in)    :: fields(:)     !< Fields to be written
type(datetime),                 intent(in)    :: vdate         !< DateTime

integer :: var, lev, nf, ncid, filei, varid
integer, target :: dimids2_tile(4), dimids3_tile(5), dimidse_tile(5), dimids4_tile(5)
integer, target :: dimids2_noti(3), dimids3_noti(4), dimidse_noti(4), dimids4_noti(4)
integer, pointer :: dimids2(:), dimids3(:), dimidse(:), dimids4(:), dimids(:)
integer, allocatable, target :: is_r3_tile(:), is_r3_noti(:)
real(kind=kind_real), allocatable :: arrayg(:,:)
integer, pointer :: istart(:), icount(:)


! Whole level of tile array
! -------------------------
allocate(arrayg(1:geom%npx-1,1:geom%npy-1))


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
    elseif (fields(var)%npz == geom%npz) then
      dimids => dimids3
    elseif (fields(var)%npz == geom%npz+1) then
      dimids => dimidse
    elseif (fields(var)%npz == 4) then
      dimids => dimids4
    else
      call abor1_ftn("write_geos: vertical dimension not supported")
    endif

    ! Define field
    call nccheck( nf90_def_var(ncid, trim(fields(var)%short_name), NF90_DOUBLE, dimids, varid), &
                  "nf90_def_var"//trim(fields(var)%short_name))

    ! Write attributes if not clobbering
    if (.not.self%clobber) then

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
      call gather_tile(geom, self%tcomm, 1, fields(var)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev), arrayg)
    else
      arrayg = fields(var)%array(geom%isc:geom%iec,geom%jsc:geom%jec,lev)
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

subroutine gather_tile(geom, comm, nlev, array_l, array_g)

implicit none

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec,1:nlev)  ! Local array
real(kind=kind_real), intent(inout) :: array_g(1:geom%npx-1,1:geom%npy-1,1:nlev)            ! Gathered array (only valid on root)

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
npx_g = geom%npx-1
npy_g = geom%npy-1
npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(geom%isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
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
  do jj = geom%jsc,geom%jec
    do ji = geom%isc,geom%iec
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

subroutine scatter_tile(geom, comm, nlev, array_g, array_l)

implicit none

type(fv3jedi_geom),   intent(in)    :: geom
integer,              intent(in)    :: comm
integer,              intent(in)    :: nlev
real(kind=kind_real), intent(in)    :: array_g(1:geom%npx-1,1:geom%npy-1,nlev)            ! Gathered array (only valid on root)
real(kind=kind_real), intent(inout) :: array_l(geom%isc:geom%iec,geom%jsc:geom%jec,nlev)  ! Local array

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
npx_g = geom%npx-1
npy_g = geom%npy-1
npx_l = geom%iec-geom%isc+1
npy_l = geom%jec-geom%jsc+1

!Gather local dimensions
allocate(isc_l(comm_size), iec_l(comm_size), jsc_l(comm_size), jec_l(comm_size))
call mpi_allgatherv(geom%isc, 1, mpi_int, isc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%iec, 1, mpi_int, iec_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jsc, 1, mpi_int, jsc_l, counts, displs, mpi_int, comm, ierr)
call mpi_allgatherv(geom%jec, 1, mpi_int, jec_l, counts, displs, mpi_int, comm, ierr)
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
  do jj = geom%jsc,geom%jec
    do ji = geom%isc,geom%iec
      n = n+1
      array_l(ji,jj,jk) = vector_l(n)
    enddo
  enddo
enddo
deallocate(vector_l)

end subroutine scatter_tile

! --------------------------------------------------------------------------------------------------

subroutine dummy_final(self)
type(fv3jedi_io_geos), intent(inout) :: self
end subroutine dummy_final

! --------------------------------------------------------------------------------------------------

subroutine string_from_conf(f_conf,varstring,var,default,memberswap)

implicit none
type(fckit_configuration),  intent(in)  :: f_conf
character(len=*),           intent(in)  :: varstring
character(len=*),           intent(out) :: var
character(len=*), optional, intent(in)  :: default
logical,          optional, intent(in)  :: memberswap

character(len=:), allocatable :: str

if (.not. f_conf%get(trim(varstring),str)) then

  if (present(default)) then
    var = trim(default)
  else
    call abor1_ftn("fv3jedi_io_geos_mod.string_from_conf: "//trim(varstring)//" not found in config&
                    and no default provided. Aborting")
  endif

else

  if (present(memberswap) .and. memberswap) call swap_name_member(f_conf, str)

  var = trim(str)

endif

if (allocated(str)) deallocate(str)

end subroutine string_from_conf

! --------------------------------------------------------------------------------------------------

subroutine str_check(str,maxlen)

implicit none
character(len=*), intent(in) :: str
integer,          intent(in) :: maxlen

character(len=maxlen) :: maxlenstr

write (maxlenstr, *) maxlen

if (len(str) > maxstring) then
  call abor1_ftn('Reading '//trim(str)//'from configuration. Too long, max length = '//trim(maxlenstr))
endif

end subroutine str_check

! --------------------------------------------------------------------------------------------------

subroutine vdate_to_datestring(vdate,datest,date,yyyy,mm,dd,hh,min,ss)

implicit none
type(datetime),              intent(in)  :: vdate
character(len=*), optional,  intent(out) :: datest
integer,          optional,  intent(out) :: date(6)
character(len=4), optional,  intent(out) :: yyyy
character(len=2), optional,  intent(out) :: mm
character(len=2), optional,  intent(out) :: dd
character(len=2), optional,  intent(out) :: hh
character(len=2), optional,  intent(out) :: min
character(len=2), optional,  intent(out) :: ss

integer :: dateloc(6)
integer(kind=c_int) :: idate, isecs

! Outputs various forms of datetime

call datetime_to_ifs(vdate, idate, isecs)
dateloc(1) = idate/10000
dateloc(2) = idate/100 - dateloc(1)*100
dateloc(3) = idate - (dateloc(1)*10000 + dateloc(2)*100)
dateloc(4) = isecs/3600
dateloc(5) = (isecs - dateloc(4)*3600)/60
dateloc(6) = isecs - (dateloc(4)*3600 + dateloc(5)*60)

if (present(datest)) &
write(datest,'(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)') dateloc(1),dateloc(2),dateloc(3),"_",&
                                                 dateloc(4),dateloc(5),dateloc(6)

!Optionally pass date back
if (present(date)) date = dateloc

! Optionally pass back individual strings of datetime
if (present(yyyy)) write(yyyy,'(I4)  ') dateloc(1)
if (present(mm  )) write(mm  ,'(I0.2)') dateloc(2)
if (present(dd  )) write(dd  ,'(I0.2)') dateloc(3)
if (present(hh  )) write(hh  ,'(I0.2)') dateloc(4)
if (present(min )) write(min ,'(I0.2)') dateloc(5)
if (present(ss  )) write(ss  ,'(I0.2)') dateloc(6)

end subroutine vdate_to_datestring

! --------------------------------------------------------------------------------------------------

subroutine set_file_names(self,var,grp)

implicit none
type(fv3jedi_io_geos), intent(inout) :: self
integer,               intent(in)    :: var
integer,               intent(in)    :: grp

self%ncid_forfield(var) = grp
self%ncid_isneeded(grp) = .true.

end subroutine set_file_names

! --------------------------------------------------------------------------------------------------

function replace_text (inputstr,search,replace)  result(outputstr)

implicit none
character(len=*), intent(in) :: inputstr
character(len=*), intent(in) :: search
character(len=*), intent(in) :: replace
character(len(inputstr)+100) :: outputstr

! Locals
integer :: i, nt, nr

outputstr = inputstr
nt = len_trim(search)
nr = len_trim(replace)

do
  i = index(outputstr,search(:nt)) ; if (i == 0) exit
  outputstr = outputstr(:i-1) // replace(:nr) // outputstr(i+nt:)
end do

end function replace_text

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_geos_mod
                                                                                          
