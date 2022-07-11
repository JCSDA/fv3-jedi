! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_io_utils_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module,   only: fckit_configuration

! oops
use datetime_mod
use string_utils, only: swap_name_member

implicit none
public

integer, parameter :: maxstring = 2048

! --------------------------------------------------------------------------------------------------

contains

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

subroutine vdate_to_datestring(vdate,datest,isodate,ufsdate,date,yyyy,mm,dd,hh,min,ss)

implicit none
type(datetime),              intent(in)  :: vdate
character(len=*), optional,  intent(out) :: datest
character(len=*), optional,  intent(out) :: isodate
character(len=*), optional,  intent(out) :: ufsdate
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
write(datest,'(I4,I0.2,I0.2,"_",I0.2,I0.2,I0.2)') dateloc(1),dateloc(2),dateloc(3),&
                                                 dateloc(4),dateloc(5),dateloc(6)
if (present(isodate)) &
write(isodate,'(I4,"-",I0.2,"-",I0.2,"T",I0.2,":",I0.2,":",I0.2,"Z")') dateloc(1),dateloc(2),dateloc(3),&
                                                 dateloc(4),dateloc(5),dateloc(6)
if (present(ufsdate)) &
write(ufsdate,'(I4,"-",I0.2,"-",I0.2,"T",I0.2,":",I0.2,":",I0.2)') dateloc(1),dateloc(2),dateloc(3),&
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

function replace_text (inputstr,search,replace) result(outputstr)

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
    call abor1_ftn("fv3jedi_io_utils_mod.string_from_conf: "//trim(varstring)// &
                    " not found in config and no default provided. Aborting")
  endif

else

  if (present(memberswap) .and. memberswap) call swap_name_member(f_conf, str)

  var = trim(str)

endif

if (allocated(str)) deallocate(str)

end subroutine string_from_conf

! --------------------------------------------------------------------------------------------------

subroutine add_iteration(f_conf,str)

implicit none
type(fckit_configuration),  intent(in)  :: f_conf
character(len=:), allocatable, intent(inout) :: str

integer :: i, j
character(len=:), allocatable :: iter

if (f_conf%has("iteration")) then
  call f_conf%get_or_die("iteration", iter)
  i = index(str, 'nc')
  j = index(str, 'res')
  if (i == 0 .AND. j == 0) then
    str = str // iter // '.'
  else
    str = iter // '.' // str
  endif
endif

end subroutine add_iteration

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_utils_mod
