module fv3jedi_fieldfail_mod

implicit none
private
public field_fail

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine field_fail(field)

character(len=*), intent(in) :: field

call abor1_ftn("Field_fail: Field "//trim(field)//" cannot be obtained from input fields.")

end subroutine field_fail

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fieldfail_mod
