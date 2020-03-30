module chebval
implicit none

contains

real recursive function chebpoly(x_val, degree) result(polval)
    real, intent(in) :: x_val
    integer, intent(in) :: degree

    if (degree == 0) then
        polval = 1
    else if (degree == 1) then
        polval = x_val
    else
        polval = 2 * x_val * chebpoly(x_val, degree-1) - chebpoly(x_val, degree-2)
    endif
end function chebpoly

real function chebserie(coeffs, tdim, x_val) result(serie)
    implicit none
    integer, intent(in) :: tdim
    real, intent(in) :: x_val
    real, dimension(tdim), intent(in) :: coeffs
    integer :: i

    do i = 1, tdim
       serie = serie + coeffs(i) * chebpoly(x_val, tdim)
    enddo 
    
end function chebserie
    
end module chebval