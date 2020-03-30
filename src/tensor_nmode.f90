module tensor_nmode 
implicit none

contains

subroutine first_mode(ndim, gdim, flattensor, tvector, newflattensor)
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(in) :: gdim
    real, dimension(:), intent(in) :: tvector
    real, dimension(product(gdim)), intent(in) :: flattensor
    real, dimension(:), intent(out) :: newflattensor
    integer, dimension(2) :: newshape
    integer :: mkappa, prod_kappa

    mkappa = gdim(1) ! Assumed product by first mode
    prod_kappa = product(gdim) / gdim(1)
    newshape = (/mkappa, prod_kappa/)

    newflattensor = matmul(reshape(flattensor, newshape), tvector)

end subroutine first_mode

    
end module tensor_nmode 