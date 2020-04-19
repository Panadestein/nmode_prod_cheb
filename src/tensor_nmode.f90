module tensor_nmode 
implicit none

contains

function n_mode(mode, gdim, newsize, core, spp) result(tenprod)
    integer, intent(in) :: mode
    integer, intent(in) :: newsize
    integer, dimension(:), intent(in) :: gdim
    real, dimension(:), intent(in) :: spp 
    real, dimension(product(gdim)), intent(in) :: core 
    real, dimension(newsize) :: tenprod

    integer, dimension(2) :: newshape

    if (gdim(mode) == 0) then
        tenprod = dot_product(core, spp)
    else
        newshape = (/newsize, gdim(mode)/)
        tenprod = matmul(reshape(core, newshape), spp)
    end if

end function n_mode

    
end module tensor_nmode 