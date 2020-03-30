program sopfbr
  use tensor_nmode
  use chebval
  implicit none
  integer :: i, j, k, jkappa, tkappa, mode, dims
  integer :: ndim = 3
  integer, dimension(3) :: gdim = (/5, 5, 5/)
  integer :: tdim = 4
  real, dimension(3) :: point
  real, dimension(15) :: u_vects
  real, dimension(60) :: coef_u_vects
  real, dimension(125) :: flattensor
  real, dimension(:), allocatable :: newflattensor_1
  real, dimension(:), allocatable :: newflattensor_2
  real :: serieval

  ! A test geometry

  point = (/-0.33, 0.15, 0.88/)

  ! A test with unit tensor and Chebyshev series

  coef_u_vects = 1.0
  flattensor = 1.0
  u_vects = 0.0

  ! Initialize factor vectors

  jkappa = 0
  tkappa = 0
  do i = 1, ndim
     do j = 1, gdim(i)
        jkappa = jkappa + 1
        serieval = 0.0
        do k = 1, tdim
           tkappa = tkappa + 1
           serieval = serieval + coef_u_vects(tkappa) * chebpoly(point(i), k - 1)
        enddo
        u_vects(jkappa) = u_vects(jkappa) + serieval
     enddo
  enddo

  dims = ndim
  allocate(newflattensor_1(125))
  do mode = 1, ndim
     allocate(newflattensor_2(100))
     call first_mode(dims, gdim(mode:), flattensor, u_vects, newflattensor)
     dims = dims - 1
     deallocate(newflattensor_1)
     deallocate(newflattensor_2)
  enddo

end program sopfbr
