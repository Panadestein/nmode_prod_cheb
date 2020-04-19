program sopfbr
  use tensor_nmode
  use chebval
  implicit none
  integer :: i, j, k, jkappa, tkappa, mode, newsize, chebslice
  integer :: ndim = 3
  integer, dimension(3) :: gdim = (/5, 5, 5/)
  integer :: tdim = 4
  real, dimension(3) :: point
  real, dimension(15) :: u_vects
  real, dimension(60) :: coef_u_vects
  real, dimension(125) :: core ! Flattened version of the core tensor
  real, dimension(:), allocatable :: tensor_holder
  real, dimension(:), allocatable :: tensor_prod
  real :: serieval

  ! A test geometry

  point = (/-0.33, 0.15, 0.88/)

  ! A test core tensor

  core = 1.

  ! A test with unit tensor and Chebyshev series

  coef_u_vects = 1.0
  core = 1.0
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

  ! Compute tensor n-mode product

  allocate(tensor_holder(125))
  tensor_holder = core

  do mode = 1, ndim

     if (allocated(tensor_prod)) deallocate(tensor_prod)

     newsize = product(gdim(mode:)) / gdim(mode)
     chebslice = 5 * (mode - 1)
     allocate(tensor_prod(newsize))

     tensor_prod = n_mode(mode, gdim(mode:), newsize, &
                   tensor_holder, u_vects(chebslice:(chebslice + 4)))

     deallocate(tensor_holder)
     allocate(tensor_holder(newsize))
     tensor_holder = tensor_prod
     
  enddo

  print *, tensor_holder

end program sopfbr
