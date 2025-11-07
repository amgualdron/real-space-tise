program main 

    !import modules
    use constants
    use potentials

    implicit none 

    !import lapack eigenvalue solver
    external :: dstev

    !Parameters
    !-----------------

    !number of N discrete points, increase for more precision
    integer, parameter :: N = 5000
    real(kind = dp), parameter :: L = 100 !in [A]
    
    !Variables
    !------------------
    real(kind = dp) :: dx, T, x_curr 
    real(kind=dp), allocatable :: x(:)
    integer :: i

    !Lapack arrays (dstev)
    real(kind = dp), allocatable :: D(:)    !main diagonal
    real(kind = dp), allocatable :: K(:)    !super diagonal(and sub)
    real(kind = dp), allocatable :: E(:)    !energy eigenvalues
    real(kind = dp), allocatable :: psi(:,:)!eigenvectors

    ! required by dstev
    real(kind=dp), allocatable :: work(:)

    ! INFO: Error flag
    integer :: INFO



end program main