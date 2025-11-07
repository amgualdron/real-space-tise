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
    real(kind = dp), parameter :: L = 100.0_dp !in [A]
    
    !Variables
    !------------------
    real(kind = dp) :: dx, T, x_curr 
    real(kind=dp), allocatable :: x(:) !discretized space x_i
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

    dx = L/ real(N-1,dp)
    T = (hbar**2)/(2.0_dp * m_e * (dx**2))

    !Allocation
    !-------------------
    allocate(x(N), D(N), K(N-1), E(N), psi(N,N), work(2*N - 2))

    !fill arrays
    !--------------------
    x(i) = real(-L/2.0_dp,dp)!set x_0 = -L/2
    do i = 0,N
        x(i) = x(1) + real(i-1,dp) * dx !fill x with N discrete points up to L/2
    end do

    D = qho_potential(x) + 2.0_dp * T
    K = -T


end program main