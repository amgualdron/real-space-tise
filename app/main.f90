program main 

    !import modules
    use constants
    use potentials

    implicit none 

    !import lapack eigenvalue solver
    external :: dstev

    !Parameters
    !-----------------

    !file numbers
    integer, parameter :: file_energies = 10
    integer, parameter :: file_wavefuncs = 11

    !number of N discrete points, increase for more precision
    integer, parameter :: N = 5000
    real(kind = dp), parameter :: L = 100.0_dp !in [A]
    
    !Variables
    !------------------
    real(kind = dp) :: dx, T
    real(kind=dp), allocatable :: x(:) !discretized space x_i
    integer :: i

    !Lapack arrays (dstev)
    real(kind = dp), allocatable :: D(:)    !main diagonal -> becomes eigenvalues after lapack
    real(kind = dp), allocatable :: K(:)    !super diagonal(and sub)
    real(kind = dp), allocatable :: psi(:,:)!eigenvectors

    ! required by dstev
    real(kind=dp), allocatable :: work(:)

    ! INFO: Error flag
    integer :: INFO

    dx = L/ real(N-1,dp)
    T = (hbar**2)/(2.0_dp * m_e * (dx**2))

    !Allocation
    !-------------------
    allocate(x(N), D(N), K(N-1), psi(N,N), work(2*N - 2))

    print *, "starting simulation, building hamiltonian"
    !fill arrays
    !--------------------
    x(1) = real(-L/2.0_dp,dp)!set x_0 = -L/2
    do i = 2,N
        x(i) = x(1) + real(i-1,dp) * dx !fill x with N discrete points up to L/2
    end do

    D = qho_potential(x) + 2.0_dp * T
    K = -T

    print *, "calling lapack"
    call dstev('V', N, D, K, psi, N, work, info)

    if(info /= 0) then 
        error stop "Lapack error"
    else 
        print*, "calculations done."
    end if

    print *, "writing data to results/"

    !write eigenvalues (energies) to file
    !---------------------------------------
    open(unit = file_energies, file = 'results/qho_energies.txt', status = 'replace')
    write(file_energies, '(a)') "# qho eigenvalues for N = 5000 L = 100A hbarw = 10meV"
    write(file_energies, '(a)') "# n, Energy(eV)"
    do i = 1,5
        write(file_energies, '(i3, 2x, f12.6)')i, D(i)
    end do

    close(file_energies)
    !-----------------------------------------

    !Write eigenvectors to file
    !------------------------------------------
    open(unit = file_wavefuncs, file = 'results/qho_wavefunctions.txt', status = 'replace')
    write(file_wavefuncs, '(a)') "# qho eigenfunctions for N = 5000 L = 100A hbarw = 10meV"
    write(file_wavefuncs, '(a)') "# x, psi_1, psi_2, psi_3, psi_4, psi_5"
    do i = 1,N
        write(file_wavefuncs, *) x(i), psi(i, 1), psi(i, 2), psi(i, 3), psi(i, 4), psi(i, 5)   
    end do

    close(file_wavefuncs)
    !-----------------------------------------

    deallocate(x, D, K, psi, work)

end program main