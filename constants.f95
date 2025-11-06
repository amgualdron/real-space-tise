module constants 
    implicit none 
    private 

    !working precision
    integer, parameter, public :: dp = SELECTED_REAL_KIND(16)

    real(kind = dp), parameter, public :: pi = 4.0_dp * atan(1.0_dp)

    !Speed of light(m/s)
    real(kind = dp), parameter, public :: c = 299792458.0_dp 

    !plancks constant and reduced plancks constant(eV * s)
    real(kind = dp), parameter, public :: h = 4.135667696E-15 
    real(kind = dp), parameter, public :: hbar = 6.582119569E-16 

    !electron rest mass energy equivalence (in eV/c^2)
    real(kind = dp), parameter :: m_e = 0.51099895069E6

end module constants 