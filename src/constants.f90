module constants 
    implicit none 
    private 

    !working precision
    integer, parameter, public :: dp = SELECTED_REAL_KIND(15)

    !pi
    real(kind = dp), parameter, public :: pi = 4.0_dp * atan(1.0_dp)

    !Working units: 
    ! distance: Angstroms [A]
    ! time: seconds [s]
    ! energy: electron Volts [eV]
    !----------------------------

    !Speed of light(A/s)
    real(kind = dp), parameter, public :: c = 2.997924580E18_dp 

    !plancks constant and reduced plancks constant(eV * s)
    real(kind = dp), parameter, public :: h = 4.135667696E-15_dp 
    real(kind = dp), parameter, public :: hbar = 6.582119569E-16_dp 
    
    !useful constant h*c (eV * A)
    real(kind = dp), parameter, public :: hc = h*c

    !hbarc (eV * A)
    real(kind = dp), parameter, public :: hbarc = hbar * c

    !electron rest mass energy (eV)
    real(kind = dp), parameter, public :: m_e_c2 = 0.51099895069E6_dp

    !electron mass(kg)
    real(kind = dp), parameter, public :: m_e = m_e_c2 / (c**2)
    
    !infinity (safe so that wont overflow variables)
    !can adjust if necessary
    real(kind = dp), parameter, public :: inf = 1E20_dp

end module constants 