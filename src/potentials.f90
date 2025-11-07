module potentials
    
    use constants 
    implicit none 

    contains 
    !----------------------------------------------
    ! 1 dimensional potentials 
    !----------------------------------------------
    real(kind = dp) pure elemental function infinite_square_well(x,x_min,x_max) result(v)
    !
    !   Infinite square well potential function
    !   Arguments: 
    !       x: independent variable at which function returns V(x)
    !       x_min, x_max: boundaries of well, xmin has to be lower than x max, the well
    !                     is non-inclusive meaning v(boudary) = infinity    
    !   Units: 
    !       x,x_min,x_max: Angstrom [A]
    !
    !   returns: 
    !       potential energy at x in [eV]
    !
    !   The value of the "infinite" boundary can be tuned in constants module
    !===================================================================================
        real(kind = dp), intent(in) :: x, x_min, x_max

        if( x > x_min .and. x < x_max) then
            v = 0.0_dp
        else 
            v = inf
        end if

    end function infinite_square_well

    real(kind = dp) pure elemental function finite_square_well(x,x_min,x_max,v_wall) result(v)
    !   
    !   Finite square well potential function 
    !   Arguments:
    !       x: independent variable at which function returns V(x) 
    !       x_min, x_max: boundaries of well, xmin has to be lower than x max, the well
    !                     is non-inclusive meaning v(boundary) = infinity
    !       v_wall: Magnitude of the potential outside of the well
    !   units: 
    !       x, x_min, x_max: Angstroms [A ]
    !       v_wall: electron volts     [eV]
    !   returns: 
    !       potential energy at x in [eV]
    !=========================================================================================
        real(kind = dp), intent(in) :: x, x_min, x_max, v_wall

        if( x > x_min .and. x < x_max) then
            v = 0.0_dp
        else 
            v = v_wall
        end if

    end function finite_square_well

    real(kind = dp) pure elemental function qho_potential(x) result(v)
    !   
    !   Parabolic well potential function (QHO)
    !   Arguments:
    !       x: independent variable at which function returns V(x) 
    !   units: 
    !       x: Angstroms               [A ]
    !       hbar_omega: electron volts [eV]
    !   Parameters: 
    !       hbar_omega: (Local parameter) The characteristic energy spacing, in [eV].
    !   returns: 
    !       potential energy at x in [eV]
    !=========================================================================================
    
        real(kind = dp), intent(in) :: x![A]

        real(kind = dp), parameter :: hbar_omega = 1.0_dp ![eV]

        v = (1.0_dp/2.0_dp) * m_e *((hbar_omega/hbar)**2) * (x**2)

    end function qho_potential

end module