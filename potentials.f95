module potentials
    use constants 
    implicit none 

    contains 
    !----------------------------------------------
    ! 1 dimensional potentials 
    !----------------------------------------------
    real pure elemental function infinite_square_well() result(v)
    end function infinite_square_well

    real pure elemental function finite_square_well() result(v)
    end function finite_square_well

    real pure elemental function qho() result(v)
    end function qho

end module