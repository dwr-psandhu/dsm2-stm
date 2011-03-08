!<license>
!    Copyright (C) 1996, 1997, 1998, 2001, 2007, 2009 State of California,
!    Department of Water Resources.
!    This file is part of DSM2.
!
!    The Delta Simulation Model 2 (DSM2) is free software: 
!    you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    DSM2 is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with DSM2.  If not, see <http://www.gnu.org/licenses>.
!</license>

!> Defines the input variables for the sediment transport sources.
!> This replaces the notion of "rough I/O" suggested by Jamie Anderson originally.
!> The module defines the basic parameters for the sediment transport sources based
!> on the table agreed on 2/11/2011.
!>@ingroup sediment
module sediment_variables
    use stm_precision
    integer :: ncell        !< Number of computational cells
    
    !> Example of variables to define, where the variable is constant.
    !> It is just a number    
    real(stm_real),save :: g_acceleration                 !< acceleration of gravity; it must be in SI units (constant)
!todo: we need to add the other variables as in the table  
  
    !> Example of variables to define, where the variable is a function of space only. 
    real(stm_real),save,allocatable :: manning_n(:)  !< Just function of space
!todo: we need to add the other variables as in the table
    
    
    
    contains
    
     
    !> Initial value is LARGEREAL
    subroutine set_sediment_constant
        use error_handling
        implicit none
        
        real(stm_real) :: g_acceleration  
       
        g_acceleration = LARGEREAL
       
        return
    end subroutine
    
    
    subroutine allocate_sediment_static(ncell)
        use error_handling
        implicit none
        integer,intent(in):: ncell    !<Number of cells
                   
       allocate(manning_n(ncell))
       
        return
    end subroutine
    
    
    !> Deallocate the sediment static variable
    !> including manning's n and width
    !> and reset ncell and nvar to zero.
    subroutine deallocate_sediment_static
        implicit none

        ncell = 0
      
        deallocate(manning_n)
      
        return
    end subroutine

end module



