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
!> on the table agreed on 2/11/11.
!>@ingroup sediment
module sediment_variables
    use stm_precision
    
    !> Example of variables to define, where the variable is constant.
    !> It is just a number    
    real(stm_real),save :: g_acceleration                 !< acceleration of gravity; it must be in SI units (constant)
!todo: we need to add the other variables as in the table  
  
    !> Example of variables to define, where the variable is a function of space only. 
    real(stm_real),save,allocatable :: manning_n(:)  !< Just function of space
!todo: we need to add the other variables as in the table

    !> Example of variables to define, where the variable is a function of both space and time.    
    real(stm_real),save,allocatable :: flow_hi(:,:)  !< Function of time and space
!todo: we need to add the other variables as in the table    
    
    contains
    
    !> Allocate the state variables consistently
    !> including concentration and hydrodynamics.
    !> Initial value is LARGEREAL
    subroutine allocate_state(a_ncell,a_nvar)
        use error_handling
        implicit none
        character(LEN=128) :: message
        integer :: istat = 0
        integer, intent(in) :: a_ncell !< Number of requested cells
        integer, intent(in) :: a_nvar  !< Number of constituents
        ncell = a_ncell
        nvar  = a_nvar
        write(message,*)"Could not allocate state variable. " //&
         "This could be due to allocating several times in " // &
         "a row without deallocating (memory leak)"
        
        allocate(conc(ncell,nvar), conc_prev(ncell,nvar), stat = istat)
        if (istat .ne. 0 )then
           call stm_fatal(message)
        end if
        conc      = LARGEREAL  ! absurd value helps expose bugs  
        conc_prev = LARGEREAL
        allocate(mass(ncell,nvar), mass_prev(ncell,nvar),stat = istat)
        if (istat .ne. 0 )then
           call stm_fatal(message)
        end if
        mass      = LARGEREAL  ! absurd value helps expose bugs  
        mass_prev = LARGEREAL       
        allocate(area(ncell), area_prev(ncell), area_lo(ncell), area_hi(ncell), &
                 area_lo_prev(ncell), area_hi_prev(ncell),stat = istat)
        if (istat .ne. 0 )then
           call stm_fatal(message)
        end if
        area      = LARGEREAL
        area_prev = LARGEREAL
        area_lo   = LARGEREAL
        area_hi   = LARGEREAL
        area_lo_prev   = LARGEREAL
        area_hi_prev   = LARGEREAL
        
        allocate(flow(ncell),flow_lo(ncell), flow_hi(ncell),stat = istat)
        if (istat .ne. 0 )then
           call stm_fatal(message)
        end if
        flow      = LARGEREAL
        flow_lo   = LARGEREAL
        flow_hi   = LARGEREAL
        return
    end subroutine
    
    !> Deallocate the state variables
    !> including concentration and hydrodynamics
    !> and reset ncell and nvar to zero.
    subroutine deallocate_state
        implicit none
        ncell = 0
        nvar  = 0
        deallocate(conc, conc_prev, mass, mass_prev)
        deallocate(area)
        deallocate(area_prev)
        deallocate(area_lo,area_hi)
        deallocate(area_lo_prev, area_hi_prev)
        deallocate(flow, flow_lo, flow_hi)
        return
    end subroutine

end module



