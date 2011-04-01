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

!> Routines provide the calculation for non-cohesive suspended sediment erosion and deposition functions.
!> The functions in this module are: 1.	Garcia and Parker (1991) 2. van Rijn (1984b)
!>  3.	Smith and McLean (1977) 4.	Zyserman and Fredsoe (1994)  for more details see ASCE manual #110 page 117.
!>@ingroup sediment 

module non_cohesive_source


contains 

subroutine source_non_cohesive(pick_up_flag)

use stm_precision 
use suspended_sediment_variable
use sediment_variables

implicit none
character(len=32), optional, intent(in) :: pick_up_flag
!---local
character :: default_function



call
call




    !subroutine entrainment_garcia_parker()
    !use suspended_utility
    !implicit none 

   ! end subroutine 

end subroutine

end module 