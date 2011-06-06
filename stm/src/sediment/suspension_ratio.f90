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

!> Routine provides the calculation for sediment driver to determine allocation coeficient
!> todo: this routine may move to another place
!>@ingroup sediment 

module suspenstion_ratio


contains
! The formula here is adopted from B. Greimann, Y. Lai and J. Huang, 2008
subroutine allocation_ratio(suspended_percent,    &
                            rouse_number,         &
                            diameter,             &
                            nclass,               &
                            ncell)
                            
use stm_precision

implicit none
integer,intent(in) :: nclass                                    !< Number of sediment classes 
integer,intent(in) :: ncell                                     !< Number of cells



return
end subroutine 


end module