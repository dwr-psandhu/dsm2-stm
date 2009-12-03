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

!> Simple source term subroutines for testing
!>@ingroup test
module example_sources

contains

!> Empty source implementation
subroutine no_source(source,conc,area,flow,a_ncell,a_nvar)
 use stm_precision
 implicit none
 !--- args
 integer,intent(in)  :: a_ncell  !< Number of cells
 integer,intent(in)  :: a_nvar   !< Number of variables
 real(STM_REAL),intent(out) :: source(a_ncell,a_nvar) !< cell centered source 
 real(STM_REAL),intent(in) :: conc(a_ncell,a_nvar)    !< Concentration
 real(STM_REAL),intent(in) :: area(a_ncell,a_nvar)    !< area at source
 !> flow at source location
 real(STM_REAL),intent(in) :: flow(a_ncell,a_nvar)
 source = zero
 return
end subroutine

end module