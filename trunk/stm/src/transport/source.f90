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

!> Source interface to be fulfilled by driver or application
!>@ingroup transport
module source_if
 !> Calculate source
 interface compute_source
   !> Generic interface for calculating source that should be fulfilled by
   !> client programs
   subroutine compute_source(source,conc,area,flow,ncell,nvar,time)
     use stm_precision
     implicit none
     !--- args
     integer,intent(in)  :: ncell  !< Number of cells
     integer,intent(in)  :: nvar   !< Number of variables
     real(stm_real),intent(out) :: source(ncell,nvar) !< cell centered source 
     real(stm_real),intent(in)  :: conc(ncell,nvar)   !< Concentration
     real(stm_real),intent(in)  :: area(ncell)        !< area at source     
     real(stm_real),intent(in)  :: flow(ncell)        !< flow at source location
     real(stm_real),intent(in)  :: time               !< flow at source location
   end subroutine
 end interface
end module