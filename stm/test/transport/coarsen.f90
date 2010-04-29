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

!> Routines for grid refinement/coarsening operations.
module grid_refinement

contains

!> Coarsen a solution at a fine level of resolution
subroutine coarsen(coarse_data,fine_data,ncell_fine,ncell_coarse, nvar)

use stm_precision
use error_handling

implicit none
!---arg
integer,intent(in) :: ncell_coarse
integer,intent(in) :: ncell_fine
integer,intent(in) :: nvar
real(stm_real), intent(in) :: fine_data(ncell_fine,nvar)
real(stm_real), intent(out):: coarse_data(ncell_coarse,nvar)

!---locals
real(stm_real) :: coarsen_factor
integer :: ivar
integer :: icell
integer :: i_coarse


if ( mod(ncell_fine , ncell_coarse) /= 0) then
    call stm_fatal("Coarsening factor is not an integer!")  
else

coarsen_factor = ncell_fine/ncell_coarse
    do ivar=1,nvar
        do icell=1,ncell_coarse
            coarse_data(icell,ivar) = zero
            i_coarse = 0
            do while (i_coarse < coarsen_factor) 
              coarse_data(icell,ivar) = coarse_data(icell,ivar)+ fine_data(icell*coarsen_factor-i_coarse,ivar)
              i_coarse= i_coarse + 1   
            end do
            coarse_data(icell,ivar)= coarse_data(icell,ivar)/dble(coarsen_factor)
        end do
    end do
    
end if


return
end subroutine coarsen

end module grid_refinement