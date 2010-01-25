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

!> todo: write tests for boundary diffusive flux subroutine
!>@ingroup test
module test_boundary_difussive_flux

use diffusion 
use fruit
use stm_precision

contains

subroutine test_boundary_dif_flux
  use diffusion 
  implicit none
  
  integer,parameter :: nvar = 2  !< Number of variables
  integer,parameter :: ncell = 10 ! todo do we need this?
  
real(stm_real) :: diffusive_flux_boundary_lo_prev(nvar)         !< Explicit diffusive boundary flux low side old time
real(stm_real) :: diffusive_flux_boundary_hi_prev(nvar)         !< Explicit diffusive boundary flux high side old time
real(stm_real) :: diffusive_flux_boundary_lo(nvar)              !< Explicit diffusive boundary flux low side new time
real(stm_real) :: diffusive_flux_boundary_hi(nvar)              !< Explicit diffusive boundary flux high side new time
real(stm_real) :: disp_coef_lo (ncell,nvar)! todo do we need this?              !< Low side constituent dispersion coef. at new time
real(stm_real) :: area_lo (ncell) ! todo do we need this?                       !< Low side area at new time

   disp_coef_lo(:,:) = LARGEREAL ! todo do we need this?
   area_lo(:) = LARGEREAL  ! todo do we need this?
   
   
call boundary_diffusive_flux(diffusive_flux_boundary_lo,            &
                                   diffusive_flux_boundary_hi,            &
                                   diffusive_flux_boundary_lo_prev,       &
                                   diffusive_flux_boundary_hi_prev,       &
                                   nvar)
                                   
                                   
                                   
  call assertEquals (diffusive_flux_boundary_lo(nvar),zero,1d-8,"Error in diffusive boundary flux low at new time")
  call assertEquals (diffusive_flux_boundary_hi(nvar),zero,1d-8,"Error in diffusive boundary flux high at new time")
  call assertEquals (diffusive_flux_boundary_lo_prev(nvar),zero,1d-8,"Error in diffusive boundary flux low at old time")
  call assertEquals (diffusive_flux_boundary_hi_prev(nvar),zero,1d-8,"Error in diffusive boundary flux high at old time")



return
end subroutine test_boundary_dif_flux

end module