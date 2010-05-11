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
use boundary_diffusion

contains

subroutine test_boundary_dif_flux

  use diffusion 
  use boundary_diffusion
  
  implicit none
  
  integer,parameter :: nvar = 1                        !< Number of variables
  integer,parameter :: ncell = 10                      !<Number of cells
  
real(stm_real) :: diffusive_flux_lo(ncell,nvar)        !< Explicit diffusive boundary flux low side old time
real(stm_real) :: diffusive_flux_hi(ncell,nvar)        !< Explicit diffusive boundary flux high side old time
real(stm_real) :: conc(ncell,nvar)                     !< Explicit diffusive boundary flux low side new time
real(stm_real) :: area_lo         (ncell)              !< Low side area centered at old time
real(stm_real) :: area_hi         (ncell)              !< High side area centered at old time
real(stm_real) :: disp_coef_lo (ncell,nvar)            !< Low side constituent dispersion coef.
real(stm_real) :: disp_coef_hi (ncell,nvar)            !< High side constituent dispersion coef.

real(stm_real) :: time = zero                          !< time 
real(stm_real) :: dx = zero                            !< dx
              
  
   
boundary_diffusion_flux => neumann_no_flow_diffusive_flux

call boundary_diffusion_flux(diffusive_flux_lo, &
                             diffusive_flux_hi, &
                             conc,              &
                             area_lo,           &
                             area_hi,           &
                             disp_coef_lo,      &  
                             disp_coef_hi,      &
                             ncell,             &
                             nvar,              &
                             time)
                                                           
  call assertEquals (diffusive_flux_lo(1,nvar),zero,1d-8,"Error in diffusive boundary flux low at new time")
  call assertEquals (diffusive_flux_hi(ncell,nvar),zero,1d-8,"Error in diffusive boundary flux high at new time")
  call add_fail("Boundary diffusion flux not really tested in non-trivial way.")


return
end subroutine test_boundary_dif_flux

end module