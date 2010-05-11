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

!> todo: write tests boundary diffusion matrix
!>@ingroup test
module test_boundary_diffusion_matrix

use boundary_diffusion
use fruit
use stm_precision

contains

subroutine sub_boundary_diffusion_matrix

use boundary_diffusion
use stm_precision

     implicit none
         !--- args
                                       
        integer, parameter :: ncell =10                               !< Number of cells
        integer, parameter :: nvar = 2                                !< Number of variables

        real(stm_real) :: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
        real(stm_real) :: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
        real(stm_real) :: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
        real(stm_real) :: area (ncell)                                !< Cell centered area at new time 
        real(stm_real) :: area_lo(ncell)                              !< Low side area at new time
        real(stm_real) :: area_hi(ncell)                              !< High side area at new time 
        real(stm_real) :: disp_coef_lo (ncell,nvar)                   !< Low side constituent dispersion coef. at new time
        real(stm_real) :: disp_coef_hi (ncell,nvar)                   !< High side constituent dispersion coef. at new time
        real(stm_real) :: time                                        !< Current time
        real(stm_real) :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
        real(stm_real) :: dx                                          !< Spatial step  
        real(stm_real) :: dt                                          !< Time step  


down_diag (:,:)= two
center_diag(:,:) = two
up_diag (:,:) = two
area (:)= one
area_lo (:) = one
area_hi(:) = one
disp_coef_lo(:,:)= 0.1
disp_coef_hi(:,:) = 0.1
time = largereal
theta_stm =0.5d0 
dx = 0.5d0
dt = 0.01d0

 
!
!
!  call assertEquals (conc(1),613.689083199382D0,1d-9,"problem in solving cell num.1")
! 
return
end subroutine sub_boundary_diffusion_matrix

end module