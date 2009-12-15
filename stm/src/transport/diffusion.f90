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

!> Module orchestrating the diffusion scheme. The main
!> routine in the module is diffuse().
!> Explicit and implicit diffusion operators are included here.
!>@ingroup transport
module diffusion
use stm_precision

contains

! ///////////////////////////////////////////////////////////////////

!> Calculates the diffusive portion of the constituent transport.
!> It contains an explicit version of the diffusion operator and a general (involving all
!> potential cases) diffusion operator as well, with a coefficient theta_stm for 
!> selecting the level of implicitness. (theta_stm=0.5 is Crank Nicolson.).
!> The matrix is solved via a tri-diagonal solver.  
subroutine diffuse(conc,     &
                  conc_prev,&
                  area,     &
                  area_prev,&
                  area_lo,  &
                  area_hi,  &
                  ks_lo,    &
                  ks_hi,    &
                  ncell,    &
                  nvar,     &
                  time,     &
                  theta_stm,&
                  dt,       &
                  dx)

use primitive_variable_conversion
implicit none

integer, intent (in) :: ncell !< Number of cells
integer, intent (in) :: nvar  !< Number of variables

real(stm_real), intent (out) :: conc(ncell,nvar)     !< Concentration at new time
real(stm_real), intent (in) :: conc_prev(ncell,nvar) !< Concentration at old time
real(stm_real), intent (in) :: area (ncell,nvar)     !< Cell-centered area at new time
real(stm_real), intent (in) :: area_prev (ncell,nvar)!< Cell-centered area at old time
real(stm_real), intent (in) :: area_lo (ncell,nvar)  !< Low side area centered intime
real(stm_real), intent (in) :: area_hi (ncell,nvar)  !< High side area centered in time 
real(stm_real), intent (in) :: ks_lo (ncell,nvar)    !< Low side sediment dispersion coef. centered in time
real(stm_real), intent (in) :: ks_hi (ncell,nvar)    !< High side sediment dispersion coef. centered in time
real(stm_real), intent (in) :: time                  !< Current time
real(stm_real), intent (in) :: theta_stm             !< Expilicitness coefficient 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in) :: dt                    !< Time step   
real(stm_real), intent (in) :: dx                    !< Spacial step  




! This routine gives the effects of diffusion fluxes on each cell
! for a single time step (ie, explicit). This is needed for the advection step.
! It is also probably part of the right hand side of the implicit diffusion solver 
! matrix calculation. 
call explicit_diffusion_operator(,     &
                  dx)

! 
call construct_diffusion_matrix()

call construct_right_hand_side()

call apply_diffusion_boundary()

call solve

end subroutine

end module diffusion
