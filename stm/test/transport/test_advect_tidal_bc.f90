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

!> Testing advection of mass which is subjected to a tidal boundary
!>@ingroup test
module test_advect_tidal_bc

use fruit

contains

!> Test the extrapolation part of the predictor step of the advection algorithm
!> There are three cells and two constituents to test symmetry in space across constituents
subroutine test_tidal_advection

use stm_precision
use advection
use example_initial_conditions

implicit none 

  integer,parameter :: ncell =  2048       !interior and two ends
  integer,parameter :: nvar = 2
 
  real(stm_real) :: mass(ncell,nvar)
  real(stm_real) :: mass_prev(ncell,nvar)
  !todo: do we need these?
  real(stm_real) :: conc(ncell,nvar)
  real(stm_real) :: conc_hi(ncell,nvar)
  real(stm_real) :: conc_lo(ncell,nvar)
  !++++++++++++++++++++++++++++++++
  !todo: is it needed?
  real(stm_real) :: source(ncell,nvar)
  real(stm_real) :: flow(ncell)
  real(stm_real) :: flow_hi(ncell)
  real(stm_real) :: flow_lo(ncell)
  real(stm_real) :: area(ncell)  
  real(stm_real) :: dx
  real(stm_real) :: dt
  real(stm_real) :: time
  integer        :: ivar,jvar
  
  !----- local
  real(stm_real),parameter :: amp = half
  real(stm_real),parameter :: gravity = 9.80d0
  real(stm_real),parameter :: big_l
  real(stm_real),parameter ::
  
  
  ! initail mass at t=0
  
  dx = big_l/ncell
  time = zero
  dt = 1  !hr
  area(:)= two
  
   
  
  
 
 
 call fill_gaussian(mass_prev,ncell,origin,dx,mean,sd,scale)
 

return
end subroutine

end module



