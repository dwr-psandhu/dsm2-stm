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



contains

!> Test the extrapolation part of the predictor step of the advection algorithm
!> There are three cells and two constituents to test symmetry in space across constituents
subroutine test_tidal_advection

use stm_precision
use advection
use example_initial_conditions
use fruit

implicit none 

  integer,parameter :: ncell =  512     !interior and two ends
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
  real(stm_real) :: area_hi(ncell)
  real(stm_real) :: area_lo(ncell)
  real(stm_real) :: area(ncell)
  real(stm_real) :: area_prev(ncell)  
  real(stm_real) :: dx
  real(stm_real) :: dt
  real(stm_real) :: time
  integer        :: ivar,jvar,kvar
  
  !----- local
  real(stm_real),parameter :: origin = zero 
  real(stm_real),parameter :: amp = half
  real(stm_real),parameter :: gravity = 9.80d0
  real(stm_real),parameter :: big_l = 409600.0d0
  real(stm_real),parameter :: depth =16.0d0
  real(stm_real),parameter :: omega= 0.506708d0 ! hr
  real(stm_real),parameter :: big_a = 0.4188704167062d0
  real(stm_real),parameter :: big_b = 0.040465522644d0
  real(stm_real),parameter :: k_0 = 10.066637844459d0
  !real(stm_real):: xpos (ncell)
  real(stm_real),parameter :: start_time = zero 
  real(stm_real),parameter :: end_time = 124d0 ! 12.4 is one exact M2 cycle of tide in hours
  real(stm_real):: sd 
  real(stm_real):: mean
  real(stm_real):: scale
  integer ,parameter :: ntime = 124
   
  ! initail mass at t=0
  
  dx = big_l/dble(ncell)
  time = zero
  dt = one  !hr
 
  area(:)= two
  area_lo(:)= two
  area_hi(:)= two
  scale = one
  sd = big_l/(two*two*two)
  mean = big_l/two
! time?
do jvar=1,nvar
 
 time=start_time
  
call fill_gaussian(mass_prev(:,jvar),ncell,origin,dx,mean,sd,scale)
  
 
do ivar=1,ntime

 time = time + dt 
  
do kvar = 1,ncell
   flow_lo =area_lo(kvar)*big_a*sin(big_b*(big_l- (kvar*800.0d0-800.0d0)))*sin(omega*time)
   flow_hi = area_hi(kvar)*big_a*sin(big_b*(big_l- (kvar*800.0d0)))*sin(omega*time)
   flow = area(kvar)*big_a*sin(big_b*(big_l- (kvar*800.0d0-400.0d0)))*sin(omega*time)
end do
  
  call advect(mass,     &
                  mass_prev,&
                  flow,     &                  
                  flow_lo,  &
                  flow_hi,  &
                  area,     &
                  area_prev,&
                  area_lo,  &
                  area_hi,  &
                  ncell,    &
                  nvar,     &
                  time,     &
                  dt,       &
                  dx)
  
  

 
  mass_prev = mass
  
  end do

 
call fill_gaussian(mass,ncell,origin,dx,mean,sd,scale) 

call assertEquals(mass(ncell/2-20,jvar) ,mass_prev(ncell/2-20,jvar),1.0d-4, "error in tidal boundary") 
call assertEquals(mass(ncell/2-2,jvar) ,mass_prev(ncell/2-2,jvar),1.0d-4, "error in tidal boundary") 
call assertEquals(mass(ncell/2-1,jvar) ,mass_prev(ncell/2-1,jvar),1.0d-4, "error in tidal boundary") 
call assertEquals(mass(ncell/2,jvar) ,mass_prev(ncell/2,jvar),1.0d-4, "error in tidal boundary")  
call assertEquals(mass(ncell/2 +1,jvar) ,mass_prev(ncell/2+1,jvar),1.0d-4, "error in tidal boundary")
call assertEquals(mass(ncell/2 +2,jvar) ,mass_prev(ncell/2+2,jvar),1.0d-4, "error in tidal boundary")   
call assertEquals(mass(ncell/2+20,jvar) ,mass_prev(ncell/2+20,jvar),1.0d-4, "error in tidal boundary")       
 
 
 
 
end do ! on nvar
 
return
end subroutine

end module



