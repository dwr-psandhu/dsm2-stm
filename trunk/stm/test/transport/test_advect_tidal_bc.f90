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
use state_variables
use stm_precision
use advection
use example_initial_conditions
use fruit

implicit none 

real(stm_real) :: dx
real(stm_real) :: dt
real(stm_real) :: time
integer        :: itime,ivar,icell

!----- local
real(stm_real),parameter :: eps = 1.d-4 ! criterion for error 
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

call allocate_state(ncell,nvar)
 
! initail mass at t=0

dx = big_l/dble(ncell)
time = zero
dt = one  !hr

! todo: Kaveh, these need to change?
area(:)= two
area_lo(:)= two
area_hi(:)= two
scale = one
sd = big_l/(two*two*two)
mean = big_l/two

time=start_time
do ivar=1,nvar
  call fill_gaussian(mass_prev(:,ivar),ncell,origin,dx,mean,sd,scale)
end do
 
do itime=1,ntime
  time = time + dt 
  
  do icell = 1,ncell
    flow_lo =area_lo(icell)*big_a*sin(big_b*(big_l- (icell*800.0d0-800.0d0)))*sin(omega*time)
    flow_hi = area_hi(icell)*big_a*sin(big_b*(big_l- (icell*800.0d0)))*sin(omega*time)
    flow = area(icell)*big_a*sin(big_b*(big_l- (icell*800.0d0-400.0d0)))*sin(omega*time)
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

do ivar=1,nvar 
  call fill_gaussian(mass(:,ivar),ncell,origin,dx,mean,sd,scale) 
end do

do ivar=1,nvar
  call assertEquals(mass(ncell/2-20,ivar),mass_prev(ncell/2-20,ivar) ,eps, "error in tidal boundary (1)") 
  call assertEquals(mass(ncell/2-2,ivar) ,mass_prev(ncell/2-2,ivar)  ,eps, "error in tidal boundary (2)") 
  call assertEquals(mass(ncell/2-1,ivar) ,mass_prev(ncell/2-1,ivar)  ,eps, "error in tidal boundary (3)") 
  call assertEquals(mass(ncell/2,ivar)   ,mass_prev(ncell/2,ivar)    ,eps, "error in tidal boundary (4)")  
  call assertEquals(mass(ncell/2 +1,ivar) ,mass_prev(ncell/2+1,ivar) ,eps, "error in tidal boundary (5)")
  call assertEquals(mass(ncell/2 +2,ivar) ,mass_prev(ncell/2+2,ivar) ,eps, "error in tidal boundary (6)")   
  call assertEquals(mass(ncell/2+20,ivar) ,mass_prev(ncell/2+20,ivar),eps, "error in tidal boundary (7)")       
end do

call deallocate_state
return
end subroutine

end module



