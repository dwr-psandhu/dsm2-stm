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
module test_advection_tidal
use stm_precision
!----- module variables
! todo: make the names more meaningful
real(stm_real),parameter :: origin = zero                  !< Left hand side of the channel
real(stm_real),parameter :: domain_length = 102400.0d0     !< Domain Length in meter
real(stm_real),parameter :: amplitude = fourth             !< Tidal amplitude in meter    
real(stm_real),parameter :: gravity = 9.80d0               !< Gravitational acceleration in m/s^2
real(stm_real),parameter :: depth = 16.0d0                 !< Channel depth in meter
real(stm_real),parameter :: sec_per_hr = 60.d0*60.d0       !< Convert factor of hour to second 
real(stm_real),parameter :: m2_period = 12.4d0*sec_per_hr  !< M2 tidal period 
real(stm_real),parameter :: freq=two*pi/m2_period          !< Frequency of tidal oscillation


contains

!> Tests the convergence of error rate in advection of mass which is exposed a tidal boundary 
subroutine test_tidal_advection_convergence(verbose)

use test_single_channel_advection
use hydro_data
use source_sink

    
implicit none
procedure(hydro_data_if),pointer :: tidal_hydro          !< The pointer points to tidal flow data

integer, parameter  :: nconc = 2                         !< Number of constituents
integer, parameter  :: nstep_base = 128                  !< Number of time steps in finer discritization
integer, parameter  :: nx_base    = 512                  !< Number of spatial discritization in finer mesh 
logical :: verbose
real(stm_real), parameter :: total_time = m2_period      !< total time of the test
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: ic_center = domain_length/two          !< Center of initial condition
real(stm_real) :: solution_center = domain_length/two    !< Center of final solution 
real(stm_real) :: ic_gaussian_sd = domain_length/64.d0   !< Standard deviation of initial values 
real(stm_real) :: solution_gaussian_sd = domain_length/64.d0 !< Standard deviation of final values

character(LEN=*),parameter :: label = "tidal_flow"      
tidal_hydro=> tidal_flow
compute_source => no_source

!> load the initial values and reference final values to feed the test routine
call initial_fine_solution_tidal(fine_initial_condition, &
                                 fine_solution,          &
                                 nx_base,                &
                                 nconc,                  &
                                 origin,                 &
                                 domain_length,          &
                                 ic_gaussian_sd,         &
                                 solution_gaussian_sd,   &
                                 ic_center,              &
                                 solution_center)


!> The general subroutine which gets the fine initial and reference values from the privious subroutine and 
!> compute the norms, after each step coarsen the values and repeat computation.
!> at the end  calculates the ratio of the norms and prints a log 
call test_advection_convergence(label,                  &
                                tidal_hydro,            &
                                domain_length,          &
                                total_time,             &
                                fine_initial_condition, &
                                fine_solution,          &            
                                nstep_base,             &
                                nx_base,                &
                                nconc,                  &
                                verbose)

return
end subroutine
!-------------------------------------------
!> Generates a fine initial and final solution of mass
subroutine initial_fine_solution_tidal(fine_initial_condition, &
                                       fine_solution,          &
                                       nx_base,                &
                                       nconc,                  &
                                       origin,                 &
                                       domain_length,          &
                                       ic_gaussian_sd,         &
                                       solution_gaussin_sd,    &
                                       ic_center,              &
                                       solution_center)
                                       


use gaussian_init_boundary_condition
use stm_precision

implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 
real(stm_real),intent(out) :: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real),intent(in)  :: ic_center                             !< Center of fine initial value
real(stm_real),intent(in)  :: solution_center                       !< Center of solution
real(stm_real),intent(in)  :: ic_gaussian_sd                        !< Standard deviation of initial value
real(stm_real),intent(in)  :: solution_gaussin_sd                   !< Standard deviation of solution 
real(stm_real),intent(in)  :: origin                                !< Left hand side of the channel
real(stm_real),intent(in)  :: domain_length                         !< Domain length
!----local
real(stm_real):: dx

dx = domain_length/nx_base

call fill_gaussian(fine_initial_condition(:,1),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd,0.1d0)
call fill_gaussian(fine_initial_condition(:,2),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd,0.1d0)

call fill_gaussian(fine_solution(:,1),nx_base,origin,dx, &
                   solution_center,solution_gaussin_sd,0.1d0)
call fill_gaussian(fine_solution(:,2),nx_base,origin,dx, &
                   solution_center,solution_gaussin_sd,0.1d0)

return
end subroutine
!///////////////////////////////////
!> generates a tidal flow to feed a pointer
subroutine tidal_flow(flow,    &
                      flow_lo, &
                      flow_hi, &
                      area,    &
                      area_lo, &
                      area_hi, &
                      ncell,   &
                      time,    &
                      dx,      &
                      dt)
                      
use stm_precision

implicit none
integer, intent(in) :: ncell                   !< number of cells
real(stm_real), intent(in)  :: time            !< time of request "old time"
real(stm_real), intent(in)  :: dx              !< spatial step 
real(stm_real), intent(in)  :: dt              !< time step 
real(stm_real), intent(out) :: flow(ncell)     !< cell and time centered flow
real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow, time centered
real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow, time centered
real(stm_real), intent(out) :: area(ncell)     !< cell center area, old time
real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face, time centered
real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face, time centered

!--- local
real(stm_real) :: big_b 
real(stm_real) :: big_a 
real(stm_real) :: vel_lo
real(stm_real) :: vel_hi
real(stm_real) :: vel
integer :: icell

big_b = freq/sqrt(gravity*depth)
big_a = amplitude* sqrt(gravity*depth)/(depth*cos(big_b*domain_length))

! width is assumed to be equal to 1 meter 
do icell = 1,ncell  
  area(icell)    = depth + amplitude * cos(big_b*(domain_length-(dble(icell)-half)*dx))/cos(big_b*domain_length)*cos(freq*time)  
  area_lo(icell) = depth + amplitude * cos(big_b*(domain_length-(dble(icell-1)*dx)))   /cos(big_b*domain_length)*cos(freq*time)  
  area_hi(icell) = depth + amplitude * cos(big_b*(domain_length-(dble(icell)*dx)))     /cos(big_b*domain_length)*cos(freq*time)  
  vel_lo = big_a*sin(big_b*(domain_length - (dble(icell-1)*dx))     )*sin(freq*time)
  vel_hi = big_a*sin(big_b*(domain_length - (dble(icell)*dx  ))     )*sin(freq*time)
  vel    = big_a*sin(big_b*(domain_length - ((dble(icell)-half)*dx)))*sin(freq*time)
  flow(icell)    = area(icell)*vel
  flow_lo(icell) = area_lo(icell)*vel_lo
  flow_hi(icell) = area_hi(icell)*vel_hi
end do  

return
end subroutine 

end module


