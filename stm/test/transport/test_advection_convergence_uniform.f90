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

!> Test of transport in uniform flow
!>@ingroup test
module test_uniform_flow
use stm_precision
integer :: istep = 0
integer, parameter  :: nstep_base = 40
integer, parameter  :: nx_base = 256
real(stm_real), parameter :: total_time = 6400.D0

contains


!> Subroutine that runs a small advective simulation
subroutine test_uniform_advection_convergence()

use test_single_channel_advection
use hydro_data

procedure(hydro_data_if),pointer :: uniform_hydro

integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base = 256
integer, parameter  :: nconc = 2
character(LEN=12),parameter :: label = "uniform flow"
real(stm_real)   ,parameter :: domain_length = 51200.d0
real(stm_real)   ,parameter :: origin =zero
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: fine_initial_area(nx_base)  !< initial area at finest resolution
real(stm_real) :: fine_final_area(nx_base)    !< final area at finest resolution
real(stm_real) :: ic_center = domain_length/two
real(stm_real) :: solution_center = domain_length/two
real(stm_real) :: ic_gaussian_sd = domain_length/sixteen
real(stm_real) :: solution_gaussian_sd = domain_length/sixteen

uniform_hydro=> uniform_flow

call initial_fine_solution_uniform(fine_initial_condition, &
                                   fine_solution,          &
                                   nx_base,                &
                                   nconc,                  &
                                   origin,                 &
                                   domain_length,          &
                                   ic_gaussian_sd,         &
                                   solution_gaussian_sd,   &
                                   ic_center,              &
                                   solution_center   )


call test_round_trip(label,                  &
                     uniform_hydro,          &
                     domain_length,          &
                     total_time,             &
                     fine_initial_condition, &
                     fine_solution,          &            
                     nstep_base,             &
                     nx_base,                &
                     nconc)

end subroutine
!=========================
!>generat constant area and constant flow foreward and backward
subroutine uniform_flow(flow,    &
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
    real(stm_real), intent(in) :: time             !< time of request
    real(stm_real), intent(in) :: dx               !< spatial step
    real(stm_real), intent(in) :: dt               !< time step 
    real(stm_real), intent(out) :: flow(ncell)     !< cell and time centered flow
    real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow, time centered
    real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow, time centered
    real(stm_real), intent(out) :: area(ncell)     !< cell center area, old time
    real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face, time centered
    real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face, time centered

    
    !> local
    real(stm_real), parameter :: constant_flow = 3.D2
    real(stm_real), parameter :: constant_area = 27.D1 


    if (time <= total_time/two) then
      flow = constant_flow
    else
      flow = minus * constant_flow
    end if
    flow_hi = flow
    flow_lo = flow
    area = constant_area
    area_lo = constant_area
    area_hi = constant_area
    return
end subroutine
! todo: ic_center and solution center must have dimension of NCONC
subroutine initial_fine_solution_uniform(fine_initial_condition,   &
                                           fine_solution,          &
                                           nx_base,                &
                                           nconc,                  &
                                           origin,                 &
                                           domain_length,          &
                                           ic_gaussian_sd,         &
                                           solution_gaussian_sd,   &
                                           ic_center,              &
                                           solution_center   )
                                   
use example_initial_conditions
use stm_precision
use grid_refinement
implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 

real(stm_real),intent(out) :: fine_initial_condition(nx_base,nconc)!< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)         !< reference solution at finest resolution

real(stm_real),intent(in)  :: ic_center
real(stm_real),intent(in)  :: solution_center
real(stm_real),intent(in)  :: ic_gaussian_sd
real(stm_real),intent(in)  :: solution_gaussian_sd
real(stm_real),intent(in)  :: origin 
real(stm_real),intent(in)  :: domain_length  
!----local
real(stm_real):: dx
real(stm_real):: area_const 
dx = domain_length/nx_base
!---initial condition
call fill_gaussian(fine_initial_condition(:,1),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd)
call fill_gaussian(fine_initial_condition(:,2),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd)                  

!---final solution
call fill_gaussian(fine_solution(:,1),nx_base,origin,dx, &
                   solution_center,solution_gaussian_sd)
call fill_gaussian(fine_solution(:,2),nx_base,origin,dx, &
                   solution_center,solution_gaussian_sd)                 

return
end subroutine

end module