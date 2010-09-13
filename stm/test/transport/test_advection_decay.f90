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

!> Test of decay of mass convergence in uniform flow
!>@ingroup test
module test_linear_decay_const_flow

use stm_precision

integer, parameter  :: nstep_base = 128
integer, parameter  :: nx_base = 512*2
integer, parameter  :: nconc = 2
real(stm_real), parameter :: total_time = 600.d0
real(stm_real), parameter :: start_time = zero           !< starts at zero
real(stm_real), parameter :: decay_rate = 0.01d0
real(stm_real), parameter :: constant_flow = sixteen
real(stm_real), parameter :: constant_area = two

contains
!> Subroutine that runs a small advective simulation
subroutine test_advection_decay_convergence(verbose)

use test_convergence_transport
use hydro_data
use source_sink
use boundary_advection
use boundary_diffusion
use logging

implicit none
procedure(hydro_data_if),pointer :: uniform_hydro

logical  :: verbose

real(stm_real)   ,parameter :: domain_length = 51200.d0
real(stm_real)   ,parameter :: origin =zero
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: ic_center = domain_length/three
real(stm_real) :: solution_center
real(stm_real) :: ic_gaussian_sd = domain_length/(sixteen*two)
real(stm_real) :: solution_gaussian_sd 
real(stm_real),dimension(nconc) :: decay_rates = (/decay_rate,decay_rate/)
real(stm_real) :: dx

character(LEN=64) :: label = "uniform_flow_linear_decay"
uniform_hydro=> uniform_flow
call set_linear_decay(decay_rates,nconc)
advection_boundary_flux  => zero_advective_flux
boundary_diffusion_flux => no_diffusion_flux
boundary_diffusion_matrix => no_diffusion_matrix
compute_source => linear_decay_source


dx = domain_length/dble(nx_base)
solution_center = ic_center + total_time*constant_flow/constant_area
solution_gaussian_sd = ic_gaussian_sd
!> Subroutine which generates fine initial values and reference values to compare with 
!> and feed the covvergence test subroutine.
call initial_fine_solution_uniform(fine_initial_condition, &
                                   fine_solution,          &
                                   nx_base,                &
                                   nconc,                  &
                                   origin,                 &
                                   domain_length,          &
                                   ic_gaussian_sd,         &
                                   solution_gaussian_sd,   &
                                   ic_center,              &
                                   solution_center)
                                   
                         
fine_solution = fine_solution* exp(- total_time*decay_rate)

call test_convergence(label,                  &
                      uniform_hydro,          &
                      zero_advective_flux,    & !todo: this is wrong
                      no_diffusion_flux,      &
                      no_diffusion_matrix,    &
                      linear_decay_source,    &
                      domain_length,          &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose)

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
    ! todo: remove
    if (time <= (total_time/two)) then
      flow = constant_flow
    else
      !flow = minus * constant_flow !todo:
       flow =constant_flow
    end if
    flow_hi = flow
    flow_lo = flow
    area = constant_area
    area_lo = constant_area
    area_hi = constant_area
    return
end subroutine
! todo: ic_center and solution center must have dimension of NCONC
!> Generates fine solution of initial condition and final values to compare for uniform flow advection 
subroutine initial_fine_solution_uniform(fine_initial_condition, &
                                         fine_solution,          &
                                         nx_base,                &
                                         nconc,                  &
                                         origin,                 &
                                         domain_length,          &
                                         ic_gaussian_sd,         &
                                         solution_gaussian_sd,   &
                                         ic_center,              &
                                         solution_center)
                                   
use gaussian_init_boundary_condition
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