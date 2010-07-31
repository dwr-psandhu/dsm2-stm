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

integer, parameter  :: nstep_base = 256*4
integer, parameter  :: nx_base = 512*4
integer, parameter  :: nconc = 2
real(stm_real), parameter :: total_time = 640.D0
real(stm_real), parameter :: decay_rate = 0.002d0

! exp(0.001*640) =0.527

contains
!> Subroutine that runs a small advective simulation
subroutine test_advection_decay_convergence(verbose)

use test_single_channel_advection
use hydro_data
use source_sink
use boundary_advection_module

implicit none
procedure(hydro_data_if),pointer :: uniform_hydro

logical  :: verbose

real(stm_real)   ,parameter :: domain_length = 51200.d0
real(stm_real)   ,parameter :: origin =zero
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: fine_initial_area(nx_base)  !< initial area at finest resolution
real(stm_real) :: fine_final_area(nx_base)    !< final area at finest resolution
real(stm_real) :: ic_center = domain_length/two
real(stm_real) :: solution_center = domain_length/two
real(stm_real) :: ic_gaussian_sd = domain_length/(ten*two)
real(stm_real) :: solution_gaussian_sd = domain_length/(ten*two)
integer :: icell
character(LEN=*),parameter :: label = "uniform flow_liner decay"
uniform_hydro=> uniform_flow
compute_source => linear_decay_source
replace_boundary_flux      => neumann_advective_flux


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

call test_advection_convergence(label,                   &
                                 uniform_hydro,          &
                                 domain_length,          &
                                 total_time,             &
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
    real(stm_real), parameter :: constant_flow = 254.0d1
    real(stm_real), parameter :: constant_area = 27.0d1 


    if (time <= (total_time/two)) then
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
!===========
subroutine linear_decay_source(source, & 
                               conc,   &
                               area,   &
                               flow,   &
                               ncell,  &
                               nvar,   &
                               time)
                                     

 use  primitive_variable_conversion
 implicit none
 
 !--- args
integer,intent(in)  :: ncell                      !< Number of cells
integer,intent(in)  :: nvar                       !< Number of variables
real(stm_real),intent(inout) :: source(ncell,nvar)!< cell centered source 
real(stm_real),intent(in)  :: conc(ncell,nvar)    !< Concentration
real(stm_real),intent(in)  :: area(ncell)         !< area at source     
real(stm_real),intent(in)  :: flow(ncell)         !< flow at source location
real(stm_real),intent(in)  :: time                !< time 
!--- local just for test
real(stm_real) :: mass (ncell,nvar)
real(stm_real) :: rate_1
real(stm_real) :: rate_2

rate_1 = decay_rate
rate_2 = rate_1

! source must be in primitive variable 
call prim2cons(mass,conc,area,ncell,nvar)
source(:,1) = -rate_1*mass(:,1)
source(:,2) = -rate_2*mass(:,2) 
 
return
end subroutine 

end module