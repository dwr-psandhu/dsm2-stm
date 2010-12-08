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

!> Test of mass transport convergence in uniform flow
!>@ingroup test
module test_uniform_flow

use stm_precision

integer, parameter  :: nstep_base = 256
integer, parameter  :: nx_base = 256
real(stm_real), parameter :: total_time = 25600.d0
real(stm_real), parameter :: start_time = zero           !< starts at zero

contains
!> Subroutine that runs a small advective simulation
subroutine test_uniform_adv_biidirectional_convergence(verbose)

use test_convergence_transport
use hydro_data
use hydro_uniform_flow
use source_sink
use diffusion
use boundary_advection
use boundary_diffusion
use gaussian_init_boundary_condition

implicit none

procedure(hydro_data_if),pointer :: uniform_hydro
integer, parameter  :: nconc = 2
logical  :: verbose

real(stm_real), parameter :: domain_length = 25600.d0
real(stm_real), parameter :: origin =zero
real(stm_real), parameter :: constant_flow = 600.d0
real(stm_real), parameter :: constant_area = 1000.d0 
real(stm_real), parameter :: reverse_time =  total_time /two
real(stm_real), parameter :: ic_center = origin + domain_length/three
real(stm_real), parameter :: ic_peak = one
real(stm_real), parameter :: ic_gaussian_sd = domain_length/32.d0
real(stm_real) :: solution_gaussian_sd = ic_gaussian_sd
real(stm_real) :: solution_center = ic_center
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution

character(LEN=*),parameter :: label = "advection_bidirectional_uniform_dirichlet"

call set_uniform_flow_area(constant_flow,constant_area,reverse_time)
! todo: force these to be set so they aren't just left over from last test
! reverse_time is optional
advection_boundary_flux => zero_advective_flux
uniform_hydro=> uniform_flow_area
compute_source => no_source



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


call test_convergence(label,                  &
                      uniform_hydro,          &
                      zero_advective_flux,    &
                      no_diffusion_flux,      &
                      no_diffusion_matrix,    &
                      no_source,              &
                      domain_length,          &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose,                &
                      detail_printout =.true.)

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