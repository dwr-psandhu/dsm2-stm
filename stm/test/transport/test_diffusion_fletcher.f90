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

!> todo: write tests for diffusion subroutine in case of Neumann BC at the left and Dirichlet BC at the right side 
!>@ingroup test
!> Test of transport diffusion convergence test for a single channel
!>@ingroup test
module test_diffusion_fletcher

contains

!> Subroutine that checks the error convergence ratio for diffusion routine 
subroutine test_diffusion_convergence_fletcher(verbose)
use test_convergence_transport
use stm_precision
use state_variables
use primitive_variable_conversion
use boundary_diffusion
use boundary_advection
use hydro_uniform_flow
use gaussian_init_boundary_condition
use diffusion
use hydro_data
use source_sink
use log_convergence
use example_initial_conditions
use error_metric
use error_handling
use fruit
use logging

implicit none

!--- Problem variables

integer, parameter  :: nstep_base = 64
integer, parameter  :: nx_base = 64

integer :: icoarse = 0
integer :: nstep
integer :: nx

integer, parameter  :: nconc = 2
real(stm_real), parameter :: domain_length = 0.9d0
real(stm_real), parameter :: origin = 0.1d0   
real(stm_real), parameter :: total_time    = one
real(stm_real), parameter :: disp_coef     = 2.1d0
real(stm_real) :: theta = half                       !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real),allocatable :: disp_coef_lo(:)       !< Low side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_hi(:)       !< High side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_lo_prev(:)   !< Low side constituent dispersion coef. at old time
real(stm_real),allocatable :: disp_coef_hi_prev(:)   !< High side constituent dispersion coef. at old time
logical, optional :: verbose

real(stm_real) :: fine_initial_condition(nx_base,nconc)    !< initial condition f concentration at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)             !< reference solution at finest resolution

real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real), parameter :: constant_area = 100.0d0
real(stm_real), parameter :: start_time = zero    !< Initial condition depends on this
real(stm_real), parameter :: end_time = start_time + total_time

real(stm_real) :: time
procedure(hydro_data_if), pointer :: uniform_hydro => null()

integer :: itime
integer :: icell ! debug only -- remove later
!------
integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
integer :: coarsening
integer, parameter :: nrefine = 3
integer :: which_cell(nrefine)
real(stm_real),allocatable :: reference(:)
real(stm_real) :: norm_error(3,nrefine)

character(LEN=64):: label = 'test_diffusion_fletcher_dirichlet'

call set_uniform_flow_area(zero,constant_area)
uniform_hydro => uniform_flow_area

boundary_diffusion_matrix  => dirichlet_test_diffusion_matrix
boundary_diffusion_flux    => dirichlet_test_diffusive_flux
advection_boundary_flux   => zero_advective_flux
compute_source            => no_source
const_dispersion = disp_coef




call initial_final_solution(fine_initial_condition,&
                            fine_solution,         &
                            start_time,            &
                            disp_coef,             &
                            total_time,            &
                            origin,                &
                            domain_length,         &
                            nx_base,               &
                            nconc)
                                  
!=====Dirichlet
label = 'test_diffusion_fletcher_dirichlet'
call test_convergence(label,                     &
                      uniform_hydro,         &
                      zero_advective_flux,    &
                      dirichlet_test_diffusive_flux,    &
                      dirichlet_test_diffusion_matrix , &
                      no_source,              &
                      domain_length,          &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose)


!=====Neumann 
label = 'test_diffusion_fletcher_neumann'
boundary_diffusion_matrix  => n_d_test_diffusion_matrix
boundary_diffusion_flux    => n_d_test_diffusive_flux

call initial_final_solution(fine_initial_condition,&
                            fine_solution,         &
                            start_time,            &
                            disp_coef,             &
                            total_time,            &
                            origin,                &
                            domain_length,         &
                            nx_base,               &
                            nconc)

call test_convergence(label,                     &
                      uniform_hydro,           &
                      zero_advective_flux,      &
                      n_d_test_diffusive_flux,    &
                      n_d_test_diffusion_matrix , &
                      no_source,              &
                      domain_length,          &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose,.true.)

return
end subroutine



!> produce fine initial condition and reference solution 
subroutine initial_final_solution(fine_initial_conc,     &
                                  fine_solution,         &
                                  init_time,             &
                                  dispersion_coef,       &
                                  total_time,            &
                                  origin,                &
                                  domain_length,         &
                                  nx_base,               &
                                  nconc)
                                  
use gaussian_init_boundary_condition

implicit none
integer, intent(in) :: nx_base
integer, intent(in) :: nconc
real(stm_real),intent(out) :: fine_initial_conc(nx_base,nconc)     !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)    !< reference solution at finest resolution
real(stm_real),intent(in)  :: init_time
real(stm_real),intent(in)  :: dispersion_coef
real(stm_real),intent(in)  :: total_time 
real(stm_real),intent(in)  :: origin
real(stm_real),intent(in)  :: domain_length
!--local
integer :: ivar
integer :: icell

real(stm_real) :: dx
real(stm_real) :: xposition(nx_base)
real(stm_real) :: current_time

dx = domain_length/nx_base

do icell = 1,nx_base
  xposition(icell) = dx*(dble(icell)-half)+ origin
end do

current_time = init_time

fine_initial_conc(:,1) = two*xposition(:) + four*cos(half*pi*xposition(:))*exp(-dispersion_coef*current_time*pi*pi/four )
fine_initial_conc(:,2) = fine_initial_conc(:,1)

current_time = init_time + total_time
fine_solution(:,1) = two*xposition(:) + four*cos(half*pi*xposition(:))*exp(-dispersion_coef*current_time*pi*pi/four )
fine_solution(:,2) = fine_solution(:,1)

return
end subroutine





end module