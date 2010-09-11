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

!> Test of transport diffusion convergence test for a single channel
!>@ingroup test
module test_diffusion_single_channel

contains

!> Subroutine that checks the error convergence ratio for diffusion routine 
subroutine test_diffusion_convergence_single_channel(verbose)

use stm_precision
use test_convergence_transport
use state_variables
use primitive_variable_conversion
use boundary_diffusion
use diffusion
use source_sink
use gaussian_init_boundary_condition
use hydro_data
use hydro_uniform_flow
use error_metric
use error_handling
use fruit
use logging
use log_convergence

implicit none

!--- Problem variables

integer, parameter  :: nstep_base = 16*16
integer, parameter  :: nx_base = 128*8

integer :: icoarse 
integer :: nstep
integer :: nx

integer, parameter  :: nconc = 2
real(stm_real), parameter :: domain_length = 51200.d0
real(stm_real), parameter :: origin = zero   
real(stm_real), parameter :: total_time    = 8000.0d0
real(stm_real), parameter :: disp_coef     = 10.5d0
real(stm_real), parameter :: theta = half                       !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real),allocatable :: disp_coef_lo(:)     !< Low side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_hi(:)     !< High side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_lo_prev(:) !< Low side constituent dispersion coef. at old time
real(stm_real),allocatable :: disp_coef_hi_prev(:) !< High side constituent dispersion coef. at old time
logical, optional :: verbose

real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real), parameter :: channel_area = 120.0d0
real(stm_real), parameter :: start_time = 2048.d0     !< Initial condition depends on this
real(stm_real), parameter :: end_time = start_time + total_time

real(stm_real) :: time
real(stm_real), allocatable :: xposition(:)

integer :: itime
integer :: icell 
integer :: ivar
!------
integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
integer :: coarsening
integer, parameter :: nrefine = 3
integer :: which_cell(nrefine)

real(stm_real),allocatable :: reference(:)
real(stm_real) norm_error(3,nrefine)
character(LEN=64) filename
procedure(hydro_data_if),pointer :: no_flow_hydro


real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: ic_center = domain_length/two          !< Center of initial condition
real(stm_real) :: solution_center = domain_length/two    !< Center of final solution 
real(stm_real) :: ic_gaussian_sd = domain_length/64.d0   !< Standard deviation of initial values 
real(stm_real) :: solution_gaussian_sd = domain_length/64.d0 !< Standard deviation of final values
real(stm_real) :: fine_dx = domain_length/dble(nx_base)

character(LEN=*),parameter :: label = "diffusion_uniform" 
compute_source => no_source
boundary_diffusion_matrix  => neumann_diffusion_matrix
boundary_diffusion_flux    => neumann_gaussian_diffusive_flux
call set_uniform_flow_area(zero, channel_area)
no_flow_hydro => uniform_flow_area
const_dispersion = disp_coef  ! todo: this is still a bit rough
use_diffusion = .true.

do ivar = 1,nconc
   call fill_gaussian(fine_initial_condition(:,ivar),nx_base,origin,fine_dx,half*domain_length, & 
                   sqrt(two*disp_coef*start_time),one)
   call fill_gaussian(fine_solution(:,ivar),nx_base,origin,fine_dx,half*domain_length, &
                       sqrt(two*disp_coef*end_time),sqrt(start_time/end_time))                   
end do                   
 


!> The general subroutine which gets the fine initial and reference values from the privious subroutine and 
!> compute the norms, after each step coarsen the values and repeat computation.
!> at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      no_flow_hydro,          &
                      domain_length,          &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose)

return
end subroutine


end module