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

!> Test of mass transport convergence no flow
!>@ingroup test
module test_linear_decay_no_flow

use stm_precision
use gaussian_init_boundary_condition

integer :: istep = 0
integer, parameter  :: nstep_base = 4*8*4
integer, parameter  :: nx_base = 4
real(stm_real), parameter :: total_time = 2048.d0
real(stm_real), parameter :: start_time = zero           !< starts at zero

contains

!> Subroutine that runs a decay in the absence of diffusion
subroutine test_linear_decay_convergence(verbose)

use test_convergence_transport
use hydro_data
use hydro_uniform_flow
use boundary_advection
use boundary_diffusion
use source_sink
use error_metric

implicit none

procedure(hydro_data_if),pointer :: no_flow_hydro
integer, parameter  :: nconc = 2
integer :: icell
character(LEN=*),parameter :: label = "linear_decay_no_flow"
logical  :: verbose

real(stm_real),parameter :: domain_length = 512.d0
real(stm_real) :: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real) :: xposition (nx_base)
real(stm_real), parameter :: decay_rate = 0.0003d0
real(stm_real),dimension(nconc) :: decay_rates = (/decay_rate,decay_rate/)

call set_uniform_flow_area(zero,120.d0)
no_flow_hydro => uniform_flow_area
call set_linear_decay(decay_rates,nconc)
boundary_diffusion_flux => no_diffusion_flux
boundary_diffusion_matrix => no_diffusion_matrix
compute_source => linear_decay_source

!> Subroutine which generates fine initial values and reference values to compare with 
!> and feed the covvergence test subroutine.

fine_initial_condition(:,1) = ten
fine_initial_condition(:,2) = fine_initial_condition(:,1)

fine_solution(:,1) = fine_initial_condition(:,1) *exp(- decay_rate*total_time)
fine_solution(:,2) = fine_initial_condition(:,2) *exp(- decay_rate*total_time)


call test_convergence(label,                  &
                      no_flow_hydro,          &
                      zero_advective_flux,      &
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
                      verbose,.true.)

end subroutine


end module