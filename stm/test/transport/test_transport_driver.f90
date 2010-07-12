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

!> Main program unit for testing teransport
!>@ingroup test
program test_transport_driver

use fruit

use test_extrapolate
use test_prim_cons_conversion
use test_uniform_flow
use test_matrix_solver
use example_initial_conditions
use test_boundary_diffusion
use test_diffusive_flux
use test_explicit_diffusion_operator
use test_interior_coef_matrix
use test_construct_r_h_s
use test_diffusion_convergence_boundaries
use source_sink
use test_diffusion_single_channel
use test_advection_tidal
use test_coarsening
use test_uniform_flow
use test_prim_increment_to_cons
use test_gradient
use test_diffusion_convergence_boundaries_dirichlet
use test_linear_decay_no_flow
use test_advect_diffuse_react
use test_linear_decay_const_flow

implicit none
logical :: verbose = .true.

call init_fruit
!/////// Advection
call test_gradient_calc
call test_limiter
call test_prim_cons_convert
call test_prim_increment2cons
call test_example_initial_conditions
call test_extrapolation
call test_tidal_advection_convergence(verbose)
call test_uniform_advection_convergence(verbose)
!    
!!/////// Diffusion
call test_tridi_solver
call test_boundary_diffusion_flux
call test_make_dif_flux_sub
call test_explicit_interior_diffusion_op
call test_interior_coef_matrix_sub
call test_construct_elemnts_rhs 
call test_coarsen
call test_diffusion_dirichlet
call test_diffusion_neumann
call test_diffusion_convergence_single_channel
!!////// reaction
call  test_linear_decay_convergence(verbose)
call  test_advection_decay_convergence(verbose)
call  test_advect_diffuse_reaction(verbose)

  
call fruit_summary
pause

end program 