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
!    along with DSM2. If not, see <http://www.gnu.org/licenses>.
!</license>

!> Main program unit for testing teransport
!>@ingroup test
program test_transport_driver

use fruit

use test_extrapolate
use test_prim_cons_conversion
use test_uniform_flow
use test_matrix_solver
use test_boundary_diffusion
use test_diffusive_flux
use test_explicit_diffusion_operator
use test_interior_coef_matrix
use test_construct_r_h_s
use source_sink

use test_hydro
use test_advection_reaction_tidal
use test_coarsening
use test_uniform_flow
use test_prim_increment_to_cons
use test_gradient
use test_diffusion_fletcher
use test_diffusion_nonlinear_decay

!&&&&&&&&&&&&&&&&&&&
use test_convergence_transport_uniform
!&&&&&&&&&&&&&&&&&&
use test_zoppou_advection_dispersion
use test_time_dependent_advection_dispersion
use test_mms_advection_dispersion

implicit none

logical :: verbose = .true.

call init_fruit

! todo: we have 6 pointers, [2 diff + 1 adv + 1 source + 1 hydro + 1 disp_coef]
! can we right nullify(the six pointer) after each test?  

!//////// Advection unit tests
call test_gradient_calc
call test_limiter
call test_prim_cons_convert
call test_prim_increment2cons
call test_extrapolation
call test_tidal_hydro

!/// Advection-diffusion-reaction convergence in uniform flow,
!    operators are layered in successively
call test_converge_transport_uniform(verbose)

!///////// Advection convergence

call test_bidirectional_advection_convergence(verbose)

!/////// Diffusion unit tests
call test_tridi_solver
call test_boundary_diffusion_flux
call test_make_dif_flux_sub
call test_explicit_interior_diffusion_op
call test_interior_coef_matrix_sub
call test_construct_elemnts_rhs 
call test_coarsen
call test_detect_wiggle
call test_mass_comparison

!////// Diffusion convergence
!call test_diffusion_convergence_fletcher(verbose)

!///// Diffusion-reaction convergence
!call test_diffusion_cubic_decay(verbose)

! Advection - reaction problems
! todo: need to set an automatic check for hitting the boundary with coarse meshes
!       this frequently causes problems that are undetected without scrutiny
!call test_tidal_advection_reaction(verbose)

!/////Advection-Diffusion tests
!call test_zoppou_flow()    ! unit test that goes with convergence test
call test_advection_diffusion_zoppou(verbose)
call test_advection_diffusion_t_dependent(verbose)
!call test_advection_diffusion_time_dependent(verbose)

!/// Advection-diffusion-reaction
!call test_advection_diffusion_mms(verbose)

call fruit_summary

pause
end program 