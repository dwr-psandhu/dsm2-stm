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

!> Testing advection and diffusion with temporal variable coefficents vs. analytical solution
!>@ingroup test
module test_time_dependent_advection_dispersion
use stm_precision
!----- module variables

! the problem here was with CFL larger than one 
integer, parameter  :: nconc = 2                      !< Number of constituents
integer, parameter  :: nstep_base = 256               !< Number of time steps in finer discritization
integer, parameter  :: nx_base    = 128               !< Number of spatial discritization in finer mesh 
real(stm_real),parameter :: origin = zero             !< Origin
real(stm_real),parameter :: x0 = 1.0d0                !< Location of the initial condition discontinuity
real(stm_real),parameter :: x_left = 1.0d0            !< Left hand side of the channel
real(stm_real),parameter :: x_right = 2.0d0           !< Right hand side of the channel
real(stm_real),parameter :: start_time = zero         !< Starts at zero sec (second)
real(stm_real),parameter :: end_time = two            !< Ends at two (second)
real(stm_real),parameter :: a0 = 100.d0               !< Constant of area A0
real(stm_real),parameter :: c0 = sixteen              !< Constant concentration
real(stm_real),parameter :: d0 = half                 !< Constant of dispersion coefficent 
real(stm_real),parameter :: u0 = one                  !< Constant of velocity 

contains

!> Tests the convergence of error rate in advection and dispersion of 
!> mass with temporally varing velocity, dispersion coefficent, and area 
subroutine test_advection_diffusion_time_dependent(verbose)

use hydro_data
use boundary_advection
use boundary_diffusion
use error_handling
use dispersion_coefficient
use source_sink
use test_convergence_transport
use test_convergence_transport_uniform
use single_channel_boundary
use dispersion_coefficient

implicit none
procedure(hydro_data_if),    pointer :: time_dependent_hydro => null() !< The pointer points to the test's flow data
procedure(diffusion_coef_if),pointer :: dispersion_coef_t    => null() !< The pointer points to the test's dispersion coefficient data   
logical :: verbose                                            !< The flag for showing the details on the screen
logical :: detail_printout=.true.                             !< The flag for printing out the details
real(stm_real) :: fine_initial_condition(nx_base,nconc)       !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)                !< reference solution at finest resolution
real(stm_real) :: test_domain_length                          !< Domain length
real(stm_real) :: total_time                                  !< Total time of testing  
character(LEN=64) :: label                                    !< Test's name label
real(stm_real) :: cfl_number                                  !< Courant number
real(stm_real) :: point_value                                 !< Point value of the analytical solution or solution on boundary

procedure(boundary_advective_flux_if),  pointer :: bc_advect_flux => null() !< Pointer for boundary advective flux to be filled by driver
procedure(boundary_diffusive_flux_if),  pointer :: bc_diff_flux   => null() !< Pointer for boundary diffusive flux to be filled by driver
procedure(boundary_diffusive_matrix_if),pointer :: bc_diff_matrix => null() !< Pointer for boundary diffusin matrix to be filled by driver
 
! this flow generator is mass conservative
! todo: use test_convergence_transport_uniform as a model. You will be using dirichlet
!       (probably). Use code similar to the stuff around line 204. You will be using the
!       existing single_channel boundary conditions, but providing data that is appropriate
!       for zoppou. This is needed for both advection and dispersion. You will also need 
!       to use the proper API for setting dispersion.

time_dependent_hydro => time_dependent_flow 
compute_source => no_source
dispersion_coef_t => time_dependent_disp_coef

label = 'advection_dispersion_time_dependent' 
test_domain_length = x_right - x_left
total_time = end_time - start_time

cfl_number = u0*x_right*total_time*nx_base/nstep_base/test_domain_length

if (cfl_number > one) then
   call stm_fatal('Courant Number Larger Than One, Time Dependent Test!') 
end if

!> load the initial values and reference final values to feed the test routine
call initial_fine_solution_time_dependent(fine_initial_condition, &
                                          fine_solution,          &
                                          nx_base,                &
                                          nstep_base,             &
                                          nconc)                  
                                    
                                      
call set_single_channel_boundary(dirichlet_advective_flux_lo, bc_data_time_dependent, &
                                 dirichlet_advective_flux_hi, bc_data_time_dependent, &
                                 dirichlet_diffusive_flux_lo, bc_data_time_dependent, &
                                 dirichlet_diffusive_flux_hi, extrapolate_hi_boundary_data )

boundary_diffusion_flux => single_channel_boundary_diffusive_flux
boundary_diffusion_matrix => single_channel_boundary_diffusive_matrix

! todo: doxygen comment remove
!> The general subroutine which gets the fine initial and reference values from the privious subroutine and 
!> compute the norms, after each step coarsen the values and repeat computation.
!> at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      time_dependent_hydro,   &
                      single_channel_boundary_advective_flux,   &
                      boundary_diffusion_flux,&
                      boundary_diffusion_matrix,&
                      no_source,              &
                      test_domain_length,     &
                      total_time,             &
                      start_time,             &
                      fine_initial_condition, &
                      fine_solution,          &            
                      nstep_base,             &
                      nx_base,                &
                      nconc,                  &
                      verbose,                &
                      detail_printout=.true.)
                      
return                      
end subroutine

subroutine time_dependent_solution(value_time_dependent, &
                                   xpos,                 &
                                   time)                
                                    
use stm_precision                                       
implicit none

real(stm_real),intent(out):: value_time_dependent !< Dirichlet initial condition at left side of channel
real(stm_real),intent(in) :: xpos                 !< Location where data is requested
real(stm_real),intent(in) :: time                 !< Time

!----local
real(stm_real):: c_term1
real(stm_real):: c_term2

c_term1 =  erfc((xpos-x0-two*u0*(time+ dsin(pi*time/two)/pi))/(two*dsqrt(d0*(two*time+two*dsin(pi*time/two)/pi))))
c_term2 =  erfc((xpos-x0+two*u0*(time+ dsin(pi*time/two)/pi))/(two*dsqrt(d0*(two*time+two*dsin(pi*time/two)/pi))))*dexp(u0*(xpos-x0)/d0)

value_time_dependent = (c0*half)*(c_term1 + c_term2)

return
end subroutine

!-------------------------------------------
!> Generates a fine initial and final solution of analytical mass distribution 
!> The cell averaging is done by the Composite Simpson's rule 
!> int (f,a,b) = 1/12 *(Fa+ 4*F2 + 2*F3 + 4*F4 + Fb)
subroutine initial_fine_solution_time_dependent(fine_initial_condition, &
                                                fine_solution,          &
                                                nx_base,                &
                                                nstep_base,             &
                                                nconc)
                                       
implicit none

integer,intent(in) :: nconc                                        !< Number of variables 
integer,intent(in) :: nx_base                                      !< Number of cells at finest grid
integer,intent(in) :: nstep_base                                   !< Number of time steps at finest grid
real(stm_real),intent(out):: fine_initial_condition(nx_base,nconc) !< Initial condition at finest resolution
real(stm_real),intent(out):: fine_solution(nx_base,nconc)          !< Reference solution at finest resolution
!----local
real(stm_real):: dx
real(stm_real):: dxby2
real(stm_real):: xpos
real(stm_real):: test_domain_length
real(stm_real):: point_value
integer :: icell

dx = (x_right - x_left)/dble(nx_base)
fine_solution = zero
fine_initial_condition = zero

do icell=1,nx_base
  ! x = x0
  xpos    = x_left +(dble(icell)-one)*dx
  call time_dependent_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = fine_initial_condition(icell,:) + (half/six)*point_value
  call time_dependent_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = fine_solution(icell,:) + (half/six)*point_value

  ! x = x0 + 1/4L
  xpos    = x_left +(dble(icell)- three/four)*dx
  call time_dependent_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = fine_initial_condition(icell,:) + (two/six)*point_value
  call time_dependent_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = fine_solution(icell,:) + (two/six)*point_value
  
  ! x = x0 + 2/4L
  xpos    = x_left +(dble(icell)- half)*dx
  call time_dependent_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = fine_initial_condition(icell,:) + (one/six)*point_value
  call time_dependent_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = fine_solution(icell,:) + (one/six)*point_value
  
  ! x = x0 + 3/4L
  xpos    = x_left +(dble(icell)- fourth)*dx
  call time_dependent_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = fine_initial_condition(icell,:) + (two/six)*point_value
  call time_dependent_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = fine_solution(icell,:) + (two/six)*point_value
  
  ! x = x0 + 4/4L = x_right
  xpos    = x_left +(dble(icell))*dx
  call time_dependent_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = fine_initial_condition(icell,:) + (half/six)*point_value
  call time_dependent_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = fine_solution(icell,:) + (half/six)*point_value
  
end do

return
end subroutine
!///////////////////////////////////
!> time dependent flow and area in the finite volume form
subroutine time_dependent_flow (flow,    &
                                flow_lo, &
                                flow_hi, &
                                area,    &
                                area_lo, &
                                area_hi, &
                                ncell,   &
                                time,    &
                                dx,      &
                                dt)
                      
implicit none
integer, intent(in) :: ncell                  !< Number of cells
real(stm_real), intent(in) :: time            !< Time of request
real(stm_real), intent(in) :: dx              !< Spatial step 
real(stm_real), intent(in) :: dt              !< Time step 
real(stm_real), intent(out):: flow(ncell)     !< Cell centered flow
real(stm_real), intent(out):: flow_lo(ncell)  !< Low face flow
real(stm_real), intent(out):: flow_hi(ncell)  !< High face flow
real(stm_real), intent(out):: area(ncell)     !< Cell center area
real(stm_real), intent(out):: area_lo(ncell)  !< Area low face
real(stm_real), intent(out):: area_hi(ncell)  !< Area high face

  area(:)    = a0
  area_lo(:) = a0
  area_hi(:) = a0

  flow(:)    = u0*(two+dcos(pi*time/two))*a0   
  flow_lo(:) = u0*(two+dcos(pi*time/two))*a0  
  flow_hi(:) = u0*(two+dcos(pi*time/two))*a0  
  
return
end subroutine

subroutine time_dependent_disp_coef(disp_coef_lo,         &
                                    disp_coef_hi,         &
                                    flow,                 &
                                    flow_lo,              &
                                    flow_hi,              &
                                    time,                 &
                                    dx,                   &
                                    dt,                   &
                                    ncell,                &
                                    nvar)  
     
     use stm_precision
         
     implicit none
      !--- args          
    real(stm_real),intent(out):: disp_coef_lo(ncell)     !< Low side constituent dispersion coef
    real(stm_real),intent(out):: disp_coef_hi(ncell)     !< High side constituent dispersion coef      
    integer,intent(in)  :: ncell                         !< Number of cells
    integer,intent(in)  :: nvar                          !< Number of variables   
    real(stm_real),intent(in) :: time                    !< Current time
    real(stm_real),intent(in) :: dx                      !< Spatial step  
    real(stm_real),intent(in) :: dt                      !< Time step 
    real(stm_real),intent(in) :: flow_lo(ncell)          !< Flow on lo side of cells centered in time
    real(stm_real),intent(in) :: flow_hi(ncell)          !< Flow on hi side of cells centered in time       
    real(stm_real),intent(in) :: flow(ncell)             !< Flow on center of cells 
    !--
       
      disp_coef_lo(:) = d0*(two+dcos(pi*time/two))
      disp_coef_hi(:) = d0*(two+dcos(pi*time/two))
                
     return
 end subroutine

subroutine bc_data_time_dependent(bc_value_t_dependent,&
                                  xloc,                &
                                  conc,                &
                                  nx_base,             &
                                  nconc,               &
                                  origin,              &
                                  time,                &
                                  dx,                  &
                                  dt)                
                                    
use stm_precision                                       
implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 
real(stm_real),intent(out):: bc_value_t_dependent(nconc)!< Dirichlet initial condition at left side of channel
real(stm_real),intent(in) :: xloc                       !< Location where data is requested
real(stm_real),intent(in) :: time                       !< Time
real(stm_real),intent(in) :: dt                         !< Time step
real(stm_real),intent(in) :: dx                         !< Spacial mesh size
real(stm_real),intent(in) :: conc(nx_base,nconc)        !< Concentration 
real(stm_real),intent(in) :: origin                     !< Space origin

!----local
real(stm_real):: c_term1
real(stm_real):: c_term2
real(stm_real):: xpos
real(stm_real):: point_value

xpos = xloc + x_left  ! value comes in relative to zero origin right now

call time_dependent_solution(point_value,xpos,time)
bc_value_t_dependent(:) = point_value

return
end subroutine

end module