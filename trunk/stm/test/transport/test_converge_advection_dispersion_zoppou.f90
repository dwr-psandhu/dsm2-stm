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

!> Testing advection and diffusion with spacially variable coefficents vs. analytical solution
!>@ingroup test

module test_zoppou_advection_dispersion
use stm_precision
!----- module variables
! todo: make the names more meaningful
! NOTE: the parameters here should not change, the have been chosen in a range to be 
! meaningful 
! the problem here was with CFL larger than one r
integer, parameter  :: nconc = 2                              !< Number of constituents
integer, parameter  :: nstep_base = 512                  !< Number of time steps in finer discritization
integer, parameter  :: nx_base    = 128                       !< Number of spatial discritization in finer mesh 
real(stm_real),parameter :: origin = zero                     !< origin
real(stm_real),parameter :: x0 = 10000.0d0                    !< left hand side of the channel
real(stm_real),parameter :: xl = 15000.0d0                    !< right hand side of the channel
real(stm_real),parameter :: start_time = 2000.0d0             !< starts at 100000 sec (second)
real(stm_real),parameter :: end_time = 8000.0d0              !< ends at 190000 (second)
real(stm_real),parameter :: a0 = 1.0d7                        !< constant of area A=A0*(x^-1)
real(stm_real),parameter :: c0 = sixteen                       !< constant concentration
real(stm_real),parameter :: d0 = 5.0d-7                       !< constant of dispersion coefficent D=D0*(x^2)
real(stm_real),parameter :: u0 = 1.0d-4                       !< constant of velocity U=u0*x

contains

!> Tests the convergence of error rate in advection and dispersion of 
!> mass with spacially varing velocity, dispersion coefficent, and area 
subroutine test_advection_diffusion_zoppou(verbose)

use hydro_data
use boundary_advection
use boundary_diffusion
use dispersion_coefficient
use source_sink
use test_convergence_transport
use test_convergence_transport_uniform
use single_channel_boundary
use dispersion_coefficient

implicit none
procedure(hydro_data_if),pointer :: zoppou_hydro          !< The pointer points to the test's flow data
logical :: verbose
logical :: detail_printout=.true.
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)           !< reference solution at finest resolution
real(stm_real) :: test_domain_length
real(stm_real) :: total_time
character(LEN=64) :: label 
procedure(boundary_advective_flux_if),pointer :: bc_advect_flux => null()
procedure(boundary_diffusive_flux_if),pointer :: bc_diff_flux => null()
procedure(boundary_diffusive_matrix_if),pointer :: bc_diff_matrix => null()
  
! this flow generator is mass conservative
! todo: use test_convergence_transport_uniform as a model. You will be using dirichlet
!       (probably). Use code similar to the stuff around line 204. You will be using the
!       existing single_channel boundary conditions, but providing data that is appropriate
!       for zoppou. This is needed for both advection and dispersion. You will also need 
!       to use the proper API for setting dispersion.

zoppou_hydro => zoppou_flow 
compute_source => no_source
dispersion_coef => zoppou_disp_coef


label = 'advection_dispersion_zoppou' 
test_domain_length = xl - x0
total_time = end_time - start_time

!> load the initial values and reference final values to feed the test routine
call initial_fine_solution_zoppou(fine_initial_condition, &
                                  fine_solution,          &
                                  nx_base,                &
                                  nstep_base,             &
                                  nconc)                  
                                    
                                      
call set_single_channel_boundary(dirichlet_advective_flux_lo, bc_data_zoppou, &
                                 dirichlet_advective_flux_hi, bc_data_zoppou, &
                                 dirichlet_diffusive_flux_lo, bc_data_zoppou, &
                                 dirichlet_diffusive_flux_hi, extrapolate_hi_boundary_data )
                                 ! todo: I do not get the extrapolate here

boundary_diffusion_flux => single_channel_boundary_diffusive_flux
boundary_diffusion_matrix => single_channel_boundary_diffusive_matrix

!> The general subroutine which gets the fine initial and reference values from the privious subroutine and 
!> compute the norms, after each step coarsen the values and repeat computation.
!> at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      zoppou_hydro ,          &
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

subroutine zoppou_solution(value_zoppou, &
                           xpos,         &
                           time)                
                                    
use stm_precision                                       
implicit none

real(stm_real),intent(out):: value_zoppou           !< Dirichlet initial condition at left side of channel
real(stm_real),intent(in) :: xpos                   !< Location where data is requested
real(stm_real),intent(in) :: time                   !< Time

!----local
real(stm_real):: c_term1
real(stm_real):: c_term2

c_term1 =  (x0/xpos)*erfc((log(xpos/x0)-u0*time)/(two*sqrt(d0*time)))
c_term2 =  erfc((log(xpos/x0)+u0*time)/(two*sqrt(d0*time)))*exp(u0*log(xpos/x0)/d0)

value_zoppou = (c0*half)*(c_term1 + c_term2)

return
end subroutine



!-------------------------------------------
!> Generates a fine initial and final solution of analytical mass distribution 
subroutine initial_fine_solution_zoppou(fine_initial_condition, &
                                        fine_solution,          &
                                        nx_base,                &
                                        nstep_base,             &
                                        nconc)
                                       


implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base
integer,intent(in) :: nstep_base 
real(stm_real),intent(out):: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real),intent(out):: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
!----local
real(stm_real):: dx
real(stm_real):: xpos
real(stm_real):: test_domain_length
real(stm_real):: point_value
integer :: icell

dx = (xl - x0)/dble(nx_base)
do icell=1,nx_base
  xpos    = x0 +(dble(icell)-half)*dx
  call zoppou_solution(point_value,xpos,start_time)
  fine_initial_condition(icell,:) = point_value
  call zoppou_solution(point_value,xpos,end_time)
  fine_solution(icell,:) = point_value
end do

return
end subroutine
!///////////////////////////////////
!> zoppou flow and area in the finite volume form
subroutine zoppou_flow(flow,    &
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
integer, intent(in) :: ncell                   !< number of cells
real(stm_real), intent(in) :: time            !< time of request
real(stm_real), intent(in) :: dx              !< spatial step 
real(stm_real), intent(in) :: dt              !< time step 
real(stm_real), intent(out):: flow(ncell)     !< cell centered flow
real(stm_real), intent(out):: flow_lo(ncell)  !< lo face flow
real(stm_real), intent(out):: flow_hi(ncell)  !< hi face flow
real(stm_real), intent(out):: area(ncell)     !< cell center area
real(stm_real), intent(out):: area_lo(ncell)  !< area lo face
real(stm_real), intent(out):: area_hi(ncell)  !< area hi face

!--- local
real(stm_real) :: xpos_lo
real(stm_real) :: xpos_hi
real(stm_real) :: xpos
integer :: icell

do icell = 1,ncell  
  xpos_lo = x0 + dble(icell-1)*dx
  xpos_hi = x0 + dble(icell)*dx
  xpos    = x0 +(dble(icell)-half)*dx
  area(icell)    = (a0/dx)*(log(xpos_hi)-log(xpos_lo)) 
  area_lo(icell) = a0/xpos_lo
  area_hi(icell) = a0/xpos_hi
end do
 
  flow(:)    = u0*a0   
  flow_lo(:) = u0*a0
  flow_hi(:) = u0*a0
  
return
end subroutine

subroutine zoppou_disp_coef(disp_coef_lo,         &
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
    real(stm_real),intent(in) :: flow_lo(ncell)          !< flow on lo side of cells centered in time
    real(stm_real),intent(in) :: flow_hi(ncell)          !< flow on hi side of cells centered in time       
    real(stm_real),intent(in) :: flow(ncell)             !< flow on center of cells 
    !--
    integer :: ivar
    integer :: icell
    real(stm_real) :: xpos
    real(stm_real) :: xpos_lo
    real(stm_real) :: xpos_hi
        
    do icell = 1,ncell
      xpos_lo = x0 + dble(icell-1)*dx
      xpos_hi = x0 + dble(icell  )*dx
      disp_coef_lo(icell) = d0*xpos_lo**two 
      disp_coef_hi(icell) = d0*xpos_hi**two 
    end do
                
     return
 end subroutine




subroutine bc_data_zoppou(bc_value_zoppou, &
                           xloc,           &
                           conc,           &
                           nx_base,        &
                           nconc,          &
                           origin,         &
                           time,           &
                           dx,             &
                           dt)                
                                    
use stm_precision                                       
implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 
real(stm_real),intent(out):: bc_value_zoppou(nconc)     !< Dirichlet initial condition at left side of channel
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

xpos = xloc + x0  ! value comes in relative to zero origin right now
call zoppou_solution(point_value,xpos,time)
bc_value_zoppou(:) = point_value

return
end subroutine

end module