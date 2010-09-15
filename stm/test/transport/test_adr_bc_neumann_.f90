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

!> Test advection diffusion and reaction in a single channel. 
!  BC for advection is Dirichlet, BC for inflow is Neumann (in) and extrapolation (out)
!  Todo: these BC are not really valid for operator splitting, but it is good to remember we
!  don't have any other Neumann test for diffusion.
!  todo: move to test framework.
!>@ingroup test
module test_adr_neumann

use stm_precision

integer, parameter  :: nstep_base = 512
integer, parameter  :: nx_base = 512/2
integer, parameter  :: nconc = 2
real(stm_real), parameter :: start_time = 256.0d0 ! sec
real(stm_real), parameter :: total_time =2048.0d0 ! sec
! todo:  since the bc is set to be zero flux, total_time and other parameters should be set 
! in the way solution does not reach the edges of channel.
real(stm_real), parameter :: domain_length = 8000.0d0 ! 4000m
real(stm_real), parameter :: origin = 12800.0d0 ! low side of channel
real(stm_real), parameter :: const_area = 100.0d0 ! m^2
real(stm_real), parameter :: dispersion_coef = 90.0d0 !do not play with this number 
real(stm_real), parameter :: const_velocity = 2.9d0 ! 2.9 m/s
real(stm_real), parameter :: decay_rate = zero
real(stm_real), parameter :: ic_center = origin + domain_length/(ten)  
real(stm_real), parameter :: ic_peak = one
real(stm_real) :: end_time = start_time + total_time

contains

!> Subroutine that tests error convergence of advection diffusion reaction 
subroutine test_advect_diffuse_reaction_neumann(verbose)

use fruit
use error_metric
use advection
use diffusion
use boundary_advection
use boundary_diffusion
use primitive_variable_conversion
use hydro_data
use source_sink
use log_convergence
use grid_refinement
use state_variables
use test_convergence_transport
use logging

implicit none

procedure(hydro_data_if),pointer:: hydro_adr                  !< This pointer, points to uniform flow data
character(LEN=64) :: label                                      !< unique label for test
logical :: verbose                                              !< whether to output convergence results

real(stm_real) :: fine_initial_condition(nx_base,nconc)              !< initial condition f concentration at finest resolution
real(stm_real) :: fine_initial_mass(nx_base,nconc)              !< initial condition of mass at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)             !< reference solution at finest resolution

!---local
integer, parameter :: nrefine = 3
integer, parameter :: coarsen_factor = 2          
integer :: itime
integer :: icell 
integer :: icoarse 
integer :: nstep
integer :: nx
integer :: coarsening

logical, parameter :: limit_slope = .false.

real(stm_real), allocatable :: reference(:,:)
real(stm_real), allocatable :: x_center(:)
real(stm_real), allocatable :: disp_coef_lo(:)       !< Low side constituent dispersion coef. at new time
real(stm_real), allocatable :: disp_coef_hi(:)       !< High side constituent dispersion coef. at new time
real(stm_real), allocatable :: disp_coef_lo_prev(:)   !< Low side constituent dispersion coef. at old time
real(stm_real) ,allocatable :: disp_coef_hi_prev(:)   !< High side constituent dispersion coef. at old time
real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real) :: time
real(stm_real) :: norm_error(3,nrefine)
integer :: which_cell(nrefine)
real(stm_real) :: theta = half  
boundary_diffusion_matrix => neumann_adr_diffusion_matrix 
boundary_diffusion_flux   => neumann_adr_diffusive_flux 
advection_boundary_flux   => neumann_adr_advective_flux
hydro_adr                 => uniform_flow_adr
compute_source            => adr_linear_decay
!------
label = 'uniform_adr_neumann'
const_dispersion = 90.d0


call initial_final_solution(fine_initial_condition,     &
                            fine_solution,    &
                            ic_center,             &
                            ic_peak,               &
                            const_velocity,        &
                            decay_rate,            &
                            total_time,            &
                            origin,                &
                            domain_length,         &
                            nx_base,               &
                            nconc)

call test_convergence(label,                     &
                      hydro_adr,                 &
                      neumann_adr_advective_flux,    &
                      neumann_adr_diffusive_flux,    &
                      neumann_adr_diffusion_matrix , &
                      adr_linear_decay,              &
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
!===========
!> produce fine initial condition and reference solution 
subroutine initial_final_solution(fine_initial_conc,     &
                                  fine_solution_conc,    &
                                  ic_center,             &
                                  ic_peak,               &
                                  const_velocity,        &
                                  decay_rate,            &
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
real(stm_real),intent(out) :: fine_solution_conc(nx_base,nconc)    !< reference solution at finest resolution
real(stm_real),intent(in)  :: ic_center
real(stm_real),intent(in)  :: ic_peak
real(stm_real),intent(in)  :: const_velocity
real(stm_real),intent(in)  :: decay_rate
real(stm_real),intent(in)  :: total_time 
real(stm_real),intent(in)  :: origin
real(stm_real),intent(in)  :: domain_length
!--local
integer :: ivar
real(stm_real) :: final_peak
real(stm_real) :: final_center
real(stm_real) :: final_stand_dev
real(stm_real) :: dx

dx = domain_length/nx_base

final_center = ic_center  + const_velocity * total_time


call fill_gaussian(fine_initial_conc(:,1),nx_base,origin,dx,ic_center,sqrt(two*dispersion_coef*start_time),ic_peak)

call fill_gaussian(fine_solution_conc(:,1),nx_base,origin,dx,final_center,sqrt(two*dispersion_coef*end_time),ic_peak*sqrt(start_time/end_time))

fine_initial_conc(:,2) = fine_initial_conc(:,1) 
fine_solution_conc(:,2)= fine_solution_conc(:,1)

fine_solution_conc = fine_solution_conc * exp(-decay_rate*total_time)

return
end subroutine
!===========
!>generat constant area and constant flow foreward and backward
subroutine uniform_flow_adr(flow,    &
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

area = const_area
area_lo = const_area
area_hi = const_area

flow = const_area*const_velocity
flow_hi = flow
flow_lo = flow

return
end subroutine
!============
subroutine adr_linear_decay(source, & 
                            conc,   &
                            area,   &
                            flow,   &
                            ncell,  &
                            nvar,   &
                            time)
                                    

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
integer ivar

do ivar = 1,nvar
    source(:,ivar) = -decay_rate*conc(:,ivar)
end do 
return
end subroutine 

! todo: kaveh this seems to have mistakes. You are altering the outflow boundary and this doesn't seem "neumann"
subroutine neumann_adr_advective_flux (flux_lo,    &
                                      flux_hi,    &
                                      conc_lo,    &
                                      conc_hi,    &
                                      flow_lo,    &
                                      flow_hi,    &
                                      ncell,      &
                                      nvar,       &
                                      time,       &
                                      dt,         &
                                      dx)
     
     use stm_precision
     use gaussian_init_boundary_condition
     
        implicit none
         !--- args          
        integer,intent(in)  :: ncell  !< Number of cells
        integer,intent(in)  :: nvar   !< Number of variables
        ! todo: check the intents
        real(stm_real),intent(inout):: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
        real(stm_real),intent(inout):: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
        real(stm_real),intent(in)   :: flow_lo(ncell)          !< flow on lo side of cells centered in time
        real(stm_real),intent(in)   :: flow_hi(ncell)          !< flow on hi side of cells centered in time
        real(stm_real),intent(in)   :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
        real(stm_real),intent(in)   :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
        real(stm_real),intent(in)   :: time                    !< current time
        real(stm_real),intent(in)   :: dx                      !< spatial step  
        real(stm_real),intent(in)   :: dt                      !< time step     
      !-----local
      
       real(stm_real) :: gaussian_bell_center
       real(stm_real) :: channel_start
       real(stm_real) :: conc_mid (nx_base+1,nvar)
       real(stm_real) :: local_time

     
      !local_time = time ! next time
      local_time = time - half*dt ! half_time
      !local_time = time - dt  ! old time
      
      gaussian_bell_center = ic_center + const_velocity*(local_time-start_time)
      
      call fill_gaussian(conc_mid(:,1),nx_base,origin,dx,gaussian_bell_center,sqrt(two*dispersion_coef*local_time),ic_peak*sqrt(start_time/local_time))
      conc_mid (:,2) = conc_mid (:,1)
       
      ! todo: time or time +1/2  
      conc_mid = conc_mid * exp(-decay_rate*(local_time -start_time))
      flux_lo(1,:)     = flow_lo(1)    * conc_mid(1,:)
      
      ! outflow boundary
     ! flux_lo(1,:)     = flow_hi(ncell)    * conc_mid(ncell,:)
      flux_hi(ncell,:) = flow_hi(ncell) *(two*flux_hi(ncell-1,:)/flow_hi(ncell-1)-flux_hi(ncell-2,:)/flow_hi(ncell-2))
             
      return
 end subroutine
 
 subroutine neumann_adr_diffusion_matrix(center_diag ,       &
                                         up_diag,            &     
                                         down_diag,          &
                                         right_hand_side,    & 
                                         conc,               &
                                         explicit_diffuse_op,&
                                         area,               &
                                         area_lo,            &
                                         area_hi,            &          
                                         disp_coef_lo,       &
                                         disp_coef_hi,       &
                                         theta_stm,          &
                                         ncell,              &
                                         time,               & 
                                         nvar,               & 
                                         dx,                 &
                                         dt)
 use stm_precision
 use gaussian_init_boundary_condition
 implicit none
     !--- args
                                   
integer, intent (in) :: ncell                                               !< Number of cells
integer, intent (in) :: nvar                                                !< Number of variables

real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of the right hand side
real(stm_real),intent  (in)  :: conc(ncell,nvar)
real(stm_real), intent (in)  :: explicit_diffuse_op(ncell,nvar)
real(stm_real), intent (in)  :: area (ncell)                                !< Cell centered area at new time 
real(stm_real), intent (in)  :: area_lo(ncell)                              !< Low side area at new time
real(stm_real), intent (in)  :: area_hi(ncell)                              !< High side area at new time 
real(stm_real), intent (in)  :: disp_coef_lo(ncell)                        !< Low side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: disp_coef_hi(ncell)                        !< High side constituent dispersion coef. at new time
real(stm_real), intent (in)  :: time                                        !< Current time
real(stm_real), intent (in)  :: theta_stm                                   !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real), intent (in)  :: dx                                          !< Spatial step  
real(stm_real), intent (in)  :: dt                                          !< Time step     
!---local
real(stm_real) :: dt_by_dxsq
real(stm_real) :: flux_start(nvar)
real(stm_real) :: flux_end(nvar) 
real(stm_real) :: grad_start(nvar)
real(stm_real) :: grad_end(nvar)  
real(stm_real) :: next_time   
real(stm_real) :: xend = origin + domain_length
real(stm_real) :: xstart = origin
real(stm_real) :: half_time
real(stm_real) :: gaussian_bell_center
real(stm_real) :: channel_start
real(stm_real) :: old_time
real(stm_real) :: local_time
real(stm_real) :: resolution

dt_by_dxsq=dt/(dx*dx)
 
next_time = time 
half_time = time - half*dt
old_time  = time - dt

local_time = half_time

gaussian_bell_center = ic_center + const_velocity*(local_time-start_time)

call derivative_gaussian(grad_start(1),xstart,gaussian_bell_center,sqrt(two*dispersion_coef*local_time),ic_peak*sqrt(start_time/local_time))
grad_start(2) = grad_start(1)
!todo: remove one
call derivative_gaussian(grad_end(1),xend,gaussian_bell_center,sqrt(two*dispersion_coef*local_time),ic_peak*sqrt(start_time/local_time))
grad_end(2) = grad_end(1)
grad_end = (conc(ncell-1,:)-conc(ncell-2,:))/dx

 
flux_start = -area_lo(1)*disp_coef_lo(1)*grad_start
flux_end = -area_hi(ncell)*disp_coef_hi(ncell)*grad_end 
! todo:
!flux_end = -area_hi(ncell)*disp_coef_hi(ncell,:)*(conc(ncell-1,:)-conc(ncell-2,:))/dx

!todo: kaveh are the lo and hi correct here?
center_diag(1,:)= area(1)+ theta_stm*dt_by_dxsq* area_hi(1)*disp_coef_hi(1)  
right_hand_side(1,:) = right_hand_side(1,:) &
                            + theta_stm*(dt/dx)*flux_start(:)
 
center_diag(ncell,:)= area(ncell)+ theta_stm*dt_by_dxsq* area_lo(ncell)*disp_coef_lo(ncell)
right_hand_side(ncell,:)= right_hand_side(ncell,:) &
                               - theta_stm*(dt/dx)*flux_end(:)


                                                             
     return
 end subroutine
 
 subroutine neumann_adr_diffusive_flux(diffusive_flux_lo, &
                                       diffusive_flux_hi, &
                                       conc,              &
                                       area_lo,           &
                                       area_hi,           &
                                       disp_coef_lo,      &  
                                       disp_coef_hi,      &
                                       ncell,             &
                                       nvar,              &
                                       time,              &
                                       dx,                &
                                       dt)
 use stm_precision
 use gaussian_init_boundary_condition
 
implicit none
!--- args
integer, intent(in)  :: ncell                                   !< number of cells
integer, intent(in)  :: nvar                                    !< number of variables
real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
real(stm_real), intent (in)   :: area_lo(ncell)        !< Low side area centered at old time
real(stm_real), intent (in)   :: area_hi(ncell)        !< High side area centered at old time
real(stm_real), intent (in)   :: time                           !< time
real(stm_real), intent (in)   :: conc(ncell,nvar)               !< concentration 
real(stm_real), intent (in)   :: disp_coef_lo(ncell)      !< Low side constituent dispersion coef.
real(stm_real), intent (in)   :: disp_coef_hi(ncell)      !< High side constituent dispersion coef.
real(stm_real), intent (in)   :: dt
real(stm_real), intent (in)   :: dx
        !----local
real(stm_real) :: dt_by_dxsq
real(stm_real) :: flux_start(nvar)
real(stm_real) :: flux_end(nvar)
real(stm_real) :: grad_start(nvar)
real(stm_real) :: grad_end(nvar)   
real(stm_real) :: next_time   
real(stm_real) :: xend = origin + domain_length
real(stm_real) :: xstart = origin
real(stm_real) :: half_time
real(stm_real) :: gaussian_bell_center
real(stm_real) :: old_time
real(stm_real) :: local_time

dt_by_dxsq=dt/(dx*dx)

! here time is time_prev 
next_time = time + dt
half_time = time + half*dt
old_time  = time 

local_time = half_time

gaussian_bell_center = ic_center + const_velocity*(local_time-start_time)

! derivative_gaussian(OUTPUT or df/dx,x,center or miu,sigma,scale or a)
call derivative_gaussian(grad_start(1),xstart,gaussian_bell_center,sqrt(two*dispersion_coef*local_time),ic_peak*sqrt(start_time/local_time))
grad_start(2) = grad_start(1)
!todo: remove one
call derivative_gaussian(grad_end(1),xend,gaussian_bell_center,sqrt(two*dispersion_coef*local_time),ic_peak*sqrt(start_time/local_time))
grad_end(2) = grad_end(1)
grad_end = (conc(ncell-1,:)-conc(ncell-2,:))/dx

diffusive_flux_lo(1,:) = -area_lo(1)*disp_coef_lo(1)*grad_start
diffusive_flux_hi(ncell,:) = -area_hi(ncell)*disp_coef_hi(ncell)*grad_end



 
    return
 end subroutine
 
end module 
