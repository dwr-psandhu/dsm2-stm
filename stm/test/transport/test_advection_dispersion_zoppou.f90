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

!> Testing advection of mass which is subjected to a tidal boundary
!>@ingroup test

!todo: this test fixes test_advection_covergence_tidal by having a mass
! conserving hydro interface. We need to do housekeeping and also compare them carefullly
! to make sure the fix is correct (things like dx and other parameters)


module test_zoppou_advection_dispersion
use stm_precision
!----- module variables
! todo: make the names more meaningful
integer, parameter  :: nconc = 2                              !< Number of constituents
integer, parameter  :: nstep_base = 128                       !< Number of time steps in finer discritization
integer, parameter  :: nx_base    = 256                       !< Number of spatial discritization in finer mesh 
real(stm_real),parameter :: origin = zero                     !< origin
real(stm_real),parameter :: x0 = 10000.0d0                    !< left hand side of the channel
real(stm_real),parameter :: xl = 12048.0d0                    !< right hand side of the channel
real(stm_real),parameter :: start_time = zero                 !< starts at zero (second)
real(stm_real),parameter :: end_time = 2000.0d0               !< ends at 2000 (second)
real(stm_real),parameter :: a0 = 1.0d7                        !< constant of area A=A0*(x^-1)
real(stm_real),parameter :: c0 = three                        !< constant concentration
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
use single_channel_boundary
    
implicit none
procedure(hydro_data_if),pointer :: zoppou_hydro          !< The pointer points to the test's flow data
logical :: verbose
logical :: detail_printout=.true.
real(stm_real) :: fine_initial_condition(nx_base,nconc)  !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)            !< reference solution at finest resolution
real(stm_real) :: domain_length
real(stm_real) :: total_time
character(LEN=64) :: label 
 
       ! this flow generator is mass conservative

! todo: use test_convergence_transport_uniform as a model. You will be using dirichlet
!       (probably). Use code similar to the stuff around line 204. You will be using the
!       existing single_channel boundary conditions, but providing data that is appropriate
!       for zoppou. This is needed for both advection and dispersion. You will also need 
!       to use the proper API for setting dispersion.

zoppou_hydro=> zoppou_flow 
advection_boundary_flux => zoppou_advective_flux 
boundary_diffusion_flux => zoppou_diffusive_flux
boundary_diffusion_matrix => zoppou_diffusion_matrix
compute_source => no_source
dispersion_coef => zoppou_disp_coef


label = 'advection_dispersion_zoppou' 
domain_length = xl - x0
total_time = end_time - start_time

!> load the initial values and reference final values to feed the test routine
call initial_fine_solution_zoppou(fine_initial_condition, &
                                  fine_solution,          &
                                  nx_base,                &
                                  nstep_base,             &
                                  nconc,                  &
                                  c0,                     &
                                  u0,                     &
                                  d0,                     &  
                                  a0,                     &
                                  x0,                     &
                                  xl,                     &
                                  start_time,             &
                                  end_time)
                                      
                                      
!                                      
!call set_single_channel_boundary(dirichlet_advective_flux_lo, left_bc_dirichlet_zoppou, &
!                                 dirichlet_advective_flux_hi, right_bc_dirichlet_zoppou, &
!                                 dirichlet_diffusive_flux_lo, left_bc_dirichlet_zoppou, &
!                                 dirichlet_diffusive_flux_hi, extrapolate_hi_boundary_data )

!> The general subroutine which gets the fine initial and reference values from the privious subroutine and 
!> compute the norms, after each step coarsen the values and repeat computation.
!> at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      zoppou_hydro ,          &
                      zoppou_advective_flux,  &
                      zoppou_diffusive_flux,  &
                      zoppou_diffusion_matrix,&
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
                      detail_printout=.true.)
                      
return                      
end subroutine

!-------------------------------------------
!> Generates a fine initial and final solution of analytical mass distribution 
subroutine initial_fine_solution_zoppou(fine_initial_condition, &
                                        fine_solution,          &
                                        nx_base,                &
                                        nstep_base,             &
                                        nconc,                  &
                                        c0,                     &
                                        u0,                     &
                                        d0,                     &  
                                        a0,                     &
                                        x0,                     &
                                        xl,                     &
                                        start_time,             &
                                        end_time)
                                       


implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base
integer,intent(in) :: nstep_base 
real(stm_real),intent(out):: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real),intent(out):: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real),intent(in) :: a0                                    !< constant of area A=A0*(x^-1)
real(stm_real),intent(in) :: c0                                    !< constant concentration
real(stm_real),intent(in) :: d0                                    !< constant of dispersion coefficent D=D0*(x^2)
real(stm_real),intent(in) :: u0                                    !< constant of velocity U=u0*x
real(stm_real),intent(in) :: x0                                    !< beginning of the channel
real(stm_real),intent(in) :: xl                                    !< end of the channel
real(stm_real),intent(in) :: start_time
real(stm_real),intent(in) :: end_time
!----local
real(stm_real):: dx
real(stm_real):: xpos
real(stm_real):: x_lo
real(stm_real):: x_hi
real(stm_real):: domain_length
real(stm_real):: c_term1
real(stm_real):: c_term2
real(stm_real):: xpos_hi
real(stm_real):: xpos_lo
real(stm_real):: time
integer :: icell

domain_length = xl - x0
dx = domain_length/nx_base

do icell=1,nx_base
  xpos_lo = x0 + dble(icell-1)*dx
  xpos_hi = x0 + dble(icell)*dx
  xpos    = x0 +(dble(icell)-half)*dx
  ! initial
  time = start_time
  
           ! f(x_hi)
  c_term1 = x0*exp((2d0-u0)*time)* erf((log(x0/xpos_hi)+(three*d0-u0)*time)/(two*sqrt(d0*time)))  &
            + xpos_hi*erfc(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time)))                    &
            ! - f(x_lo)
            -x0*exp((2d0-u0)*time)* erf((log(x0/xpos_lo)+(three*d0-u0)*time)/(two*sqrt(d0*time))) &
            - xpos_lo*erfc(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time)))
             ! f(x_hi) 
  c_term2 = (one/(d0-u0))*d0*x0*erf(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time)))                   &
            + d0*xpos_hi*((x0/xpos_hi)**(d0/u0))*erfc(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time))) &
            ! -f(x_lo)
            -(one/(d0-u0))*d0*x0*erf(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time)))                  &
            - d0*xpos_lo*((x0/xpos_lo)**(d0/u0))*erfc(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time))) 
 
  fine_initial_condition(icell,1) = (c0*half/dx)*(c_term1 + c_term2)
  
  ! final
  time = end_time
            ! f(x_hi)
  c_term1 = x0*exp((2d0-u0)*time)* erf((log(x0/xpos_hi)+(three*d0-u0)*time)/(two*sqrt(d0*time)))  &
            + xpos_hi*erfc(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time)))                    &
            ! - f(x_lo)
            -x0*exp((2d0-u0)*time)* erf((log(x0/xpos_lo)+(three*d0-u0)*time)/(two*sqrt(d0*time))) &
            - xpos_lo*erfc(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time)))
             ! f(x_hi) 
  c_term2 = (one/(d0-u0))*d0*x0*erf(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time)))   &
            + d0*xpos_hi*((x0/xpos_hi)**(d0/u0))*erfc(((d0-u0)*time+log(x0/xpos_hi))/(two*sqrt(d0*time))) &
            ! -f(x_lo)
            -(one/(d0-u0))*d0*x0*erf(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time)))   &
            - d0*xpos_lo*((x0/xpos_lo)**(d0/u0))*erfc(((d0-u0)*time+log(x0/xpos_lo))/(two*sqrt(d0*time))) 
 
  fine_solution(icell,1) = (c0*half/dx)*(c_term1 + c_term2)
end do

fine_initial_condition(:,2)=fine_initial_condition(:,1)
fine_solution(:,2) = fine_solution(:,1)

return
end subroutine
!///////////////////////////////////
!> tidal flow and area for a rectangular basin with periodic forcing in the finite volume form
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
real(stm_real) :: x0
real(stm_real) :: xl
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

subroutine zoppou_advective_flux(flux_lo,    &
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
     use error_handling
     implicit none
      !--- args          
     integer,intent(in)  :: ncell  !< Number of cells
     integer,intent(in)  :: nvar   !< Number of variables
     ! todo: check the intents
     real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
     real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
     real(stm_real),intent(in)    :: flow_lo(ncell)          !< flow on lo side of cells centered in time
     real(stm_real),intent(in)    :: flow_hi(ncell)          !< flow on hi side of cells centered in time
     real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
     real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
     real(stm_real),intent(in)    :: time                    !< Current time
     real(stm_real),intent(in)    :: dx                      !< Spatial step  
     real(stm_real),intent(in)    :: dt                      !< Time step    
     
     flux_lo(1,:) = zero
     flux_hi(ncell,:) = zero
     return
 end subroutine
 
 subroutine zoppou_diffusive_flux(diffusive_flux_lo, &
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
    implicit none
    !--- args
    integer, intent(in)  :: ncell                                   !< number of cells
    integer, intent(in)  :: nvar                                    !< number of variables
    real(stm_real), intent (inout):: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
    real(stm_real), intent (inout):: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
    real(stm_real), intent (in)   :: area_lo(ncell)        !< Low side area centered at time
    real(stm_real), intent (in)   :: area_hi(ncell)        !< High side area centered at time
    real(stm_real), intent (in)   :: time                           !< time
    real(stm_real), intent (in)   :: conc(ncell,nvar)               !< concentration 
    real(stm_real), intent (in)   :: disp_coef_lo(ncell)      !< Low side constituent dispersion coef.
    real(stm_real), intent (in)   :: disp_coef_hi(ncell)      !< High side constituent dispersion coef.
    real(stm_real), intent (in)   :: dt
    real(stm_real), intent (in)   :: dx
    
    diffusive_flux_lo(1,:) = zero
    diffusive_flux_hi(ncell,:) = zero
    
       
    return
 end subroutine
   
subroutine zoppou_diffusion_matrix(center_diag ,       &
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
     implicit none
         !--- args
                                       
     integer, intent (in) :: ncell                                               !< Number of cells
     integer, intent (in) :: nvar                                                !< Number of variables
     real(stm_real),intent (inout):: down_diag(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
     real(stm_real),intent (inout):: center_diag(ncell,nvar)                     !< Values of the coefficients at the diagonal in matrix
     real(stm_real),intent (inout):: up_diag(ncell,nvar)                         !< Values of the coefficients above the diagonal in matrix
     real(stm_real),intent (inout):: right_hand_side(ncell,nvar)                 !< Values of the coefficients of right hand side vector
     real(stm_real), intent (in)  :: conc(ncell,nvar)
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
     real(stm_real) :: xstart
     real(stm_real) :: xend  
     real(stm_real) :: flux_start(nvar)
     real(stm_real) :: flux_end (nvar)
    
    dt_by_dxsq = dt/(dx*dx)
    xstart = 0.1d0
    xend = one 
    ! area is new?!? 
    !todo: call flux
    !todo: which is dirchlet?
    flux_start(:) = - area_lo(1)*disp_coef_lo(1)*(two-two*pi*sin(pi*xstart/two)*exp(-disp_coef_lo(1)*pi*pi*time/four)) 
    flux_end(:) = - area_hi(ncell)*disp_coef_hi(ncell)*(two-two*pi*sin(pi*xend/two)*exp(-disp_coef_hi(ncell)*pi*pi*time/four))
         
    center_diag(1,:)= area(1)+ theta_stm*dt_by_dxsq* area_hi(1)*disp_coef_hi(1)  
    right_hand_side(1,:) = right_hand_side(1,:) &
                                + theta_stm*(dt/dx)*flux_start(:)
    
    !todo: disp_coef_lo???
    center_diag(ncell,:)= area(ncell)+ theta_stm*dt_by_dxsq* area_lo(ncell)*disp_coef_lo(1)
    right_hand_side(ncell,:)= right_hand_side(ncell,:) &
                                   - theta_stm*(dt/dx)*flux_end(:)
      
    return
end subroutine

subroutine zoppou_disp_coef(disp_coef,            &
                            disp_coef_lo,         &
                            disp_coef_hi,         &
                            flow,                 &
                            flow_lo,              &
                            flow_hi,              &
                            time,                 &
                            dx,                   &
                            dt,                   &
                            origin,               &
                            ncell,                &
                            nvar,                 &
                            const_disp_coef)  
     
     use stm_precision
         
     implicit none
      !--- args          
    integer,intent(in)  :: ncell                         !< Number of cells
    integer,intent(in)  :: nvar                          !< Number of variables   
    real(stm_real),intent(in) :: time                    !< Current time
    real(stm_real),intent(in) :: dx                      !< Spatial step  
    real(stm_real),intent(in) :: dt                      !< Time step 
    real(stm_real),intent(in) :: origin                  !< Left side of the channel
    real(stm_real),intent(in) :: flow_lo(ncell)          !< flow on lo side of cells centered in time
    real(stm_real),intent(in) :: flow_hi(ncell)          !< flow on hi side of cells centered in time       
    real(stm_real),intent(in) :: flow(ncell)             !< flow on center of cells 
    real(stm_real),intent(out):: disp_coef(ncell,nvar)   !< center constituent dispersion coef. at new time
    real(stm_real),intent(out):: disp_coef_lo(ncell,nvar)!< Low side constituent dispersion coef. at new time
    real(stm_real),intent(out):: disp_coef_hi(ncell,nvar)!< High side constituent dispersion coef. at new time
    real(stm_real),intent(in),optional :: const_disp_coef(nvar)!< Constant value of dispersion coef. 
    !--
    integer :: ivar
    integer :: icell
    real(stm_real) :: xpos
    real(stm_real) :: xpos_lo
    real(stm_real) :: xpos_hi
    
    do ivar = 1,nconc
        do icell = 1,ncell
        xpos_lo = x0 + dble(icell-1)*dx
        xpos_hi = x0 + dble(icell)*dx
        xpos    = x0 +(dble(icell)-half)*dx
        
        end do
    end do
                
    disp_coef = LARGEREAL
    disp_coef_lo = LARGEREAL
    disp_coef_hi = LARGEREAL
    
        
     return
 end subroutine


subroutine left_bc_dirichlet_zoppou(left_bc_value_zoppou, &
                                     nconc,                  &
                                     c0,                     &
                                     u0,                     &
                                     d0,                     &  
                                     a0,                     &
                                     x0,                     &
                                     xl,                     &
                                     time,                   &
                                     dt)
                                       


implicit none

integer,intent(in) :: nconc 
real(stm_real),intent(out):: left_bc_value_zoppou(nconc)           !< Dirichlet initial condition at left side of channel
real(stm_real),intent(in) :: a0                                    !< constant of area A=A0*(x^-1)
real(stm_real),intent(in) :: c0                                    !< constant concentration
real(stm_real),intent(in) :: d0                                    !< constant of dispersion coefficent D=D0*(x^2)
real(stm_real),intent(in) :: u0                                    !< constant of velocity U=u0*x
real(stm_real),intent(in) :: x0                                    !< beginning of the channel
real(stm_real),intent(in) :: xl                                    !< end of the channel
real(stm_real),intent(in) :: time
real(stm_real),intent(in) :: dt  
!----local
real(stm_real):: c_term1
real(stm_real):: c_term2

  c_term1 =  erfc((log(x0/x0)-(u0-d0)*time)/(two*sqrt(d0*time)))
  c_term2 =  erfc((log(x0/x0)+(u0-d0)*time)/(two*sqrt(d0*time)))*exp(u0*log(x0/x0)/d0)
  left_bc_value_zoppou(1) = (c0*half)*(c_term1 + c_term2)
  
  left_bc_value_zoppou(2) = left_bc_value_zoppou(1)
  
return
end subroutine

subroutine right_bc_dirichlet_zoppou(right_bc_value_zoppou, &
                                     nconc,                  &
                                     c0,                     &
                                     u0,                     &
                                     d0,                     &  
                                     a0,                     &
                                     x0,                     &
                                     xl,                     &
                                     time,                   &
                                     dt)
                                       


implicit none

integer,intent(in) :: nconc 
real(stm_real),intent(out):: right_bc_value_zoppou(nconc)          !< Dirichlet initial condition at right side of channel
real(stm_real),intent(in) :: a0                                    !< constant of area A=A0*(x^-1)
real(stm_real),intent(in) :: c0                                    !< constant concentration
real(stm_real),intent(in) :: d0                                    !< constant of dispersion coefficent D=D0*(x^2)
real(stm_real),intent(in) :: u0                                    !< constant of velocity U=u0*x
real(stm_real),intent(in) :: x0                                    !< beginning of the channel
real(stm_real),intent(in) :: xl                                    !< end of the channel
real(stm_real),intent(in) :: time
real(stm_real),intent(in) :: dt  
!----local
real(stm_real):: c_term1
real(stm_real):: c_term2

  c_term1 =  erfc((log(x0/xl)-(u0-d0)*time)/(two*sqrt(d0*time)))
  c_term2 =  erfc((log(x0/xl)+(u0-d0)*time)/(two*sqrt(d0*time)))*exp(u0*log(x0/xl)/d0)
  right_bc_value_zoppou(1) = (c0*half)*(c_term1 + c_term2)
  
  right_bc_value_zoppou(2) = right_bc_value_zoppou(1)
  
return
end subroutine


end module