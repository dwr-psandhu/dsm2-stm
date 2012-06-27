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


module test_advection_reaction_tidal
use stm_precision
!----- module variables
! todo: make the names more meaningful
integer, parameter  :: nconc = 2                              !< Number of constituents
integer, parameter  :: nstep_base = 128*2                       !< Number of time steps in finer discritization
integer, parameter  :: nx_base    = 256*2                       !< Number of spatial discritization in finer mesh 
real(stm_real),parameter :: origin = zero                     !< Left hand side of the channel
real(stm_real),parameter :: domain_length = 204800.d0/two     !< Domain Length in meter
real(stm_real),parameter :: amplitude = half                  !< Tidal amplitude in meter    
real(stm_real),parameter :: gravity = 9.80d0                  !< Gravitational acceleration in m/s^2
real(stm_real),parameter :: depth = 16.0d0                    !< Channel depth in meter
real(stm_real),parameter :: sec_per_hr = 60.d0*60.d0          !< Convert factor of hour to second 
real(stm_real),parameter :: m2_period = 12.4d0*sec_per_hr     !< M2 tidal period 
real(stm_real),parameter :: freq=two*pi/m2_period             !< Frequency of tidal oscillation
real(stm_real),parameter :: dye_length = domain_length/six    !< Length of cosine distribution of mass
real(stm_real),parameter :: dye_center = domain_length/two    !< Center of cosine distribution of mass
real(stm_real),parameter :: ic_gaussian_sd = domain_length/32.d0   !< Standard deviation of initial values 
real(stm_real),parameter :: ic_center = domain_length/two     !< Center of initial condition
real(stm_real),parameter :: total_time = m2_period            !< Total time of the test
real(stm_real),parameter :: start_time = zero                 !< Starts at zero
real(stm_real),parameter :: const_tidal_decay_rate = 1.552749d-5 !< Decay rate (set in the way it reduces half of the mass after one tidal period)
contains

!todo: this test could be much less weird-looking if the tidal excursion moved toward the domain
!      rather than into the domain.
!      the tidal range is very short, so the discretization of the plume is fairly coarse despite the apparently
!      coarse discretization
!todo: to have an approximation in mind the tide hight is about one meter in Antioch where the water depth is ?? 


!> Tests the convergence of error rate in advection of mass which is exposed a tidal boundary 
subroutine test_tidal_advection_reaction(verbose)

use hydro_data
use boundary_advection
use boundary_diffusion
use gaussian_init_boundary_condition
use source_sink
use test_convergence_transport
use diffusion
use dispersion_coefficient
  
implicit none
procedure(hydro_data_if),pointer :: tidal_hydro           !< The pointer points to tidal flow data
logical :: verbose                                        !< Falg for showing the detail on the screen
logical :: detail_printout=.true.                         !< Flag for printing out the details
real(stm_real) :: fine_initial_condition(nx_base,nconc)   !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)            !< reference solution at finest resolution
real(stm_real) :: solution_center = ic_center             !< Center of final solution 
real(stm_real) :: solution_gaussian_sd = ic_gaussian_sd   !< Standard deviation of final values
real(stm_real) :: tidal_ar_decay_rate                     !< Tidal decay rate
character(LEN=64) :: label                                !< Test name label
 
tidal_hydro=> tidal_flow_modified ! this flow generator is mass conservative
! do not remove it 
!tidal_hydro=> tidal_flow_cell_average ! this flow generator is NOT mass conservative but it is cell averaged

advection_boundary_flux => zero_advective_flux !todo: move this so it isn't hardwired
boundary_diffusion_flux => no_diffusion_flux
boundary_diffusion_matrix => no_diffusion_matrix

call set_constant_dispersion(zero)

tidal_ar_decay_rate = zero
compute_source => no_source

label = 'advection_tidal_gaussian' 

! load the initial values and reference final values to feed the test routine
call initial_fine_solution_tidal_gaussian(fine_initial_condition, &
                                          fine_solution,          &
                                          nx_base,                &
                                          nconc,                  &
                                          origin,                 &
                                          domain_length,          &
                                          ic_gaussian_sd,         &
                                          solution_gaussian_sd,    &
                                          ic_center,              &
                                          solution_center,        &        
                                          tidal_ar_decay_rate)

! The general subroutine which gets the fine initial and reference values from the privious subroutine and 
! compute the norms, after each step coarsen the values and repeat computation.
! at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      tidal_hydro ,           &
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
                      detail_printout=.true.)
                      
label = "advection_tidal_sinusoidal" 
! load the initial values and reference final values to feed the test routine
call initial_fine_solution_tidal_sinusoidal(fine_initial_condition, &
                                            fine_solution,          &
                                            nx_base,                &
                                            nconc,                  &
                                            domain_length,          &
                                            dye_center,             &
                                            dye_length,             &
                                            tidal_ar_decay_rate)

! The general subroutine which gets the fine initial and reference values from the privious subroutine and 
! compute the norms, after each step coarsen the values and repeat computation.
! at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      tidal_hydro ,           &
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
                      detail_printout=.true.)

!!!!!!!!!!!!!!!!!!!!!!
tidal_ar_decay_rate = const_tidal_decay_rate
compute_source => tidal_reaction_source


label = 'advection_reaction_tidal_gaussian' 

! load the initial values and reference final values to feed the test routine
call initial_fine_solution_tidal_gaussian(fine_initial_condition, &
                                          fine_solution,          &
                                          nx_base,                &
                                          nconc,                  &
                                          origin,                 &
                                          domain_length,          &
                                          ic_gaussian_sd,         &
                                          solution_gaussian_sd,    &
                                          ic_center,              &
                                          solution_center,        &        
                                          tidal_ar_decay_rate)

! The general subroutine which gets the fine initial and reference values from the privious subroutine and 
! compute the norms, after each step coarsen the values and repeat computation.
! at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      tidal_hydro ,           &
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
                      detail_printout=.true.)
                      
label = "advection_reaction_tidal_sinusoidal" 
! load the initial values and reference final values to feed the test routine
call initial_fine_solution_tidal_sinusoidal(fine_initial_condition, &
                                            fine_solution,          &
                                            nx_base,                &
                                            nconc,                  &
                                            domain_length,          &
                                            dye_center,             &
                                            dye_length,             &
                                            tidal_ar_decay_rate)

! The general subroutine which gets the fine initial and reference values from the privious subroutine and 
! compute the norms, after each step coarsen the values and repeat computation.
! at the end  calculates the ratio of the norms and prints a log 
call test_convergence(label,                  &
                      tidal_hydro ,           &
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
                      detail_printout=.true.)

end subroutine
!-------------------------------------------
!> Generates a fine initial and final solution of a Gaussian  mass distribution 
subroutine initial_fine_solution_tidal_gaussian(fine_initial_condition, &
                                                fine_solution,          &
                                                nx_base,                &
                                                nconc,                  &
                                                origin,                 &
                                                domain_length,          &
                                                ic_gaussian_sd,         &
                                                solution_gaussin_sd,    &
                                                ic_center,              &
                                                solution_center,        &        
                                                tidal_ar_decay_rate)

use gaussian_init_boundary_condition
use stm_precision

implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 
real(stm_real),intent(out) :: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real),intent(in)  :: ic_center                             !< Center of fine initial value
real(stm_real),intent(in)  :: solution_center                       !< Center of solution
real(stm_real),intent(in)  :: ic_gaussian_sd                        !< Standard deviation of initial value
real(stm_real),intent(in)  :: solution_gaussin_sd                   !< Standard deviation of solution 
real(stm_real),intent(in)  :: origin                                !< Left hand side of the channel
real(stm_real),intent(in)  :: domain_length                         !< Domain length
real(stm_real),intent(in)  ::  tidal_ar_decay_rate                  !< Decay rate
!----local
real(stm_real):: dx
real(stm_real):: xposition(nx_base)
integer :: icell

dx = domain_length/nx_base

!todo: remove these gaussian
call fill_gaussian(fine_initial_condition(:,1),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd,one)
call fill_gaussian(fine_initial_condition(:,2),nx_base,origin,dx, &
                   ic_center,ic_gaussian_sd,one)

call fill_gaussian(fine_solution(:,1),nx_base,origin,dx, &
                   solution_center,solution_gaussin_sd,one)
call fill_gaussian(fine_solution(:,2),nx_base,origin,dx, &
                   solution_center,solution_gaussin_sd,one)

fine_initial_condition(:,2)=fine_initial_condition(:,1)
fine_solution = exp(-total_time*tidal_ar_decay_rate)*fine_initial_condition

return
end subroutine

!-------------------------------------------
!> Generates a fine initial and final solution of a sinusoidal  mass distribution 
subroutine initial_fine_solution_tidal_sinusoidal(fine_initial_condition, &
                                                  fine_solution,          &
                                                  nx_base,                &
                                                  nconc,                  &
                                                  domain_length,          &
                                                  dye_center,             &
                                                  dye_length,             &
                                                  tidal_ar_decay_rate)
                                       


use gaussian_init_boundary_condition
use stm_precision

implicit none

integer,intent(in) :: nconc 
integer,intent(in) :: nx_base 
real(stm_real),intent(out) :: fine_initial_condition(nx_base,nconc) !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)          !< reference solution at finest resolution
real(stm_real),intent(in)  :: domain_length                         !< Domain length
real(stm_real),intent(in)  :: dye_center                            !< center of sinusoidal mass
real(stm_real),intent(in)  :: dye_length                            !< length of mass at the middle of the domain
real(stm_real),intent(in)  ::  tidal_ar_decay_rate                  !< Decay rate
!----local
real(stm_real):: dx
real(stm_real):: xposition
real(stm_real):: x_lo
real(stm_real):: x_hi
integer :: icell

dx = domain_length/nx_base

do icell=1,nx_base
  x_lo     = (dble(icell)-one)*dx
  xposition= (dble(icell)-half)*dx 
  x_hi     = (dble(icell))*dx
  if (( x_lo > (dye_center + dye_length*half)) .or. & 
      ( x_hi < (dye_center - dye_length*half))  ) then
        fine_initial_condition(icell,1) = zero
  else
        fine_initial_condition(icell,1)= one + &
            (dye_length*half/pi)*(sin((x_hi - dye_center)*two*pi/dye_length) - &
                                  sin((x_lo-dye_center)*two*pi/dye_length))/dx
  end if 
end do

fine_initial_condition(:,2)=fine_initial_condition(:,1)
fine_solution = exp(-total_time*tidal_ar_decay_rate)*fine_initial_condition

return
end subroutine
!///////////////////////////////////
!> Tidal flow and area for a rectangular basin with periodic forcing in the cell averaged form
subroutine tidal_flow_cell_average(flow,    &
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
real(stm_real), intent(in)  :: time            !< time of request
real(stm_real), intent(in)  :: dx              !< spatial step 
real(stm_real), intent(in)  :: dt              !< time step 
real(stm_real), intent(out) :: flow(ncell)     !< cell centered flow
real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow
real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow
real(stm_real), intent(out) :: area(ncell)     !< cell center area
real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face
real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face

!--- local
real(stm_real) :: half_time
real(stm_real) :: old_time
real(stm_real) :: big_b 
real(stm_real) :: big_a 
real(stm_real) :: vel_lo
real(stm_real) :: vel_hi
real(stm_real) :: vel
real(stm_real) :: xpos_lo
real(stm_real) :: xpos_hi
real(stm_real) :: xpos
real(stm_real) :: flow_term1
real(stm_real) :: flow_term2

integer :: icell

half_time = time - half*dt
old_time = time - dt
big_b = freq/sqrt(gravity*depth)
big_a = amplitude* sqrt(gravity*depth)/(depth*cos(big_b*domain_length))

! width is assumed to be equal to 1 meter 
do icell = 1,ncell  
  
  ! this is L-x
  xpos_lo = domain_length- dble(icell-1)*dx
  xpos_hi = domain_length- dble(icell)*dx
  xpos    = domain_length-(dble(icell)-half)*dx
  
 
  area(icell)    = depth + (amplitude*cos(freq*time))/(dx*big_b*cos(big_b*domain_length)) &
                    *(sin(big_b*xpos_hi) - sin(big_b*xpos_lo))
  
  area_lo(icell) = depth + amplitude * cos(big_b*xpos_lo)/cos(big_b*domain_length)*cos(freq*half_time)  
  area_hi(icell) = depth + amplitude * cos(big_b*xpos_hi)/cos(big_b*domain_length)*cos(freq*half_time)  
 

   flow(icell)    = big_a*depth*sin(freq*time)*(cos(big_b*xpos_hi)-cos(big_b*xpos_lo))/(dx*big_b) &
                   + big_a*amplitude*sin(two*freq*time)*(one/(four*dx*big_b*cos(big_b*domain_length))) * &
                   (cos(big_b*xpos_hi)**two - cos(big_b*xpos_lo)**two)

   flow_term1 = - big_a*depth*sin(big_b*xpos_lo)*(cos(freq*time)-cos(freq*old_time))/(dt*freq)
   flow_term2 = - big_a*amplitude*sin(two*big_b*xpos_lo)*(cos(freq*time)**two-cos(freq*old_time)**two)&
                   /(four*freq*cos(big_b*domain_length)*dt)
   flow_lo(icell) = flow_term1+flow_term2
          

   flow_term1 = - big_a*depth*sin(big_b*xpos_hi)*(cos(freq*time)-cos(freq*old_time))/(dt*freq)
   flow_term2 = - big_a*amplitude*sin(two*big_b*xpos_hi)*(cos(freq*time)**two-cos(freq*old_time)**two)&
                   /(four*freq*cos(big_b*domain_length)*dt)
   flow_hi(icell) = flow_term1+flow_term2
   
   
end do  
return
end subroutine


!> Tidal flow for a rectangular basin with periodic forcing in the cell centerd 
!> the area here retrived from dQ/dx+dA/dt=0 --> A =width*(depth+int(-dQ/dx,t)) 
subroutine tidal_flow_modified(flow,    &
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
real(stm_real), intent(in)  :: time            !< time of request
real(stm_real), intent(in)  :: dx              !< spatial step 
real(stm_real), intent(in)  :: dt              !< time step 
real(stm_real), intent(out) :: flow(ncell)     !< cell centered flow
real(stm_real), intent(out) :: flow_lo(ncell)  !< lo face flow
real(stm_real), intent(out) :: flow_hi(ncell)  !< hi face flow
real(stm_real), intent(out) :: area(ncell)     !< cell center area
real(stm_real), intent(out) :: area_lo(ncell)  !< area lo face
real(stm_real), intent(out) :: area_hi(ncell)  !< area hi face

!--- local
real(stm_real) :: half_time
real(stm_real) :: old_time
real(stm_real) :: big_b 
real(stm_real) :: big_a 
real(stm_real) :: vel_lo
real(stm_real) :: vel_hi
real(stm_real) :: vel
real(stm_real) :: xpos_lo
real(stm_real) :: xpos_hi
real(stm_real) :: xpos
real(stm_real) :: flow_term1
real(stm_real) :: flow_term2

integer :: icell

half_time = time - half*dt
old_time = time - dt
big_b = freq/sqrt(gravity*depth)
big_a = amplitude* sqrt(gravity*depth)/(depth*cos(big_b*domain_length))

! width is assumed to be equal to 1 meter 
do icell = 1,ncell  
  
  ! this is L-x
  xpos_lo = domain_length- dble(icell-1)*dx
  xpos_hi = domain_length- dble(icell)*dx
  xpos    = domain_length-(dble(icell)-half)*dx
  
   area(icell) = (big_a/(freq*dx))*(depth*cos(freq*time)*(sin(big_b*xpos_hi)-sin(big_b*xpos_lo))+&
                                  amplitude*cos(two*freq*time)*(sin(two*big_b*xpos_hi)-sin(two*big_b*xpos_lo))/(eight*cos(big_b*domain_length)))
 ! todo: is it a correction factor for constant in integeration?
 ! check this
   area(icell) = depth + area(icell)
   area_lo(icell) = depth + amplitude * cos(big_b*xpos_lo)/cos(big_b*domain_length)*cos(freq*half_time)  

   area_lo(icell) = big_a*big_b*(-depth*cos(big_b*xpos_lo)*cos(freq*half_time)/freq - &
                                amplitude*cos(two*freq*half_time)*cos(two*big_b*xpos_lo)&
                               /(four*freq*cos(big_b*domain_length)))
  ! todo:
   ! i just match them base on old values ???                             
  area_lo(icell) = depth + area_lo(icell)
 
  area_hi(icell) = depth + amplitude * cos(big_b*xpos_hi)/cos(big_b*domain_length)*cos(freq*half_time) 
   
  area_hi(icell) = big_a*big_b*(-depth*cos(big_b*xpos_hi)*cos(freq*half_time)/freq - &
                                amplitude*cos(two*freq*half_time)*cos(two*big_b*xpos_hi)&
                               /(four*freq*cos(big_b*domain_length)))
   
   area_hi(icell) = depth + area_hi(icell)
                         
    
   flow(icell)    = big_a*depth*sin(freq*time)*(cos(big_b*xpos_hi)-cos(big_b*xpos_lo))/(dx*big_b) &
                   + big_a*amplitude*sin(two*freq*time)*(one/(four*dx*big_b*cos(big_b*domain_length))) * &
                   (cos(big_b*xpos_hi)**two - cos(big_b*xpos_lo)**two)
  
   flow_term1 = - big_a*depth*sin(big_b*xpos_lo)*(cos(freq*time)-cos(freq*old_time))/(dt*freq)
   flow_term2 = - big_a*amplitude*sin(two*big_b*xpos_lo)*(cos(freq*time)**two-cos(freq*old_time)**two)&
                   /(four*freq*cos(big_b*domain_length)*dt)
   flow_lo(icell) = flow_term1+flow_term2
          

   flow_term1 = - big_a*depth*sin(big_b*xpos_hi)*(cos(freq*time)-cos(freq*old_time))/(dt*freq)
   flow_term2 = - big_a*amplitude*sin(two*big_b*xpos_hi)*(cos(freq*time)**two-cos(freq*old_time)**two)&
                   /(four*freq*cos(big_b*domain_length)*dt)
   flow_hi(icell) = flow_term1+flow_term2
   
   
end do  

return
end subroutine 
!> Subroutine provides source term (constant decay) to the tidal test
subroutine tidal_reaction_source(source,             & 
                                 conc,               &
                                 area,               &
                                 flow,               &
                                 ncell,              &
                                 nvar,               &
                                 time)
                               
                                     

 use  primitive_variable_conversion
 implicit none
 
 !--- args
integer,intent(in)  :: ncell                      !< Number of cells
integer,intent(in)  :: nvar                       !< Number of variables
real(stm_real),intent(inout):: source(ncell,nvar) !< cell centered source 
real(stm_real),intent(in)   :: conc(ncell,nvar)   !< Concentration
real(stm_real),intent(in)   :: area(ncell)        !< area at source     
real(stm_real),intent(in)   :: flow(ncell)        !< flow at source location
real(stm_real),intent(in)   :: time               !< time 

! source must be in primitive variable 
source = -const_tidal_decay_rate*conc
  
return
end subroutine 

end module


