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

!> Test advection diffusion and reaction in a single channel with uniform flow
!>@ingroup test
module test_advect_diffuse_react

use stm_precision

integer, parameter  :: nstep_base = 64
integer, parameter  :: nx_base = 2048
integer, parameter  :: nconc = 2
real(stm_real), parameter :: start_time = 256.0d0 ! sec
real(stm_real), parameter :: total_time = 512.0d0 ! sec
! todo:  since the bc is set to be zero flux, total_time and other parameters should be set 
! in the way solution does not reach the edges of channel.
real(stm_real), parameter :: domain_length = 51200.0d0 ! m
real(stm_real), parameter :: origin = zero ! low side of channel
real(stm_real), parameter :: const_area = 110.0d0 ! m^2
real(stm_real), parameter :: const_disp_coef = 150.0d0 !todo: is it in a correct range? 
real(stm_real), parameter :: const_velocity = 2.9d0 ! m/s
real(stm_real), parameter :: decay_rate = 0.0005d0 
real(stm_real), parameter :: ic_center = domain_length/(four)
real(stm_real), parameter :: ic_peak = one
real(stm_real) :: end_time = start_time + total_time

contains

!> Subroutine that tests error convergence of advection diffusion reaction 
subroutine test_advect_diffuse_reaction(verbose)

use fruit
use error_metric
use advection
use diffusion
use boundary_advection_module
use boundary_diffusion
use primitive_variable_conversion
use hydro_data
use source_sink
use log_convergence
use grid_refinement
use state_variables
use logging

implicit none

procedure(hydro_data_if),pointer:: hydro_adr                    !< This pointer, points to uniform flow data
character(LEN=64) :: label                                      !< unique label for test
logical :: verbose                                              !< whether to output convergence results

real(stm_real) :: fine_initial_conc(nx_base,nconc)              !< initial condition f concentration at finest resolution
real(stm_real) :: fine_initial_mass(nx_base,nconc)              !< initial condition of mass at finest resolution
real(stm_real) :: fine_solution_conc(nx_base,nconc)             !< reference solution at finest resolution

!---local
integer, parameter :: nrefine = 3
integer, parameter :: coarsen_factor = 2          
integer :: itime
integer :: icell 
integer :: icoarse 
integer :: nstep
integer :: nx
integer :: coarsening
integer :: which_cell
logical, parameter :: limit_slope = .false.

real(stm_real), allocatable :: solution_mass(:,:)
real(stm_real), allocatable :: reference(:,:)
real(stm_real), allocatable :: x_center(:)
real(stm_real), allocatable :: disp_coef_lo (:,:)     !< Low side constituent dispersion coef. at new time
real(stm_real), allocatable :: disp_coef_hi (:,:)     !< High side constituent dispersion coef. at new time
real(stm_real), allocatable :: disp_coef_lo_prev(:,:) !< Low side constituent dispersion coef. at old time
real(stm_real) ,allocatable :: disp_coef_hi_prev(:,:) !< High side constituent dispersion coef. at old time
real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real) :: time
real(stm_real) :: norm_error(3,nrefine)
real(stm_real) :: theta = half  

boundary_diffusion_impose  => neumann_diffusion_matrix 
boundary_diffusion_flux    => neumann_no_flow_diffusive_flux 
replace_boundary_flux      => neumann_advective_flux 
hydro_adr                  => uniform_flow_adr  
compute_source             => adr_linear_decay 
!------
label = 'ADR uniform flow, const A & Ks'

call initial_final_solution(fine_initial_conc,     &
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

open (unit = 11, file= label//'_fine_solution_ic_adr.txt')

write (11,*) 'v',const_velocity, 'dx' ,domain_length/nx_base
do icell = 1,nx_base
  write (11,*) (icell-half)*domain_length/nx_base, fine_initial_conc(icell,2)
end do
write (11,*) '=============================='
write (11,*) 'Fine solution'
do icell = 1,nx_base
  write (11,*) (icell-half)*domain_length/nx_base, fine_solution_conc(icell,2)
end do

close (11)

do icoarse = 1,nrefine

    coarsening = coarsen_factor**(icoarse - 1)
    nx = nx_base/(coarsening)
    nstep = nstep_base/(coarsening)
    call allocate_state(nx,nconc)
    area = const_area
    area_prev = const_area
    area_lo_prev = const_area
    area_hi_prev = const_area
    area_lo = const_area
    area_hi = const_area

    allocate(disp_coef_lo(nx,nconc),disp_coef_hi(nx,nconc), &
             disp_coef_lo_prev(nx,nconc),disp_coef_hi_prev(nx,nconc))
    disp_coef_lo = const_disp_coef
    disp_coef_hi = const_disp_coef
    disp_coef_lo_prev = const_disp_coef
    disp_coef_hi_prev = const_disp_coef

    allocate(reference(nx,nconc))
    call coarsen(reference,fine_solution_conc,nx_base,nx,nvar)

    allocate(x_center(nx))
    
     ! discretization parameters
    dx = domain_length/dble(nx)
    dt = total_time/dble(nstep)
          
    do icell = 1,nx
        x_center(icell) = dx*(dble(icell)-half)+origin
    end do   
    
   ! todo: we need a satbility check for Advection Diffusion splitting
        
    time = zero
    call hydro_adr(flow,    &
                   flow_lo, &
                   flow_hi, &
                   area,    &
                   area_lo, &
                   area_hi, &
                   nx,      &
                   time,    &
                   dx,      &                  
                   dt)
    area_prev = area
    area_lo_prev = area_lo
    area_hi_prev = area_hi
           
    if (icoarse == 1)then
        call prim2cons(fine_initial_mass,fine_initial_conc,area,nx,nconc)
    end if
        
    call coarsen(mass,fine_initial_mass,nx_base,nx, nconc)
    mass_prev = mass
    
    do itime = 1,nstep
       time = start_time + itime*dt
       call hydro_adr(flow,    &
                      flow_lo, &
                      flow_hi, &
                      area,    &
                      area_lo, &
                      area_hi, &
                      nx,      &
                      time,    &
                      dx,      &                  
                      dt)
   
      ! call transport using linear the callback 

      call advect(mass,     &
                  mass_prev,&  
                  flow,     &
                  flow_lo,  &
                  flow_hi,  &
                  area,     &
                  area_prev,&
                  area_lo,  &
                  area_hi,  &
                  nx,       &
                  nconc,    &
                  time,     &
                  dt,       &
                  dx,       &
                  limit_slope)

      
      area_prev = area
      call cons2prim(conc,mass,area,nx,nconc) 
      conc_prev = conc
      
      call diffuse(conc,              &
                   conc_prev,         &
                   area,              &
                   area_prev,         &
                   area_lo,           &
                   area_hi,           &
                   area_lo_prev,      &
                   area_hi_prev,      &
                   disp_coef_lo,      &  
                   disp_coef_hi,      &
                   disp_coef_lo_prev, &  
                   disp_coef_hi_prev, &
                   nx,                &
                   nconc,             &
                   time,              &
                   theta,             &
                   dt,                &
                   dx)
                   
    call prim2cons(mass,conc,area,nx,nconc)
    mass_prev = mass
    
    end do 
   
     call error_norm(norm_error(1,icoarse), &
                     norm_error(2,icoarse), &
                     norm_error(3,icoarse), &
                     which_cell,            &
                     conc(:,1),             &
                     reference(:,1),        &
                     nx,                    &
                     dx)
                     
   write(label, "(a\i4\'.txt')"), "result_adr_", nx 
    call printout(conc(:,2),x_center,label)
  

    call deallocate_state
    deallocate(disp_coef_lo,disp_coef_hi, &
               disp_coef_lo_prev,disp_coef_hi_prev)
    deallocate(x_center)
    deallocate(reference)
   
end do !icoarse

call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on " // trim(label))
call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on " // trim(label))
call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on " // trim(label))

if (verbose == .true.) then
   call log_convergence_results(norm_error,nrefine,dx,dt,const_velocity,label,which_cell,nx_base)
   print *, "Grid Peclet Number : " , const_disp_coef *dt/dx/dx
   print *, "Peclet number L=dx : " , dx*const_velocity/const_disp_coef
   print *, "maximum error in : ", which_cell
   print *, "decay rate is : " , decay_rate
   print *, 'flux_limiter :' , limit_slope
end if

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
                                  
use example_initial_conditions

implicit none
integer, intent(in) :: nx_base
integer, intent(in) :: nconc
real(stm_real),intent(out) :: fine_initial_conc(nx_base,nconc)              !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution_conc(nx_base,nconc)                  !< reference solution at finest resolution
real(stm_real),intent(in)  :: ic_center
real(stm_real),intent(in)  :: ic_peak
real(stm_real),intent(in)  :: const_velocity
real(stm_real),intent(in)  :: decay_rate
real(stm_real),intent(in)  :: total_time 
real(stm_real),intent(in)  :: origin
real(stm_real),intent(in)  :: domain_length
!--local
integer :: ivar
!real(stm_real) :: origin
real(stm_real) :: final_peak
real(stm_real) :: final_center
real(stm_real) :: final_stand_dev
real(stm_real) :: dx

dx = domain_length/nx_base

final_center = ic_center  + const_velocity * total_time


call fill_gaussian(fine_initial_conc(:,1),nx_base,origin,dx,ic_center,sqrt(two*const_disp_coef*start_time),ic_peak)

call fill_gaussian(fine_solution_conc(:,1),nx_base,origin,dx,final_center,sqrt(two*const_disp_coef*end_time),ic_peak*sqrt(start_time/end_time))

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
                                     

 use  primitive_variable_conversion
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
real(stm_real) :: mass (ncell,nvar)

! source must be in primitive variable 
call prim2cons(mass,conc,area,ncell,nvar)

source = -decay_rate*mass

return
end subroutine 

end module 
