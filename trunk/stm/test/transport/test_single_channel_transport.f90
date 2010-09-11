!!<license>
!!    Copyright (C) 1996, 1997, 1998, 2001, 2007, 2009 State of California,
!!    Department of Water Resources.
!!    This file is part of DSM2.
!!
!!    The Delta Simulation Model 2 (DSM2) is free software: 
!!    you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    DSM2 is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with DSM2.  If not, see <http://www.gnu.org/licenses>.
!!</license>
!
!!> Test advection in a single channel for a flow regime that
!!> goes back and forth and ends up in the same spot periodically
!!>@ingroup test
!module test_single_channel_advection
!
!contains
!
!!> Subroutine that tests transport convergence
!!> A fine grid initial condition and solution at time total_time
!!> must be provided. 
!subroutine test_transport_convergence(label,                &
!                                      hydro,                &
!                                      domain_length,        &
!                                      total_time,           &
!                                      fine_initial_conc,    &
!                                      fine_solution,        &           
!                                      nstep_base,           &
!                                      nx_base,              &
!                                      nconc,                &
!                                      verbose)
!                             
!use hydro_data
!use boundary_advection
!use stm_precision
!use state_variables
!use primitive_variable_conversion 
!use advection
!use gaussian_init_boundary_condition
!use error_metric
!use fruit
!use logging
!use grid_refinement
!use log_convergence
!use source_sink
!
!implicit none
!
!!--- Problem variables
!procedure(hydro_data_if), pointer :: hydro                      !< This pointer, points to uniform flow data
!character(LEN=*),intent(in) :: label                            !< unique label for test
!logical,intent(in) :: verbose                                   !< whether to output convergence results
!integer, intent(in) :: nconc                                    !< number of constituents
!integer, intent(in) :: nstep_base                               !< number of steps at finest resolution
!integer, intent(in) :: nx_base                                  !< number of cells at finest resolution
!real(stm_real), intent(in) :: fine_initial_conc(nx_base,nconc)  !< initial condition at finest resolution
!real(stm_real), intent(in) :: fine_solution(nx_base,nconc)      !< reference solution at finest resolution
!real(stm_real), intent(in) :: total_time                        !< total time of simulation
!real(stm_real), intent(in) :: domain_length                     !< length of domain
!
!!---local
!integer, parameter :: nrefine = 3
!integer, parameter :: coarsen_factor = 2                 ! coarsening factor used for convergence test
!integer :: itime
!integer :: icell 
!integer :: icoarse 
!integer :: nstep
!integer :: nx
!integer :: which_cell(nrefine)
!integer :: coarsening
!character(LEN=64)  ::  filename 
!logical, parameter :: limit_slope = .true.
!real(stm_real), allocatable :: solution_mass(:,:)
!real(stm_real), allocatable :: reference(:,:)
!real(stm_real), allocatable :: x_center(:)
!real(stm_real), allocatable :: velocity (:)
!real(stm_real) :: max_velocity 
!
!real(stm_real) :: fine_initial_mass(nx_base,nconc)   !< initial condition at finest resolution
!real(stm_real) :: fine_solution_mass(nx_base,nconc)  !< reference solution at finest resolution
!
!real(stm_real) :: dt              ! seconds
!real(stm_real) :: dx              ! meters
!real(stm_real) :: time
!real(stm_real) :: norm_error(3,nrefine)
!
!
!!todo: this is really "no flux"
!replace_advection_boundary_flux => neumann_advective_flux
!filename=label
!! coarsening factor in convergence test
!do icoarse = 1,nrefine
!    coarsening = coarsen_factor**(icoarse-1)
!    nx = nx_base/(coarsening)
!    nstep = nstep_base/(coarsening)
!    call allocate_state(nx,nconc)
!    allocate(x_center(nx))
!    allocate(reference(nx,nconc))
!    allocate(solution_mass(nx,nconc))
!    allocate(velocity (nx))
!    
!    dx = domain_length/dble(nx)  
!    dt = total_time/dble(nstep)
!
!    do icell = 1,nx
!        x_center(icell) = (dble(icell)-half)*dx
!    end do
!
!    time = zero
!    call hydro(flow,    &
!               flow_lo, &
!               flow_hi, &
!               area,    &
!               area_lo, &
!               area_hi, &
!               nx,      &
!               time,    &
!               dx,      &                  
!               dt)
!    area_prev = area
!            
!    if (icoarse == 1)then
!        call prim2cons(fine_initial_mass,fine_initial_conc,area,nx,nconc)
!    end if
!        
!    call coarsen(mass,fine_initial_mass,nx_base,nx, nconc)
!    mass_prev = mass
!    call cons2prim(conc,mass,area,nx,nconc)
!    conc_prev = conc
!    
!    max_velocity =zero
!    do itime = 1,nstep
!       time = time + dt
!       call hydro(flow,    &
!                  flow_lo, &
!                  flow_hi, &
!                  area,    &
!                  area_lo, &
!                  area_hi, &
!                  nx,      &
!                  time,    &
!                  dx,      &                  
!                  dt)
!                  
!      
!      if (maxval(abs(flow)/area) >=  max_velocity) then
!          max_velocity = maxval(abs(flow)/area)
!      end if
!       
!
!      ! call transport using no_source as the callback 
!      call advect(mass,     &
!                  mass_prev,&  
!                  flow,     &
!                  flow_lo,  &
!                  flow_hi,  &
!                  area,     &
!                  area_prev,&
!                  area_lo,  &
!                  area_hi,  &
!                  nx,       &
!                  nconc,     &
!                  time,     &
!                  dt,       &
!                  dx,       &
!                  limit_slope)
!
!      mass_prev = mass
!      area_prev = area
!      call cons2prim(conc,mass,area,nx,nconc) 
!    
!    if (label == 'linear decay no flow') then
!        if (minval (conc) < zero)then
!            print *,'Negative concentration !!!!!','Conc =',minval(conc)             
!        end if       
!    end if
!    
!    end do
!
!    ! Now take fine solution (provided in concentration) and coarsen it to
!    ! a reference solution at the current level of refinement. This needs to 
!    ! be done by converting it to mass, coarsening, then converting back to
!    ! a reference concentration
!    if (icoarse == 1) then
!        call prim2cons(fine_solution_mass,fine_solution,area,nx,nvar)
!    end if
!
!    call coarsen(solution_mass,fine_solution_mass,nx_base,nx, nvar)
!    call cons2prim(reference,solution_mass,area,nx,nconc)
!    
!    !todo: look up how to remove spaces in filename
!    write(filename, "(a\i4\'.txt')"), "uniform_gaussian_start_", nx        
!    call printout(reference(:,2),x_center(:),filename)
!    write(filename, "(a\i\'.txt')"), "uniform_gaussian_end_", nx 
!    call printout(conc(:,2),x_center(:),filename)
!    ! test error norm over part of domain
!    call error_norm(norm_error(1,icoarse), &
!                    norm_error(2,icoarse), &
!                    norm_error(3,icoarse), &
!                    which_cell(icoarse),   &
!                    conc(:,2),reference(:,2),nx,dx) 
!    deallocate(solution_mass)
!    deallocate(reference)
!    deallocate(x_center)
!    deallocate(velocity)
!    call deallocate_state
!end do
!
!call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on " // trim(label))
!call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on " // trim(label))
!call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on " // trim(label))
!
!call assert_true(norm_error(1,3)/norm_error(1,2) > four,"L-1 second order convergence on " // trim(label))
!call assert_true(norm_error(2,3)/norm_error(2,2) > four,"L-2 second order convergence on " // trim(label))
!call assert_true(norm_error(3,3)/norm_error(3,2) > four,"L-inf second order convergence on " // trim(label))
!
!if (verbose == .true.) then
!   call log_convergence_results(norm_error ,                   &
!                                nrefine,                       &
!                                dx,                            &
!                                dt,                            &
!                                max_velocity= max_velocity,    &
!                                label = label,                 &
!                                which_cell=which_cell,         &
!                                ncell_base = nx_base,          &
!                                ntime_base = nstep_base,       &
!                                reaction_rate = zero,          &
!                                scheme_order = two,            &
!                                length_scale = dx,             &
!                                limiter_switch = limit_slope)
!                                 
!end if
!
!return
!end subroutine
!
!end module