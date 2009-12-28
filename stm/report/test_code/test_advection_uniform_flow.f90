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

!> Test of advection in uniform flow
!>@ingroup test
module test_advection_uniform_flow

    contains
    !> Subroutine that runs a small advective simulation
    subroutine test_uniform_flow_advection
        use stm_precision
        use state_variables
        use primitive_variable_conversion
        use advection
        use example_initial_conditions
        use logging

        !use example_hydro_data
        !use example_sources
        !use error_metric
        !use fruit
        implicit none

        !--- Problem variables

        integer, parameter  :: nstep  = 5
        integer, parameter  :: nx = 20
        real(STM_REAL), parameter :: cfl = 0.8


        integer, parameter  :: nconc = 2
        real(STM_REAL), parameter :: origin = zero   ! meters
        real(STM_REAL), parameter :: domain_length = 100
        real(STM_REAL) :: dt              ! seconds
        real(STM_REAL) :: dx              ! meters
        real(STM_REAL), parameter :: ic_center      = three*fourth*domain_length
        real(STM_REAL), parameter :: ic_gaussian_sd = domain_length/sixteen
        real(STM_REAL), parameter :: constant_flow = 100
        real(STM_REAL), parameter :: constant_area = 100
        real(STM_REAL) :: vel
        real(STM_REAL) :: time
        integer :: itime = 0
        !integer :: icell ! debug only -- remove later
        !------
        real(STM_REAL),allocatable :: reference(:)
        character(LEN=64) :: filename

        call allocate_state(nx,nconc)
        
        area = constant_area
        area_prev = area
        area_lo = area     ! todo: used?
        area_hi = area     ! todo: used? remove from advect?

        flow = constant_flow 
        flow_hi = flow
        flow_lo = flow
        vel = constant_flow/constant_area
        dx = domain_length/dble(nx)
        dt = cfl*dx/vel

        call fill_gaussian(conc(:,1),nx,origin,dx,0.75*domain_length,ic_gaussian_sd)
        call fill_gaussian(conc(:,2),nx,origin,dx,0.25*domain_length,ic_gaussian_sd)
        
        call fill_rectangular(conc(:,1),nx,1,nx,5.0D0,2.0D0)
        call fill_rectangular(conc(:,2),nx,1,nx,5.0D0,2.0D0)
        
        call prim2cons( mass_prev,conc,area,nx,nconc)
        mass = mass_prev
        allocate(reference(ncell))  ! reference copy of initial state
        reference = conc(:,2)

        time = zero
        ! forwards
        do itime = 1,nstep
            time = time + dt
           
              ! call transport using no_source as the callback 
              call advect(mass,     &
                          mass_prev,&  
                          flow,     &
                          flow_lo,  &
                          flow_hi,  &
                          area,     &
                          area_prev,&
                          area_lo,  &
                          area_hi,  &
                          ncell,    &
                          nvar,     &
                          time,     &
                          dt,       &
                          dx)

              mass_prev = mass
              call cons2prim(conc,mass,area,nx,nconc) 

        end do

        write(filename, "(a\i3\'.txt')"), "uniform_rectangular_", ncell 
        call printout(conc(:,2),filename)




        deallocate(reference)
        call deallocate_state


      return
      end subroutine 



      end module test_advection_uniform_flow