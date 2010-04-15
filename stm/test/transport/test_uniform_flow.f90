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

!> Test of transport in uniform flow
!>@ingroup test
module test_uniform_flow

contains

!> Subroutine that runs a small advective simulation
subroutine test_uniform_flow_advection()

use stm_precision
use state_variables
use primitive_variable_conversion
use advection
use example_initial_conditions
use example_hydro_data
use example_sources
use error_metric
use fruit
use logging

implicit none

!--- Problem variables

integer, parameter  :: nstep_base = 40 
integer, parameter  :: nx_base = 256
real(stm_real), parameter :: cfl = 0.8 

integer :: icoarse = 0
integer :: nstep
integer  :: nx

integer, parameter  :: nconc = 2
real(stm_real), parameter :: origin = zero   ! meters
real(stm_real), parameter :: domain_length = 51200
real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real), parameter :: ic_center      = three*fourth*domain_length
real(stm_real), parameter :: ic_gaussian_sd = domain_length/sixteen
real(stm_real), parameter :: constant_flow = 1.D2
real(stm_real), parameter :: constant_area = 1.D2
real(stm_real) :: vel
real(stm_real) :: time
integer :: itime = 0
integer :: icell ! debug only -- remove later
!------
integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
integer :: coarsening
integer, parameter :: nrefine = 3
real(stm_real),allocatable :: reference(:)
real(stm_real) norm_error(3,nrefine)
character(LEN=64) filename

! coarsening factor in convergence test
do icoarse = 1,nrefine
    coarsening = coarsen_factor**icoarse
    nx = nx_base/(coarsening)
    nstep = nstep_base/(coarsening)
    call allocate_state(nx,nconc)
    area = constant_area
    area_prev = area
!    area_lo = area     ! todo: used?
!    area_hi = area     ! todo: used? remove from advect?

    flow = constant_flow 
    flow_hi = flow
    flow_lo = flow
    vel = constant_area/constant_flow
    dx = domain_length/dble(nx)
    dt = cfl*dx/vel

    call fill_gaussian(conc(:,1),nx,origin,dx,three*fourth*domain_length,ic_gaussian_sd)
    call fill_gaussian(conc(:,2),nx,origin,dx,one*fourth*domain_length,ic_gaussian_sd)
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

    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_mid_", ncell 
    !call printout(conc(:,2),filename)


    ! back
    flow = -flow
    flow_hi=-flow_hi
    flow_lo=-flow_lo
    do itime = 1, nstep
      time = time + dt

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
    
    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_start_", ncell 
    !call printout(reference,filename)
    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_end_", ncell 
    !call printout(conc(:,2),filename)
    call error_norm(norm_error(1,icoarse), &
                    norm_error(2,icoarse), &
                    norm_error(3,icoarse), &
                    conc(:,2),reference,ncell,dx)

    deallocate(reference)
    call deallocate_state
end do

! todo: four?
call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergemce on uniform flow")
call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergemce on uniform flow")
! This is known not to pass for second order convergence

! todo: why 2.0d5 why not 4 ? Because the L-inf norm often fails locally with a limiter,
! causing lower L-inf convergence than L-1 or L-2. If you get rid of the limiter, this 
! would be a 4
call assert_true(norm_error(3,2)/norm_error(3,1) > 2.D5,"L-inf second order convergemce on uniform flow")

return
end subroutine

end module