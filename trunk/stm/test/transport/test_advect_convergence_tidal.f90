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
!!> Test of transport in tidal flow
!!>@ingroup test
!module test_uniform_flow_convergence
!
!contains
!
!!> Subroutine that runs a tidal advective simulation
!subroutine test_uniform_flow_advection_convergence()
!
!use stm_precision
!use state_variables
!use primitive_variable_conversion
!use advection
!use example_initial_conditions
!use example_hydro_data
!use example_sources
!use error_metric
!use fruit
!use logging
!
!implicit none
!
!!--- Problem variables
!
!integer, parameter  :: nstep_base = 40 
!integer, parameter  :: nx_base = 512
!real(stm_real), parameter :: cfl = 0.8 
!
!integer :: icoarse = 0
!integer :: nstep
!integer  :: nx
!
!integer, parameter  :: nconc = 2
!real(stm_real), parameter :: origin = zero   ! meters
!real(stm_real), parameter :: domain_length = 51200
!real(stm_real) :: dt              ! seconds
!real(stm_real) :: dx              ! meters
!real(stm_real), parameter :: ic_center      = three*fourth*domain_length
!real(stm_real), parameter :: ic_gaussian_sd = domain_length/sixteen
!real(stm_real), parameter :: constant_area = 1.D2
!real(stm_real) :: vel
!real(stm_real) :: time
!integer :: itime = 0
!integer :: icell ! debug only -- remove later
!!------
!integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
!integer :: coarsening
!integer, parameter :: nrefine = 3
!real(stm_real),allocatable :: reference(:)
!real(stm_real) norm_error(3,nrefine)
!character(LEN=64) filename
!
!! coarsening factor in convergence test
!do icoarse = 1,nrefine
!    coarsening = coarsen_factor**icoarse
!    nx = nx_base/(coarsening)
!    nstep = nstep_base/(coarsening)
!    call allocate_state(nx,nconc)
!    area = constant_area
!    area_prev = area
!!    area_lo = area     ! todo: used?
!!    area_hi = area     ! todo: used? remove from advect?
!
!    
!    vel = constant_area/constant_flow
!    dx = domain_length/dble(nx)
!    dt = cfl*dx/vel
!
!    call fill_gaussian(conc(:,1),nx,origin,dx,three*fourth*domain_length,ic_gaussian_sd)
!    call fill_gaussian(conc(:,2),nx,origin,dx,one*fourth*domain_length,ic_gaussian_sd)
!    call prim2cons( mass_prev,conc,area,nx,nconc)
!    mass = mass_prev
!    allocate(reference(ncell))  ! reference copy of initial state
!    reference = conc(:,2)
!
!    time = zero
!    ! forwards
!    do itime = 1,nstep
!       time = time + dt
!       
!!     flow = constant_flow 
!!     flow_hi = flow
!!     flow_lo = flow
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
!                  ncell,    &
!                  nvar,     &
!                  time,     &
!                  dt,       &
!                  dx)
!
!      mass_prev = mass
!      call cons2prim(conc,mass,area,nx,nconc) 
!
!    end do
!
!    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_mid_", ncell 
!    !call printout(conc(:,2),filename)
!
!
!    ! back
!    flow = -flow
!    flow_hi=-flow_hi
!    flow_lo=-flow_lo
!    do itime = 1, nstep
!      time = time + dt
!
!      call advect(mass,     &
!                  mass_prev,&  
!                  flow,     &
!                  flow_lo,  &
!                  flow_hi,  &
!                  area,     &
!                  area_prev,&
!                  area_lo,  &
!                  area_hi,  &
!                  ncell,    &
!                  nvar,     &
!                  time,     &
!                  dt,       &
!                  dx)
!
!      mass_prev = mass
!      call cons2prim(conc,mass,area,nx,nconc) 
!    end do
!    
!    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_start_", ncell 
!    !call printout(reference,filename)
!    !write(filename, "(a\i3\'.txt')"), "uniform_gaussian_end_", ncell 
!    !call printout(conc(:,2),filename)
!    call error_norm(norm_error(1,icoarse), &
!                    norm_error(2,icoarse), &
!                    norm_error(3,icoarse), &
!                    conc(:,2),reference,ncell,dx)
!
!    deallocate(reference)
!    call deallocate_state
!end do
!
!! todo: four?
!call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergemce on uniform flow")
!call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergemce on uniform flow")
!! This is known not to pass for second order convergence
!
!! todo: why 2.0d5 why not 4 ? Because the L-inf norm often fails locally with a limiter,
!! causing lower L-inf convergence than L-1 or L-2. If you get rid of the limiter, this 
!! would be a 4
!call assert_true(norm_error(3,2)/norm_error(3,1) > 2.D5,"L-inf second order convergemce on uniform flow")
!
!return
!end subroutine
!
!end module


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
module test_tidal_flow_convergence



contains

!> Test the extrapolation part of the predictor step of the advection algorithm
!> There are three cells and two constituents to test symmetry in space across constituents
subroutine test_tidal_flow_advection_convergence()

use stm_precision
use advection
use example_initial_conditions
use fruit
use error_metric

implicit none 

  integer :: ncell     !interior and two ends
  integer,parameter :: nvar = 2
   integer,parameter :: num_fine = 3
 
  real(stm_real),allocatable :: mass(:,:)
  real(stm_real),allocatable :: mass_prev(:,:)
  !todo: do we need these?
!  real(stm_real) :: conc(ncell,nvar)
!  real(stm_real) :: conc_hi(ncell,nvar)
!  real(stm_real) :: conc_lo(ncell,nvar)
  !++++++++++++++++++++++++++++++++
  !todo: is it needed?
!  real(stm_real) :: source(ncell,nvar)
  real(stm_real),allocatable :: flow(:)
  real(stm_real),allocatable :: flow_hi(:)
  real(stm_real),allocatable :: flow_lo(:)
  real(stm_real),allocatable :: area_hi(:)
  real(stm_real),allocatable :: area_lo(:)
  real(stm_real),allocatable :: area(:)
  real(stm_real),allocatable :: area_prev(:)  
  real(stm_real) :: dx
  real(stm_real) :: dt
  real(stm_real) :: time
  integer        :: ivar,jvar,kvar,pvar
  
  !----- local
  real(stm_real),parameter :: origin = zero 
  real(stm_real),parameter :: amp = half
  real(stm_real),parameter :: gravity = 9.80d0
  real(stm_real),parameter :: big_l = 409600.0d0
  real(stm_real),parameter :: depth =16.0d0
  real(stm_real),parameter :: omega= 0.506708d0 ! hr
  real(stm_real),parameter :: big_a = 0.4188704167062d0
  real(stm_real),parameter :: big_b = 0.040465522644d0
  real(stm_real),parameter :: k_0 = 10.066637844459d0
  !real(stm_real):: xpos (ncell)
  real(stm_real),parameter :: start_time = zero 
  real(stm_real),parameter :: end_time = 124d0 ! 12.4 is one exact M2 cycle of tide in hours
  real(stm_real):: sd 
  real(stm_real):: norm_1,norm_2,norm_inf
  real(stm_real):: mean
  real(stm_real):: scale
  real(stm_real):: norms (3,num_fine,nvar)
  integer  :: ntime 
   
  ! initail mass at t=0
  
  
  time = zero
  dt = 4.0d0  !hr
  ntime = 12*2
  area(:)= two
  area_lo(:)= two
  area_hi(:)= two
  scale = one
  sd = big_l/(two*two*two)
  mean = big_l/two
  ncell = 1024*4
! time?

do pvar=1,num_fine

ncell = ncell/2
ntime = ntime/2

dx = big_l/dble(ncell)
dt = 40.0d0 /dble(ntime) 

allocate(area(ncell), area_prev(ncell), area_lo(ncell), area_hi(ncell))
allocate (flow(ncell),flow_lo(ncell),flow_hi(ncell))
allocate ( mass(ncell,nvar),mass_prev(ncell,nvar))

do jvar=1,nvar
 
 time=start_time
  
call fill_gaussian(mass_prev(:,jvar),ncell,origin,dx,mean,sd,scale)
  
 
do ivar=1,ntime

     time = time + dt 
      
                do kvar = 1,ncell
                   flow_lo =area_lo(kvar)*big_a*sin(big_b*(big_l- (dble(kvar)*dx-dx)))*sin(omega*time)
                   flow_hi = area_hi(kvar)*big_a*sin(big_b*(big_l- (dble(kvar)*dx)))*sin(omega*time)
                   flow = area(kvar)*big_a*sin(big_b*(big_l- (dble(kvar)*dx-dx/two)))*sin(omega*time)
                end do
      
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
      
      end do 

 
call fill_gaussian(mass,ncell,origin,dx,mean,sd,scale) 

call error_norm(norm_1,norm_2,norm_inf,mass_prev,mass,ncell,dx)

norms(1,pvar,jvar) = norm_1
norms(2,pvar,jvar) = norm_2
norms(3,pvar,jvar) = norm_inf

end do ! on nvar

deallocate( mass,mass_prev)
deallocate(area, area_prev, area_lo, area_hi)
deallocate(flow, flow_lo, flow_hi)
 
end do ! on pvar
print *, '======= var 1========='
 print *, norms(1,1,1)/norms(1,2,1) , norms(1,2,1)/norms(1,3,1),'L-1'
 print *, norms(2,1,1)/norms(2,2,1) , norms(2,2,1)/norms(2,3,1),'L-2'
 print *, norms(3,1,1)/norms(3,2,1) , norms(3,2,1)/norms(3,3,1) ,'L-inf'
print *, '======= var 2========='
 print *, norms(1,1,2)/norms(1,2,2) , norms(1,2,2)/norms(1,3,2),'L-1'
 print *, norms(2,1,2)/norms(2,2,2) , norms(2,2,2)/norms(2,3,2),'L-2'
 print *, norms(3,1,2)/norms(3,2,2) , norms(3,2,2)/norms(3,3,2) ,'L-inf'
 
 ! todo change the numbers 3.48
call assert_true(norms(1,1,2)/norms(1,2,2) > 3.48d0,"L-1 second order convergemce on tidal advection")
call assert_true(norms(2,1,2)/norms(2,2,2) > 3.48d0,"L-2 second order convergemce on tidal advection")
call assert_true(norms(3,1,2)/norms(3,2,2) > 3.48d0,"L-inf second order convergemce on tidal advection")

call assert_true(norms(1,2,2)/norms(1,3,2) > 3.48d0,"L-1 second order convergemce on tidal advection")
call assert_true(norms(2,2,2)/norms(2,3,2) > 3.48d0,"L-2 second order convergemce on tidal advection")
call assert_true(norms(3,2,2)/norms(3,3,2) > 3.48d0,"L-inf second order convergemce on tidal advection")

 
 
return
end subroutine

end module



