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

!> Test of transport diffusion convergence test for a single channel
!>@ingroup test
module test_diffusion_single_channel

contains

!> Subroutine that checks the error convergence ratio for diffusion routine 
subroutine test_diffusion_convergence_single_channel

use stm_precision
use state_variables
use primitive_variable_conversion
use boundary_diffusion
use diffusion
use example_initial_conditions
use error_metric
use error_handling
use fruit
use logging

implicit none

!--- Problem variables

integer, parameter  :: nstep_base = 128
integer, parameter  :: nx_base = 256

integer :: icoarse = 0
integer :: nstep
integer :: nx

integer, parameter  :: nconc = 2
real(stm_real), parameter :: domain_length = 51200.d0
real(stm_real), parameter :: origin = zero   
real(stm_real), parameter :: total_time    = 2048.d0
real(stm_real), parameter :: disp_coef     = 512.d0
real(stm_real) :: theta = half                       !< Explicitness coefficient; 0 is explicit, 0.5 Crank-Nicolson, 1 full implicit  
real(stm_real),allocatable :: disp_coef_lo (:,:)     !< Low side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_hi (:,:)     !< High side constituent dispersion coef. at new time
real(stm_real),allocatable :: disp_coef_lo_prev(:,:) !< Low side constituent dispersion coef. at old time
real(stm_real),allocatable :: disp_coef_hi_prev(:,:) !< High side constituent dispersion coef. at old time


real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real), parameter :: ic_center      = three*fourth*domain_length
real(stm_real), parameter :: ic_gaussian_sd = domain_length/sixteen
real(stm_real), parameter :: ic_peak = one
real(stm_real), parameter :: constant_area = 110.0d0
real(stm_real) :: vel
real(stm_real), parameter :: start_time = 256.d0
real(stm_real), parameter :: end_time = start_time + total_time

real(stm_real) :: time
real(stm_real), allocatable :: xposition(:)


integer :: itime
integer :: icell ! debug only -- remove later
!------
integer, parameter :: coarsen_factor = 2      ! coarsening factor used for convergence test
integer :: coarsening
integer, parameter :: nrefine = 3
real(stm_real),allocatable :: reference(:)
real(stm_real) norm_error(3,nrefine)
character(LEN=64) filename

boundary_diffusion_impose  => neumann_diffusion_matrix
boundary_diffusion_flux    => neumann_no_flow_diffusive_flux

! coarsening factor in convergence test
do icoarse = 1,nrefine
    coarsening = coarsen_factor**(icoarse - 1)
    nx = nx_base/(coarsening)
    nstep = nstep_base/(coarsening)
    call allocate_state(nx,nconc)
    area = constant_area
    area_prev = constant_area
    area_lo_prev = constant_area
    area_hi_prev = constant_area
    area_lo = constant_area
    area_hi = constant_area
    allocate(disp_coef_lo(ncell,nvar),disp_coef_hi(ncell,nvar), &
             disp_coef_lo_prev(ncell,nvar),disp_coef_hi_prev(ncell,nvar))
    disp_coef_lo = disp_coef
    disp_coef_hi = disp_coef
    disp_coef_lo_prev = disp_coef
    disp_coef_hi_prev = disp_coef
    
    ! discretization parameters
    dx = domain_length/dble(nx)
    dt = total_time/dble(nstep)
    
    allocate(xposition(nx))
    do icell = 1,nx
        xposition(icell) = dx*(dble(icell)-half)+origin
    end do      

    call fill_gaussian(conc(:,1),nx,origin,dx,half*domain_length, & 
                       sqrt(two*disp_coef*start_time),one)
    conc(:,2) = conc(:,1)
    
    allocate(reference(ncell))  ! reference solution
    call fill_gaussian(reference,nx,origin,dx,half*domain_length, &
                       sqrt(two*disp_coef*end_time),sqrt(start_time/end_time))

    time = zero
    ! forwards

    timemarch: do itime = 1,nstep
        conc_prev=conc
        time = start_time + itime*dt
        call diffuse(conc,             &
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
                     ncell,             &
                     nvar,              &
                     time,              &
                     theta,             &
                     dt,                &
                     dx                 ) 
    end do timemarch
    write(filename, "(a\i3\'.txt')"), "diffuse_gaussian_reference_", ncell 
!    call printout(reference,xposition,filename)
    write(filename, "(a\i3\'.txt')"), "diffuse_gaussian_solution_", ncell 
 !todo:   call printout(conc(:,2),xposition,filename)
    call error_norm(norm_error(1,icoarse), &
                    norm_error(2,icoarse), &
                    norm_error(3,icoarse), &
                    conc(:,2),reference,ncell,dx)

    deallocate(reference)
    deallocate(xposition)
    deallocate(disp_coef_lo,disp_coef_hi,disp_coef_lo_prev,disp_coef_hi_prev)    
    call deallocate_state
end do

call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on diffusion")
call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on diffusion")
call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on diffusion")

call assert_true(norm_error(1,3)/norm_error(1,2) > four,"L-1 second order convergence on diffusion")
call assert_true(norm_error(2,3)/norm_error(2,2) > four,"L-2 second order convergence on diffusion")
call assert_true(norm_error(3,3)/norm_error(3,2) > four,"L-inf second order convergence on diffusion")

!todo: remove priints
!print *,norm_error(1,3)/norm_error(1,2)
!print *,norm_error(2,3)/norm_error(2,2)
!print *,norm_error(3,3)/norm_error(3,2)
!
!print *,norm_error(1,2)/norm_error(1,1)
!print *,norm_error(2,2)/norm_error(2,1)
!print *,norm_error(3,2)/norm_error(3,1)
!
!print *, norm_error
!pause


return
end subroutine

end module