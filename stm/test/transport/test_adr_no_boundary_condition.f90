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

integer, parameter  :: nstep_base = 40
integer, parameter  :: nx_base = 256
integer, parameter  :: nconc = 2
real(stm_real), parameter :: start_time = 1000.0d0 ! sec
real(stm_real), parameter :: total_time = 1280.0d0 ! sec
! todo:  since the bc is set to be zero flux, total_time and other parameters should be set 
! in the way solution does not reach the edges of channel.
real(stm_real), parameter :: domain_length = 51200.0d0 ! m
real(stm_real), parameter :: origin = zero ! low side of channel
real(stm_real), parameter :: const_area = 500.0d0 ! m^2
real(stm_real), parameter :: const_disp_coef = 0.01d0 !todo: is it in a correct range? 
real(stm_real), parameter :: const_velocity = 0.8d0 ! m/s
real(stm_real), parameter :: decay_rate = 0.001d0
real(stm_real), parameter :: ic_center = domain_length/four
real(stm_real), parameter :: ic_stand_dev = domain_length/(four*four)
real(stm_real), parameter :: ic_peak = ten

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

implicit none

procedure(hydro_data_if), pointer :: hydro                      !< This pointer, points to uniform flow data
character(LEN=64) :: label                                       !< unique label for test
logical :: verbose                                              !< whether to output convergence results
real(stm_real) :: fine_initial_conc(nx_base,nconc)              !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)                  !< reference solution at finest resolution
!---local
integer, parameter :: nrefine = 3
integer, parameter :: coarsen_factor = 2                 ! coarsening factor used for convergence test
integer :: itime
integer :: icell 
integer :: icoarse 
integer :: nstep
integer :: nx
integer :: coarsening
logical, parameter :: limit_slope = .false.
real(stm_real), allocatable :: solution_mass(:,:)
real(stm_real), allocatable :: reference(:,:)
real(stm_real), allocatable :: x_center(:)
real(stm_real), allocatable :: velocity (:)
real(stm_real) :: dt              ! seconds
real(stm_real) :: dx              ! meters
real(stm_real) :: time
real(stm_real) :: norm_error(3,nrefine)
!------
label = 'ADR uniform flow, const A & Ks'

call initial_final_solution(fine_initial_conc,fine_solution,ic_center,ic_stand_dev,ic_peak,const_velocity,decay_rate,total_time,origin,domain_length,nx_base,nconc)

do icoarse = 1,nrefine

  !  call allocate_state
!    call coarsen
    
    do itime = 1,nstep
    
 !       call hydro 
 !       call advect
  !      call prim2cons
  !      call diffuse
    
    end do ! itime
    
  !  call printout
 !   call error_norm
 !   call deallocate_state

end do !icoarse

!call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on " // trim(label))
!call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on " // trim(label))
!call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on " // trim(label))

if (verbose == .true.) then
!   call log_convergence_results(norm_error,nrefine,dx,dt,max_velocity,label)
end if

return
end subroutine 
!===========
!> produce fine initial condition and reference solution 
subroutine initial_final_solution(fine_initial_conc,fine_solution,ic_center,ic_stand_dev,ic_peak,const_velocity,decay_rate,total_time,origin,domain_length,nx_base,nconc)

use example_initial_conditions

implicit none
integer, intent(in) :: nx_base
integer, intent(in) :: nconc
real(stm_real),intent(out) :: fine_initial_conc(nx_base,nconc)              !< initial condition at finest resolution
real(stm_real),intent(out) :: fine_solution(nx_base,nconc)                  !< reference solution at finest resolution
real(stm_real),intent(in)  :: ic_center
real(stm_real),intent(in)  :: ic_stand_dev
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

final_peak = ic_peak * exp(- decay_rate * total_time)
final_center = ic_center + const_velocity * total_time
final_stand_dev = ic_stand_dev * (total_time + start_time)/start_time

do ivar = 1,nconc
    call fill_gaussian(fine_initial_conc(1,ivar),nx_base,origin,dx,ic_center,ic_stand_dev,ic_peak)
    call fill_gaussian(fine_solution(1,ivar),nx_base,origin,dx,final_center,final_stand_dev,final_peak)
end do

fine_initial_conc(2,:) = fine_initial_conc(1,:) 
fine_solution(2,:)     = fine_solution(1,:)

return
end subroutine

end module 
