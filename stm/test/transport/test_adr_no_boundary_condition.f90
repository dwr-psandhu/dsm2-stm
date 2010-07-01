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

contains

use stm_precision

integer, parameter  :: nstep_base = 40
integer, parameter  :: nx_base = 256
integer, parameter  :: nconc = 2
real(stm_real), parameter :: total_time = 6400.0d0 ! sec
real(stm_real), parameter :: domain_length = 51200.0d0 ! m
real(stm_real), parameter :: const_area = 500.0d0 ! m^2
real(stm_real), parameter :: const_disp_coef = 0.01d0 !todo: is it in a correct range? 
real(stm_real), parameter :: const_velocity = 0.8d0 ! m/s
real(stm_real), parameter :: decay_rate = 0.005d0
real(stm_real), parameter :: ic_center = domain_length/four
real(stm_real), parameter :: ic_stand_dev = domain_length/(four*four)
real(stm_real), parameter :: ic_peak = ten

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
character(LEN=*),intent(in) :: label                            !< unique label for test
logical :: verbose                                              !< whether to output convergence results
real(stm_real) :: fine_initial_conc(nx_base,nconc)              !< initial condition at finest resolution
real(stm_real) :: fine_solution(nx_base,nconc)                  !< reference solution at finest resolution

lable = 'ADR uniform flow, const A & Ks'

call initial_final_solution(fine_initial_conc,fine_solution,ic_center,ic_stand_dev,ic_peak,const_velocity,decay_rate,nx_base,nconc)

do icoarse = 1,nrefine

    call allocate_state
    call coarsen
    
    do itime = 1,nstep
    
        call hydro 
        call advect
        call prim2cons
        call diffuse
    
    end do ! itime
    
    call printout
    call error_norm
    call deallocate_state

end do !icoarse

call assert_true(norm_error(1,2)/norm_error(1,1) > four,"L-1 second order convergence on " // trim(label))
call assert_true(norm_error(2,2)/norm_error(2,1) > four,"L-2 second order convergence on " // trim(label))
call assert_true(norm_error(3,2)/norm_error(3,1) > four,"L-inf second order convergence on " // trim(label))

if (verbose == .true.) then
   call log_convergence_results(norm_error,nrefine,dx,dt,max_velocity,label)
end if

return
end subroutine 

end module 
