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

!> Test of mass transport convergence in uniform flow
!>@ingroup test
module test_rk3

use stm_precision

integer, parameter  :: nstep_base = 40
integer, parameter  :: nvar = 2
real(stm_real), parameter :: total_time = 160.d0
real(stm_real), parameter :: decay_coef = 0.01d0

contains
!> Subroutine that runs a small reaction simulation
subroutine test_reaction_decay_convergence_rk3(verbose)

use source_sink
use reaction
use error_metric
use fruit
use log_convergence

implicit none

integer, parameter :: ncell = 4              !< Number of cells
integer, parameter :: nrefine = 3
real(stm_real), parameter :: start_time = zero

! todo: what do we need for reaction? 
real(stm_real):: reference(ncell,nvar)
real(stm_real):: conc(ncell,nvar)              !< Concentration at new time
real(stm_real):: conc_prev(ncell,nvar)         !< Concentration at old time
real(stm_real):: area (ncell)                  !< Cell-centered area at new time
real(stm_real):: area_prev (ncell)             !< Cell-centered area at old time
real(stm_real):: flow(ncell)                   !< cell-centered flow 
real(stm_real):: time                          !< current time
real(stm_real):: dt                            !< Time step   
real(stm_real):: nstep
real(stm_real):: dx
real(stm_real), parameter :: area_const = 17.0d0
real(stm_real), parameter :: flow_const = ten
real(stm_real), parameter :: conc_start = ten
real(stm_real) :: norm_error(3,nrefine)
integer :: itime
integer :: which_cell(nrefine)
integer :: icoarse
logical, optional :: verbose
character(LEN =64) :: test_label 

! for print
 real(stm_real) :: scheme_order
 real(stm_real) :: max_velocity
 real(stm_real) :: reaction_rate
 real(stm_real) :: dispersion
 real(stm_real) :: length_scale

compute_source => linear_decay_source_rk3


test_label = 'Runge-Kutta reaction 3rd order'
area = area_const
area_prev = area_const
flow = flow_const

reference = conc_start * exp (-decay_coef*total_time)

do icoarse= 1,nrefine
     
    nstep =nstep_base /2**(icoarse-1)
    conc_prev = conc_start
    time =zero
    dt = total_time/nstep
    
    do itime =1,nstep

        time = time + dt ! todo:
        call  react(conc,       &
                    conc_prev,  &
                    area,       &
                    area_prev,  &
                    flow,       &
                    ncell,      &
                    nvar,       &
                    time,       &
                    dt)

        conc_prev = conc
        
        call error_norm(norm_error(1,icoarse),  &
                    norm_error(2,icoarse),  &
                    norm_error(3,icoarse),  &
                    which_cell(icoarse),    &
                    conc,                   &
                    reference,              &
                    ncell,                  &
                    dt)
          
                
    end do
    
    
end do

!if (verbose == .true.) then
call log_convergence_results(norm_error ,              &
                             nrefine,                  &
                             dx,                       &
                             dt,                       &
                             max_velocity= zero,       &
                             dispersion = zero,        &
                             label = test_label,       &
                             reaction_rate=decay_coef, &
                             which_cell=which_cell,    &
                             ncell_base = ncell,       &
                             ntime_base = nstep_base,  &
                             scheme_order = three,     &
                             length_scale = dx,        &
                             limiter_switch = .false.)

!end if
 
call assert_true(norm_error(1,2)/norm_error(1,1) > eight,"L-1 second order convergence on " // trim(test_label))
call assert_true(norm_error(2,2)/norm_error(2,1) > eight, "L-2 second order convergence on " // trim(test_label))
call assert_true(norm_error(3,2)/norm_error(3,1) > eight,"L-inf second order convergence on " // trim(test_label))


                                   
end subroutine

!===========
subroutine linear_decay_source_rk3(source, & 
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


! source must be in primitive variable 
source = -decay_coef*conc
  
return
end subroutine 

end module