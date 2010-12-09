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

!> Routines for grid refinement/coarsening operations.
!>@ingroup test
module test_utility

contains

!> Coarsen a solution at a fine level of resolution
subroutine coarsen(coarse_data, &
                   fine_data,   &
                   ncell_fine,  &
                   ncell_coarse,&
                   nvar)

use stm_precision
use error_handling

implicit none
!---arg
integer,intent(in) :: ncell_coarse                              !< Number of coarsened array cells 
integer,intent(in) :: ncell_fine                                !< Number of fine initial array cells
integer,intent(in) :: nvar                                      !< Number of constituents
real(stm_real), intent(in) :: fine_data(ncell_fine,nvar)        !< Fine initial data  (input)
real(stm_real), intent(out):: coarse_data(ncell_coarse,nvar)    !< Coarsened finial data (output)

!---locals
real(stm_real) :: coarsen_factor                                !< Coarsening factor (must be an integer)
integer :: ivar                                                 !< Counter on constituents
integer :: icell                                                !< Counter
integer :: i_coarse                                             !< Counter

!> Check if the coarsening factor is an integer and if not it bails.
if ( mod(ncell_fine , ncell_coarse) /= 0) then
    call stm_fatal("Coarsening factor is not an integer!")  
else
    coarsen_factor = ncell_fine/ncell_coarse
!> Computes coarsened array base on the coarsening factor from fine input array.
    do ivar=1,nvar
        do icell=1,ncell_coarse
            coarse_data(icell,ivar) = zero
            i_coarse = 0
            do while (i_coarse < coarsen_factor) 
              coarse_data(icell,ivar) = coarse_data(icell,ivar)+ fine_data(icell*coarsen_factor-i_coarse,ivar)
              i_coarse= i_coarse + 1   
            end do
            coarse_data(icell,ivar)= coarse_data(icell,ivar)/dble(coarsen_factor)
        end do
    end do
    
end if

return
end subroutine 


!================================

!> Calculate the L-1, L-2 and L-inf error norms for calculated values and reference solution
subroutine error_norm(norm_1,    &
                      norm_2,    &
                      norm_inf,  &
                      which_cell,&
                      vals,      &
                      reference, &
                      ncell,     &
                      dx)

use stm_precision

implicit none

integer, intent(in) :: ncell                     !< Number of cells
integer, intent(out):: which_cell                !< The cell in which largest error occurs
real(stm_real), intent(out) :: norm_1            !< L-1   error norm
real(stm_real), intent(out) :: norm_2            !< L-2   error norm
real(stm_real), intent(out) :: norm_inf          !< L-inf error norm

real(stm_real), intent(in) :: vals(ncell)        !< Calculated values
real(stm_real), intent(in) :: reference(ncell)   !< Reference or 'other' values
real(stm_real), optional   :: dx                 !< Spatial step !todo: do we use this????

!------ locals                                    
integer :: icell                                
real(stm_real) :: err
real(stm_real) :: sq_error
real(stm_real) :: abs_error
!> initial value of all norms is zero
norm_1=zero
norm_2=zero
norm_inf=zero
!> sum up the L-1 and L-2, and fid the largest error in the domain (L-inf) 
do icell=1,ncell
   err = vals(icell) - reference(icell)
   abs_error = abs(err)
   sq_error = err*err
   if (abs_error > norm_inf) then
       norm_inf = abs_error
       which_cell = icell
   end if
   norm_1 = norm_1 + abs_error
   norm_2 = norm_2 + sq_error
end do
norm_1 = norm_1/dble(ncell)
norm_2 = sqrt(norm_2)/dble(ncell)

return
end subroutine

!> claculates the total mass and checks the oscillations
subroutine mass_calculator(total_mass,  &
                           vals,        &
                           num_cell,    &
                           dx,          &
                           mass_alarm)
use stm_precision
              
implicit none

real(stm_real), intent(out) :: total_mass        !< mass
integer, intent(in) :: num_cell                  !< number of cells
real(stm_real), intent(in) :: vals(num_cell)     !< input values
real(stm_real), intent(in) :: dx                 !< space discretization
logical, intent(in), optional :: mass_alarm      !< negative mass alarm
!---local
integer :: icell

total_mass=zero

total_mass = sum(vals)
total_mass = total_mass*dx

if (present(mass_alarm)) then
    if (mass_alarm ==.true.) then
        do icell=2,num_cell
          if (vals(icell)*vals(icell-1) < zero) then
           ! todo: call error report
           ! todo: how we should find oscillations in case we shift the concentrations up
           !  should we check gradients? grad_left*grad_right >0  
          end if
        end do        
    end if
end if

return
end subroutine

!======================================

!> Create a assert message based on an error ratio
subroutine create_converge_message(converge_message,norm_name,label,ratio)
use stm_precision
implicit none
character(LEN=*),intent(out) :: converge_message !< message
character(LEN=*),intent(in)  :: norm_name        !< name of norm (something like 'L-2 (fine)'
character(LEN=*),intent(in)  :: label            !< label identifying problem
real(stm_real),intent(in)    :: ratio            !< error ratio

write(converge_message,"(a,' 2nd order on ',a, ' (error=', f7.4,')')")norm_name,label,ratio
return
end subroutine
 
 !> Logs convergence results to a file
 !> Outputs the norm-p errors, maximum velocity, discretization parameters and CFL
 ! todo: add peclet number and grid peclet number and other needed terms here
 ! todo: this subroutine assumes a lot, like a scalar reaction rate. Makes it less general
 subroutine log_convergence_results(norm_error,    &
                                    nrefine,       &
                                    dx,            &
                                    dt,            &
                                    max_velocity,  &
                                    label,         &
                                    which_cell,    &
                                    ncell_base,    &
                                    ntime_base,    &
                                    reaction_rate, &
                                    dispersion,    &
                                    scheme_order,  &
                                    length_scale,  &
                                    limiter_switch )
 use stm_precision
 implicit none
 
 integer, parameter  :: log_unit = 91
 integer, intent(in) :: nrefine
 integer, intent(in) :: which_cell(nrefine)
 integer, intent(in) :: ncell_base
 integer, intent(in) :: ntime_base
 real(stm_real),  intent(in) :: norm_error(3,nrefine)
 real(stm_real),  intent(in) :: dx
 real(stm_real),  intent(in) :: dt
 character(LEN=*),intent(in) :: label
 
 real(stm_real),intent(in),optional :: scheme_order
 real(stm_real),intent(in),optional :: max_velocity
 real(stm_real),intent(in),optional :: reaction_rate
 real(stm_real),intent(in),optional :: dispersion
 real(stm_real),intent(in),optional :: length_scale
 real(stm_real) :: refine_rate = two
 logical,intent(in),optional :: limiter_switch
  ! local
 real(stm_real) :: order

if (present (scheme_order)) then
    order = scheme_order
else
    order = two
end if

!todo : which one do we want
open(unit = log_unit, file= trim(label)//'_convergence_log.txt', &
      status='unknown')
!open (unit = log_unit, file= 'log_of_run.txt', status='keep')
    
write(log_unit,*)"==== Convergence test results "// label,' ====' 
write(log_unit,*)
write(log_unit,*)'(from finer to coarser)'
write(log_unit,*)'L-inf error ratio '
write(log_unit,*) norm_error(3,2)/norm_error(3,1),norm_error(3,3)/norm_error(3,2)
write(log_unit,*)'L-2 error ratio '
write(log_unit,*) norm_error(2,2)/norm_error(2,1),norm_error(2,3)/norm_error(2,2)
write(log_unit,*)'L-1 error ratio '
write(log_unit,*) norm_error(1,2)/norm_error(1,1),norm_error(1,3)/norm_error(1,2)
write(log_unit,*)
write(log_unit,*)'L-inf convergence rate estimate'
write(log_unit,*)'fine :',log(norm_error(3,2)/norm_error(3,1))/log(refine_rate),' coarse :',log(norm_error(3,3)/norm_error(3,2))/log(refine_rate)  
write(log_unit,*)'L-2 convergence rate estimate'
write(log_unit,*)'fine :',log(norm_error(2,2)/norm_error(2,1))/log(refine_rate),' coarse :',log(norm_error(2,3)/norm_error(2,2))/log(refine_rate)  
write(log_unit,*)'L-1 convergence rate estimate'
write(log_unit,*)'fine :',log(norm_error(1,2)/norm_error(1,1))/log(refine_rate),' coarse :',log(norm_error(1,3)/norm_error(1,2))/log(refine_rate)  
write(log_unit,*)
write(log_unit,*)'number of cells : ',ncell_base,ncell_base/2,ncell_base/4
write(log_unit,*)'L-inf occures at :',which_cell(1),which_cell(2),which_cell(3)
write(log_unit,*) 'dx :', dx/four,dx/two,dx 
write(log_unit,*)
write(log_unit,*)'number of steps : ',ntime_base,ntime_base/2,ntime_base/4
write(log_unit,*)'dt :',dt/four,dt/two,dt
write(log_unit,*)
write(log_unit,*)'Error L-inf '//label//' : '
write(log_unit,*) norm_error (3,:)
write(log_unit,*)'Error L-2 '//label//' : '
write(log_unit,*) norm_error (2,:)
write(log_unit,*)'Error L-1 '//label//' : '
write(log_unit,*) norm_error (1,:)
write(log_unit,*)

if (present(reaction_rate)) then
    write(log_unit,*) 'Decay rate : ', reaction_rate
    write(log_unit,*) 'Kdt : '
    write(log_unit,*) reaction_rate*dt/four,reaction_rate*dt/two,reaction_rate*dt
    write(log_unit,*)
end if 

if (present(dispersion)) then
    write(log_unit,*) 'Dispersion coefficient : ', dispersion
    write(log_unit,*) 'Diffusion Number (Ddt/dx2) : '
    write(log_unit,*) four*dispersion*dt/dx/dx,two*dispersion*dt/dx/dx,dispersion*dt/dx/dx
    write(log_unit,*)
end if 

if (present(max_velocity)) then
    write(log_unit,*) 'CFL : (<1)' , max_velocity*dt/dx 
    write(log_unit,*) 'Max Velocity', max_velocity
    if (present(limiter_switch)) then
        if (limiter_switch == .true.)then
            write(log_unit,*) 'Flux Limiter : ON '
        else
            write(log_unit,*) 'Flux Limiter : OFF'
        end if
        write(log_unit,*)
    end if  
    
    if (present(dispersion)) then
        write(log_unit,*) 'Mesh Peclet Number(Vdx/D) :'
        write(log_unit,*) max_velocity*dx/dispersion/four,max_velocity*dx/dispersion/two,max_velocity*dx/dispersion
        write(log_unit,*)
        ! todo: do we also need Peclet number? I don't think
    end if
    
    if (present(reaction_rate)) then
        ! todo: this part is the scale of Damkohler
        write(log_unit,*) 'Da : Advection Time Scale/ Reaction Time Scale'
        write(log_unit,*) 'Da =', reaction_rate*dx/max_velocity   
    end if 
    
    write(log_unit,*)
end if
write (log_unit,*) '====================================' 
close (log_unit)
   
return
end subroutine


end module 