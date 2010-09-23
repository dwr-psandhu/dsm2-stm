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

!> Log results of a convergence test
!>@ingroup tes
module log_convergence
 contains

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
write(log_unit,*)'fine :',log(norm_error(3,2)/norm_error(3,1))/log(order),' coarse :',log(norm_error(3,3)/norm_error(3,2))/log(refine_rate)  
write(log_unit,*)'L-2 convergence rate estimate'
write(log_unit,*)'fine :',log(norm_error(2,2)/norm_error(2,1))/log(order),' coarse :',log(norm_error(2,3)/norm_error(2,2))/log(refine_rate)  
write(log_unit,*)'L-1 convergence rate estimate'
write(log_unit,*)'fine :',log(norm_error(1,2)/norm_error(1,1))/log(order),' coarse :',log(norm_error(1,3)/norm_error(1,2))/log(refine_rate)  
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