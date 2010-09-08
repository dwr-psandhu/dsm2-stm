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
 
 !> Logs convergence results to a file
 !> Outputs the norm-p errors, maximum velocity, discretization parameters and CFL
 ! todo: add peclet number and grid peclet number and other needed terms here
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
 logical,intent(in),optional :: limiter_switch
  ! local
 integer :: icell
 real(stm_real) :: order

if (present (scheme_order)) then
    order = scheme_order
else
    order = two
end if

open (unit = 3, file= label//'_error_log.txt', status='unknown')
    
write(3,*)"==== Log of connvergence test "// label,' ====' 
write(3,*)
write(3,*)'L-inf convergence ratio : '
write(3,*)'fine :',log(norm_error(3,2)/norm_error(3,1))/log(order),' coarse :',log(norm_error(3,3)/norm_error(3,2))/log(order)  
write(3,*)'L-2 convergence ratio : '
write(3,*)'fine :',log(norm_error(2,2)/norm_error(2,1))/log(order),' coarse :',log(norm_error(2,3)/norm_error(2,2))/log(order)  
write(3,*)'L-1 convergence ratio : '
write(3,*)'fine :',log(norm_error(1,2)/norm_error(1,1))/log(order),' coarse :',log(norm_error(1,3)/norm_error(1,2))/log(order)  
write(3,*)
write(3,*)'number of cells : ',ncell_base,ncell_base/2,ncell_base/4
write(3,*)'L-inf occures at :',which_cell(1),which_cell(2),which_cell(3)
write(3,*) 'dx :', dx/four,dx/two,dx 
write(3,*)
write(3,*)'number of steps : ',ntime_base,ntime_base/2,ntime_base/4
write(3,*)'dt :',dt/four,dt/two,dt
write(3,*)
write(3,*)'Error L-inf '//label//' : '
write(3,*) norm_error (3,:)
write(3,*)'Error L-2 '//label//' : '
write(3,*) norm_error (2,:)
write(3,*)'Error L-1 '//label//' : '
write(3,*) norm_error (1,:)
write(3,*)

if (present(reaction_rate)) then
    write(3,*) 'Decay rate : ', reaction_rate
    write(3,*) 'Kdt : '
    write(3,*) reaction_rate*dt/four,reaction_rate*dt/two,reaction_rate*dt
    write(3,*)
end if 

if (present(dispersion)) then
    write(3,*) 'Dispersion coefficient : ', dispersion
    write(3,*) 'Diffusion Number (Ddt/dx2<1) : '
    write(3,*) four*dispersion*dt/dx/dx,two*dispersion*dt/dx/dx,dispersion*dt/dx/dx
    write(3,*)
end if 

if (present(max_velocity)) then
    
    write (3,*) 'CFL : (<1)' , max_velocity*dt/dx 
    
    if (present(limiter_switch)) then
        if (limiter_switch == .true.)then
            write (3,*) 'Flux Limiter : ON '
        else
            write(3,*) 'Flux Limiter : OFF'
        end if
        write(3,*)
    end if  
    
    if (present(dispersion)) then
        write(3,*) 'Mesh Peclet Number(Vdx/D) :'
        write(3,*) max_velocity*dx/dispersion/four,max_velocity*dx/dispersion/two,max_velocity*dx/dispersion
        write(3,*)
        ! todo: do we also need Peclet number? I don't think
    end if
    
    if (present(reaction_rate)) then
        ! todo: this part is the scale of Damkohler
        write(3,*) 'Da : Advection Time Scale/ Reaction Time Scale'
        write(3,*) 'Da =', reaction_rate*dx/max_velocity   
    end if 
    
    write(3,*)
end if
write (3,*) '====================================' 
close (3)
   
return
end subroutine

end module