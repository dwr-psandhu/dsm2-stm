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

!> write the log of a convergence test in additio to fruit's output
!>@ingroup tes
module log_convergence

    use stm_precision
    contains
    
    !> Gets the norms, maximum velocity, time step size and spacial step size and 
    !> Compute CFL. Then the subroutine print all of the above 
    subroutine log_convergence_results(norm_error,nrefine,dx,dt,max_velocity,label)

        implicit none

        integer, intent(in) :: nrefine
        real(stm_real),  intent(in) :: norm_error(3,nrefine)
        real(stm_real),  intent(in) :: dx
        real(stm_real),  intent(in) :: dt
        real(stm_real),  intent(in) :: max_velocity
        character(LEN=*),intent(in) :: label

        print *, '========'
        print *, label
        print *,'L-inf convergence ratio : ', norm_error(3,2)/norm_error(3,1)
        print *,'L-2 convergence ratio : ',norm_error(2,2)/norm_error(2,1)
        print *,'L-1 convergence ratio : ',norm_error(1,2)/norm_error(1,1)
        print *, 'dt:',dt,'dx:',dx
        print *, ' CFL = ', max_velocity*dt/dx, 'Max_Velocity', max_velocity 
        print *, 'Error norms '//label//' : ',norm_error
        print *, '========'
!        open (unit=4, file='convergence_log.txt',status = 'unknown')
!        write(filename, "(a\i3\'.txt')"), "uniform_gaussian_end_", nx 

    return
    end subroutine


end module