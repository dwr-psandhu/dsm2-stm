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

!> Routines containing error metrics for assessing convergence or accuracy
!>@ingroup test
module error_handling

 use stm_precision

    contains
        
    !> Logs the input message as an error and stops execution with error code
    subroutine stm_fatal(message)
   
    use logging, only: stm_log, ERROR
    
    implicit none
    character(LEN=*), intent(in)  :: message  !< Message concerning error
    
    call stm_log(ERROR,message)
    stop 1
    
    return
    end subroutine
    
    !> Prints an array to file
    subroutine printout(arr,x,filename)

        implicit none
        real(stm_real),intent(in) :: arr(:)         !< array values
        real(stm_real),intent(in) :: x(:)           !< x values
        character(LEN=*)          :: filename       !< name of file to write
        
        !--local
        integer                   :: icell
        integer                   :: nx
        
        nx = size(arr)
        
        open(unit = 11, file = filename)
        
        do icell = 1,nx
          write(11,'(f10.5, f20.10)') x(icell), arr(icell)
        end do
        close(11)
    end subroutine

    !> Prints an array to file
    subroutine printout_append(arr,x,time,funit)

        implicit none
        real(stm_real),   intent(in)   :: arr(:)         !< array values
        real(stm_real),   intent(in)   :: x(:)           !< x coordinate
        real(stm_real),   intent(in)   :: time           !< time
        integer,          intent(in)   :: funit          !< file unit
                
        !--local
        integer           :: nx
        integer           :: icell       

        
        nx = size(arr)        

        write(funit,*) 'variables = "x", "conc"'
        write(funit,*) "zone i = ", nx, ', t="', time, '"'
        do icell = 1,nx
          write(funit,'(f10.5, f20.10)') x(icell), arr(icell)
        end do

    end subroutine


end module