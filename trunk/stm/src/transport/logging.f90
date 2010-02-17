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
module logging

    use stm_precision
    ! todo: an initialization expression is required when using the PARAMETER attribute
    ! i added one
    integer,parameter :: ERROR = 8
    integer,parameter :: WARNING = 6
    integer,parameter :: FINE = 4
    integer,parameter :: INFO =2
    
    
    interface printout
        module procedure printout
        module procedure printout_append
    end interface
    
    
    contains
    
    subroutine stm_log(level,message)
    implicit none
    ! todo: This name does not have a type, and must have an explicit type
    ! I added something just to compile
    integer :: level
    character(LEN=*) :: message
    !todo: do real logging
    print*,message
    return
    end subroutine
    
    !< Prints an array to file
    subroutine printout(arr,x,filename)

        implicit none
        real(STM_REAL),intent(in) :: arr(:)         !< array values
        real(STM_REAL),intent(in) :: x(:)           !< x values
        character(LEN=*)          :: filename       !< name of file to write
        integer                   :: icell
        
        !--local
        integer                   :: nx
        
        nx = size(arr)
        
        open(unit = 11, file = filename)
        
        do icell = 1,nx
          write(11,'(f10.5, f20.10)') x(icell), arr(icell)
        end do
        close(11)
    end subroutine

    !< Prints an array to file
    subroutine printout_append(arr,x,time,funit)

        implicit none
        real(STM_REAL),   intent(in)   :: arr(:)         !< array values
        real(STM_REAL),   intent(in)   :: x(:)           !< x coordinate
        real(STM_REAL),   intent(in)   :: time           !< time
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