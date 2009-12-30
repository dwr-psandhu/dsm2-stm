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

contains

!> Calculate the L1, L2 and Linf error between calculated values and a reference

!< Prints an array to file
subroutine printout(arr,x,filename)

    use stm_precision
    implicit none
    real(STM_REAL),intent(in) :: arr(:)         !< array values
    real(STM_REAL),intent(in) :: x(:)           !< x values
    character(LEN=*)          :: filename       !< name of file to write
    integer                   :: icell
    
    !--local
    integer                   :: nx
    
    nx = size(arr)
    
    open(unit = 11, file = filename)
    !write(11,'(a,i5)')    "nx    ", nx
    
    do icell = 1,nx
      write(11,'(f10.5, f20.10)') x(icell), arr(icell)
    end do
    close(11)
end subroutine

end module