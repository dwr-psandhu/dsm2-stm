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

module Multiply

contains

  subroutine multiplyFiveTimes(a,b)
    implicit none
    integer,intent(IN) :: a
    integer,intent(OUT) :: b

    b=5*a
    continue
  end subroutine

  subroutine multiplyTenTimes(a,b) 
    implicit none
    integer,intent(IN) :: a
    integer,intent(OUT) :: b
    
    b=10*a
  end subroutine

  subroutine mainSub(a, multiply)

    integer, intent(IN) :: a
    integer :: b
    integer :: i
    interface
      subroutine multiply(c,d)
         integer :: c,d
      end subroutine
    end interface
    
    do i = 1, 7 
      call multiply(a,b)
      print *,b
    enddo
  
  end subroutine

end module Multiply

program callback_simple
use Multiply
implicit none

 call mainSub(35, multiplyTenTimes)

  pause
end program



