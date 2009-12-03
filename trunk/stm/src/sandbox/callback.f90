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

program callback_example
implicit none
include 'callback.fi'
integer :: intarg = 1
external callback1
external callback2
procedure(callback),pointer :: callback_ptr1 
callback_ptr1 => callback1
call sub1(intarg,callback_ptr1)
callback_ptr1 => callback2
call sub1(intarg,callback_ptr1)
end program

subroutine callback1(intarg)
implicit none
integer,intent(in) :: intarg
print*,"Here is the argument",intarg
end subroutine

subroutine callback2(intarg)
implicit none
integer,intent(in) :: intarg
print*,"Here is the argument times 2",2*intarg
end subroutine

subroutine sub1(intarg,callback)
implicit none
include 'callback.fi'
integer,intent(in) :: intarg
call callback(intarg)
end subroutine



