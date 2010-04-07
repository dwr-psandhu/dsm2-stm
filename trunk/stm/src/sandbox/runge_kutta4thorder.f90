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

!> Module orchestrating the diffusion scheme. The main
!> routine in the module is diffuse().
!> Explicit and implicit diffusion operators are included here.
!>@ingroup sand box
module rk4th
  
  contains
!> ODE solver fourth order Runge-Kutta
subroutine rk4(yout,t_start,t_end,y_initial,nstep)

implicit none
use stm precision 
!---arg
integer, intent(in)  :: nstep
real(stm_real),intent(in)  :: t_start
real(stm_real),intent(in)  :: t_end
real(stm_real),intent(in)  :: y_initial(nstep)
real(stm_real),intent(out)  :: yout(nstep)

!----local
real(stm_real)::dh
real(stm_real):: t_position (nstep)
real(stm_real) :: k1,k2,k3,k4

dh = (t_end - t_start)/nstep
t_position =  (/ (t_start + dble(ivar-1)*dh, ivar = 1,nstep+1) /)

yout=y_initial

do ivar=1,nstep
    call func_y(yprim,t_position(ivar),yout(ivar))
    k1=dh*yprim;
    call func_y(yprim,(t_position(ivar)+ dh*half),(yout(ivar) + k1*half))
    k2=dh*yprim;
    call func_y(yprim,(t_position(ivar)+ dh*half),(yout(ivar) + k2*half))
    k3=dh*yprim;
    call func_y(yprim,(t_position(ivar)+ dh*half),(yout(ivar) + k3))
    k4=dh*yprim;
    y(ivar+1)=y(ivar)+ k1/6.0d0+ k2/3.0d0+k3/3.0d0+k4/6.0d0;
    
end do

return
end subroutine rk4

pure subroutine func_y(yprim,t,y)

implicit none

use stm_real
!integer, intent(in)  :: nstep
real(stm_real),intent(out) :: yprim
real(stm_real),intent(in)  :: t
real(stm_real),intent(in)  :: y
!---local
integer :: ivar
!todo:chnage this

yprim = y*sin(t)
!do ivar=1,nstep
!yprim(ivar) = y(ivar)*sin(t(ivar))
!end do

return
end subroutine derivs

end module
