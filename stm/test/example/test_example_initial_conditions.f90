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

!> Sample initial conditions for tests
!>@ingroup test
subroutine test_example_initial_conditions

use stm_precision
use example_initial_conditions
use gaussian_init_boundary_condition
use fruit

implicit none

integer, parameter :: nloc = 100
real(stm_real) :: vals(nloc,2)
real(stm_real), parameter  :: dx = ten
real(stm_real), parameter  :: origin = zero
real(stm_real), parameter  :: center1 = 405
real(stm_real), parameter  :: center2 = 605.  ! middle of cell 61
real(stm_real), parameter  :: sd = dx*4.0d0
real(stm_real) :: offline_calc
real(stm_real) :: cell_calc41
real(stm_real) :: cell_calc61
real(stm_real) :: miu
real(stm_real) :: sigma
real(stm_real) :: velocity
real(stm_real) :: xpos
real(stm_real) :: dfbydx1
real(stm_real) :: dfbydx2

character(LEN=32) :: message

integer :: icell

! Check symmetry and two constituents
call fill_gaussian(vals(:,1),nloc,origin,dx,center1,sd,one)

call fill_gaussian(vals(:,2),nloc,origin,dx,center2,sd,one)

! test the center cell for each plume. 
! The cell edges lo/hi are half a dx = 1/8 of sd from center
! so use tabulated values of the cdf at the mean +/- 1/8*sigma
!offline_calc = 0.9947645
offline_calc = half + half*erf((410.d0-center1)/(sqrt(two)*sd)) - (half + half*erf((400.0d0-center1)/(sqrt(two)*sd)))
offline_calc = offline_calc*sqrt(two*pi*sd*sd)/dx ! cell average

cell_calc41 = vals(41,1)

call assertEquals(cell_calc41,offline_calc,weak_eps,"Integral gaussian in cell 41")

cell_calc61 = vals(61,2)

call assertEquals(cell_calc41,cell_calc61,eps,"Symmetry of integral gaussian in cells 41, 61")

call fill_gaussian(vals(:,1),nloc,origin,dx/two,center1,sd)

call fill_discontinuity(vals(:,1),nloc,origin,dx,center1,zero,one)
call fill_discontinuity(vals(:,2),nloc,origin,dx,center2,zero,one)
call assertEquals(vals(41,1),half,"Discontinuity IC halfway in cell (41)")
call assertEquals(vals(41,2),zero,"Discontinuity IC halfway in cell (41) -- non-discontinuous constituent")
call assertEquals(vals(61,2),half,"Discontinuity IC halfway in cell (61)")
call assertEquals(vals(61,1),one,"Discontinuity IC halfway in cell (61) -- non-discontinuous constituent")

do icell = 1,40
    message = "Discontinuity IC, lo cells (constituent 1)"
    call assertEquals(vals(icell,1),zero, message)
    message = "Discontinuity IC, lo cells (constituent 2)"    
    call assertEquals(vals(icell,2),zero, message)
end do

do icell = 42,60
    message = "Discontinuity IC, middle (constituent 1)"
    call assertEquals(vals(icell,1),one, message)
    message = "Discontinuity IC, middle (constituent 2)"    
    call assertEquals(vals(icell,2),zero, message)
end do

do icell = 62,100
    message = "Discontinuity IC, hi cells (constituent 1)"
    call assertEquals(vals(icell,1),one, message)
    message = "Discontinuity IC, hi cells (constituent 2)"
    call assertEquals(vals(icell,2),one, message)
end do

! test derivative_gaussian(val,xposition,center,sd,scale)
miu = ten  
sigma = LARGEREAL
velocity = two
xpos = ten

call derivative_gaussian(dfbydx1,xpos,miu,sigma,sixteen)

message = "Error in derivative_gaussian"    
call assertEquals(dfbydx1,zero, message)
  
sigma = five    
xpos = eight
call derivative_gaussian(dfbydx1,xpos,miu,sigma,one)    

xpos = ten+two
call derivative_gaussian(dfbydx2,xpos,miu,sigma,one)  

message = "Error in derivative_gaussian"    
call assertEquals(dfbydx1,-dfbydx2, message)  

sigma=sqrt(two)
miu = one
xpos= one 
call derivative_gaussian(dfbydx1,xpos,miu,sigma,one)
message = "Error in derivative_gaussian"    
call assertEquals(dfbydx1,zero, message)

xpos=two
call derivative_gaussian(dfbydx1,xpos,miu,sigma,one)
message = "Error in derivative_gaussian"    
call assertEquals(dfbydx1,-half*exp(-half*half), weak_eps,message)



return
end subroutine