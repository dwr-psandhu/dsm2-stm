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

!> todo: write tests for matrix solver
!>@ingroup test
module test_matrix_solver

use matrix_solver
use fruit
use stm_precision

contains

subroutine test_tridi_solver

use matrix_solver

implicit none
  
  integer,parameter :: ncell = 11             !< Number of volumes 
    
    real(stm_real)  :: down_diag(ncell)       !< Values of the coefficients below diagonal in matrix
    real(stm_real)  :: center_diag(ncell)     !< Values of the coefficients at the diagonal in matrix
    real(stm_real)  :: up_diag(ncell)         !< Values of the coefficients above the diagonal in matrix
    real(stm_real)  :: right_hand_side(ncell) !< Values of the right hand side vector
    real(stm_real)  :: conc(ncell)            !< Values of the computed solution
   
 !--- Small numbers on center diag large numbers on up and down diag
  
 center_diag = (/0.17D0,0.18D0,0.19D0,0.2D0,0.21D0,0.22D0,0.23D0,0.24D0,0.25D0,0.26D0,0.27D0/)
 down_diag = (/LARGEREAL,26.D0,31.D0,36.D0,41.D0,46.D0,51.D0,56.D0,61.D0,66.D0,71.D0/)
 up_diag = (/36.D0,37.D0,38.D0,39.D0,40.D0,41.D0,42.D0,43.D0,44.D0,45.D0,LARGEREAL/)
 right_hand_side = (/0.01D0,0.02D0,0.03D0,0.04D0,500.0D0,400.0D0,0.07D0,0.08D0,0.09D0,0.1D0,0.11D0/)

call tridi_solver (center_diag,up_diag,down_diag,right_hand_side,conc,ncell)


  call assertEquals (conc(1),613.689083199382D0,1d-9,"problem in solving cell num.1")
  call assertEquals (conc(4),4.520833065491D0,1d-9,"problem in solving cell num.4")
  call assertEquals (conc(11),-834.471026100105D0,1d-9,"problem in solving cell num.11")

!--- Large numbers on center diag small numbers on up and down diag

    center_diag = (/17D0,18D0,19D0,20D0,21D0,22D0,23D0,24D0,25D0,26D0,27D0/)
    down_diag = (/LARGEREAL,0.05d0,0.1d0,0.15d0,0.2d0,0.25d0,0.3d0,0.35d0,0.4d0,0.45d0,0.5d0/)
    up_diag = (/0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1D0,LARGEREAL/)
    right_hand_side = (/0.01D0,0.02D0,0.03D0,0.04D0,0.05d0,0.06D0,0.07D0,0.08D0,0.09D0,0.1D0,0.11D0/)

call tridi_solver (center_diag,up_diag,down_diag,right_hand_side,conc,ncell)


  call assertEquals (conc(1),0.000581809672D0,1d-9,"problem in solving cell num.1")
  call assertEquals (conc(4),0.001942430407D0,1d-9,"problem in solving cell num.4")
  call assertEquals (conc(11),0.004006798484D0,1d-9,"problem in solving cell num.11")

return
end subroutine test_tridi_solver

end module