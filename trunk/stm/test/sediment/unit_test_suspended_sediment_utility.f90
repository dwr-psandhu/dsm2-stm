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

!> Tests the suspended sedimentutility subroutine
!>@ingroup test 
!todo: is it a correct place in doxygen?

module unit_test_suspend_sed_utility

contains

!> Tests the coarsening subroutine
subroutine test_settling_velocity

use fruit
use suspended_utility
use stm_precision

implicit none
!---arg
real(stm_real) :: w_s                            !< Settling velocity
real(stm_real), parameter :: nu =1.0d-6          !< Kinematic viscosity 
real(stm_real), parameter :: specific_g = 2.65d0 !< Specific gravity of particle (~2.65)
real(stm_real) :: diameter                       !< Particle diameter in meter
real(stm_real), parameter :: g_accel = 9.80d0    !< Gravitational acceleration
!---- local
real(stm_real) :: hand_calc_value                !< Value of the function which is known

! Small value
diameter = 0.1d-6
hand_calc_value = minus * LARGEREAL
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel) 
                       
call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, small diameter!")
             
! Zero value          
diameter = zero
hand_calc_value = minus * LARGEREAL
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel) 
                       
call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, small diameter!")

! Medium size          
diameter = 5.0d-4  ! meter
hand_calc_value = 0.072114059730314d0
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel) 
                       
call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, medium diameter!")

! Large size
diameter = 2.0d-3  ! meter
hand_calc_value = 0.19781658171144d0
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel) 
                       
call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, large diameter!")

return
end subroutine 


end module
