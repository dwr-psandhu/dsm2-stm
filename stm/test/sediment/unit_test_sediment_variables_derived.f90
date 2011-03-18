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

!> Tests the suspended sediment utility subroutine
!>@ingroup test sediment
module unit_test_suspend_sed_utility
!todo: is it a correct place in doxygen?
! todo: this test must be in the test_sediment project but temporarily placed here
contains

!> Tests the coarsening subroutine
!subroutine test_settling_velocity
!
!use fruit
!use suspended_utility
!use stm_precision
!
!implicit none
!!---arg
!real(stm_real) :: w_s                            !< Settling velocity
!real(stm_real), parameter :: nu =1.0d-6          !< Kinematic viscosity 
!real(stm_real), parameter :: specific_g = 2.65d0 !< Specific gravity of particle (~2.65)
!real(stm_real) :: diameter                       !< Particle diameter in meter
!real(stm_real), parameter :: g_accel = 9.80d0    !< Gravitational acceleration
!!---- local
!real(stm_real) :: hand_calc_value                !< Value of the function which is known
!
!! Small value
!diameter = 0.8d-7
!hand_calc_value = minus * LARGEREAL
!call settling_velocity(w_s,              &
!                       nu,               &
!                       specific_g,       &
!                       diameter,         &
!                       g_accel) 
!                       
!call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, small diameter!")
!             
!! Zero value          
!diameter = zero
!hand_calc_value = minus * LARGEREAL
!call settling_velocity(w_s,              &
!                       nu,               &
!                       specific_g,       &
!                       diameter,         &
!                       g_accel) 
!                       
!call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, zero diameter!")
!
!! Medium size          
!diameter = 5.0d-4  ! meter
!hand_calc_value = 0.072114059730314d0
!call settling_velocity(w_s,              &
!                       nu,               &
!                       specific_g,       &
!                       diameter,         &
!                       g_accel) 
!                       
!call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, medium diameter!")
!
!! Large size
!diameter = 2.0d-3  ! meter
!hand_calc_value = 0.19781658171144d0
!call settling_velocity(w_s,              &
!                       nu,               &
!                       specific_g,       &
!                       diameter,         &
!                       g_accel) 
!                       
!call assertEquals(w_s,hand_calc_value,weak_eps,"Error in settling velocity, large diameter!")
!
!return
!end subroutine 

subroutine test_submerged_specific_gravity

use fruit
use suspended_utility
use stm_precision

implicit none
!---args
real(stm_real) :: big_r
real(stm_real) :: rho_w
real(stm_real) :: rho_s
real(stm_real) :: hand_calc_value

rho_w = 1000d0
rho_s = 2650d0
hand_calc_value = 1.65d0

call submerged_specific_gravity(big_r,       &
                                rho_w,       &
                                rho_s)
                                
call assertEquals(big_r,hand_calc_value,weak_eps,"Error in submerged_specific_gravity subroutine!")

return
end subroutine

subroutine test_explicit_particle_reynolds_number

use fruit
use suspended_utility
use stm_precision

implicit none
!---args
integer,parameter  :: nclas = 3            !< Number of sediment diameter classes
real(stm_real)  :: exp_re_p(nclas)         !< Explicit particle reynolds number
real(stm_real)  :: diameter(nclas)         !< Particle diameter
real(stm_real)  :: capital_r               !< Submerged specific gravity of sediment particles  
real(stm_real)  :: g_acceleration          !< Gravitational acceleration 
real(stm_real)  :: kinematic_viscosity     !< Kinematic viscosity (m2/sec)
real(stm_real)  :: hand_calc_value(nclas)

diameter =  (/1d-2,2d-2,1d-3/)
g_acceleration = 9.81d0
capital_r = 1.65d0
kinematic_viscosity = 1.0d-6 

hand_calc_value =  (/4023.2449589852d0,11379.455171492d0,127.22617655184d0/)
 
call explicit_particle_reynolds_number(exp_re_p,           &
                                       diameter,           &
                                       capital_r,          &
                                       g_acceleration,     &
                                       kinematic_viscosity,&
                                       nclas)
                                       
call assertEquals(hand_calc_value(1),exp_re_p(1),weak_eps,"Error in subroutine explicit_particle_reynolds_number!")
call assertEquals(hand_calc_value(2),exp_re_p(2),weak_eps,"Error in subroutine explicit_particle_reynolds_number!")
call assertEquals(hand_calc_value(3),exp_re_p(3),weak_eps,"Error in subroutine explicit_particle_reynolds_number!")

return
end subroutine




end module
