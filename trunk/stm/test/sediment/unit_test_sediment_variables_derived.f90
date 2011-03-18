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

subroutine test_particle_reynolds_number

use fruit
use suspended_utility
use stm_precision

implicit none
!---args
integer, parameter :: nclas = 3       !< Number of sediment diameter classes
real(stm_real):: re_p(nclas)          !< Particle Reynolds number
real(stm_real):: settling_v(nclas)    !< Settling velocity
real(stm_real):: diameter(nclas)      !< Particle diameter
real(stm_real):: kinematic_viscosity  !< Kinematic viscosity (m2/sec)
real(stm_real):: hand_calc_value(nclas)


diameter =  (/2d-3,0.25d-3,0.031d-3/) ! coarse silt medium sand and sand
kinematic_viscosity = 1.0d-6 
settling_v = (/162d-3,25.7d-3,0.49d-3/)

hand_calc_value =  (/324.0d0,6.425d0,0.01519d0/)

call particle_reynolds_number(re_p,                &
                              settling_v,          &
                              diameter,            &
                              kinematic_viscosity, &
                              nclas)
                              
                              

call assertEquals(hand_calc_value(1),re_p(1),weak_eps,"Error in subroutine particle_reynolds_number!")
call assertEquals(hand_calc_value(2),re_p(2),weak_eps,"Error in subroutine particle_reynolds_number!")
call assertEquals(hand_calc_value(3),re_p(3),weak_eps,"Error in subroutine particle_reynolds_number!")
                                         
return
end subroutine

subroutine test_dimless_particle_diameter

use fruit
use suspended_utility
use stm_precision

implicit none

integer, parameter :: nclas = 2        !< Number of sediment diameter classes
real(stm_real):: d_star(nclas)          
real(stm_real):: capital_r             !< Submerged specific gravity of sediment particles  
real(stm_real):: g_accel               !< Gravitational acceleration 
real(stm_real):: diameter(nclas)       !< Particle diameter
real(stm_real):: kinematic_viscosity   !< Kinematic viscosity (m2/sec)
real(stm_real):: hand_calc_value(nclas)

diameter =  (/2d-3,0.25d-3/) ! coarse silt and medium sand 
kinematic_viscosity = 1.0d-6
g_accel = 9.81d0
capital_r = 1.65d0

hand_calc_value =  (/50.591898800422d0,6.3239873500d0/)

call dimless_particle_diameter(d_star,                 &
                               g_accel,                &
                               diameter,               &
                               kinematic_viscosity,    &
                               capital_r,              &
                               nclas)                        

call assertEquals(hand_calc_value(1),d_star(1),weak_eps,"Error in subroutine dimensionless particle number!")
call assertEquals(hand_calc_value(2),d_star(2),weak_eps,"Error in subroutine dimensionless particle number!")

return 
end subroutine

subroutine test_critical_shields_parameter

use fruit
use suspended_utility
use stm_precision

implicit none

integer, parameter :: nclas = 7          !< Number of sediment diameter classes
real(stm_real):: d_star(nclas)          
real(stm_real):: cr_shields_prmtr(nclas) !< Critical Shields parameter                                      
real(stm_real):: hand_calc_value(nclas)
integer :: iclas

d_star =  (/160d0,21d0,15d0,10d0,2d0,1d0,-4d0/) ! coarse silt and medium sand 

hand_calc_value =  (/0.055d0,             &
                     0.031433080718165d0, &
                     0.030510608231307d0, &
                     0.0320721471387d0,   &
                     0.12d0,              &
                     LARGEREAL,           &
                     LARGEREAL/)

call critical_shields_parameter(cr_shields_prmtr,   &
                                d_star,             &
                                nclas)                        

do iclas=1,nclas
    call assertEquals(hand_calc_value(iclas),cr_shields_prmtr(iclas),weak_eps,"Error in subroutine dcritical_shields_parameter!")
end do


return 
end subroutine



end module
