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
subroutine test_settling_velocity

use fruit
use suspended_utility
use stm_precision

implicit none
!---arg
integer,parameter :: nclas = 5
real(stm_real) :: w_s(nclas)                     !< Settling velocity
real(stm_real), parameter :: nu =1.0d-6          !< Kinematic viscosity 
real(stm_real), parameter :: specific_g = 2.65d0 !< Specific gravity of particle (~2.65)
real(stm_real) :: diameter(nclas)                !< Particle diameter in meter
real(stm_real), parameter :: g_accel = 9.80d0    !< Gravitational acceleration
real(stm_real) :: hand_calc_value(nclas)         !< Value of the function which is known
logical :: pick_up_function
integer :: iclas
! Small value
diameter = [0.8d-7,zero,5.0d-4,5.0d-5,2.0d-3]
hand_calc_value = [-LARGEREAL,-LARGEREAL,0.072114059730314d0,0.0022458333333d0,0.19781658171144d0] 
! van Rijn 
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel,          &
                       nclas) 
do iclas=1,nclas                                                                   
    call assertEquals(w_s(iclas),hand_calc_value(iclas),weak_eps,"Error in settling velocity, van Rijn, no optional input!")
end do

pick_up_function =.true.
! agian van Rijn
call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel,          &
                       nclas,            &
                       pick_up_function) 
                       
do iclas=1,nclas                                                                   
    call assertEquals(w_s(iclas),hand_calc_value(iclas),weak_eps,"Error in settling velocity, van Rijn optional input=.true.!")
end do

!Dietrich 
pick_up_function =.false.

diameter = [100d-3,10d-3,1d-3,0.1d-3,0.01d-3]

hand_calc_value = [1.9697877833755d0,0.741517868347728d0,0.15497120869d0,0.0074779137192d0,7.9999d-05]


call settling_velocity(w_s,              &
                       nu,               &
                       specific_g,       &
                       diameter,         &
                       g_accel,          &
                       nclas,            &
                       pick_up_function)


do iclas=1,nclas                                                                   
    call assertEquals(w_s(iclas),hand_calc_value(iclas),weak_eps,"Error in settling velocity, Dietrich optional input=.false.!")
end do

return
end subroutine 

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

diameter =  [1d-2,2d-2,1d-3]
g_acceleration = 9.81d0
capital_r = 1.65d0
kinematic_viscosity = 1.0d-6 

hand_calc_value =  [4023.2449589852d0,11379.455171492d0,127.22617655184d0]
 
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


diameter =  [2d-3,0.25d-3,0.031d-3] ! coarse silt medium sand and sand
kinematic_viscosity = 1.0d-6 
settling_v = [162d-3,25.7d-3,0.49d-3]

hand_calc_value =  [324.0d0,6.425d0,0.01519d0]

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

diameter =  [2d-3,0.25d-3] ! coarse silt and medium sand 
kinematic_viscosity = 1.0d-6
g_accel = 9.81d0
capital_r = 1.65d0

hand_calc_value =  [50.591898800422d0,6.3239873500d0]

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

d_star =  [160d0,21d0,15d0,10d0,2d0,1d0,-4d0] ! coarse silt and medium sand 

hand_calc_value =  [0.055d0,              &
                     0.031433080718165d0, &
                     0.030510608231307d0, &
                     0.0320721471387d0,   &
                     0.12d0,              &
                     LARGEREAL,           &
                     LARGEREAL]

call critical_shields_parameter(cr_shields_prmtr,   &
                                d_star,             &
                                nclas)                        

do iclas=1,nclas
    call assertEquals(hand_calc_value(iclas),cr_shields_prmtr(iclas),weak_eps,"Error in subroutine critical_shields_parameter!")
end do

return 
end subroutine

subroutine test_shear_velocity

use fruit
use suspended_utility
use stm_precision

implicit none

integer, parameter :: ncell = 3          !< Number cells 
real(stm_real):: vel(ncell)              !< Velocity          
real(stm_real):: manning_n(ncell)        !< Manning's n                                     
real(stm_real):: hand_calc_value(ncell)  !< The sought output 
real(stm_real):: big_r(ncell)            !< Hydraulic radius 
real(stm_real):: gravity                !< Gravity
real(stm_real):: shear_v(ncell)          !< Shear velocity 
logical :: si_br                        !< SI and British unit switch
integer :: iclas

vel =  [1.1d0,.7d0,-1.5d0]    ! values for a river
manning_n = [0.02d0,0.03d0,0.045d0]
gravity = 9.8d0
big_r = [three,five,seven]

hand_calc_value =  [ 0.057347634619921d0,   0.050273292832295d0,   0.152780222207618d0]

call shear_velocity_calculator(shear_v,   &
                               vel,       &
                               manning_n, &
                               gravity,   &
                               big_r,     &
                               ncell)                      

do iclas=1,ncell
    call assertEquals(hand_calc_value(iclas),shear_v(iclas),weak_eps,"Error in subroutine Shear Velocity!")
end do

si_br =.true.

call shear_velocity_calculator(shear_v,    &
                                vel,       &
                                manning_n, &
                                gravity,   &
                                big_r,     &
                                ncell,      &
                                si_br)                      

do iclas=1,ncell
    call assertEquals(hand_calc_value(iclas),shear_v(iclas),weak_eps,"Error in subroutine Shear Velocity, SI unit!")
end do

gravity =32.2d0
si_br = .false.
hand_calc_value =  [ 0.069953846244591d0,   0.061324415911969d0,   0.186364516066955d0]

call shear_velocity_calculator(shear_v,    &
                                vel,       &
                                manning_n, &
                                gravity,   &
                                big_r,     &
                                ncell,      &
                                si_br)                      

do iclas=1,ncell
    call assertEquals(hand_calc_value(iclas),shear_v(iclas),weak_eps,"Error in subroutine Shear Velocity, British unit!")
end do


return 
end subroutine

subroutine test_rouse_number()

use fruit
use suspended_utility
use stm_precision

implicit none

integer, parameter :: nclas   = 2
integer, parameter :: nvolume = 3

real(stm_real) :: rouse_num(nvolume,nclas)   !< Rouse dimensionless number  
real(stm_real) :: fall_vel(nclas)            !< Settling velocity
real(stm_real) :: shear_vel(nvolume)         !< Shear velocity 
real(stm_real) :: hand_value(nvolume,nclas)  !< Calculated values 
real(stm_real) :: von_karman                 !< Von karman constant, Kappa = 0.41
!---local
integer:: iclas,ivol

fall_vel  = [0.001d0, 0.1d0]
shear_vel = [one,two,five]/ten

hand_value = reshape ([0.024390244d0,	0.012195122d0,	0.004878049d0, &
                       2.439024390d0,	1.219512195d0,	0.487804878d0 ],[3,2])

call rouse_dimensionless_number(rouse_num,   &
                                fall_vel,    &
                                shear_vel,   &
                                nvolume,     &
                                nclas)
                                
do iclas=1,nclas
    do ivol =1, nvolume
        call assertEquals(hand_value(ivol,iclas),rouse_num(ivol,iclas),weak_eps,"Error in subroutine Rouse number!")
    end do
end do

hand_value = hand_value/two
von_karman = 0.82d0

call rouse_dimensionless_number(rouse_num,   &
                                fall_vel,    &
                                shear_vel,   &
                                nvolume,     &
                                nclas,       &
                                von_karman)
                                
do iclas=1,nclas
    do ivol =1, nvolume
        call assertEquals(hand_value(ivol,iclas),rouse_num(ivol,iclas),weak_eps,"Error in subroutine Rouse number!")
    end do
end do
                                


return
end subroutine


subroutine test_allocation_ratio()

use fruit
use suspended_utility
use stm_precision

implicit none

integer, parameter :: nclas   = 2
integer, parameter :: nvolume = 3

real(stm_real) :: rouse_num(nvolume,nclas)    !< Rouse dimensionless number  
real(stm_real) :: susp_percent(nvolume,nclas) !< Percentage in suspension  
real(stm_real) :: bed_percent(nvolume,nclas)  !< Percentage in bedload
real(stm_real) :: hand_value(nvolume,nclas)   !< Calculated value
!---local
integer:: iclas,ivol

rouse_num  = reshape ([0.5d0,	one ,	1.1d0, &
                       2d0,	 5.5d0,	8.5d0 ],[3,2])
                                             

hand_value = reshape ([1.000000000000000d0,   0.919698602928606d0,   0.832177709245199d0, &
                       0.338338208091532d0,   0.010216928596160d0,   0.000508670922527d0],[3,2])

call allocation_ratio(susp_percent,    &
                      bed_percent,     &
                      rouse_num,       &
                      nclas,           &
                      nvolume)
                                
do iclas=1,nclas
    do ivol =1, nvolume
        call assertEquals(hand_value(ivol,iclas),susp_percent(ivol,iclas),weak_eps,"Error in subroutine bedload allocation ratio!")
    end do
end do

hand_value = one - hand_value

do iclas=1,nclas
    do ivol =1, nvolume
        call assertEquals(hand_value(ivol,iclas),bed_percent(ivol,iclas),weak_eps,"Error in subroutine bedload allocation ratio!")
    end do
end do


return
end subroutine




end module
