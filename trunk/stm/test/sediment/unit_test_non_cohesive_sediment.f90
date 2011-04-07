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

!> Tests the suspended non cohesive sediment subroutines
!>@ingroup test sediment
module test_non_cohesive

contains

!> Tests Einsstein's first integral 
! todo: incase the main subroutine replaced somewhere else this counterpart should place in the correct test package
subroutine test_first_einstein_integral

use fruit
use non_cohesive_source
use stm_precision

implicit none
!---args

real(stm_real)  :: exp_re_p(nclas)         !< Explicit particle reynolds number
real(stm_real)  :: diameter(nclas)         !< Particle diameter
real(stm_real)  :: capital_r               !< Submerged specific gravity of sediment particles  
!--- local
real(stm_real)  :: hand_calc_value

diameter =  (/1d-2,2d-2,1d-3/)
g_acceleration = 9.81d0
capital_r = 1.65d0
kinematic_viscosity = 1.0d-6 

hand_calc_value = 
 
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
