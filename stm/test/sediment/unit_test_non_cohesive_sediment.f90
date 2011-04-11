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

real(stm_real)  :: rouse         !< Rouse dimenssionless number  
real(stm_real)  :: delta         !< Relative bed layer thickness = b/H 
real(stm_real)  :: J_1           !< First Einstein integral value     
!--- local
real(stm_real)  :: hand_calc_value

delta = 0.01d0
rouse = 0.1d0
hand_calc_value = 0.630990839362793d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral!")

rouse = 0.7d0
hand_calc_value = 0.075646372654714d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral!")

rouse = 1.7d0
hand_calc_value = 0.011612330444738d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral!")

rouse = 2.7d0
hand_calc_value = 0.005925241451994d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral!")


rouse = one
hand_calc_value = 0.03651687056d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral integer=1!")


rouse = two
hand_calc_value =   0.009262285443120d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral integer=2!")

rouse = three
hand_calc_value =  0.004859662341771d0 !MATLAB calculation
 
call first_einstein_integral(J_1,    &
                             delta,  &
                             rouse)
                                       
call assertEquals(hand_calc_value,J_1,weak_eps,"Error in subroutine first Einstein integral integer=3!")


return
end subroutine

subroutine test_es_garcia_parker()

use fruit
use non_cohesive_source
use stm_precision

integer,parameter :: nvol = 3                  !< Number of computational volumes in a channel
integer,parameter :: nclass = 2                !< Number of non-cohesive sediment grain classes
real(stm_real) :: e_s(nvol,nclass)             !< Dimenssionless rate of entrainment of bed sediment into suspension 
real(stm_real) :: shear_v(nvol)                !< Shear Velocity
real(stm_real) :: exp_re_p(nclass)             !< Explicit particle Reynolds number
real(stm_real) :: settling_v(nclass)           !< Settling velocity
!---local
real(stm_real) :: hand_calc_value(nvol,nclass)
integer :: ivol

shear_v =[0.1d0,0.4d0,one]
exp_re_p =[two,ten]
settling_v = [0.001d0,0.1d0]

hand_calc_value = reshape ([0.2999514d0,	0.3000000d0,	0.3000000d0, &
                            0.0001299d0,	0.0922054d0,	0.2932331d0],[3,2])


call es_garcia_parker(e_s,         &
                      shear_v,     &
                      exp_re_p,    &
                      settling_v,  & 
                      nclass,      &
                      nvol)
                      
do ivol=1,nvol
  call assertEquals(hand_calc_value(ivol,1),e_s(ivol,1),weak_eps,"Error in subroutine es_garcia_parker")
  call assertEquals(hand_calc_value(ivol,2),e_s(ivol,2),weak_eps,"Error in subroutine es_garcia_parker")
end do 




return
end subroutine


end module
