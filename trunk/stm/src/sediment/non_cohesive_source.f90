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

!> Routines provide the calculation for non-cohesive suspended sediment erosion and deposition functions.
!> The functions in this module are: 1.	Garcia and Parker (1991) 2. van Rijn (1984b)
!>  3.	Smith and McLean (1977) 4.	Zyserman and Fredsoe (1994)  for more details see ASCE manual #110 page 117.
!>@ingroup sediment 

module non_cohesive_source

contains 
subroutine source_non_cohesive(nvol,nclass,pick_up_flag,dx,dt,time)

use stm_precision 
use suspended_sediment_variable
use sediment_variables
use suspended_utility

implicit none
real(stm_real),intent(in) :: dx
real(stm_real),intent(in) :: dt
real(stm_real),intent(in) :: time
integer, intent(in) :: nvol
integer, intent(in) :: nclass
character(len=32), optional, intent(in) :: pick_up_flag
!---local
character :: pick_up_function 
procedure(sediment_hydro_if),pointer :: velocity_non_cohesive 

velocity_non_cohesive => sediment_velocity_width

pick_up_function = 'garcia_parker'

if (present(pick_up_flag)) then
     pick_up_function = pick_up_flag
end if

!- initialization to LARGEREAL
call set_sediment_constants
call allocate_sediment_parameters(nvol,nclass)
! gives manning n + width and non_cohesive diameters
!- getting the values
call set_sediment_values(gravity,                 &                 
                         water_density,           &           
                         sediment_density,        &        
                         kappa,                   &                   
                         kinematic_viscosity,     &     
                         floc_density,            &            
                         cohesive_diameter,       &       
                         crit_stress_full_dep,    &    
                         density_wet_bulk,        &        
                         crit_stress_partial_dep, & 
                         crit_stress_surf_erosion,&
                         density_dry_bulk,        &        
                         ta_floc)  
                         
call shear_velocity 

select case (pick_up_function)
   
    case('garcia_parker')
    
   call entrainment_garcia_parker()
   
   call deposition()
   
      
!   case('zyserman_fredsoe')
!   
!    subroutine entrainment_zyserman_fredsoe()
!      implicit none 
!
!    end subroutine entrainment_zyserman_fredsoe ! todo: remove 
!   
!   case('van_rijn')
!   
!   subroutine entrainment_van_rijn()
!      implicit none 
!
!    end subroutine  entrainment_van_rijn
!    
!   case('smith_mclean')
!    call entrainment_smith_mclean()
  
end select
   
call deallocate_sediment_static()

end subroutine 

    subroutine entrainment_garcia_parker()
      implicit none 

    end subroutine 

!> Calculates the first Einstein integral values
!> This subroutine is developed based on analtycal solution of Guo and Julien (2004)
!> To avoid disambiguation: C_bar = c_b_bar * first_einstein_integral
!> the out put of the subroutine is equal to J_1 in the page 116 of ASCE sediment manual  
!> To avoid singularities here an analytical solution used for integers    
! todo: Should we place this subroutine here? another separate file? or sediment derived variable?
! I think we will use it again in the bedload

subroutine first_einstein_integral(I_1,      &
                                   delta_b,  &
                                   rouse_num) 
                                   
use stm_precision
implicit none

real(stm_real),intent(in) :: rouse_num        !< Rouse dimenssionless number  
real(stm_real),intent(in) :: delta_b          !< Relative bed layer thickness = b/H 
real(stm_real),intent(out):: I_1              !< First Einstein integral value


if (rouse_num > 3.98d0) then
!todo: I am not sure if we need this subroutine in bed load or not 
    call stm_fatal("This is not a Rouse number value for suspended sediment!")
elseif (abs(rouse_num - three)< 0.02d0) then
    I_1 = -three*log(delta_b)* + one/(two*delta_b*delta_b) - three/delta_b + three/two + delta_b
elseif (abs(rouse_num - two)< 0.02d0) then
    I_1 = 1/delta_b + half + two * log(delta_b) - delta_b
elseif(abs(rouse_num - one)< 0.02d0) then  
    I_1 = -log(delta_b) + delta_b - one
else
I_1   = (rouse_num*pi/sin(rouse_num*pi) - ((1-delta_b)**rouse_num)/(delta_b**(rouse_num-1)) &
         - rouse_num*(((delta_b/(1-delta_b))**(1-rouse_num))/(1-rouse_num))                 & 
         + rouse_num*(((delta_b/(1-delta_b))**(2-rouse_num))/(1-rouse_num))                 &
         - rouse_num*(((delta_b/(1-delta_b))**(3-rouse_num))/(1-rouse_num))                 &
         + rouse_num*(((delta_b/(1-delta_b))**(4-rouse_num))/(1-rouse_num))                 &
         - rouse_num*(((delta_b/(1-delta_b))**(5-rouse_num))/(1-rouse_num))                 &
         + rouse_num*(((delta_b/(1-delta_b))**(6-rouse_num))/(1-rouse_num))                 &
         - rouse_num*(((delta_b/(1-delta_b))**(7-rouse_num))/(1-rouse_num))                 &
         + rouse_num*(((delta_b/(1-delta_b))**(8-rouse_num))/(1-rouse_num)))                &
         * (delta_b**(rouse_num)/((1-delta_b)**rouse_num))
         
end if
                                   

end subroutine

end module 