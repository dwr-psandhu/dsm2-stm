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
subroutine source_non_cohesive(vertical_flux,    &
                               conc,             &
                               nvol,             &
                               nclass,           &
                               pick_up_flag,     &
                               dx,               &
                               dt,               &
                               time)

use stm_precision 
use suspended_sediment_variable
use sediment_variables
use suspended_utility

implicit none
real(stm_real),intent(out):: vertical_flux(nvol,nclass)    !< Vertical sediment net flux into the water column
real(stm_real),intent(in) :: conc(nvol)                    !< Concentration at new time
real(stm_real),intent(in) :: dx                            !< Grid size in space
real(stm_real),intent(in) :: dt                            !< Step size in time
real(stm_real),intent(in) :: time                          !< Current time
integer, intent(in)       :: nvol                          !< Number of cells 
integer, intent(in)       :: nclass                        !< Number of classes for non-cohesive sediment in suspension
character(len=32), optional, intent(in) :: pick_up_flag    !< Switch for sediment pickup function

!---local
real(stm_real) :: c_bar_bed(nvol,nclass)        !< Near bed vaule of mean volumetric sediment concentration
real(stm_real) :: fall_vel(nclass)              !< Settling velocity                       
character :: pick_up_function 
procedure(sediment_hydro_if),pointer :: velocity_non_cohesive 

! set velocity and width
velocity_non_cohesive => sediment_velocity_width

pick_up_function = 'garcia_parker'

if (present(pick_up_flag)) then
     pick_up_function = pick_up_flag
end if

!- initialization to LARGEREAL
call set_sediment_constants

! allocate manning n + width and non_cohesive diameters
call allocate_sediment_parameters(nvol,nclass)

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
!------ set the values of manning's n, width and diameters of grains                     
call set_manning_width_diameter(manning_n,    &
                                width,        &
                                diameter,     &
                                nclass,       &
                                nvol)
                                
! here verfical_net_sediment_flux = settling_vel * (Es - c_bar_sub_b)




   
call deallocate_sediment_static()

end subroutine 

!-----------------------------------------------------------------
 
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
use error_handling
implicit none

real(stm_real),intent(in) :: rouse_num        !< Rouse dimenssionless number  
real(stm_real),intent(in) :: delta_b          !< Relative bed layer thickness = b/H 
real(stm_real),intent(out):: I_1              !< First Einstein integral value

!-- local
real(stm_real) :: ro_l   
real(stm_real) :: ro_r    !right
real(stm_real) :: i_1_l
real(stm_real) :: i_1_r   !right



if (rouse_num > 3.98d0) then
!todo: I am not sure if we need this subroutine in bed load or not 
    call stm_fatal("This is not a Rouse number value for suspended sediment!")
elseif (abs(rouse_num - three)< 0.01d0) then

ro_l = three - 0.05d0
ro_r = three + 0.05d0 

call inside_i_1(i_1_l,delta_b,ro_l)
call inside_i_1(i_1_r,delta_b,ro_r)

I_1 = (i_1_r + i_1_l) / two

         
elseif (abs(rouse_num - two)< 0.01d0) then

ro_l = two - 0.05d0
ro_r = two + 0.05d0 

call inside_i_1(i_1_l,delta_b,ro_l)
call inside_i_1(i_1_r,delta_b,ro_r)

I_1 = (i_1_r + i_1_l) / two
  
     
elseif(abs(rouse_num - one)< 0.01d0) then  

ro_l = one - 0.05d0
ro_r = one + 0.05d0 

call inside_i_1(i_1_l,delta_b,ro_l)
call inside_i_1(i_1_r,delta_b,ro_r)

I_1 = (i_1_r + i_1_l) / two
    
else
   call inside_i_1(I_1,      &
                   delta_b,  &
                   rouse_num)
         
end if

end subroutine

pure subroutine inside_i_1(I_1,      &
                           delta_b,  &
                           rouse_num)   
                            
use stm_precision
implicit none
real(stm_real),intent(in) :: rouse_num        !< Rouse dimenssionless number  
real(stm_real),intent(in) :: delta_b          !< Relative bed layer thickness = b/H 
real(stm_real),intent(out):: I_1              !< First Einstein integral value

I_1   = (rouse_num*pi/sin(rouse_num*pi) - ((one-delta_b)**rouse_num)/(delta_b**(rouse_num-one))     &
         - rouse_num*(((delta_b/(one-delta_b))**(one-rouse_num))  /(one-rouse_num))                 & 
         + rouse_num*(((delta_b/(one-delta_b))**(two-rouse_num))  /(one-rouse_num))                 &
         - rouse_num*(((delta_b/(one-delta_b))**(three-rouse_num))/(one-rouse_num))                 &
         + rouse_num*(((delta_b/(one-delta_b))**(four-rouse_num)) /(one-rouse_num))                 &
         - rouse_num*(((delta_b/(one-delta_b))**(five-rouse_num)) /(one-rouse_num))                 &
         + rouse_num*(((delta_b/(one-delta_b))**(six-rouse_num))  /(one-rouse_num))                 &
         - rouse_num*(((delta_b/(one-delta_b))**(seven-rouse_num))/(one-rouse_num))                 &
         + rouse_num*(((delta_b/(one-delta_b))**(eight-rouse_num))/(one-rouse_num))                 &
         - rouse_num*(((delta_b/(one-delta_b))**(nine-rouse_num)) /(one-rouse_num))                 &
         + rouse_num*(((delta_b/(one-delta_b))**(ten -rouse_num)) /(one-rouse_num)))                &
         * (delta_b**(rouse_num)/((one-delta_b)**rouse_num))
                               

end subroutine 

!------------------------------------------------------------------------

subroutine es_garcia_parker(big_e_sub_s,       &
                            shear_v,           &
                            exp_re_p,          &
                            settling_v,        & 
                            nclass,            &
                            nvol)
use stm_precision
implicit none

!-- arg
integer, intent(in):: nvol                                !< Number of computational volumes in a channel
integer, intent(in):: nclass                              !< Number of non-cohesive sediment grain classes
real(stm_real),intent(out):: big_e_sub_s(nvol,nclass)     !< Dimenssionless rate of entrainment of bed sediment into suspension (i.e., vol entrained sediment/unit bed area/time)                                       
!  big_e_sub_s is in a range of 0.0002 ~ 0.06
real(stm_real),intent(in) :: shear_v(nvol)                !< Shear Velocity
real(stm_real),intent(in) :: exp_re_p(nclass)             !< Explicit particle Reynolds number
real(stm_real),intent(in) :: settling_v(nclass)           !< Settling velocity
!---local
real(stm_real) :: z_u(nvol,nclass)                        !< Captial z sub u a measure for strength of shear stress but it also takes into account the particle size in Garcia notation
real(stm_real), parameter :: cap_a = 1.3d-7               ! Constant value (see ASCE sediment manual no. 110 page 118)
integer :: ivol

do ivol=1,nvol

    z_u(ivol,:) = shear_v(ivol)*(exp_re_p**0.6d0)/settling_v
    where (exp_re_p < 3.5d0)
       z_u(ivol,:) = 0.708d0*shear_v(ivol)*(exp_re_p**0.6d0)/settling_v
    end where 

end do

big_e_sub_s  = cap_a*(z_u**five)/(one + (z_u**five)*cap_a/0.3d0)                                  
                                     
end subroutine

end module 