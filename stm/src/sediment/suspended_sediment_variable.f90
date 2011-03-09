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

!> Routines provide the general calculation for suspended sediment sink/source subroutines.
!> All the constant based drived variables are here
!>@ingroup sediment !todo: test or test sediment

module suspended_utility

contains

!> Calculating particle's settling velocity. NOTE: the subroutine works with SI units.
!> Settling velocity formula based on Leo van Rijn (1984b).
!> The subroutine does not consider particles smaller than 10 microns (fine clay).
!> The smaller particles are assumed to be either part of wash load or to take part in flocs. 
!> The subroutine is for non-cohesive particles.
subroutine settling_velocity(settling_v,&
                             kinematic_viscosity,               &
                             specific_gravity, &
                             diameter,         &
                             g_acceleration,   &
                             function_van_rijn) 
               
use stm_precision
implicit none
!--- arg
real(stm_real),intent(out) :: settling_v !< Settling velocity (m/s)
real(stm_real),intent(in)  :: kinematic_viscosity                !< Kinematic viscosity (m2/sec)
real(stm_real),intent(in)  :: specific_gravity  !< Specific gravity of particle (~2.65)
real(stm_real),intent(in)  :: diameter          !< Particle diameter in meter
real(stm_real),intent(in)  :: g_acceleration    !< Gravitational acceleration (m/sec2)
logical, optional          :: function_van_rijn !< Flag for using van Rijn (1984) formula or Dietrich (1982) the default is Dietrich
!--local
 logical :: van_rijn_flag
 ! todo: I checked the Journal article by Dietrich and these numbers are not the same
 ! todo: I am not sure if the log10 or log e 
 real(stm_real) :: b_1 = 2.891394d0
 real(stm_real) :: b_2 = 0.95296d0
 real(stm_real) :: b_3 = 0.056835d0
 real(stm_real) :: b_4 = 0.002892d0
 real(stm_real) :: b_5 = 0.000245d0

 real(stm_real) :: dimless_fall_velocity
 real(stm_real) :: exp_re_p        !< Explicit Reynols particle number 
 real(stm_real) :: capital_r       !< Submerged specific gravity of sediment particles 
  



 van_rijn_flag = .true.
if (present(function_van_rijn)) then
     van_rijn_flag = function_van_rijn
end if
 
select case (van_rijn_flag)
 
   case (.true.)
        if (diameter > 1.0d-3)    then
            settling_v = 1.1d0*sqrt((specific_gravity - one)*g_acceleration*diameter)
        elseif (diameter > 1.0d-4)  then
            settling_v = (ten*kinematic_viscosity/diameter)*(sqrt(one + (0.01d0*(specific_gravity - one)*g_acceleration*diameter**three)/kinematic_viscosity**two)- one)
        elseif (diameter > 0.9d-7) then
            ! Stokes low
            ! todo: what is the lower limit? for diameter
            settling_v = ((specific_gravity - one)*g_acceleration*diameter**two) /(18.0d0*kinematic_viscosity)
        else
           settling_v = minus * LARGEREAL
           ! todo: the stm_fatal can not be called here because settling velocity is a pure subroutine
        end if 
   case(.false.)
   
        capital_r = specific_gravity - one
        ! Stokes fall velocity
        if ( diameter < 1.0d-4 )    then
            settling_v = (capital_r*g_acceleration*diameter**two) /(18.0d0*kinematic_viscosity)
        else
            call explicit_particle_reynolds_number(exp_re_p,       &
                                                   diameter,       &
                                                   capital_r,      &
                                                   g_acceleration, &
                                                   kinematic_viscosity)
            
            dimless_fall_velocity = exp(- b_1 + b_2*log(exp_re_p)**two - b_3*log(exp_re_p)**three &
                                    - b_4*log(exp_re_p)**four + b_5*log(exp_re_p)**five)
                                    
            settling_v = dimless_fall_velocity * sqrt(capital_r*g_acceleration*diameter)
                 
       end if   
    
end select 

return
end subroutine

!> Calculates the submereged specific gravity
pure subroutine submerged_specific_gravity(capital_r,  &
                                           rho,        &
                                           rho_s)
 use stm_precision
 implicit none
 !-- arguments
real(stm_real),intent(out) :: capital_r      !< Submerged specific gravity of sediment particles     
real(stm_real),intent(in)  :: rho            !< Water density  
real(stm_real),intent(in)  :: rho_s          !< Solid particle density

capital_r = rho_s/rho  - one                                     

return 
end subroutine

!> Calculates the explicit particle Reynolds number
pure subroutine explicit_particle_reynolds_number(exp_re_p,      &
                                                  diameter,      &
                                                  capital_r,     &
                                                  g_acceleration, &
                                                  kinematic_viscosity)
use stm_precision
implicit none
!--- arguments 
real(stm_real),intent(out) :: exp_re_p       !< Explicit particle reynolds number
real(stm_real),intent(in)  :: diameter       !< Particle diameter
real(stm_real),intent(in)  :: capital_r      !< Submerged specific gravity of sediment particles  
real(stm_real),intent(in)  :: g_acceleration !< Gravitational acceleration 
real(stm_real),intent(in)  :: kinematic_viscosity             !< Kinematic viscosity (m2/sec)

exp_re_p = diameter*sqrt(diameter*capital_r*diameter)/kinematic_viscosity

return
end subroutine

!> Calculates particle Reynolds number
pure subroutine particle_reynolds_number(re_p,      &
                                         settling_v,       &
                                         diameter,  &
                                         kinematic_viscosity)

use stm_precision
implicit none
!--- arguments 
real(stm_real),intent(out) :: re_p               !< Particle reynolds number
real(stm_real),intent(in)  :: settling_v  !< Settling velocity
! todo: Eli; do we need settling_velocity and diameter as settling_velocity(nvar)?
real(stm_real),intent(in)  :: diameter           !< Particle diameter
real(stm_real),intent(in)  :: kinematic_viscosity                !< Kinematic viscosity (m2/sec)                            
 
 re_p = settling_v*diameter/kinematic_viscosity
 
return
end subroutine

!> Calculates dimensionless particle diameter
pure subroutine dimless_particle_diameter(d_star,                 &
                                          g_acceleration,         &
                                          diameter,               &
                                          kinematic_viscosity,    &
                                          capital_r)

use stm_precision
implicit none
!--- arguments 
real(stm_real),intent(out) :: d_star              !< Dimensionless particle diameter
real(stm_real),intent(in)  :: g_acceleration      !< Gravitational acceleration 
real(stm_real),intent(in)  :: diameter            !< Particle diameter
real(stm_real),intent(in)  :: kinematic_viscosity !< Kinematic viscosity of water sediment mixture (m2/sec)                            
real(stm_real),intent(in)  :: capital_r           !< Submerged specific gravity of sediment particles     

d_star = diameter*(capital_r*g_acceleration/(kinematic_viscosity**two))**third
 
return
end subroutine

!> Calculates critical shields parameter based on Yalin (1972) formula
!> See van Rijn book equation (4.1.11)
! todo: add Parker formula here
pure subroutine critical_shields_parameter(cr_shields_prmtr,   &
                                           d_star)
                                           
    use stm_precision
    implicit none
    !--- arguments  
    real(stm_real),intent(out):: cr_shields_prmtr !< Critical Shields parameter                                      
    real(stm_real),intent(in) :: d_star           !< Dimensionless particle diameter

    if    (d_star > 150.0d0) then
        cr_shields_prmtr = 0.055d0   
    elseif (d_star > 20.0d0) then
        cr_shields_prmtr = 0.013d0*d_star**0.29d0 
    elseif (d_star > 20.0d0) then
        cr_shields_prmtr = 0.04d0*d_star**(-0.1d0)
    elseif (d_star > 4.0d0)  then
        cr_shields_prmtr = 0.14d0*d_star**(-0.64d0)
    elseif (d_star > one)    then
        cr_shields_prmtr = 0.24d0/d_star
    else
        cr_shields_prmtr = minus*LARGEREAL
        ! the number set here to prevent bad input (stm_fatal can not be called)
    end if                                                      
                                           
return
end subroutine


end module
