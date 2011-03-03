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
!>@ingroup sediment !todo: test or test sediment

module suspended_utility

contains

!> Calculating particle's settling velocity. NOTE: the subroutine works with SI units.
!> Settling velocity formula based on Leo van Rijn (1984b).
!> The subroutine does not consider particles smaller than 0.9 microns (fine clay).
!> The smaller particles are assumed to be either part of wash load or to take part in flocs. 
pure subroutine settling_velocity_van_rijn(w_s,              &
                                           nu,               &
                                           specific_gravity, &
                                           diameter,         &
                                           g_acceleration) 
               
use stm_precision
implicit none
!--- arg
real(stm_real),intent(out) :: w_s              !< Settling velocity (m/s)
real(stm_real),intent(in)  :: nu               !< Kinematic viscosity (m2/sec)
real(stm_real),intent(in)  :: specific_gravity !< Specific gravity of particle (~2.65)
real(stm_real),intent(in)  :: diameter         !< Particle diameter in meter
real(stm_real),intent(in)  :: g_acceleration   !< Gravitational acceleration (m/sec2)

if ( diameter > 1.0d-3 )    then
    w_s = 1.1d0*sqrt((specific_gravity - one)*g_acceleration*diameter)
elseif (diameter > 1.0d-4)  then
    w_s = (ten*nu/diameter)*(sqrt(one + (0.01d0*(specific_gravity - one)*g_acceleration*diameter**three)/nu**two)- one)
elseif (diameter > 0.9d-6 ) then
    ! Stokes low
    w_s = ((specific_gravity - one)*g_acceleration*diameter**two) /(18.0d0*nu)
else
   w_s = minus * LARGEREAL
   ! todo: the stm_fatal can not be called here because settling velocity is a pure subroutine
end if 

return
end subroutine

!> Calculates the submereged specific gravity
pure subroutine submerged_specific_gravity(capital_r,  &
                                           rho,        &
                                           rho_s)
 use stm_precision
 implicit none
 !-- argument
real(stm_real),intent(out) :: capital_r      !< Submerged specific gravity of sediment particles     
real(stm_real),intent(in)  :: rho            !< Water density  
real(stm_real),intent(in)  :: rho_s          !< Solid particle density

capital_r = rho_s/rho  - one                                     

return 
end subroutine


end module
