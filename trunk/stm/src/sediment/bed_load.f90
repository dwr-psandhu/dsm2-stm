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
!> Contains different relations of bedload transport
module bed_load


contains 
!> The subroutine provides unit bedload per volume
!> q_sub_b has the dimension L2/T, Q_sub_b = integral over width q_sub_b, [L3/T]
subroutine bedload(q_sub_b,       &
                   velocity,      &
                   area,          &
                   nclass,        &
                   nvol,          &
                   !!!!!!! todo: fill other ones
                   bedload_func)                   

use stm_precision
use suspended_utility

implicit none

integer,intent(in) :: nclass                  !< Number of grain classes in bedload transport
integer,intent(in) :: nvol                    !< Number of volumes in a channel
real(stm_real),intent(in) :: velocity(nvol)                !< Velocity
real(stm_real),intent(in) :: area(nvol)                    !< Area
real(stm_real),intent(out)::q_sub_b(nvol,nclass)           !< Volumetric bedload transport rate
character, intent(in),optional :: bedload_func             !< Bedload relation 
!---
character,parameter :: default ='meyer_peter_muller'
character :: relation
 

if(present(bedload_func))then
  relation = bedload_func
end if 
 




return
end subroutine 


!> Subroutine for calculating the q_b based on q_* (see Garcia, 2008, page 70)
subroutine volumetric_bedload_transport_rate(q_sub_b,                &
                                             einstein_bedload_num,   &
                                             diameter,               &
                                             capital_r,              &
                                             gravity,                &
                                             nclass,                 &
                                             nvol)

use stm_precision
implicit none
                                             
integer,intent(in) :: nclass                                  !< Number of sediment classes 
integer,intent(in) :: nvol                                    !< Number of volume
real(stm_real),intent(in) ::gravity                           !< Gravity
real(stm_real),intent(in) ::capital_r                         !< Submerged specific gravity of sediment particles  
real(stm_real),intent(in) ::diameter(nclass)                  !< Particle diameter
real(stm_real),intent(in) ::einstein_bedload_num(nvol,nclass) !< Dimensionless belload transport rate 
real(stm_real),intent(out)::q_sub_b(nvol,nclass)              !< Volumetric bedload transport rate
!--- local
integer:: iclass

                
do iclass=1,nclass                             
    q_sub_b(:,iclass) = diameter * einstein_bedload_num(:,iclass) * sqrt(gravity*capital_r*diameter)
end do 


return
end subroutine 

end module 