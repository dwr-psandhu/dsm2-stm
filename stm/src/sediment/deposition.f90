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

!> Sediment Sink module
!>@ingroup sediment
module entrainment

contains
 !> Example deposition function
 subroutine uninitialized_deposition_if( deposition, &
       !todo: check the needed arguments and locals
                                             ncell,   &
                                             nvar,    &
                                             time)
                                         
                                           
         use stm_precision
         use error_handling
         implicit none
         !--- args
         real(stm_real), intent (out):: deposition(ncell,nvar)       !< face flux, lo side
         real(stm_real), intent (in)  ::  time                             !< time
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         !-----locals
         !todo: fill the others
         real(stm_real),parameter :: g_accel =9.80d0 
        !todo: fill this
         deposition =LARGEREAL                                     
                                                 
      
     call stm_fatal("deposition function implemented!")
     
     return
 end subroutine
 
end module
