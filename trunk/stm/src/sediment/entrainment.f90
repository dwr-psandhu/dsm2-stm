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

!> Sediment Source  module
!>@ingroup sediment
module entrainment

contains
!> Garcia_Parker Entrainment function that prints an error and bails

 subroutine garcia_parker_source(entrainment, &
       !todo: check the needed arguments and locals
                                             ncell,   &
                                             nvar,    &
                                             time)
        
         
        
         use stm_precision
         use error_handling
         implicit none
         !--- args
         real(stm_real), intent (out):: entrainment(ncell,nvar)       !< face flux, lo side
         real(stm_real), intent (in)  ::  time                             !< time
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         !-----locals
         !todo: fill the others
         real(stm_real),parameter :: g_accel =9.80d0 
        !todo: fill this
         entrainment =LARGEREAL
        
     call stm_fatal("garcia_parker entrainment was not included")
     
     return
 end subroutine 
 
 !> Smith McLean entrainment function
 subroutine smith_mclean_source(entrainment, &
       !todo: check the needed arguments and locals
                                             ncell,   &
                                             nvar,    &
                                             time)
        
         
        
         use stm_precision
         use error_handling
         implicit none
         !--- args
         real(stm_real), intent (out):: entrainment(ncell,nvar)       !< face flux, lo side
         real(stm_real), intent (in)  ::  time                             !< time
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         !-----locals
         !todo: fill the others
         !todo: we should move this
         real(stm_real),parameter :: gravity =9.80d0 
        !todo: fill this
         entrainment =LARGEREAL
        
     call stm_fatal("smith_mclean entrainment was not included")
     
     ! todo: implement and test
     return
 end subroutine
 !> van Rijn entrainment function
 subroutine van_rijn_source(entrainment, &
       !todo: check the needed arguments and locals
                                             ncell,   &
                                             nvar,    &
                                             time)
        
         
        
         use stm_precision
         use error_handling
         implicit none
         !--- args
         real(stm_real), intent (out):: entrainment(ncell,nvar)       !< face flux, lo side
         real(stm_real), intent (in)  ::  time                             !< time
         integer, intent(in)  :: ncell                                   !< number of cells
         integer, intent(in)  :: nvar                                    !< number of variables
         !-----locals
         !todo: fill the others
         real(stm_real),parameter :: g_accel =9.80d0 
        !todo: fill this
         entrainment =LARGEREAL
        
     call stm_fatal("van Rijn entrainment was not included")
     
     ! todo: implement and test
     return
 end subroutine
 
   subroutine cohesive_source(source, & 
                                conc,   &
                                area,   &
                                flow,   &
                                ncell,  &
                                nvar,   &
                                time)
                                
     use stm_precision
    
          implicit none
     !--- args
     integer,intent(in):: ncell                        !< Number of cells
     integer,intent(in):: nvar                         !< Number of variables
     real(stm_real),intent(inout) :: source(ncell,nvar)!< cell centered source 
     real(stm_real),intent(in)  :: conc(ncell,nvar)    !< Concentration
     real(stm_real),intent(in)  :: area(ncell)         !< area at source     
     real(stm_real),intent(in)  :: flow(ncell)         !< flow at source location
     real(stm_real),intent(in)  :: time                !< flow at source location
 
     ! todo: implement and test
     ! todo: add other kind of sources
     return
 end subroutine
 
 end module
