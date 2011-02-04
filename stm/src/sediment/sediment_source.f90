!!<license>
!!    Copyright (C) 1996, 1997, 1998, 2001, 2007, 2009 State of California,
!!    Department of Water Resources.
!!    This file is part of DSM2.
!!
!!    The Delta Simulation Model 2 (DSM2) is free software: 
!!    you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    DSM2 is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with DSM2.  If not, see <http://www.gnu.org/licenses>.
!!</license>
!
!!> Sediment Sink/Source interface to be fulfilled by driver or application
!!>@ingroup transport
!module sediment_source_sink
! 
! 
! 
! 
! abstract interface
!       !> Generic interface for calculating entrainment that should be fulfilled by
!       !> client programs
!       subroutine entrainment_if(entrainment, &
!       !todo: check the needed arguments and locals
!                                             ncell,   &
!                                             nvar,    &
!                                             time)
!        
!         
!        
!         use stm_precision
!         implicit none
!         !--- args
!         real(stm_real), intent (out):: entrainment(ncell,nvar)       !< face flux, lo side
!         real(stm_real), intent (in)  ::  time                             !< time
!         integer, intent(in)  :: ncell                                   !< number of cells
!         integer, intent(in)  :: nvar                                    !< number of variables
!         !-----locals
!         !todo: fill the others
!         real(stm_real),parameter :: g_accel =9.80d0 
!
!
!       end subroutine entrainment_if
! end interface
!
! !> This pointer should be set by the driver or client code to specify the 
! !> treatment of entrainment function
! procedure(entrainment_if),pointer :: entrainment_source  => null()
!
!!=========================================================
!
! abstract interface
!       !> Generic interface for calculating deposition that should be fulfilled by
!       !> client programs
!       subroutine deposition_if(deposition, &
!       !todo: check the needed arguments and locals
!                                             ncell,   &
!                                             nvar,    &
!                                             time)
!                                                      
!                                                      
!        
!         use stm_precision
!         implicit none
!         !--- args
!                                       
!        integer, intent (in) :: ncell                                               !< Number of cells
!        integer, intent (in) :: nvar                                                !< Number of variables
!        real(stm_real),intent (out):: deposition(ncell,nvar)                       !< Values of the coefficients below diagonal in matrix
!        real(stm_real), intent (in)  :: time                                        !< Current time
!         !todo: fill the others
!         real(stm_real),parameter :: g_accel =9.80d0 
!       
!
!       end subroutine deposition_if
! end interface
!
! !> This pointer should be set by the driver or client code to specify the 
! !> treatment of deposition
! ! todo: 
! procedure(deposition_if),pointer :: deposition_sink  => null()
!!===============================================================
!abstract interface 
!  subroutine cohesive_source_if(cohesive_source, & 
!                                conc,   &
!                                area,   &
!                                flow,   &
!                                ncell,  &
!                                nvar,   &
!                                time)
!                                
!     use stm_precision
!    
!          implicit none
!     !--- args
!     integer,intent(in):: ncell                        !< Number of cells
!     integer,intent(in):: nvar                         !< Number of variables
!     real(stm_real),intent(inout) :: cohesive_source(ncell,nvar)!< cell centered source 
!     real(stm_real),intent(in)  :: conc(ncell,nvar)    !< Concentration
!     real(stm_real),intent(in)  :: area(ncell)         !< area at source     
!     real(stm_real),intent(in)  :: flow(ncell)         !< flow at source location
!     real(stm_real),intent(in)  :: time                !< flow at source location
! 
!     ! todo: implement and test
!     ! todo: add other kind of sources
!    
! end subroutine cohesive_source_if
! 
! end interface
! 
! procedure(cohesive_source_if),pointer :: cohesive_source  => null()
!!================================================= 
!
! contains 
! 
! !> Compute sediment tranport source
! !> This subroutine orchestrates the calculation of
! !> entrainement, deposition and cohesion by calling
! !> functions for each and adding them up. This is just a
! !> skeleton of the sediment source...the actual
! !> choice of entrainment (etc) functions can be set with a pointer
! 
! ! todo
! !subroutine sediment_source(sed_source,conc,ncell,nvar)
!
! use stm_precision
! use entrainment
!! todo :
!! use deposition
!          
! implicit none
! ! ---arg
!  
! integer, intent(in)  :: ncell      !< number of cells
! integer, intent(in)  :: nvar       !< number of constituents  
! real(stm_real),intent(out) :: sed_source(ncell,nvar) !< source of sediment 
! 
! ! todo: this will work, but you probably want to do it elsewhere
!! entrainment => van_rijn_source        
!! deposition_sink => uninitialized_deposition_if      
!! cohesive_source => cohesive_source
!
!! ! method #1 separate data structures 
!! call entrainment(entrain,...)  ! entrain is intent out
!! call deposition_sink(deposit )
!! call cohesive_source(cohere...)
!! sed_source = entrain + deposit + cohere
!! 
!! ! method #2 each function increments rather than sets
!! call entrainment(sed_source,...)  ! sed_source is intent inout
!! call deposition_sink(sed_source,...)
!! call cohesive_source(sed_source,...)
! 
! !todo
! !end subroutine sediment_source
! 
!end module
