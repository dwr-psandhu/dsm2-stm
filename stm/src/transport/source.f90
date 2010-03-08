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

!> Source interface to be fulfilled by driver or application
!>@ingroup transport
module source_module
 !> Calculate source
 interface compute_source
   !> Generic interface for calculating source that should be fulfilled by
   !> client programs
   
   ! todo: check if the source must be inout or just out
   subroutine compute_source_if(source, & 
                                conc,   &
                                area,   &
                                flow,   &
                                ncell,  &
                                nvar,   &
                                time)
     use stm_precision
     implicit none
     !--- args
     integer,intent(in)  :: ncell                       !< Number of cells
     integer,intent(in)  :: nvar                        !< Number of variables
     real(stm_real),intent(inout) :: source(ncell,nvar) !< cell centered source 
     real(stm_real),intent(in)  :: conc(ncell,nvar)     !< Concentration
     real(stm_real),intent(in)  :: area(ncell)          !< area at source     
     real(stm_real),intent(in)  :: flow(ncell)          !< flow at source location
     real(stm_real),intent(in)  :: time                 !< flow at source location
     
   end subroutine compute_source_if
 end interface
 !> This pointer should be set by the driver or client code to specify the 
 !> source term 
 procedure(compute_source_if),pointer :: source => null() !no_source

  contains
 
 subroutine no_source(source,   & 
                        conc,   &
                        area,   &
                        flow,   &
                        ncell,  &
                        nvar,   &
                        time)
                                         
     use stm_precision 
     use error_handling
     
     
         implicit none
     !--- args
     integer,intent(in)  :: ncell                      !< Number of cells
     integer,intent(in)  :: nvar                       !< Number of variables
     real(stm_real),intent(inout) :: source(ncell,nvar)!< cell centered source 
     real(stm_real),intent(in)  :: conc(ncell,nvar)    !< Concentration
     real(stm_real),intent(in)  :: area(ncell)         !< area at source     
     real(stm_real),intent(in)  :: flow(ncell)         !< flow at source location
     real(stm_real),intent(in)  :: time                !< flow at source location
     
     call stm_fatal("No Source!")
     
     return
 end subroutine no_source
 
  subroutine linear_decay_source(source,    & 
                                    conc,   &
                                    area,   &
                                    flow,   &
                                    ncell,  &
                                    nvar,   &
                                    time)
                                         
     use stm_precision 
     use error_handling
     
     
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
     return
 end subroutine linear_decay_source
 
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
 
 subroutine non_cohesive_source(source, & 
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
     return
 end subroutine
 
end module
 
