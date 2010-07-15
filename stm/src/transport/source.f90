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
!> Specifically, the user must set the "source" pointer (a variable in this module)
!> to a subroutine that meets the compute_source interface.
!>@ingroup transport
module source_sink
 !> Calculate source
 interface
   !> Generic interface for calculating source that should be fulfilled by
   !> client programs. The source is calculated using concentration,
   !> and produces a source increment appropriate for primitives.
   !> If you want to increment a conservative variable with this source,
   !> call prim_increment_to_cons()
   
   ! todo: check if the source must be inout or just out
   ! todo: create prim_increment_to_cons()
   subroutine source_if(source, & 
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
     real(stm_real),intent(in)  :: area(ncell)          !< Cell centered area at source     
     real(stm_real),intent(in)  :: flow(ncell)          !< flow at source location
     real(stm_real),intent(in)  :: time                 !< time
   end subroutine source_if
 end interface
 
 !>\file
 !>\var source_sink::source 
 !>\brief Pointer to source term
 !>This pointer should be set by the driver or client code
 procedure(source_if), pointer :: compute_source => null() 

 contains
 
 !> Zero source implementation.
 !> This source term adds nothing (e.g., for conservative constituents)
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
     real(stm_real),intent(in)  :: time                !< time
     ! todo: source must be primitive variable
    source = zero
     
     return
 end subroutine 
 !==================================

 
end module
 
