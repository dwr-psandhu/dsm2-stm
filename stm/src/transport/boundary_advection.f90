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

!> boundary advection interface to be fulfilled by driver or application
!>@ingroup transport
module boundary_advection_module
 !> Calculate boundary advection
 interface
       !> Generic interface for calculating BC advection routine that should be fulfilled by
       !> client programs
       subroutine boundary_advective_flux(flux_lo,  &
                                        flux_hi,    &
                                        conc_lo,    &
                                        conc_hi,    &
                                        flow_lo,    &
                                        flow_hi,    &
                                        ncell,      &
                                        nvar,       &
                                        time,       &
                                        dt,         &
                                        dx)
        
         use stm_precision
         
         implicit none
         !--- args          
        integer,intent(in)  :: ncell  !< Number of cells
        integer,intent(in)  :: nvar   !< Number of variables
        ! todo: check the intents
        real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
        real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
        real(stm_real),intent(in)    :: flow_lo(ncell,nvar)     !< flow on lo side of cells centered in time
        real(stm_real),intent(in)    :: flow_hi(ncell,nvar)     !< flow on hi side of cells centered in time
        real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
        real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
        real(stm_real), intent (in)  :: time                    !< Current time
        real(stm_real), intent (in)  :: dx                      !< Spatial step  
        real(stm_real), intent (in)  :: dt                      !< Time step     
      
      
       end subroutine boundary_advective_flux
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the advection boundary condition 
 ! todo: check here
 procedure(boundary_advective_flux),pointer :: boundary_advection  => null()


 contains
 
 !> Example uninitialize that prints an error and bails
 subroutine uninitialized_advection_bc(flux_lo,     &
                                        flux_hi,    &
                                        conc_lo,    &
                                        conc_hi,    &
                                        flow_lo,    &
                                        flow_hi,    &
                                        ncell,      &
                                        nvar,       &
                                        time,       &
                                        dt,         &
                                        dx)
                                         
     use stm_precision 
     use error_handling
     
        implicit none
         !--- args          
        integer,intent(in)  :: ncell  !< Number of cells
        integer,intent(in)  :: nvar   !< Number of variables
        ! todo: check the intents
        real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
        real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
        real(stm_real),intent(in)    :: flow_lo(ncell,nvar)     !< flow on lo side of cells centered in time
        real(stm_real),intent(in)    :: flow_hi(ncell,nvar)     !< flow on hi side of cells centered in time
        real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
        real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
        real(stm_real), intent (in)  :: time                    !< Current time
        real(stm_real), intent (in)  :: dx                      !< Spatial step  
        real(stm_real), intent (in)  :: dt                      !< Time step     
      
     call stm_fatal("Boundary not implemented in advection!")
     
     return
 end subroutine 
 
 !> Example advective flux that imposes Neumann boundaries with zero flux at
 !> both ends of the channel.
 subroutine neumann_advective_flux(flux_lo,     &
                                    flux_hi,    &
                                    conc_lo,    &
                                    conc_hi,    &
                                    flow_lo,    &
                                    flow_hi,    &
                                    ncell,      &
                                    nvar,       &
                                    time,       &
                                    dt,         &
                                    dx)
     
     use stm_precision
        implicit none
         !--- args          
        integer,intent(in)  :: ncell  !< Number of cells
        integer,intent(in)  :: nvar   !< Number of variables
        ! todo: check the intents
        real(stm_real),intent(inout) :: flux_lo(ncell,nvar)     !< flux on lo side of cell, time centered
        real(stm_real),intent(inout) :: flux_hi(ncell,nvar)     !< flux on hi side of cell, time centered
        real(stm_real),intent(in)    :: flow_lo(ncell,nvar)     !< flow on lo side of cells centered in time
        real(stm_real),intent(in)    :: flow_hi(ncell,nvar)     !< flow on hi side of cells centered in time
        real(stm_real),intent(in)    :: conc_lo(ncell,nvar)     !< concentration extrapolated to lo face
        real(stm_real),intent(in)    :: conc_hi(ncell,nvar)     !< concentration extrapolated to hi face
        real(stm_real), intent (in)  :: time                    !< Current time
        real(stm_real), intent (in)  :: dx                      !< Spatial step  
        real(stm_real), intent (in)  :: dt                      !< Time step     
      
      flux_lo(1,nvar) = zero
      flux_hi(ncell,nvar) = zero
        
     ! todo: implement and test
     ! todo: add non trivial cases
     return
 end subroutine
 
end module