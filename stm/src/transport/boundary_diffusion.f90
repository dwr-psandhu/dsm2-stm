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

!> boundary diffusive flux interface to be fulfilled by driver or application
!>@ingroup transport
module boundary_diffusion
 !> Calculate boundary diffusion flux
 interface
       !> Generic interface for calculating BC that should be fulfilled by
       !> client programs
       subroutine boundary_diffusive_flux_if(diffusive_flux_lo, &
                                             diffusive_flux_hi, &
                                             conc,              &
                                             ncell,             &
                                             nvar,              &
                                             time)
        
         use stm_precision
         implicit none
         !--- args
         integer,intent(in)  :: ncell                                   !< number of cells
         integer,intent(in)  :: nvar                                    !< number of variables
         real(stm_real),intent(inout) :: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real),intent(inout) :: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real),intent(in)  :: time                             !< time
         real(stm_real),intent(in)  :: conc(ncell,nvar)                 !< concentration 
       
       end subroutine boundary_diffusive_flux_if
 end interface

 !> This pointer should be set by the driver or client code to specify the 
 !> treatment at the boundaries
 procedure(boundary_diffusive_flux_if),pointer :: boundary_diffusion_flux  => null()


 contains
 
 !> Example diffusive flux that prints an error and bails
 subroutine uninitialized_diffusive_flux(diffusive_flux_lo, &
                                         diffusive_flux_hi, &
                                         conc,              &
                                         ncell,             &
                                         nvar,              &
                                         time)
                                         
     use stm_precision 
     use error_handling
     
     implicit none
     !--- args
         integer,intent(in)  :: ncell                                   !< Number of cells
         integer,intent(in)  :: nvar                                    !< Number of variables
         real(stm_real),intent(inout) :: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real),intent(inout) :: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real),intent(in)  :: time                             !< time
         real(stm_real),intent(in)  :: conc(ncell,nvar)                 !< concentration 
     call stm_fatal("boundary not implemented")
     
     return
 end subroutine 
 
 !> Example diffusive flux that imposes Neumann boundaries with zero flux at
 !> both ends of the channel.
 subroutine neumann_no_flow_diffusive_flux(diffusive_flux_lo, &
                                             diffusive_flux_hi, &
                                             conc,              &
                                             ncell,             &
                                             nvar,              &
                                             time)
     use stm_precision
     implicit none
     !--- args
         integer,intent(in)  :: ncell                                   !< Number of cells
         integer,intent(in)  :: nvar                                    !< Number of variables
         real(stm_real),intent(inout) :: diffusive_flux_lo(ncell,nvar)  !< face flux, lo side
         real(stm_real),intent(inout) :: diffusive_flux_hi(ncell,nvar)  !< face flux, hi side
         real(stm_real),intent(in)  :: time                             !< time
         real(stm_real),intent(in)  :: conc(ncell,nvar)                 !< concentration 
     
     ! todo: add other BC 
     ! neumann default
     diffusive_flux_lo(1,:) = zero
     diffusive_flux_hi(ncell,:) =zero
     
     
     ! todo: implement and test
     return
 end subroutine
 
end module
